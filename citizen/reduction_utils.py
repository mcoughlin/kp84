#!/usr/bin/env python3
import math
import os, copy
import numpy as np
import datetime
import astropy.io.ascii as asci
from astropy.table import Table
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from copy import deepcopy
from scipy.ndimage import median_filter

import image_registration
from skimage.feature import register_translation
from photutils import DAOStarFinder

import PythonPhot as pp

import matplotlib
import matplotlib.pyplot as plt
fs = 14
matplotlib.rcParams.update({'font.size': fs})


def filter2filtstr(myfilter):
    if myfilter == "SDSS g'":
        filtstr = "sg"
    elif myfilter == "SDSS r'":
        filtstr = "sr"
    else:
        print ("unknown filter! modify your code")
        exit(0)
    return filtstr


def mylinear_fit(x, y, yerr, npar = 2):
    '''
    Ref: 
        1. Numerical Recipes, 3rd Edition, p745, 781 - 782
        2. http://web.ipac.caltech.edu/staff/fmasci/ztf/ztf_pipelines_deliverables.pdf, p38
    '''
    assert len(x) == len(y)
    assert len(y) == len(yerr)
    
    Sx = np.sum(x)
    Sy = np.sum(y)
    Sxy = np.sum(x * y)
    Sxx = np.sum(x**2)
    N = len(x)
    
    Sx_sigma = np.sum(x * yerr**2)
    Sxx_sigma = np.sum(x**2 * yerr**2)
    S_sigma = np.sum(yerr**2)
    
    if npar==1:
        Fpsf = Sxy / Sxx
        e_Fpsf = np.sqrt(Sxx_sigma) / Sxx
        a = 0
    elif npar==2:
        Fpsf = (N * Sxy - Sx * Sy) / (N * Sxx - Sx**2)
        a = (Sxx * Sy - Sx * Sxy) / (N * Sxx - Sx**2)
        e_Fpsf = np.sqrt(N**2*Sxx_sigma - 2*N*Sx*Sx_sigma + Sx**2*S_sigma) / (N * Sxx - Sx**2)
    # x_mean = np.mean(x)
    # y_mean = np.mean(y)
    # pearson_r = np.sum( (x - x_mean) * (y - y_mean) ) / np.sqrt(np.sum( (x - x_mean)**2 )) / np.sqrt(np.sum( (y - y_mean)**2 ))
    return Fpsf, e_Fpsf, a


def get_predictioon_xy(ind, xobjs, cs):
    nchuck = int(np.ceil(len(cs)/150))
    neach = int(np.ceil(len(cs)/nchuck))
    xpredict = np.zeros(len(cs))
    for i in range(nchuck):
        csub = cs[i*neach:(i+1)*neach]
        xsub = xobjs[i*neach:(i+1)*neach]
        indsub = ind[i*neach:(i+1)*neach]
        #print (np.sum(indsub))
        if np.sum(indsub)>100:
            xslope, tmp1x, xintercept = mylinear_fit(csub[indsub==1], xsub[indsub==1], np.ones(np.sum(indsub==1)), npar = 2)
        else:
            # leave one out at a time
            xslope, tmp1x, xintercept = mylinear_fit(csub, xsub, np.ones(len(csub)), npar = 2)
            xpict = csub * xslope + xintercept
            resi = xsub - xpict
            std = np.std(resi)
            #plt.plot(csub, xsub, '.')
            csub_ = deepcopy(csub)
            xsub_ = deepcopy(xsub)
            #count = 0
            while std > 1.2:
                aresi = abs(resi)
                ind_ = np.argsort(aresi)[:-1]
                csub_ = csub_[ind_]
                xsub_ = xsub_[ind_]
                #if count>20:
                #    plt.plot(csub_, xsub_, '.')
                xslope, tmp1x, xintercept = mylinear_fit(csub_, xsub_, np.ones(len(csub_)), npar = 2)
                xpict = csub_ * xslope + xintercept
                #if count>20:
                #    plt.plot(csub_, xpict)
                resi = xsub_ - xpict
                std = np.std(resi)
                #count+=1
        xpredict[i*neach:(i+1)*neach] = csub * xslope + xintercept
    return xpredict
            

def get_reference(scienceimage, cat):
    # plt.hist(cat[:,20])
    cat = cat[cat[:,20]<20] # cut sources with too high fwhm
    # plt.hist(cat[:,14])
    cat = cat[cat[:,14]<4] # Profile RMS along major axis 
    # plt.hist(cat[:,15])
    cat = cat[cat[:,15]<3] # Profile RMS along minor axis   
    
    exts = np.unique(cat[:,22])
    
    sciHDU = fits.open(scienceimage)
    for i in range(len(sciHDU)):
        if i==0:
            continue
        hdu = sciHDU[i]
        header = hdu.header
        
        catind = cat[:,22]==i
            
        plt.figure()
        plt.imshow(subdata, origin = "lower")
        plt.plot(xsub_guess-1, ysub_guess-1, 'rx')
        plt.plot(xsub_tuned-1, ysub_tuned-1, 'r+')
        plt.savefig('/Users/yuhanyao/Desktop/kped_tmp/20191117/ZTFJ01395245/figures/%d.png'%i)
        plt.close()


def register_images(fitsfile, shiftfile, xyframe, x, y, maxdist=10.):
    tbshift = asci.read(shiftfile)

    x1 = x + tbshift["xshift"][xyframe] # get the (x, y) in extension 1
    y1 = y + tbshift["yshift"][xyframe]
    
    print ("  Initial registration...")
    procHDU = fits.open(fitsfile)
    regis1HDU = fits.HDUList()
    ndata = procHDU[1].data.shape[0]
    for i in range(len(procHDU)):
        hdu = procHDU[i]
        if i==0:
            regis1HDU.append(hdu)
        else:
            x_frame = x1 - tbshift["xshift"][i]
            y_frame = y1 - tbshift["yshift"][i]
            if x_frame<1 or x_frame>(ndata) or y_frame<1 or y_frame>(ndata):
                print ("    Target outside of image! Remove extention %d!"%i)
            else:
                hdu.header["x_guess"] = (x_frame, "Initial guess of the object x Pixel Coordinate")
                hdu.header["y_guess"] = (y_frame, "Initial guess of the object y Pixel Coordinate")
                regis1HDU.append(hdu)
    
    print ("  Tuning registration...")
    
    regis2HDU = fits.HDUList()
    for i in range(len(regis1HDU)):
        hdu = regis1HDU[i]
        if i==0:
            regis2HDU.append(hdu)
        else:
            x_guess = hdu.header["x_guess"]
            y_guess = hdu.header["y_guess"]
            data = hdu.data
            # print (i, x_guess, y_guess)
            xstart = int(np.floor(x_guess))-30
            xend = int(np.ceil(x_guess))+30
            ystart = int(np.floor(y_guess))-30
            yend = int(np.ceil(y_guess))+30
            if xstart < 0:
                xstart = 0
            if xend > ndata:
                xend = ndata
            if ystart < 0:
                ystart= 0
            if yend > ndata:
                yend = ndata
            findsourceflag = 1
            subdata = data[ystart:yend, xstart:xend]
            xsub_guess = x_guess - xstart # in pixel convention (start from 1)
            ysub_guess = y_guess - ystart # in pixel convention (start from 1)
            mean, median, std = sigma_clipped_stats(subdata, sigma=3.0)  
            daofind = DAOStarFinder(fwhm=3.0, threshold=3.*std)  
            sources = daofind(subdata)   # in python convention (start from 0)
            if sources is None:
                findsourceflag = 0
            else:
                xpixel = sources["xcentroid"].data+1 # convert to pixel convention
                ypixel = sources["ycentroid"].data+1 # convert to pixel convention
                distances = np.sqrt((xpixel - xsub_guess)**2 + (ypixel - ysub_guess)**2)
                ind = distances<maxdist
                if np.sum(ind)>0:
                    distances = distances[ind]
                    xpixel = xpixel[ind]
                    ypixel = ypixel[ind]
                    ind = np.argsort(distances)[0]
                    xsub_tuned = xpixel[ind]
                    ysub_tuned = ypixel[ind]
                else:
                    findsourceflag = 0
                # print (i, "did not find 3 sigma")
            if findsourceflag==0:
                if i!=1:
                    lastheader = regis2HDU[-1].header
                    xoff = lastheader["X_DAOOBJ"] - lastheader["X_GUESS"]
                    yoff = lastheader["Y_DAOOBJ"] - lastheader["Y_GUESS"]
                    xsub_tuned = xsub_guess + xoff
                    ysub_tuned = ysub_guess + yoff
                else:
                    xsub_tuned = xsub_guess
                    ysub_tuned = ysub_guess
            """
            plt.figure()
            plt.imshow(subdata, origin = "lower")
            plt.plot(xsub_guess-1, ysub_guess-1, 'rx')
            plt.plot(xsub_tuned-1, ysub_tuned-1, 'r+')
            plt.savefig('/Users/yuhanyao/Desktop/kped_tmp/20191117/ZTFJ01395245/figures/%d.png'%i)
            plt.close()
            """
            x_tuned = xsub_tuned + xstart
            y_tuned = ysub_tuned + ystart
            hdu.header["DAOFINDFLAG"] = (findsourceflag, "DAOStarFinder flag")
            hdu.header["X_DAOOBJ"] = (x_tuned, "object x Pixel Coordinate found by DAOStarFinder")
            hdu.header["Y_DAOOBJ"] = (y_tuned, "object y Pixel Coordinate found by DAOStarFinder")
            regis2HDU.append(hdu)
            
    ind = np.array([hdu.header["DAOFINDFLAG"] for hdu in regis2HDU[1:]])
    xobjs = np.array([hdu.header["X_DAOOBJ"] for hdu in regis2HDU[1:]])
    yobjs = np.array([hdu.header["Y_DAOOBJ"] for hdu in regis2HDU[1:]])
    xinits = np.array([hdu.header["X_GUESS"] for hdu in regis2HDU[1:]])
    yinits = np.array([hdu.header["Y_GUESS"] for hdu in regis2HDU[1:]])
    cs = np.arange(len(regis2HDU[1:]))
    
    xpredict = get_predictioon_xy(ind, xobjs, cs)
    ypredict = get_predictioon_xy(ind, yobjs, cs)
    
    # print (np.std(xpredict - xobjs), np.std(ypredict - yobjs))
    
    indx = abs(xpredict - xobjs) < 1.5
    indy = abs(ypredict - yobjs) < 1.5
    kk_x, ekk, bb_x = mylinear_fit(xinits[indx], xobjs[indx], np.ones(np.sum(indx)), npar = 2)
    kk_y, ekk, bb_y = mylinear_fit(yinits[indy], yobjs[indy], np.ones(np.sum(indy)), npar = 2)
    xfinals = kk_x * xinits + bb_x
    yfinals = kk_y * yinits + bb_y
    
    plt.figure(figsize=(12, 12))
    ax1 = plt.subplot(211)
    ax1.plot(cs[ind==0], xobjs[ind==0], 'r.', alpha = 0.5)
    ax1.plot(cs[ind==1], xobjs[ind==1], 'b.', alpha = 0.5)
    ax1.plot(cs, xpredict, 'k--', alpha = 0.5)
    ax1.plot(cs, xinits, label="Initial")
    # plt.plot(cs[indx], xinits[indx])
    ax1.plot(cs[indx], xobjs[indx], label = "DAOFind")
    ax1.plot(cs, kk_x * xinits + bb_x, label = "Final")
    ax1.legend()
    ax2 = plt.subplot(212)
    ax2.plot(cs[ind==1], yobjs[ind==1], 'b.', alpha = 0.5)
    ax2.plot(cs[ind==0], yobjs[ind==0], 'r.', alpha = 0.5)
    ax2.plot(cs, ypredict, 'k--', alpha = 0.5)
    ax2.plot(cs, yinits, label="Initial")
    #plt.plot(cs[indy], yinits[indy])
    ax2.plot(cs[indy], yobjs[indy], label = "DAOFind")
    ax2.plot(cs, kk_y * yinits + bb_y, label = "Final")
    ax2.legend()
    figfile = fitsfile.replace("/processing/", "/registration/")
    figfile = figfile[:-5]+'_xypos.png'
    ax1.set_ylabel("x")
    ax2.set_ylabel("y")
    plt.tight_layout()
    plt.savefig(figfile)
    plt.close()
    
    for i in range(len(regis2HDU)):
        if i==0:
            continue
        hdu = regis2HDU[i]
        xfinal = xfinals[i-1]
        yfinal = yfinals[i-1]
        hdu.header["X_OBJ"] = (xfinal, "Final decision of the object's x Pixel Coordinate")
        hdu.header["Y_OBJ"] = (yfinal, "Final decision of the object's y Pixel Coordinate")
    
    regisfile = fitsfile.replace("/processing/", "/registration/")
    regisfile = regisfile[:-5]+'_regis.fits'
    print (" Writing to %s"%regisfile)
    regis2HDU.writeto(regisfile, overwrite=True)
            

def get_wcs_xy(ra, dec, wcsfile, fitsfile, get_distance = True):
    wcs_header = fits.open(wcsfile)[0].header
    w = WCS(wcs_header)
    x, y = w.wcs_world2pix(ra, dec, 1)
    wcsfileSplit = wcsfile.split("/")[-1]
    if "stack" in wcsfileSplit:
        xyframe = 1
    else:
        xyframe = int(wcsfileSplit.split("_")[-2])
    
    if get_distance == True:
        wcs_sci_header = fits.open(fitsfile)[0].header
        primary_header = fits.open(fitsfile)[0].header
        w_sci = WCS(wcs_sci_header)
        x_sci, y_sci = w_sci.wcs_world2pix(ra, dec, 1)
        distance = np.sqrt(((x - x_sci)*primary_header["PIXSCALX"])**2 +\
                    ((y - y_sci)*primary_header["PIXSCALY"])**2)
        print ("original astrometry is off by %.2f arcmin"%(distance/60))
    return float(x), float(y), xyframe
    

def forcedphotometry_kp(imagefile, x0, y0, fwhm=10., gain=1.0, zp=0.0):
    """
    zp: zero point for converting flux (in ADU) to magnitudes
    fwhm: scalar or 1D array of photometry aperture radii in pixel units.
    
    """
    hdulist = fits.open(imagefile)
    header = fits.getheader(imagefile)

    forcedfile = imagefile.replace(".fits",".forced")
    fid = open(forcedfile,'w')

    mjds, mags, magerrs, fluxes, fluxerrs = [], [], [], [], []
    for ii in range(len(hdulist)):
        if ii == 0: 
            continue
        header = hdulist[ii].header
        image = hdulist[ii].data
        if "GPS_TIME" in header:
            timeobs = Time(header["GPS_TIME"])
        elif "DATE" in header:
            timeobs = Time(header["DATE"])
        else:
            print("Warning: 'both GPS_TIME and DATE are missing from %s hdu %d/%d"%(imagefile,ii,len(hdulist)))
            continue
        
        
        plt.imshow(image, origin = "lower")
        
        mag, magerr, flux, fluxerr, sky, skyerr, badflag, outstr = \
            pp.aper.aper(image, x0, y0, phpadu=gain, apr=fwhm, zeropoint=zp, 
                         skyrad=[3*fwhm,5*fwhm],exact=False)

        mjds.append(timeobs.mjd)
        mags.append(mag)
        magerrs.append(magerr)
        fluxes.append(flux)
        fluxerrs.append(fluxerr)

        fid.write('%.5f %.5f %.5f %.5f %.5f\n'%(dateobs.mjd,mag,magerr,flux,fluxerr))
    fid.close()

    return np.array(mjds), np.array(mags), np.array(magerrs), np.array(fluxes), np.array(fluxerrs)

    


def stack_shifted_frames(fitsfile):
    """
    shift the image to account for drift issues, with extension = 1 as reference
    """
    median_size = 30
    
    hdulist = fits.open(fitsfile)
    reference = np.array(hdulist[1].data, dtype = float)
    blur = np.asfarray(median_filter(reference, size=(median_size,median_size)))
    reference = reference-blur
    last_image = deepcopy(reference)

    nimg = len(hdulist)-1
    nx = reference.shape[1]
    ny = nx
    
    shifty = 0 # This records the total shift in x of last image
    shiftx = 0 # This records the total shift in y of last image
    xshifts = np.ones(nimg+1)*(-1)
    yshifts = np.ones(nimg+1)*(-1)
    datas = np.zeros((ny, nx, nimg))
    masks = np.zeros((ny, nx, nimg), dtype = bool)
    for j in range(nimg):
        jj = j+1
        imagedata = hdulist[jj].data
        if jj%50==0:
            # do not need to calculate this for many times since cloud do not move much...?
            blur = np.asfarray(median_filter(imagedata, size=(median_size,median_size)))
        imagedata = imagedata-blur
        
        # shift by last images's total shift
        shifted = image_registration.fft_tools.shiftnd(imagedata, (shifty, shiftx))
            
        shift, error, diffphase = register_translation(last_image, shifted, upsample_factor=100) # compare with the last image
        shiftmax = 5
        shifttot = np.sqrt(shift[0]**2 + shift[1]**2) # the shift in this step compared with last image
        if shifttot <= shiftmax:
            shifted = image_registration.fft_tools.shiftnd(shifted, (shift[0], shift[1]))
            shifty = shifty + shift[0]
            shiftx = shiftx + shift[1]
                    
        last_image = deepcopy(shifted)
        
        if j%50==0:
            print ("    frame %d, shifty = %.2f, shiftx = %.2f"%(jj, shifty, shiftx))
        xshifts[jj] = shiftx
        yshifts[jj] = shifty
        
        mymask = np.zeros((shifted.shape[0], shifted.shape[1]))
        if shiftx>0 and shifty<0:
            mymask[:math.floor(shifty), math.ceil(shiftx):]= 1
        elif shiftx>0 and shifty>0:
            mymask[math.ceil(shifty):, math.ceil(shiftx):]= 1
        elif shiftx<0 and shifty>0:
            mymask[math.ceil(shifty):, :math.floor(shiftx)]= 1
        else:
            mymask[:math.floor(shifty), :math.floor(shiftx)]= 1
        if j==0:
            mymask[:,:] = 1
        masks[:,:,j] = mymask
        datas[:,:,j] = shifted
        """
        cmap_name = "gray"
        plt.figure(figsize=(6,6))
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(last_image, cmap = cmap_name, origin='lower', norm=norm)
        plt.gca().set_axis_off()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
                            hspace = 0, wspace = 0)
        plt.margins(0,0)
        plt.gca().xaxis.set_major_locator(NullLocator())
        plt.gca().yaxis.set_major_locator(NullLocator())
        plt.savefig('/Users/yuhanyao/Desktop/kped_tmp/20191117/ZTFJ01395245/figures/%d.png'%jj)
        plt.close()
        """
    data_stacked = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):
            values = datas[i,j,:][masks[i,j,:]]
            data_stacked[i,j] = np.median(values)
   
    tb_shift = Table(data = [xshifts, yshifts], names = ["xshift", "yshift"],
                     dtype=('f2', 'f2'))
    return tb_shift, data_stacked


def get_sortedfiles_from_fits(fitsfiles):
    nums = []
    for fitsfile in fitsfiles:
        fitsfileSplit = fitsfile.replace(".fits.fz","").replace(".fits","").replace("_proc","").split("_")  
        try:
            num = int(fitsfileSplit[-1])
        except:
            num = -1
        nums.append(num)
    fitsfiles = [fitsfiles for _,fitsfiles in sorted(zip(nums,fitsfiles))]
    return fitsfiles


def stack_images(stackDir, fitsfiles, nimages,
                 doRegistration=False, registration_size=-1, x=None, y=None): 
    
    fitsfiles = get_sortedfiles_from_fits(fitsfiles)
    shiftx, shifty = 0, 0

    cnt = 0
    for ii in range(len(fitsfiles)):
        scienceimage = os.path.join(stackDir,fitsfiles[ii].split("/")[-1])
        if os.path.isfile(scienceimage): 
            continue
        hdulist = fits.open(fitsfiles[ii])

        if cnt == 0:
            if doRegistration:
                reference = hdulist[1].data

        hdulist2 = []
        cnt = 1
        for jj in range(len(hdulist)):
            if jj == 0:
                hdulist2.append(hdulist[jj])
            else:
                if cnt == 1:
                    hdulist_hold = copy.copy(hdulist[jj])
                    xshape, yshape = hdulist_hold.data.shape
                    data = np.empty([xshape,yshape,0])

                if doRegistration:
                    if registration_size > 0:
                        #image_shape = reference.shape
                        xlow = int(x - registration_size/2)
                        xhigh = int(x + registration_size/2)
                        ylow = int(y - registration_size/2)
                        yhigh = int(y + registration_size/2)

                        shift, error, diffphase = register_translation(reference[xlow:xhigh,ylow:yhigh], hdulist[jj].data[xlow:xhigh,ylow:yhigh], upsample_factor=1)
                        shiftx, shifty = shift[0], shift[1]
                    else:
                        # shift by last images's shift
                        shifted = image_registration.fft_tools.shiftnd(hdulist[jj].data, (shiftx, shifty))                        

                        shift, error, diffphase = register_translation(reference, shifted, upsample_factor=1)

                        shiftmax = 5
                        shifttot = np.sqrt(shift[0]**2 + shift[1]**2)
                                  
                        if shifttot <= shiftmax:
                            shiftx = shiftx + shift[0]
                            shifty = shifty + shift[1]

                    shifted = image_registration.fft_tools.shiftnd(hdulist[jj].data, (shiftx,shifty))
                    data = np.append(data,np.expand_dims(shifted,axis=2),axis=2)
                else:
                    data = np.append(data,np.expand_dims(hdulist[jj].data,axis=2),axis=2)
                cnt = cnt + 1

                if cnt == nimages+1:
                    hdulist_hold.data = np.mean(data,axis=2)
                    hdulist2.append(hdulist_hold)
                    cnt = 1 

        hdulist2 = fits.HDUList(hdus=hdulist2)
        hdulist2.writeto(scienceimage,output_verify='fix',overwrite=True)
    
    
def utc2date(utcheader):
    timestr = utcheader.split("_")
    yy, mm, dd = int(timestr[0][:4]), int(timestr[0][4:6]), int(timestr[0][6:8])
    hh, mmm, ss, ums = int(timestr[1][:2]), int(timestr[1][2:4]), int(timestr[1][4:6]), int(timestr[1][7:])
    d = datetime.datetime(yy, mm, dd, hh, mmm, ss, ums)
    return Time(d, format="datetime")
    

        