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
from astropy.coordinates import SkyCoord, EarthLocation
from copy import deepcopy
from scipy.ndimage import median_filter

import image_registration
from skimage.feature import register_translation
from photutils import DAOStarFinder

import PythonPhot as pp
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from citizen.photometry_utils import ps1_query

import matplotlib
matplotlib.use('Agg')
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
            

def get_reference_pos(scienceimage, cat, objName, zp=0, passband='sg', xoff=0, yoff=0, refmag=0):
    sciHDU = fits.open(scienceimage)
    ndim = sciHDU[1].data.shape[0]
    nmulti = int(1024/ndim)
    
    if objName=="WS35":
        xoff = 86
        yoff = -8
        rafield = 13.280726
        decfield = 67.50018
        refmag = 12.8 # saturated in ps1 (otherwise return 18.67 mag one close by)
        print ("    Fixing reference star position: xoff = %.2f, yoff = %.2f"%(xoff, yoff))
        print ("    Fixing reference star: ra_ref = %.5f, dec_ref = %.5f"%(rafield, decfield))
    elif xoff!=0 and yoff!=0 and refmag!=0:
        rafield = -99
        print ("    Fixing reference star position: xoff = %.2f, yoff = %.2f"%(xoff, yoff))
    else:
        rafield = None
        # plt.hist(cat[:,20])
        cat = cat[cat[:,20]<15] # cut sources with too high fwhm
        # plt.hist(cat[:,14])
        cat = cat[cat[:,14]/cat[:,15]<2.5] # Profile RMS along major axis 
        # plt.hist(cat[:,15])
        # cat = cat[cat[:,15]<3] # Profile RMS along minor axis   
        magthreshold = zp - 8
        cat = cat[cat[:,4]<magthreshold] # mag
    
        xobjs = np.array([hdu.header["X_OBJ"] for hdu in sciHDU[1:]])
        yobjs = np.array([hdu.header["Y_OBJ"] for hdu in sciHDU[1:]])
        xtolerance = 10#max(xobjs) - min(xobjs)
        ytolerance = 10#max(yobjs) - min(yobjs)
        ndata = sciHDU[1].data.shape[0]
        ind = (cat[:,0]>xtolerance)&(cat[:,0]<(ndata-xtolerance))&(cat[:,1]>ytolerance)&(cat[:,1]<(ndata-ytolerance))
        cat = cat[ind]
        
        nrefs = np.array([np.sum(cat[:,22]==i) for i in range(len(sciHDU))])
        framenum = np.where(nrefs == np.max(nrefs))[0][0]
        xobj_frame = xobjs[framenum-1]
        yobj_frame = yobjs[framenum-1]
        xind = (cat[:,0]<(ndata-(max(xobjs)-xobj_frame)))&(cat[:,0]>(xobj_frame-min(xobjs)))
        yind = (cat[:,1]<(ndata-(max(yobjs)-yobj_frame)))&(cat[:,1]>(yobj_frame-min(yobjs)))
        ind = (cat[:,22]==framenum)&xind&yind
    
        catsub = cat[ind] # in pixel coordinate
        hdu = sciHDU[framenum]
        header = hdu.header
        #xguess = header["X_GUESS"]
        #yguess = header["Y_GUESS"]    
        xobj = header["X_OBJ"]
        yobj = header["Y_OBJ"]
        data = hdu.data
        dist = np.sqrt((catsub[:,0] - xobj)**2+(catsub[:,1] - yobj)**2)
        catsub = np.vstack((catsub.T, dist.T)).T
        catsub = catsub[catsub[:, -1]>30]
        # If there is a star that is bright (and have good psf shape)
        indbright = (catsub[:,4]<(zp-10))&(catsub[:,20]<5)
        if np.sum(indbright)!=0:
            index = np.arange(len(catsub))
            catbright = catsub[indbright]
            indexbright = index[indbright]
            ind_field = indexbright[np.argsort(catbright[:, -1])][0]
        else:
            ind_field = np.argsort(catsub[:, -1])[0]
        xfield = catsub[:,0][ind_field]
        yfield = catsub[:,1][ind_field]
    
        xoff = xfield - xobj_frame
        yoff = yfield - yobj_frame
        print ("  Field reference star: xoff = %.2f, yoff = %.2f"%(xoff, yoff))
    
        # plot the figure
        print ("  Plotting figures...")
        path_out_dir = "/".join(scienceimage.split("/")[:-1])
        # mag
        figfile1 = os.path.join(path_out_dir, "image_mag.pdf")
        plt.figure(figsize=(9,8))
        cmap_name = "gray"
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(data, origin = "lower", norm=norm, cmap = cmap_name,
                   vmax = np.percentile(data.ravel(), 99.9))
        plt.scatter(catsub[:,0]-1, catsub[:,1]-1, c=catsub[:,4], s=10)
        plt.colorbar()
        plt.plot(xobj-1, yobj-1, 'rx')
        plt.plot(xfield-1, yfield-1, 'mx')
        plt.tight_layout()
        plt.savefig(figfile1)
        plt.close()
        # FWHM
        figfile2 = os.path.join(path_out_dir, "image_fwhm.pdf")
        plt.figure(figsize=(9,8))
        cmap_name = "gray"
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(data, origin = "lower", norm=norm, cmap = cmap_name,
                   vmax = np.percentile(data.ravel(), 99.9))
        plt.scatter(catsub[:,0]-1, catsub[:,1]-1, c=catsub[:,20], s=10)
        plt.colorbar()
        plt.plot(xobj-1, yobj-1, 'rx')
        plt.plot(xfield-1, yfield-1, 'mx')
        plt.tight_layout()
        plt.savefig(figfile2)
        plt.close()
        # ELONGATION
        figfile3 = os.path.join(path_out_dir, "image_AoverB.pdf")
        plt.figure(figsize=(9,8))
        cmap_name = "gray"
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(data, origin = "lower", norm=norm, cmap = cmap_name,
                   vmax = np.percentile(data.ravel(), 99.9))
        plt.scatter(catsub[:,0]-1, catsub[:,1]-1, c=catsub[:,14]/catsub[:,15], s=10)
        plt.colorbar()
        plt.plot(xobj-1, yobj-1, 'rx')
        plt.plot(xfield-1, yfield-1, 'mx')
        plt.tight_layout()
        plt.savefig(figfile3)
        plt.close()
    
    for i in range(len(sciHDU)):
        if i==0:
            continue
        hdu = sciHDU[i]
        header = hdu.header
        header["X_FIELD"] = (header["X_OBJ"] + xoff, "Field reference x Pixel Coordinate")
        header["Y_FIELD"] = (header["Y_OBJ"] + yoff, "Field reference y Pixel Coordinate")
    
    if rafield==None:
        ra = sciHDU[0].header["RA_OBJ"]
        dec = sciHDU[0].header["DEC_OBJ"]
        wcsauto = sciHDU[0].header["WCSAUTO"]
        framenum = sciHDU[0].header["WCSFRAME"]
        if wcsauto==1:
            w0 = WCS(sciHDU[framenum].header)
            x0, y0 = w0.wcs_world2pix(ra, dec, 1)
            x0field = x0 + xoff
            y0field = y0 + yoff
            rafield, decfield = w0.wcs_pix2world(x0field, y0field, 1)
        else:
            rafield = ra - sciHDU[0].header["PIXSCALX"]*nmulti*xoff/3600
            decfield = dec + sciHDU[0].header["PIXSCALY"]*nmulti*yoff/3600
            sciHDU[0].header["RA_REF"] = float(rafield)
            sciHDU[0].header["DEC_REF"] = float(decfield)
        print ("    Find reference star: ra_ref = %.5f, dec_ref = %.5f"%(rafield, decfield))
        radius_deg = 20/60/60. # 20 arcsec
        result = ps1_query(rafield, decfield, radius_deg, maxmag=20,
                           maxsources=100)
        if len(result)==0:
            refmag = 0
            print ("    Reference ps1 query failed! Setting refmag = 0")
        else:
            index = np.argsort(result["_r"].data.data)[0]
            if passband == "sg":
                refmag = result["gmag"][index]
            elif passband == "sr":
                refmag = result["rmag"][index]

    print ("    Reference star: %.2f mag in %s band"%(refmag, passband))
    sciHDU[0].header["MAGREF"] = refmag
    return sciHDU
    

def register_images(fitsfile, shiftfile, xyframe, x, y, path_out_dir,
                    maxdist=10., aper_size=10, refit = True, offtune = False):
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
            if x_frame<aper_size or x_frame>(ndata-aper_size) or y_frame<aper_size or y_frame>(ndata-aper_size):
                print ("    Target outside of image (aperture size=%.1f). Remove extention %d!"%(aper_size, i))
            else:
                hdu.header["x_guess"] = (x_frame, "Initial guess of the object x Pixel Coordinate")
                hdu.header["y_guess"] = (y_frame, "Initial guess of the object y Pixel Coordinate")
                regis1HDU.append(hdu)
    
    if not offtune:
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
                plt.savefig('/Users/yuhanyao/Desktop/kped_tmp/20191116/ZTFJ11514412/figures/%d.png'%i)
                plt.close()
                """
                x_tuned = xsub_tuned + xstart
                y_tuned = ysub_tuned + ystart
                hdu.header["DAOFLAG"] = (findsourceflag, "DAOStarFinder flag")
                hdu.header["X_DAOOBJ"] = (x_tuned, "object x Pixel Coordinate found by DAOStarFinder")
                hdu.header["Y_DAOOBJ"] = (y_tuned, "object y Pixel Coordinate found by DAOStarFinder")
                regis2HDU.append(hdu)
            
        cs = np.arange(len(regis2HDU[1:]))
        ind = np.array([hdu.header["DAOFLAG"] for hdu in regis2HDU[1:]])
        xobjs = np.array([hdu.header["X_DAOOBJ"] for hdu in regis2HDU[1:]])
        yobjs = np.array([hdu.header["Y_DAOOBJ"] for hdu in regis2HDU[1:]])
        xinits = np.array([hdu.header["X_GUESS"] for hdu in regis2HDU[1:]])
        yinits = np.array([hdu.header["Y_GUESS"] for hdu in regis2HDU[1:]])
    
        if refit==True:
            xpredict = get_predictioon_xy(ind, xobjs, cs)
            ypredict = get_predictioon_xy(ind, yobjs, cs)
    
            indx = abs(xpredict - xobjs) < 1.5
            indy = abs(ypredict - yobjs) < 1.5
            if len(np.unique(xinits[indx])) < 2:
                kk_x = 1
                bb_x = 0
            else:
                kk_x, ekk, bb_x = mylinear_fit(xinits[indx], xobjs[indx], np.ones(np.sum(indx)), npar = 2)
            if len(np.unique(yinits[indy])) < 2:
                kk_y = 1
                bb_y = 0
            else:
                kk_y, ekk, bb_y = mylinear_fit(yinits[indy], yobjs[indy], np.ones(np.sum(indy)), npar = 2)
            xfinals = kk_x * xinits + bb_x
            yfinals = kk_y * yinits + bb_y
    
        plt.figure(figsize=(12, 12))
        ax1 = plt.subplot(211)
        ax1.plot(cs[ind==0], xobjs[ind==0], 'r.', alpha = 0.5)
        ax1.plot(cs[ind==1], xobjs[ind==1], 'b.', alpha = 0.5)
        ax1.plot(cs, xinits, label="Initial")
        # plt.plot(cs[indx], xinits[indx])
        ax2 = plt.subplot(212)
        ax2.plot(cs[ind==1], yobjs[ind==1], 'b.', alpha = 0.5)
        ax2.plot(cs[ind==0], yobjs[ind==0], 'r.', alpha = 0.5)
        ax2.plot(cs, yinits, label="Initial")
        #plt.plot(cs[indy], yinits[indy])
        if refit==True:
            ax1.plot(cs[indx], xobjs[indx], label = "DAOFind")
            ax1.plot(cs, xpredict, 'k--', alpha = 0.5)
            ax1.plot(cs, kk_x * xinits + bb_x, label = "Final")
            ax2.plot(cs[indy], yobjs[indy], label = "DAOFind")
            ax2.plot(cs, ypredict, 'k--', alpha = 0.5)
            ax2.plot(cs, kk_y * yinits + bb_y, label = "Final")
        figfile = os.path.join(path_out_dir, 'findxypos.png')
        ax1.set_ylabel("x")
        ax2.set_ylabel("y")
        ax1.legend()
        ax2.legend()
        plt.tight_layout()
        plt.savefig(figfile)
        plt.close()
    
        for i in range(len(regis2HDU)):
            if i==0:
                continue
            hdu = regis2HDU[i]
            if refit==True:
                xfinal = xfinals[i-1]
                yfinal = yfinals[i-1]
            else:
                xfinal = hdu.header["X_DAOOBJ"]
                yfinal = hdu.header["Y_DAOOBJ"]
            hdu.header["X_OBJ"] = (xfinal, "Final decision of the object's x Pixel Coordinate")
            hdu.header["Y_OBJ"] = (yfinal, "Final decision of the object's y Pixel Coordinate")
        return regis2HDU
    else:
        for i in range(len(regis1HDU)):
            hdu = regis1HDU[i]
            if i==0:
                continue
            xfinal = hdu.header["x_guess"]
            yfinal = hdu.header["y_guess"]
            hdu.header["X_OBJ"] = (xfinal, "Final decision of the object's x Pixel Coordinate")
            hdu.header["Y_OBJ"] = (yfinal, "Final decision of the object's y Pixel Coordinate")
        return regis1HDU
            

def get_wcs_xy(ra, dec, wcsfile, fitsfile, get_distance = True):
    # PIXSCALX/Y: arcsec per pixel in UNBINNED image
    
    hdus = fits.open(fitsfile)
    ndim = hdus[1].data.shape[0]
    nmulti = int(1024/ndim)
    
    wcs_header = fits.open(wcsfile)[0].header
    w = WCS(wcs_header)
    x, y = w.wcs_world2pix(ra, dec, 1)
    wcsfileSplit = wcsfile.split("/")[-1]
    if "stack" in wcsfileSplit:
        xyframe = 1
    else:
        xyframe = int(wcsfileSplit.split("_")[-2])
    
    if get_distance == True:
        wcs_sci_header = fits.open(fitsfile)[1].header
        primary_header = fits.open(fitsfile)[0].header
        w_sci = WCS(wcs_sci_header)
        x_sci, y_sci = w_sci.wcs_world2pix(ra, dec, 1)
        distance = np.sqrt(((x - x_sci)*primary_header["PIXSCALX"]*nmulti)**2 +\
                    ((y - y_sci)*primary_header["PIXSCALY"]*nmulti)**2)
        print ("  FYI: Original astrometry is off by %.2f arcmin"%(distance/60))
    
    hdus[0].header["WCSFRAME"] = (xyframe, "The frame extension used to find astrometry.")
    hdus[0].header["WCSAUTO"] = (1, "Indicate whether wcs is automatically found by astrometry.net")
    header = hdus[xyframe].header
    
    keys = ["CTYPE1", "CTYPE2", "EQUINOX", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2",
            "CUNIT1", "CUNIT2", "CD1_1", "CD1_2", "CD2_1", "CD2_2"]  
    for key in keys:
        header[key] = wcs_header[key]
    print ("  Saving correct wcs to frame %d's header of %s"%(xyframe, fitsfile))
    hdus.writeto(fitsfile, overwrite=True)
    return float(x), float(y), xyframe


def update_wcsstatus(fitsfile, xyframe):
    hdus = fits.open(fitsfile)
    hdus[0].header["WCSFRAME"] = (xyframe, "The frame extension used to find astrometry.")
    hdus[0].header["WCSAUTO"] = (0, "Indicate whether wcs is automatically found by astrometry.net")
    hdus.writeto(fitsfile, overwrite=True)
    

def forcedphotometry_kp(scienceimage, aper_size=10., gain=1.0, zp=0.0,
                        sky_in = 11, sky_out = 30,
                        xkey = "X_OBJ", ykey = "Y_OBJ"):
    """
    https://github.com/djones1040/PythonPhot/blob/master/PythonPhot/aper.py
    zp: zero point for converting flux (in ADU) to magnitudes
    fwhm: scalar or 1D array of photometry aperture radii in pixel units.
    """
    hdulist = fits.open(scienceimage)

    mjds, mags, magerrs, fluxes, fluxerrs = [], [], [], [], []
    
    exptime = hdulist[0].header["FRAMETIM"] # in second
    for ii in range(len(hdulist)):
        if ii == 0: 
            continue
        header = hdulist[ii].header
        image = hdulist[ii].data
        if "GPS_TIME" in header:
            timeobsend = Time(header["GPS_TIME"])
        elif "DATE" in header:
            timeobsend = Time(header["DATE"])
        else:
            print("Warning: 'both GPS_TIME and DATE are missing from %s hdu %d/%d"%(scienceimage,ii,len(hdulist)))
            continue
        x = header[xkey]
        y = header[ykey]
        
        mag, magerr, flux, fluxerr, sky, skyerr, badflag, outstr = \
            pp.aper.aper(image, x, y, phpadu=gain, apr=aper_size, zeropoint=zp, 
                         skyrad=[sky_in, sky_out], exact=False)

        mjds.append(timeobsend.mjd - exptime/2./86400)
        mags.append(mag)
        magerrs.append(magerr)
        fluxes.append(flux)
        fluxerrs.append(fluxerr)
    #plt.errorbar(mjds, mags, magerrs, fmt='.k')
    return np.array(mjds), np.array(mags), np.array(magerrs), np.array(fluxes), np.array(fluxerrs)


def stack_shifted_frames(fitsfile, fix_ref = False):
    """
    shift the image to account for drift issues, with extension = 1 as reference
    """
    median_size = 30
    shiftmax = 5
    
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
        if fix_ref==False:
            shifted = image_registration.fft_tools.shiftnd(imagedata, (shifty, shiftx))
            shift, error, diffphase = register_translation(last_image, shifted, upsample_factor=100) # compare with the last imag
            shifttot = np.sqrt(shift[0]**2 + shift[1]**2) # the shift in this step compared with last image
            if shifttot <= shiftmax:
                shifted = image_registration.fft_tools.shiftnd(shifted, (shift[0], shift[1]))
                shifty = shifty + shift[0]
                shiftx = shiftx + shift[1]
            last_image = deepcopy(shifted)
        else:
            shifted = imagedata
            shift, error, diffphase = register_translation(reference, imagedata, upsample_factor=100) # compare with the last image
            shifttot = np.sqrt(shift[0]**2 + shift[1]**2) 
            if shifttot <= shiftmax:
                shifty = shift[0]
                shiftx = shift[1]
        if j%50==0:
            print ("    frame %d, shifty = %.2f, shiftx = %.2f"%(jj, shifty, shiftx))
        xshifts[jj] = shiftx
        yshifts[jj] = shifty
        
        mymask = np.zeros((imagedata.shape[0], imagedata.shape[1]))
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
    

def jd2bjd(jd, RA, Dec, obs_site = "Kitt Peak"):
    """
    Convert jd into bjd
    https://en.wikipedia.org/wiki/Barycentric_Julian_Date
    https://docs.astropy.org/en/stable/time/
    """
    times = Time(jd, format='jd', scale='utc')
    c = SkyCoord(RA, Dec, unit="deg") # defaults to ICRS frame
    observatory = EarthLocation.of_site(obs_site)
    ltt_bary = times.light_travel_time(c, kind='barycentric', location=observatory)
    time_barycentre = times.tdb + ltt_bary
    # tdb:  # Barycentric Dynamical Time
    return time_barycentre


def jd2hjd(jd, RA, Dec, obs_site = "Kitt Peak"):
    """
    Convert jd into hjd
    https://en.wikipedia.org/wiki/Heliocentric_Julian_Day
    """
    times = Time(jd, format='jd', scale='utc')
    c = SkyCoord(RA, Dec, unit="deg")
    observatory = EarthLocation.of_site(obs_site)
    ltt_helio = times.light_travel_time(c, kind='heliocentric', location=observatory)
    times_heliocentre = times.utc + ltt_helio
    # tdb:  # Barycentric Dynamical Time
    return times_heliocentre
        