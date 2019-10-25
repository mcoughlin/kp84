## Code to re-run Sextractor on the LCO images 
## Change the folder that contains the unpacked LCO data

import optparse
import warnings
import glob
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import os

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from astroML.crossmatch import crossmatch_angular
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import Angle
from astropy.io import fits
from astropy.wcs import WCS
from astropy import wcs
import astropy.coordinates as coord
from astropy.time import Time
import aplpy
from astropy.nddata import Cutout2D

import PythonPhot as pp

import ztfsub.utils, ztfsub.surveys
import ztfsub.plotting

import urllib
import urllib.request

from photometry_utils import panstarrs_query,get_ztf_cand,get_ztf_object,do_crossmatch,do_sextractor,sdss_query, ps1_query, gaia_query
from penquins import Kowalski

__version__ = 1.0

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("--tmpDir",default="/tmp")
    parser.add_option("-o", "--outputDir", help="output file",default="../output/Jacquinery/")
    parser.add_option("-l","--transients", help="Transient list.", default='../data/Jacquinery/transients.dat')
    parser.add_option("-t","--transient", help="Transient name.", default='CS0123-5848')
    parser.add_option("-d","--dataDir", help="Data directory.", default='../data/Jacquinery/00h36m27_-27d47m12/') 
    parser.add_option("-s","--sextractorDir", help="sextractor directory.", default='../defaults/')
    parser.add_option("--defaultsDir",default="../defaults")
    parser.add_option("-f","--filter", help="filter.", default='g')
    parser.add_option("-p","--pixel_scale", help="pixel scale (1m=0.389, 2m=0.3).", default=0.485419262282388)

    #parser.add_option("--ra_diff",default=13.40264868710,type=float)
    #parser.add_option("--dec_diff",default=-23.83615397420,type=float)
    parser.add_option("--ra_diff",default=13.932779,type=float)
    parser.add_option("--dec_diff",default=-22.974694,type=float)

    parser.add_option("--doSubtraction",  action="store_true", default=False)
    parser.add_option("--subtractionDir",default="../subtraction")
    parser.add_option("--subtractionSource",default="decals")
    parser.add_option("--image_size",default=200,type=int)

    parser.add_option("--doForcedPhotometry",  action="store_true", default=False)
    parser.add_option("--doAstrometryNet",  action="store_true", default=False)
    parser.add_option("--doKowalski",  action="store_true", default=False)
    parser.add_option("--doPlots",  action="store_true", default=False)
    parser.add_option("--doDifferential",  action="store_true", default=False)

    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Run verbosely. (Default: False)")

    opts, args = parser.parse_args()

    # show parameters
    if opts.verbose:
        print >> sys.stderr, ""
        print >> sys.stderr, "running gwemopt_run..."
        print >> sys.stderr, "version: %s"%__version__
        print >> sys.stderr, ""
        print >> sys.stderr, "***************** PARAMETERS ********************"
        for o in opts.__dict__.items():
          print >> sys.stderr, o[0]+":"
          print >> sys.stderr, o[1]
        print >> sys.stderr, ""

    return opts

# =============================================================================
#
#                                    MAIN
#
# =============================================================================

warnings.filterwarnings("ignore")

# Parse command line
opts = parse_commandline()

data_path = opts.dataDir
outputDir = os.path.join(opts.outputDir,opts.transient)

# sextractor parameters & options
sextractor = True # change to True to run sextractor
control = True

if not os.path.isdir(outputDir):
    os.makedirs(outputDir)

rm_command = "rm %s/image*"%data_path
os.system(rm_command)

sextractor_path = opts.sextractorDir # config, params and other basic files here
workdir_sex_out = outputDir # folder to drop sextractor output
# workdir_cat = './cat/' # folder to save catalogs
workdir_cat = data_path # folder to save catalogs

#### PATHS & CONFIG ####

# output file
kped_saved = data_path
workdir_fits_red = kped_saved # folder with KPED data

if opts.doKowalski:
    # read files
    username = 'kped'
    password = 'queryitEDdy!'

    ko = Kowalski(protocol='https', host='kowalski.caltech.edu', port=443,
                  verbose=True,
                  username=username, password=password)
else:
    names = ('name', 'ra', 'dec')
    transients = Table.read(opts.transients, format='ascii',
                            names=names, data_start=0)

ra_all_tran,ra_header_=[],[]
dec_all_tran,dec_header_=[],[]
name_ =[]
filt_=[]
date_,date_jd_=[],[]

files_old = glob.glob(workdir_fits_red+'*.fit')
print(files_old)#,workdir_fits_red

for i in range(len(files_old)):

    fitsfile = files_old[i]
    print('Checking: %s' %fitsfile)

    hdu = fits.open(fitsfile)
    #ID = hdu[1].header['OBJECT']
    ra_header_.append(hdu[0].header['CRVAL1'])
    dec_header_.append(hdu[0].header['CRVAL2'])
    #filt_.append(hdu[0].header['filter'])
    filt_.append(opts.filter)
    date = Time(hdu[0].header['DATE-OBS'], format='isot', scale='utc')
    date_.append(date)
    date_jd_.append(date.jd)
    name_.append(opts.transient)

    if opts.doKowalski:
        q = {"query_type": "general_search", "query": "db['ZTF_alerts'].find({'objectId': {'$eq': '"+ID+"'}})" }
        r = ko.query(query=q,timeout=30)

        if len(r['result_data']['query_result']) >0:
            # getting metadata
            candidate = r['result_data']['query_result'][0]
            ra,dec = candidate['candidate']['ra'],candidate['candidate']['dec']
            ra_all_tran.append(ra)
            dec_all_tran.append(dec)
#           print("query worked")
#           print('solve-field '+files[i][:-5]+'_red.fits --ra '+str(ra)+' --dec '+str(dec)+' --dir /home/roboao/Tomas/output --scale-units arcsecperpix --scale-low 0.255 --scale-high 0.26 --radius 0.04 --overwrite')
        else:
            ra_all_tran.append(np.nan)
            dec_all_tran.append(np.nan)
    else:
        idx = np.where(opts.transient == transients["name"])[0]
        ra_hex = transients["ra"][idx]
        dec_hex = transients["dec"][idx]

        ra  = Angle(ra_hex, unit=u.hour).deg
        dec = Angle(dec_hex, unit=u.deg).deg

        ra_all_tran.append(ra)
        dec_all_tran.append(dec)

    hdu0 = fits.PrimaryHDU(header=hdu[0].header)
    hdu1 = fits.ImageHDU(header=hdu[0].header,data=hdu[0].data)

    #del hdu[0].header["EXTEND"]
    hdulist = fits.HDUList(hdus=[hdu0,hdu1])

    fitsfilenew = data_path+'/image_'+str(i)+'.fits'
    hdulist.writeto(fitsfilenew,overwrite=True)

    if opts.doAstrometryNet:
        print('Running Astrometry.net')
        ztfsub.utils.astrometrynet(fitsfilenew,pixel_scale=opts.pixel_scale,ra=ra,dec=dec,radius=5.0,depth=100,cutedges=-1,ext=1)

    if opts.doSubtraction:
        print('Running image differencing')
        passband = opts.filter
        refimage = os.path.join(outputDir,'ref_%s.fits' % passband)
        if os.path.isfile(refimage):
            rm_command = "rm %s"%refimage
            os.system(rm_command)

        if passband == "U":
            refband = "u"
        elif passband == "B":
            refband = "g"
        elif passband == "V":
            refband = "r"
        elif passband == "R":
            refband = "i"
        elif passband == "I":
            refband = "z"
        elif passband == "g":
            refband = "g"
        elif passband == "r":
            refband = "r"
        elif passband == "i":
            refband = "i"

        if (opts.subtractionSource == "decals") and (passband == "i"):
            refband = "i"
 
        if not os.path.isfile(refimage):
            if opts.subtractionSource == "sdss":
                refgood = ztfsub.surveys.get_sdss(opts,refimage,ra,dec,refband)
            elif opts.subtractionSource == "ps1":
                refgood = ztfsub.surveys.get_ps1(opts,refimage,ra,dec,refband)
            elif opts.subtractionSource == "decals":
                if refband == "i":
                    refgood = ztfsub.surveys.get_decals(opts,refimage,ra,dec,"r")
                    hdulist1=fits.open(refimage)
                    rm_command = "rm %s" % refimage
                    os.system(rm_command) 
                    refgood = ztfsub.surveys.get_decals(opts,refimage,ra,dec,"z")
                    hdulist=fits.open(refimage)
                    hdulist[0].data = hdulist[0].data + hdulist1[0].data
                    hdulist.writeto(refimage,overwrite=True)
                else:
                    refgood = ztfsub.surveys.get_decals(opts,refimage,ra,dec,refband)
            elif opts.subtractionSource == "dr1":
                refgood = ztfsub.surveys.get_dr1(opts,refimage,ra,dec,refband)
            else:
                print("Only PS1 and SDSS supported.")
                exit(0)
        else:
            refgood = True

        tmpdir='%s/i'%(outputDir)
        if not os.path.isdir(tmpdir):
            os.makedirs(tmpdir)

        ztfsub.utils.p60sdsssub(opts, fitsfilenew, refimage, [ra,dec],
                    distortdeg=1, scthresh1=3.0,
                    scthresh2=10.0, tu=60000, iu=60000, ig=2.3, tg=1.0,
                    stamps=None, nsx=8, nsy=8, ko=0, bgo=0, radius=10,
                    tlow=-5000.0, ilow=-5000.0, sthresh=3.0, ng=None,
                    aperture=10.0,
                    defaultsDir=opts.defaultsDir)
       
# #name_,ra_all_tran,dec_all_tran=get_ztf_cand(url_cand, username, password)
name_ = np.asarray(name_)
ra_all_tran =np.asarray(ra_all_tran)
dec_all_tran = np.asarray(dec_all_tran)

files_entire = glob.glob(workdir_fits_red+'image*.fits')
files_entire = [x for x in files_entire if (not "sub" in x) and (not "shift" in x)]

print('Analyzing: %s' % ",".join(files_entire))

if opts.doDifferential:
    ps1_table_diff = ps1_query(opts.ra_diff, opts.dec_diff, 1.5/3600.0,maxsources=1, maxmag=30.0)

# do reduction
KPED = []
for i in range(len(files_entire)):
    imagefile = files_entire[i]

    print('Extracting photometry: %s' %imagefile)

    imagefilesub = imagefile.replace(".fits",".sub.fits")
    hdu = fits.open(imagefile)

    ra_header = ra_header_[i]
    dec_header = dec_header_[i]
    ra_transient = ra_all_tran[i]
    dec_transient = dec_all_tran[i]
    filt = filt_[i]
    print('Object: ',files_entire[i])
    print('RA:' ,ra_header,'DEC:', dec_header)

    radius_deg = 0.25
    if filt == 'up':
        radius_deg = 0.18
        ps1_table = sdss_query(ra_header, dec_header, radius_deg)
        gaia_table = gaia_query(ra_header, dec_header, radius_deg)
    else:
        ps1_table = ps1_query(ra_header, dec_header, radius_deg,maxsources=3000, maxmag=30.0)
        gaia_table = gaia_query(ra_header, dec_header, radius_deg,maxsources=3000, maxmag=30.0)

    ps1_table = ps1_table.filled()
    gaia_table = gaia_table.filled()

    # sextract data from fits file
    do_sextractor(imagefile.split("/")[-1],workdir_fits_red,workdir_sex_out,workdir_cat,sextractor_path,control)
    # get photometry from sextractor
    try:
        cat = Table.read(imagefile[:-5]+'.cat',format='ascii.sextractor')
    except:
        raise ValueError('catalog not created, run the program again using sextractor = True')

    if opts.doSubtraction:
        filesub = imagefile.split("/")[-1].replace(".fits",".sub.fits")
        do_sextractor(filesub,workdir_fits_red,workdir_sex_out,workdir_cat,sextractor_path,control)      

        # get photometry from sextractor
        try:
            catsub = Table.read(workdir_cat+filesub[:-5]+'.cat',format='ascii.sextractor')
        except:
            raise ValueError('catalog not created, run the program again using sextractor = True')

    ra,dec,mag_psf,mag_psf_err,x,y = cat['ALPHA_J2000'] , cat['DELTA_J2000'],cat['MAG_BEST'],cat['MAGERR_BEST'],cat['X_IMAGE'],cat['Y_IMAGE']
    if opts.doSubtraction:
        ra_sub,dec_sub,mag_psf_sub,mag_psf_err_sub,x_sub,y_sub = catsub['ALPHA_J2000'] , catsub['DELTA_J2000'],catsub['MAG_BEST'],catsub['MAGERR_BEST'],catsub['X_IMAGE'],catsub['Y_IMAGE'] 
 
    if opts.doForcedPhotometry:
        hdulist=fits.open(imagefile)
        image = hdulist[1].data
        header = hdulist[1].header
        fwhm, zp, gain = 8.0, 0.0, 1.0
        mag_forced, mag_err_forced = [], []
        for x0, y0 in zip(x, y):
            mag,magerr,flux,fluxerr,sky,skyerr,badflag,outstr = pp.aper.aper(image,x0,y0,phpadu=gain,apr=fwhm,zeropoint=zp,skyrad=[3*fwhm,4*fwhm],exact=False)
            mag_forced.append(mag[0])
            mag_err_forced.append(magerr[0])
        mag_psf, mag_psf_err = np.array(mag_forced), np.array(mag_err_forced)

        if opts.doDifferential:
            w = WCS(header)
            x0,y0 = w.wcs_world2pix(opts.ra_diff,opts.dec_diff,1)

            xnew, ynew = pp.cntrd.cntrd(image,np.array([x0]),np.array([y0]),3*fwhm,keepcenter=True)
            if not xnew == -1:
                xnew, ynew = xnew[0], ynew[0]
                x0, y0 = xnew, ynew
            mag_psf_diff,mag_psf_err_diff,flux,fluxerr,sky,skyerr,badflag,outstr = pp.aper.aper(image,x0,y0,phpadu=gain,apr=fwhm,zeropoint=zp,skyrad=[3*fwhm,4*fwhm],exact=False)
 
        if opts.doSubtraction:
            hdulist=fits.open(imagefilesub)
            image = hdulist[0].data
            header = hdulist[0].header

            skys, skyerrs = [], []
            for x0, y0 in zip(x, y):
                mag,magerr,flux,fluxerr,sky,skyerr,badflag,outstr = pp.aper.aper(image,x0,y0,phpadu=gain,apr=fwhm,zeropoint=zp,skyrad=[3*fwhm,4*fwhm],exact=False)
                skys.append(sky)
                skyerrs.append(skyerr)
            skys, skyerrs = np.array(skys), np.array(skyerrs)
            idx = np.where(skyerrs>0)[0]
            skymed, skyerrmed = np.median(skys[idx]), np.median(skyerrs[idx])

            w = WCS(header)
            x0,y0 = w.wcs_world2pix(ra_transient,dec_transient,1)

            xnew, ynew = pp.cntrd.cntrd(image,np.array([x0]),np.array([y0]),3*fwhm,keepcenter=True)
            if not xnew == -1:
                xnew, ynew = xnew[0], ynew[0]
                x0, y0 = xnew, ynew
            mag_psf_sub,mag_psf_err_sub,flux,fluxerr,sky,skyerr,badflag,outstr = pp.aper.aper(image,x0,y0,phpadu=gain,apr=fwhm,zeropoint=zp,skyrad=[3*fwhm,4*fwhm],exact=False,badpix=[0,np.inf])

    # crossmatch to PS1
    ind_im,match_im = do_crossmatch(np.array(ra),np.array(dec),ps1_table['RAJ2000'],ps1_table['DEJ2000'])
 
    if opts.doDifferential: 
        f = filt[0]+'mag'
        ZP_median = -mag_psf_diff+ps1_table_diff[f][0]
        ZP_err = mag_psf_err_diff
    else:
        plotName = os.path.join(outputDir,'zp_%s_%d.pdf' % (filt, i)) 
        fig = plt.figure(figsize=(12,10))
        # getting ZP
        ZP_mag, ps1_mag = [], []
        for ii in range(len(ind_im)):
            if ind_im[ii] < len(ps1_table) and np.abs(mag_psf_err[ii])<0.1 :
                f = filt[0]+'mag'
                thismag = -mag_psf[ii]+ps1_table[f][ind_im[ii]]
                if (thismag > 1e10) or np.isnan(thismag): continue
                if (ps1_table[f][ind_im[ii]] < 15): continue
                ZP_mag.append(thismag)
                ps1_mag.append(ps1_table[f][ind_im[ii]])    
    
                plt.plot(mag_psf[ii], ps1_table[f][ind_im[ii]], 'kx')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig(plotName)
        plt.close()
    
        ZP_low, ZP_median, ZP_high = np.percentile(ZP_mag, [16, 50, 84])
        ZP_err = np.mean([ZP_high-ZP_median, ZP_median-ZP_low])      
        ZP_median = np.min(ZP_median) 
    
        print('number of standards used:', len(ZP_mag))
        print('ZP = ',np.round(ZP_median,3),'+-',np.round(ZP_err,3))#,np.mean(ZP_mag)-np.median(ZP_mag))

    # Get Transient and save it in KPED list
    ind_tra,match_tra = do_crossmatch(ra_transient,dec_transient,ra,dec)
    if opts.doSubtraction:
        ra_sub[0] = ra_transient
        dec_sub[0] = dec_transient
        ind_tra_sub,match_tra_sub = do_crossmatch(ra_transient,dec_transient,ra_sub,dec_sub)

    maglim = ZP_median+np.max(mag_psf[mag_psf_err<0.2])
    if opts.doSubtraction:
        if not opts.doForcedPhotometry:
            mag_psf_sub = mag_psf_sub[ind_tra_sub][0]
            mag_psf_err_sub = mag_psf_err_sub[ind_tra_sub][0]

        if match_tra_sub.any():
            mag_psf_err_tot = np.sqrt(mag_psf_err_sub**2 + ZP_err**2)
            print('LCO '+filt+' = '+str(np.round((mag_psf_sub+ZP_median),3))+' +- '+str(np.round(mag_psf_err_tot,2)))
            KPED.append([imagefile.split("/")[-1],1,filt,np.round(mag_psf_sub+ZP_median,3),np.round(mag_psf_err_tot,2)])
            print(name_[i],date_jd_[i],ra_all_tran[i],dec_all_tran[i],filt[0],'=',np.round(mag_psf_sub+ZP_median,3),'+-',str(np.round(mag_psf_err_tot,2)),np.round(maglim,2))
        else: 
            print('No transient detected',filt,'  = ',np.round(maglim,3))
            KPED.append([imagefile.split("/")[-1],0,filt,np.round(maglim,3),-99])
            print(name_[i],date_jd_[i],ra_transient,dec_transient,filt[0],'>',np.round(maglim,2),'+-',str(np.round(ZP_err,2)))
    else:
        if match_tra.any():
            print('LCO '+filt+' = '+str(np.round((mag_psf[ind_tra]+ZP_median)[0],3))+' +- '+str(np.round(mag_psf_err[ind_tra],2)))
            KPED.append([imagefile.split("/")[-1],1,filt,np.round((mag_psf[ind_tra]+ZP_median)[0],3),np.round(mag_psf_err[ind_tra],2)])
            print(name_[i],date_jd_[i],ra_all_tran[i],dec_all_tran[i],filt[0],'=',np.round((mag_psf[ind_tra][0]+ZP_median),3),'+-',str(np.round(mag_psf_err[ind_tra][0],2)),np.round(maglim,2))
        else: 
            print('No transient detected',filt,'  = ',np.round(maglim,3))
            KPED.append([imagefile.split("/")[-1],0,filt,np.round(maglim,3),-99])
            print(name_[i],date_jd_[i],ra_transient,dec_transient,filt[0],'>',np.round(maglim,2),'+-',str(np.round(ZP_err,2)))
  
    plotName = os.path.join(outputDir,'image_%s_%d.pdf' % (filt, i))
    fig = plt.figure(figsize=(12,10))
    f1 = aplpy.FITSFigure(imagefile,figure=fig)
    #f1.show_grayscale(invert=True)
    f1.show_grayscale(stretch='power',pmin=5,pmax=95)
    #idx = np.where(mag_psf_err<0.05)[0]
    f1.show_circles(ra,dec,1.5/3600.0,zorder=99,linestyle='dashed', edgecolor='red')
    f1.show_circles(ps1_table['RAJ2000'],ps1_table['DEJ2000'],3.0/3600.0,zorder=99,linestyle='dashed', edgecolor='green')
    f1.show_circles(ra_transient,dec_transient,1.5/3600.0,zorder=99,linestyle='dashed', edgecolor='blue')
    f1.show_circles(gaia_table['RA_ICRS'],gaia_table['DE_ICRS'],1.5/3600.0,zorder=99,linestyle='dashed', edgecolor='cyan')
    f1.axis_labels.set_xtext('Right Ascension')
    f1.axis_labels.set_ytext('Declination')
    f1._ax1.set_xlim(900,1100)
    f1._ax1.set_ylim(900,1100)
    f1.show_grayscale(stretch='power',pmin=5,pmax=95)
    f1.add_colorbar()
    f1.colorbar.set_location('right')
    fig.canvas.draw()
    plt.savefig(plotName)
    plt.close()

    if opts.doSubtraction:

        refimage = os.path.join(outputDir,'ref_%s.fits' % filt[0])

        plotName = os.path.join(outputDir,'image_ref_%s.pdf' % (filt))
        fig = plt.figure(figsize=(12,10))
        f1 = aplpy.FITSFigure(refimage,figure=fig)
        f1.show_grayscale(stretch='power',pmin=5,pmax=95)
        f1.show_circles(ra,dec,1.5/3600.0,zorder=99,linestyle='dashed', edgecolor='red')
        f1.show_circles(ps1_table['RAJ2000'],ps1_table['DEJ2000'],3.0/3600.0,zorder=99,linestyle='dashed', edgecolor='green')
        f1.show_circles(ra_transient,dec_transient,1.5/3600.0,zorder=99,linestyle='dashed', edgecolor='blue')
        f1.show_circles(gaia_table['RA_ICRS'],gaia_table['DE_ICRS'],1.5/3600.0,zorder=99,linestyle='dashed', edgecolor='cyan')
        f1.axis_labels.set_xtext('Right Ascension')
        f1.axis_labels.set_ytext('Declination')
        xpix, ypix = f1.world2pixel(ra_transient,dec_transient)
        f1._ax1.set_xlim(xpix-100,xpix+100)
        f1._ax1.set_ylim(ypix-100,ypix+100)
        f1.show_grayscale(stretch='power',pmin=5,pmax=95)
        f1.add_colorbar()
        f1.colorbar.set_location('right')
        fig.canvas.draw()
        plt.savefig(plotName)
        plt.close()

        plotName = os.path.join(outputDir,'image_sub_%s_%d.pdf' % (filt, i)) 
        fig = plt.figure(figsize=(12,10))
        f1 = aplpy.FITSFigure(imagefilesub,figure=fig)
        f1.show_grayscale(stretch='power',pmin=5,pmax=95)
        f1.show_circles(ra_sub,dec_sub,1.5/3600.0,zorder=99,linestyle='dashed', edgecolor='red')
        f1.show_circles(ps1_table['RAJ2000'],ps1_table['DEJ2000'],3.0/3600.0,zorder=99,linestyle='dashed', edgecolor='green')
        f1.show_circles(ra_transient,dec_transient,1.5/3600.0,zorder=99,linestyle='dashed', edgecolor='blue')
        f1.axis_labels.set_xtext('Right Ascension')
        f1.axis_labels.set_ytext('Declination')
        xpix, ypix = f1.world2pixel(ra_transient,dec_transient)
        f1._ax1.set_xlim(xpix-100,xpix+100)
        f1._ax1.set_ylim(ypix-100,ypix+100)
        f1.show_grayscale(stretch='power',pmin=5,pmax=95)
        f1.add_colorbar()
        f1.colorbar.set_location('right')
        fig.canvas.draw()
        plt.savefig(plotName)
        plt.close()
   
        hdu = fits.open(imagefilesub)
        data = hdu[0].data
        header = hdu[0].header
        w = WCS(header)
        # Make the cutout, including the WCS
        size = (36, 36)
        sizex, sizey = int(size[0]/2.0), int(size[1]/2.0)
        x, y = np.arange(-sizex, sizex), np.arange(-sizey, sizey)
        angles = np.linspace(0,2*np.pi,500)

        x0, y0 = w.wcs_world2pix(ra_transient,dec_transient,1)
        cutout = Cutout2D(data, position=(x0, y0), size=size, wcs=w)
        X, Y = np.meshgrid(x,y)
        plotName = os.path.join(outputDir,'cutout_%s_%d.pdf' % (filt, i))
        fig = plt.figure(figsize=(10,10))
        c = plt.pcolor(X,Y,cutout.data)
        fig.colorbar(c)
        radii = [2*fwhm,3*fwhm,4*fwhm]
        for radius in radii:
            plt.plot(radius*np.cos(angles),radius*np.sin(angles),'k')
        plt.savefig(plotName)
        plt.close()

        if opts.doDifferential:
            hdu = fits.open(imagefile)
            data = hdu[1].data
            header = hdu[1].header

            w = WCS(header)
            x0, y0 = w.wcs_world2pix(opts.ra_diff,opts.dec_diff,1)
            cutout = Cutout2D(data, position=(x0, y0), size=size, wcs=w)
            cutout = cutout.data - np.median(cutout.data)
            sizex, sizey = int(cutout.shape[0]/2.0), int(cutout.shape[1]/2.0)
            x, y = np.arange(-sizex, sizex), np.arange(-sizey, sizey)
            Y, X = np.meshgrid(x,y)

            plotName = os.path.join(outputDir,'cutout_diff_%s_%d.pdf' % (filt, i))
            fig = plt.figure(figsize=(10,10))
            c = plt.pcolor(X,Y,cutout.data,vmin=-100,vmax=100)
            fig.colorbar(c)
            radii = [2*fwhm,3*fwhm,4*fwhm]
            for radius in radii:
                plt.plot(radius*np.cos(angles),radius*np.sin(angles),'k')
            plt.savefig(plotName)
            plt.close()
 
    print(' \n')
    
    filename = os.path.join(outputDir,'phot.dat')
    fid = open(filename, 'w')
    for image, det, filt, mag, magerr in KPED:
        fid.write('%s %d %s %.5f %.5f\n' % (image, det, filt, mag, magerr))
    fid.close()
    
    
    # np.save(kped_saved,KPED)
    
    # print('Results saved to ',kped_saved)
