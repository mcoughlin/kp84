#!/usr/bin/env python
import os, sys, optparse, shutil, glob
import operator
import numpy as np
import astropy.io.ascii as asci
from astropy.io import fits
# sys.path.append("/Users/yuhanyao/Documents/GitHub/kp84/")
# sys.path.append("/Users/yuhanyao/Documents/GitHub/PythonPhot/")
from astropy.table import vstack
sys.path.append("/home/roboao/Michael/kp84/")
from citizen.reduction_utils import stack_images, register_images
from citizen.reduction_utils import get_wcs_xy, update_wcsstatus, forcedphotometry_kp
from citizen.reduction_utils import get_reference_pos, filter2filtstr, jd2hjd, jd2bjd
from citizen.visualize_utils import makemovie
#from astroML.crossmatch import crossmatch_angular

import matplotlib
matplotlib.use('agg')
matplotlib.pyplot.switch_backend('agg')
fs = 14
matplotlib.rcParams.update({'font.size': fs})
from matplotlib import pyplot as plt


def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser()

    parser.add_option("--day",default="20191116", help="date of observation")
    parser.add_option("--objName",default="ZTFJ19015309", help="object name")
    
    parser.add_option("-x","--xstar", default="", help="object's x pixel of frame 1, only give if astrometry.net solution failed")
    parser.add_option("-y","--ystar", default="", help="object's y pixel of frame 1, only give if astrometry.net solution failed")
    parser.add_option("--xyext", default="", help="frame number of the Multi-extension Cube that xstar and ystar are specified")
    parser.add_option("--xyfile", default="", 
                      help="Name of the Multi-extension Cube file, need to be specified if there are multiple files for a single object")
    
    parser.add_option("--xoffref",default=0, type=float, help="x_ref - x_sci")
    parser.add_option("--yoffref",default=0, type=float, help="y_ref - y_sci")
    parser.add_option("--refmag",default=0, type=float, help="AB magnitude of the reference star")
    
    parser.add_option("--doSubtractBackground",  action="store_true", default=False)
    parser.add_option("--doSkipRegis",  action="store_true", default=False)
    parser.add_option("--doSkipFindRef",  action="store_true", default=False)
    
    # stack options
    parser.add_option("-n","--nimages",default=1,type=int, help="see --doStack")
    parser.add_option("--doStack",  action="store_true", default=False, help="stack each --nimages together")
    
    parser.add_option("--doOffRefit",  action="store_true", default=False, help='If true, do not refit after tuning') 
    parser.add_option("--doOffTune",  action="store_true", default=False, help='If true, do not tune the registration') 
    parser.add_option("--maxdist",default=10.0, type=float, help="during tuning, maximum distance in pixel numbers")
    
    parser.add_option("--aper_size",default=10.0, type=float, help="aperture size in pixel numbers for sci object")
    parser.add_option("--sky_in",default=30.0, type=float, help="inner sky annulus radius in pixel numbers for sci object")
    parser.add_option("--sky_out",default=50.0, type=float, help="outer sky annulus radius in pixel numbers for sci object")
    parser.add_option("--aper_size_ref",default=10.0, type=float, help="aperture size in pixel numbers for ref object")
    parser.add_option("--sky_in_ref",default=30.0, type=float, help="inner sky annulus radius in pixel numbers for ref object")
    parser.add_option("--sky_out_ref",default=50.0, type=float, help="outer sky annulus radius in pixel numbers for ref object")
    
    parser.add_option("--moviemode", default=1,type=int, help="if not 0, then make a movie")
    
    # Transient options
    #parser.add_option("--doSubtraction",  action="store_true", default=False, help="subtract source extractor background from science image before force photometry")
    #parser.add_option("--subtractionSource", default="ps1", 
    #                  choices=('sdss', 'ps1'), help="source of the reference image")
    
    opts, args = parser.parse_args()
    return opts


"""
day = "20191222"
objName = "ZTFJ0538+1953"
setupDir = "/Users/yuhanyao/Desktop/kped_tmp/"
#xstar = "382*360"
#ystar = "297*286"
#xyext = "174*1"
#xyfile = "110815_ZTFJ11514412_cl_o*121311_ZTFJ11514412_cl_o"

day = "20200105"
objName = "ZTFJ0538+1953"
xstar = "267*256*234*191*134"
ystar = "341*341*343*351*353"
xyext = "499*86*519*548*84"
xyfile = "033949_ZTFJ0538+1953_4_cl_o*033949_ZTFJ0538+1953_4_cl_o_0000*054524_ZTFJ0538+1953_5_cl_o*074541_ZTFJ0538+1953_5_cl_o*084543_ZTFJ0538+1953_5_cl_o_0000"
xoffref = 7.23
yoffref = -48
refmag = 14.838	# SDSS DR12
# ra and dec of ref star: 84.5102826845 19.8772905537
"""
# Parse command line
opts = parse_commandline()
setupDir = "/Data3/archive_kped/data/reductions/"
day = opts.day
objName = opts.objName

inputDir = "../input"
dataDir = os.path.join(setupDir, day, objName) # the setup directory of this object
outputDir = os.path.join("../output", day, objName) # the output directory of this object
outputProDir = os.path.join(outputDir, "product") 

doSubtractBackground = opts.doSubtractBackground
doSkipRegis = opts.doSkipRegis
doSkipFindRef = opts.doSkipFindRef
doStack = opts.doStack
nimages = opts.nimages
doOffTune = opts.doOffTune
doOffRefit = opts.doOffRefit

aper_size = opts.aper_size
sky_in = opts.sky_in
sky_out = opts.sky_out
aper_size_ref = opts.aper_size_ref
sky_in_ref = opts.sky_in_ref
sky_out_ref = opts.sky_out_ref

moviemode = opts.moviemode

maxdist = opts.maxdist
xstar = opts.xstar
ystar = opts.ystar
xyext = opts.xyext
xyfile = opts.xyfile
xoffref = opts.xoffref
yoffref = opts.yoffref
refmag = opts.refmag

xyfiles_ = np.array(xyfile.split("*"))
xyexts_ = np.array(xyext.split("*"))
xs_ = np.array(xstar.split("*"))
ys_ = np.array(ystar.split("*"))

if not os.path.isdir(outputDir):
    os.makedirs(outputDir)
if not os.path.isdir(outputProDir):
    os.makedirs(outputProDir)

print ("")
print ("=================================")
print ("Getting the Coordinate of Object!")
print ("=================================")
print ("")
coofile = os.path.join(outputProDir, "coo.reg")
observedFile = "%s/observed.dat"%inputDir
lines = [line.rstrip('\n') for line in open(observedFile)]
lines = np.array(lines)
objs = []
ras = []
decs = []
for ii,line in enumerate(lines):
    lineSplit = list(filter(None,line.split(" ")))
    objs.append(lineSplit[0])
    ras.append(float(lineSplit[1]))
    decs.append(float(lineSplit[2]))
    
objs = np.array(objs)
ras = np.array(ras)
decs = np.array(decs)
if objName not in objs:
    print("%s missing from observed list, please add."%objName)
    exit(0)
else:
    ind = np.where(objs==objName)[0][0]
    ra = ras[ind]
    dec = decs[ind]
print ("%s, ra=%.5f, dec=%.5f"%(objName, ra, dec))
print ("Saving coordinate to %s"%coofile)
np.savetxt(coofile, [ra, dec])

fitsfiles = sorted(glob.glob(os.path.join(dataDir,'processing','*.fits'))) 
if day == "20191222" and objName == "ZTFJ0538+1953":
    fitsfiles = sorted(glob.glob(os.path.join(dataDir,'processing','*_073459_*0001_proc.fits'))) 
tmpheader = fits.open(fitsfiles[0])[0].header
tmpfilter = tmpheader["FILTER"]
passband = filter2filtstr(tmpfilter)
bandfile = os.path.join(outputProDir, "filter.txt")
np.savetxt(bandfile, [passband], fmt="%s")

if not doSkipRegis:
    print ("")
    print ("=================================================")
    print ("Finding the object on each frame -- Registration!")
    print ("=================================================")
    print ("")
    nfiles = len(fitsfiles)

    for i in range(nfiles):
        flag = 1
        fitsfile = fitsfiles[i]
        fitsfileSplit = fitsfile.split("/")[-1].replace(".fits","").replace("_proc","")
        path_out_dir='%s/%s'%(outputDir,fitsfileSplit)
        if not os.path.isdir(path_out_dir):
            os.makedirs(path_out_dir)
        print("%d/%d: %s"%(i+1, nfiles, fitsfileSplit))
        if fitsfileSplit[14:] in xyfiles_:
            index = np.where(xyfiles_ == fitsfileSplit[14:])[0][0]
            xyframe = int(xyexts_[index])
            x = int(xs_[index])
            y = int(ys_[index])
            update_wcsstatus(fitsfile, xyframe)
            print ("  Using the specified location: x = %d, y = %d, xyext = %d"%(x, y, xyframe))
            flag = 2
        else:
            wcsfiles = glob.glob(os.path.join(dataDir,'wcs', '%s*wcs.fits'%fitsfileSplit))
            if len(wcsfiles)!=0:
                print ("  Finding it on the wcs file extension...")
                # The astrometry is successfully found:
                wcsfile = wcsfiles[0]
                x, y, xyframe = get_wcs_xy(ra, dec, wcsfile, fitsfile, get_distance = True)
                print ("  Found the location with wcs solution: x = %d, y=%d, xyext = %d"%(x, y, xyframe))
            else: # astrometry.net solution failed:
                flag = 0
                print ("  Astrometry failed and no wcs specified -- discard %s!"%fitsfileSplit)
        if flag in [1,2]:
            shiftfiles = glob.glob(os.path.join(dataDir,'registration', '%s_proc_shift.dat'%fitsfileSplit))
            #print (shiftfiles)
            shiftfile = shiftfiles[0]
            refit = operator.not_(doOffRefit)
            regis2HDU = register_images(fitsfile, shiftfile, xyframe, x, y, path_out_dir, 
                                        maxdist = maxdist, aper_size = aper_size, 
                                        refit = refit, offtune = doOffTune)
            regisfile = fitsfile.replace("/processing/", "/registration/")
            regisfile = regisfile[:-5]+'_regis.fits'
            regis2HDU[0].header["OBJNAME"] = objName
            regis2HDU[0].header["RA_OBJ"] = ra
            regis2HDU[0].header["DEC_OBJ"] = dec
            print ("  Writing to %s"%regisfile)
            regis2HDU.writeto(regisfile, overwrite=True)
            print ("")

fitsfiles = sorted(glob.glob(os.path.join(dataDir,'registration','*_regis.fits'))) 
if day == "20191222" and objName == "ZTFJ0538+1953":
    fitsfiles = sorted(glob.glob(os.path.join(dataDir,'registration','*_073459_*_regis.fits'))) 
nfiles = len(fitsfiles)

# yyao: This call is to be checked
if doStack:
    stackDir = os.path.join(path_out_dir,'stack')
    if not os.path.isdir(stackDir):
        os.makedirs(stackDir)
    if nimages > 1:
        stack_images(stackDir,fitsfiles,opts.nimages,doRegistration=doRegistration,
                     registration_size=registration_size,x=x0,y=y0) 
        fitsfiles = sorted(glob.glob(os.path.join(stackDir,'*.fits*')))
    else:
        print("You asked to stack but with --nimages 1... passing.")


if not doSkipFindRef:
    print ("")
    print ("========================")
    print ("Find the Reference Star!")
    print ("========================")
    print ("")
    for ii in range(nfiles):
        fitsfile = fitsfiles[ii]
        _fitsfileSplit = fitsfile.split("/")[-1].replace("_regis","")
        fitsfileSplit = fitsfile.split("/")[-1].replace("_proc_regis","").replace(".fits","")
        print("%d/%d: %s"%(ii+1, nfiles, fitsfileSplit))
        path_out_dir='%s/%s'%(outputDir, fitsfileSplit)
        folderName_sextract = "%s/sextract"%(dataDir)
        catfile = os.path.join(folderName_sextract, _fitsfileSplit.replace(".fits",".cat"))
        backfile = os.path.join(folderName_sextract, _fitsfileSplit.replace(".fits",".background.fits"))
        
        cat = np.loadtxt(catfile)
        if not cat.size: 
            continue
        
        sciHDU = get_reference_pos(fitsfile, cat, objName, zp=0, passband=passband, 
                                   xoff=xoffref, yoff=yoffref, refmag=refmag)
        sciHDU.writeto(fitsfile, overwrite=True)
    
        scienceimage = '%s/science.fits'%(path_out_dir)
        sciencebkgimage = '%s/science.background.fits'%(path_out_dir)
        shutil.copy(fitsfile, scienceimage)
        shutil.copy(backfile, sciencebkgimage)
    
        if doSubtractBackground:
            print ("   Subtracting Background Before Forced Photometry...")
            scisubimage = '%s/science.subtracted.fits'%(path_out_dir)
            sciHDU = fits.open(scienceimage)
            bkgHDU = fits.open(sciencebkgimage)
            for k in range(len(sciHDU)-1):
                kk = k+1
                data0 = sciHDU[kk].data
                data1 = bkgHDU[kk].data
                sciHDU[k].data = data0-data1
                #plt.imshow(data0)
                #plt.imshow(data1)
                #plt.imshow(sciHDU[k].data)
            sciHDU[0].header.add_history('Background subtracted -- using the Source Extractor bkg')
            sciHDU.writeto(scisubimage, overwrite=True)
        

print ("")
print ("==========================")
print ("Perform Forced Photometry!")
print ("==========================")
print ("")

def differential_phot(scienceimage, aper_size, sky_in, sky_out,
                      aper_size_ref, sky_in_ref, sky_out_ref, refmag,
                      forcedfile, figname, profigname):
    mjd_forced, mag_forced, magerr_forced, flux_forced, fluxerr_forced = \
            forcedphotometry_kp(scienceimage, aper_size=aper_size, 
                                sky_in = sky_in, sky_out = sky_out,
                                xkey = "X_OBJ", ykey = "Y_OBJ")
    mjd_forced, mag_forced_field, magerr_forced_field, flux_forced_field, fluxerr_forced_field = \
            forcedphotometry_kp(scienceimage, aper_size=aper_size_ref, 
                                    sky_in = sky_in_ref, sky_out = sky_out_ref,
                                    xkey = "X_FIELD", ykey = "Y_FIELD")
    mag = mag_forced - mag_forced_field + refmag
    magerr = np.sqrt(magerr_forced**2 + magerr_forced_field**2)
    flux = flux_forced/flux_forced_field
    fluxerr = flux*np.sqrt((fluxerr_forced/flux_forced)**2 + (fluxerr_forced_field/flux_forced_field)**2)
        
    fid = open(forcedfile,'w')
    for ii in range(len(mjd_forced)):
        fid.write('%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n'%(mjd_forced[ii],
                                                     mag[ii],magerr[ii],
                                                     flux[ii],fluxerr[ii],
                                                     mag_forced_field[ii], magerr_forced_field[ii], 
                                                     flux_forced_field[ii], fluxerr_forced_field[ii],
                                                     mag_forced[ii], magerr_forced[ii], 
                                                     flux_forced[ii], fluxerr_forced[ii]))
    fid.close()
    
    plt.figure(figsize=(12,6))
    mjd0 = mjd_forced[0]
    # plt.errorbar(mjd_forced-mjd0, mag_forced, magerr_forced, fmt='.k', label = "science")
    refmtemp = mag_forced_field- np.median(mag_forced_field) + refmag
    timetmp = (mjd_forced-mjd0)*24*60
    plt.errorbar(timetmp, refmtemp, magerr_forced_field, fmt='.b', label="reference")
    plt.errorbar(timetmp, mag, magerr, fmt=".r", label="object")
    plt.xlim(min(timetmp), max(timetmp))
    plt.legend()
    plt.title(fitsfileSplit, fontsize=fs)
    plt.gca().invert_yaxis()
    plt.xlabel("time in min")
    plt.tight_layout()
    plt.savefig(figname)
    plt.close()
    shutil.copy(figname, profigname)
            
    
for ii in range(nfiles):
    fitsfile = fitsfiles[ii]
    fitsfileSplit = fitsfile.split("/")[-1].replace("_proc_regis","").replace(".fits","")
    print("%d/%d: %s"%(ii+1, nfiles, fitsfileSplit))
    path_out_dir='%s/%s'%(outputDir, fitsfileSplit)
    scienceimage = '%s/science.fits'%(path_out_dir)
    scisubimage = '%s/science.subtracted.fits'%(path_out_dir)
    refmag = fits.open(scienceimage)[0].header["MAGREF"]
    
    print ("  do forced photometry on calibrated image...")
    forcedfile = '%s/science.forced'%(path_out_dir)
    figname = os.path.join(path_out_dir, "diff_phot.pdf")
    profigname = os.path.join(outputProDir, fitsfileSplit+"_diff_phot.pdf")
    differential_phot(scienceimage, aper_size, sky_in, sky_out,
                      aper_size_ref, sky_in_ref, sky_out_ref, refmag,
                      forcedfile, figname, profigname)
    
    if doSubtractBackground:
        print ("  do forced photometry on sky subtracted image...")
        forcedfile = '%s/science.subtracted.forced'%(path_out_dir)
        figname = os.path.join(path_out_dir, "subtracted_diff_phot.pdf")
        profigname = os.path.join(outputProDir, fitsfileSplit+"_subtracted_diff_phot.pdf")
        differential_phot(scisubimage, aper_size, sky_in, sky_out,
                      aper_size_ref, sky_in_ref, sky_out_ref, refmag,
                      forcedfile, figname, profigname)
    
print ("")
print ("=======================")
print ("Make Final Light Curve!")
print ("=======================")
print ("")
for ii in range(nfiles):
    fitsfile = fitsfiles[ii]    
    fitsfileSplit = fitsfile.split("/")[-1].replace("_proc_regis","").replace(".fits","")
    print("%d/%d: %s"%(ii+1, nfiles, fitsfileSplit))
    path_out_dir='%s/%s'%(outputDir, fitsfileSplit)
    forcedfile = '%s/science.forced'%(path_out_dir)
    
    if ii == 0:
        tblforced = asci.read(forcedfile,names=['MJD','mag','magerr','flux','fluxerr',
                                                "mag_forced_field", "magerr_forced_field",
                                                "flux_forced_field", "fluxerr_forced_field",
                                                "mag_forced", "magerr_forced",
                                                "flux_forced", "fluxerr_forced"])
    else:
        tbltemp = asci.read(forcedfile,names=['MJD','mag','magerr','flux','fluxerr',
                                                "mag_forced_field", "magerr_forced_field",
                                                "flux_forced_field", "fluxerr_forced_field",
                                                "mag_forced", "magerr_forced",
                                                "flux_forced", "fluxerr_forced"]) 
        tblforced = vstack([tblforced,tbltemp])
    
    if doSubtractBackground:
        forcedsubfile = '%s/science.subtracted.forced'%(path_out_dir)
        if ii == 0:
            tblsubforced = asci.read(forcedsubfile,names=['MJD','mag','magerr','flux','fluxerr',
                                                "mag_forced_field", "magerr_forced_field",
                                                "flux_forced_field", "fluxerr_forced_field",
                                                "mag_forced", "magerr_forced",
                                                "flux_forced", "fluxerr_forced"])
        else:
            tbltemp = asci.read(forcedsubfile,names=['MJD','mag','magerr','flux','fluxerr',
                                                "mag_forced_field", "magerr_forced_field",
                                                "flux_forced_field", "fluxerr_forced_field",
                                                "mag_forced", "magerr_forced",
                                                "flux_forced", "fluxerr_forced"]) 
            tblsubforced = vstack([tblsubforced,tbltemp])
 
def save_forced_tb(tblforced, finalforcefile, plotName):
    hjd = jd2hjd(tblforced["MJD"].data+2400000.5, ra, dec)
    bjd = jd2bjd(tblforced["MJD"].data+2400000.5, ra, dec)
    mbjd = bjd-2400000.5
    tblforced["HJD"] = hjd.value
    print ("Writing to %s"%finalforcefile)
    asci.write(tblforced, finalforcefile, overwrite=True)
    
    tblforced = asci.read(finalforcefile)
    mjd_forced = tblforced['MJD'].data
    mag_forced, magerr_forced = tblforced['mag'].data, tblforced['magerr'].data
    flux_forced, fluxerr_forced = tblforced['flux'].data, tblforced['fluxerr'].data
    
    print ("Plotting mag photometry...")
    try:
        timetmp = (mjd_forced-mjd_forced[0])*24
        plt.figure(figsize=(20,8))
        plt.errorbar(timetmp,mag_forced,magerr_forced,fmt='ko')
        plt.xlabel('Time [hrs]')
        plt.ylabel('Magnitude [ab]')
        idx = np.where(np.isfinite(mag_forced))[0]
        ymed = np.nanmedian(mag_forced)
        y10, y90 = np.nanpercentile(mag_forced[idx],10), np.nanpercentile(mag_forced[idx],90)
        ystd = np.nanmedian(magerr_forced[idx])
        ymin = y10 - 3*ystd
        ymax = y90 + 3*ystd
        plt.ylim([ymin,ymax])
        plt.xlim(min(timetmp), max(timetmp))
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(plotName)
        plt.close()
    except:
        print ("plot failed")

finalforcefile = os.path.join(outputProDir,"lightcurve.forced")
plotName = os.path.join(outputProDir,'mag_forced.pdf')
save_forced_tb(tblforced, finalforcefile, plotName)

if doSubtractBackground:
    finalsubforcefile = os.path.join(outputProDir,"lightcurve.subtracted.forced")
    plotsubName = os.path.join(outputProDir,'subtracted_mag_forced.pdf')
    save_forced_tb(tblsubforced, finalsubforcefile, plotsubName)

if moviemode!=0:
    print ("")
    print ("==================")
    print ("Making the Movie!")
    print ("==================")
    print ("")
    movieDir = os.path.join(outputDir,'movie')
    if not os.path.isdir(movieDir):
        os.makedirs(movieDir)
    makemovie(movieDir, fitsfiles, moviemode=moviemode, aper_size = aper_size, sky_in = sky_in, sky_out = sky_out)
    shutil.copy(movieDir+"/movie.mpg", outputProDir+"/movie.mpg")

