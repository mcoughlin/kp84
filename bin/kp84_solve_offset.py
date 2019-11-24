#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 08:53:19 2019

@author: yuhanyao
"""
import os
import time
import pytz
import datetime
import optparse
import numpy as np
from astral import Location
from scipy.ndimage import median_filter

import astropy.io.fits as fits
from astropy.time import Time
from astropy.stats import sigma_clipped_stats
from astropy.visualization import SqrtStretch, LogStretch, LinearStretch
from astropy.visualization.mpl_normalize import ImageNormalize

import warnings
from photutils import DAOStarFinder

import matplotlib
import matplotlib.pyplot as plt
fs = 14
matplotlib.rcParams.update({'font.size': fs})


def parse_commandline():   
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser()

    parser.add_option("--doUseLatest",  action="store_true", default=False)
    parser.add_option("--doApplyRates",  action="store_true", default=False)
    parser.add_option("--doLoop",  action="store_true", default=False)

    parser.add_option("-f", "--filename",default='/Users/yuhanyao/Desktop/ZTFJ05381953/ZTFJ05381953_9_g/registration/kped_20191112_074018_ZTFJ05381953_cl_o.fits.fz')
    parser.add_option("--idstart",default=1,type=int)
    parser.add_option("--idend",default=-1,type=int)

    opts, args = parser.parse_args()

    return opts

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
    

class KPEDtube(object):
    
    def __init__(self, tubelocalpath, idstart = 1, idend = 100,
                 figdir = '/Users/yuhanyao/Desktop/ZTFJ05381953/ZTFJ05381953_9_g/figures/'):
        allhdu = fits.open(tubelocalpath)
        self.tubelocalpath = tubelocalpath
        self.filter = allhdu[0].header["FILTER"]
        self.filename = allhdu[0].header["FILENAME"]
        self.figdir = figdir
        #if idstart==-1 or idend==-1:
        #    idstart = len(allhdu)-1-100
        #    idend = len(allhdu)-1  
        self.idstart = idstart
        self.idend = idend
        self.hdus = allhdu[idstart:idend]
        self.nimg = len(self.hdus)
        
        
    def identify_bright_source(self):
        """
        hdus = kptube.hdus
        nimg = kptube.nimg
        """
        hdus = self.hdus
        nimg = self.nimg
        
        starpositions = {}
        utcs = {}
        median_size = 40
        for i in range(nimg):
            hdu = hdus[i]
            header = hdu.header
            utc = header["GPS_TIME"]
            data = hdu.data
            # only calculate global median once
            if i==0:
                data_median = np.asfarray(median_filter(data, size=(median_size, median_size)))
            data1 = data - data_median
            
            mean, median, std = sigma_clipped_stats(data1, sigma=3.0)  
            daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)  
            sources = daofind(data1 - median)  
            arg = np.argsort(sources["peak"])[::-1]
            sources = sources[arg]
            
            positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
            """
            figdir = self.figdir
            cmap_name = 'viridis'
            if i%50==0:
                print (i)
                fig, ax = plt.subplots(2, 2, figsize=(12, 12))
                norm = ImageNormalize(stretch=LinearStretch())
                ax[0,0].imshow(data_median, cmap = cmap_name, origin='lower', norm=norm)
                ax[0,0].set_title('Median', fontsize=fs+1)
            
                ax[0,1].imshow(data, cmap = cmap_name, origin='lower', norm=norm)
                ax[0,1].set_title('Data', fontsize=fs+1)
            
                norm = ImageNormalize(stretch=LinearStretch())
                ax[1,0].imshow(data1, cmap = cmap_name, origin='lower', zorder = 1, norm =norm)
                ax[1,0].set_title('Data - Median', fontsize=fs+1)
                for j in range(len(positions)):
                    xpos = positions[j,0]
                    ypos = positions[j, 1]
                    width = 7
                    lw = 1
                    if j<5:
                        alpha = 1
                    else:
                        alpha = 0.5
                    ax[1,0].plot([xpos-width, xpos-width], [ypos-width, ypos+width], 
                              'r-', linewidth = lw, alpha=alpha)
                    ax[1,0].plot([xpos+width, xpos+width], [ypos-width, ypos+width], 
                              'r-', linewidth = lw, alpha=alpha)
                    ax[1,0].plot([xpos-width, xpos+width], [ypos-width, ypos-width], 
                              'r-', linewidth = lw, alpha=alpha)
                    ax[1,0].plot([xpos-width, xpos+width], [ypos+width, ypos+width], 
                              'r-', linewidth = lw, alpha=alpha)
                
                ax[1,1].axis('off')
                plt.savefig(figdir + "%s"%(i+1)+".pdf")
                plt.close()
            """
            starpositions[i+1] = positions
            utcs[i+1] = utc
        self.starpositions = starpositions
        self.utcs = utcs
    
    
    def get_offset(self):
        # get time in the unit of second
        starpositions = self.starpositions
        utcs = self.utcs
        xthre = 5
        ythre = 3
        # it seems easier to move in the x-direction
        
        keys = np.array(list(utcs.keys()))
        times = np.zeros(len(keys))
        for (i,key) in enumerate(keys):
            times[i] = Time(utcs[key], format="isot").mjd
        times -= times[0]
        self.times = times * 86400.0
        # get position
        # choose the image with the greatest number of stars identified
        nstars = np.array([len(starpositions[key]) for key in keys])
        nstar_max = np.max(nstars)
        id_best = np.where(nstars == nstar_max)[0][0]
        keybest = keys[id_best]
        posbest = starpositions[keybest]
        # register the top 10 brightest star if nstar_max >= 10
        # else, register nstar_max stars
        
        nregister = min(nstar_max, 10)
        self.nregister = nregister
        eachstarpos = {}
        for j in range(nregister):
            eachstarpos[j] = {}
            eachstarpos[j]["x"] = {}
            eachstarpos[j]["x"][keybest] = posbest[j][0]
            eachstarpos[j]["y"] = {}
            eachstarpos[j]["y"][keybest] = posbest[j][1]
        starkeys = list(eachstarpos.keys())
        
        keys_good = keys[nstars > nstar_max//5]
        myid = np.where(keys_good==keybest)[0][0]
        keys_arranged = np.hstack([keys_good[myid:], keys_good[:myid]])
        for key in keys_arranged:
            if key==keybest:
                continue
            posnow = starpositions[key]
            for starkey in starkeys:
                starpos = eachstarpos[starkey]
                starx = starpos["x"]
                stary = starpos["y"]
                
                pastkeys = np.array(list(starx.keys()))
                nearbykey = pastkeys[np.argsort(abs(pastkeys - key))[0]]
                xx = starx[nearbykey]
                yy = stary[nearbykey]
                
                xdiff = abs(posnow[:,0] - xx)
                ydiff = abs(posnow[:,1] - yy)
                xselect = np.where(xdiff < xthre)[0]
                yselect = np.where(ydiff < ythre)[0]
                select = -99
                for elem in xselect:
                    if elem in yselect:
                        select = elem
                if select == -99:
                    # this bright star is not in the currect image
                    continue
                xnow = posnow[select, 0]
                ynow = posnow[select, 1]
                starx[key] = xnow
                stary[key] = ynow
        self.eachstarpos = eachstarpos
        
    def cal_offset(self, plot_figure = False):
        pixelscale = 0.259 # arcsec per pixel
        
        nregister = self.nregister
        eachstarpos = self.eachstarpos
        times = self.times
        
        xdots = np.zeros(nregister)
        ydots = np.zeros(nregister)
        
        if plot_figure == True:
            plt.figure(figsize= (12,6))
            ax1 = plt.subplot(121)
            ax2 = plt.subplot(122)
        
        for i in range(nregister):
            mystar = eachstarpos[i]
            starxs = mystar["x"]
            starys = mystar["y"]
            ind = np.array(list(starxs.keys()))-1
            mytime = times[ind]
            xs = np.array(list(starxs.values()))
            ys = np.array(list(starys.values()))
            ind = np.argsort(mytime)
            mytime = mytime[ind]
            xs = xs[ind]
            ys = ys[ind]
            dxdt, tmp1x, tmp2x = mylinear_fit(mytime, xs, np.ones(len(xs)), npar = 2)
            dydt, tmp1y, tmp2y = mylinear_fit(mytime, ys, np.ones(len(xs)), npar = 2)
            
            xdots[i] = dxdt
            ydots[i] = dydt
            if plot_figure == True:
                ax1.plot(mytime, xs, '.')
                ax2.plot(mytime, ys, '.')
        
        xdot = np.median(xdots) * pixelscale
        ydot = np.median(ydots) * pixelscale
        if plot_figure == True:
            ax1.set_title("dx/dt = %.2f"%(xdot*1e+3)+r"$\times 10^{-3}$"+"arcsec / s", fontsize= fs)
            ax2.set_title("dy/dt = %.2f"%(ydot*1e+3)+r"$\times 10^{-3}$"+"arcsec / s", fontsize= fs)
            ax1.set_xlabel("time (s)")
            ax1.set_ylabel("x")
            ax2.set_ylabel("y")
            
        self.xdot = xdot
        self.ydot = ydot

def all_files_under(path):
    """Iterates through all files that are under the given path."""
    for cur_path, dirnames, filenames in os.walk(path):
        if not "20" in cur_path:
            continue 
        for filename in filenames:
            if not "fits.fz" in filename:
                continue
            yield os.path.join(cur_path, filename)
            

def get_kp84_sunrise_time(tnow):
    """give current time"""
    l = Location()
    l.name = 'KP84'
    l.region = 'Kitt Peak Observatory'
    l.latitude = 31.9599 # N
    l.longitude = -111.5997 # in east  (111.5997 w)
    l.timezone = "US/Arizona"
    l.elevation = 2099 # in meters
    if tnow.hour>12:
        sun = l.sun(date = datetime.date(tnow.year, tnow.month, tnow.day)+ datetime.timedelta(days=1))
    else:
        sun = l.sun(date = datetime.date(tnow.year, tnow.month, tnow.day))
    tsunrise = sun['sunrise']
    return tsunrise


opts = parse_commandline()
idstart = opts.idstart
idend = opts.idend

if opts.doLoop is not True:

    if opts.doUseLatest:
        # filename convention: kped_yyyymmdd_hhmmss_objectname_cl_o.fits.fz
        # hhmmss is ut time
        filename = max(all_files_under('/Data/'), key=os.path.getmtime)
    else:
        filename = opts.filename 
 
    kptube = KPEDtube(tubelocalpath = filename, idstart = idstart, idend = idend)
    warnings.filterwarnings('ignore', category=UserWarning, append=True)
    kptube.identify_bright_source()
    kptube.get_offset()
    kptube.cal_offset()

    xdot = kptube.xdot
    ydot = kptube.ydot

    print ("dra/dt [arcsec/s] = %.5f, ddec/dt [arcsec/s] = %.5f"% (xdot, ydot))

    if opts.doApplyRates:
        system_command = "ssh tcs@sells.kpno.noao.edu 'tx track ra=%.5f dec=%.5f'" % (xdot, ydot)
        os.system(system_command)

else:
    tz = pytz.timezone("US/Arizona")
    tnow = datetime.datetime.now(tz)
    tsunrise = get_kp84_sunrise_time(tnow)

    while tnow < tsunrise:
        
        filename = max(all_files_under('/Data/'), key=os.path.getmtime)
        
        kptube = KPEDtube(tubelocalpath = filename, idstart = idstart, idend = idend)
        warnings.filterwarnings('ignore', category=UserWarning, append=True)
        kptube.identify_bright_source()
        kptube.get_offset()
        kptube.cal_offset()
        
        xdot = kptube.xdot
        ydot = kptube.ydot

        print ("dra/dt [arcsec/s] = %.5f, ddec/dt [arcsec/s] = %.5f"% (xdot, ydot))

        if opts.doApplyRates:
            system_command = "ssh tcs@sells.kpno.noao.edu 'tx track ra=%.5f dec=%.5f'" % (xdot, ydot)
            os.system(system_command)
            
        #last_filename = filename
        time.sleep(60)
        tnow = datetime.datetime.now(tz)
    