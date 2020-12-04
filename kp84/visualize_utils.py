#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:20:41 2019

@author: yuhanyao
"""
import os
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats, mad_std

import matplotlib
fs = 14
matplotlib.rcParams.update({'font.size': fs})
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from matplotlib.colors import LogNorm

from kp84.scheduler import load_targets


def make_triplet(fitsfiles, ra, dec):
    """
        Feed in alert packet
    """

    cutout_dict = dict()
    normalize = True
    registration_size = 64
    triplet = np.zeros((registration_size, registration_size, 3))

    for ii, fitsfile in enumerate(fitsfiles):
        hdul = fits.open(fitsfile)
        data = hdul[0].data
        header = hdul[0].header

        if "X_OBJ" in header:
            x0, y0 = header["X_OBJ"], header["Y_OBJ"]
        else:
            w0 = WCS(header)
            x0, y0 = w0.wcs_world2pix(ra, dec, 1)

        if normalize:
            data /= np.linalg.norm(data)

        #shape = data.shape
        #x, y = shape[0]/2, shape[1]/2

        xlow = int(x0 - registration_size/2)
        xhigh = int(x0 + registration_size/2)
        ylow = int(y0 - registration_size/2)
        yhigh = int(y0 + registration_size/2)        

        triplet[:, :, ii] = data[xlow:xhigh,ylow:yhigh]

    return triplet

def plot_triplet(tr, show_fig=True):
    fig = plt.figure(figsize=(8, 2), dpi=120)
    ax1 = fig.add_subplot(131)
    ax1.axis('off')
    mean, median, std = sigma_clipped_stats(tr[:, :, 0])
    ax1.imshow(tr[:, :, 0], vmin = median - 2*std, vmax = median + 3*std)
    #ax1.imshow(tr[:, :, 0], origin='upper', cmap=plt.cm.bone, norm=LogNorm())
    ax1.title.set_text('Science')
    ax2 = fig.add_subplot(132)
    ax2.axis('off')
    mean, median, std = sigma_clipped_stats(tr[:, :, 1])
    ax2.imshow(tr[:, :, 1], vmin = median - 2*std, vmax = median + 3*std)
    #ax2.imshow(tr[:, :, 1], origin='upper', cmap=plt.cm.bone, norm=LogNorm())
    ax2.title.set_text('Reference')
    ax3 = fig.add_subplot(133)
    ax3.axis('off')
    mean, median, std = sigma_clipped_stats(tr[:, :, 2])
    ax3.imshow(tr[:, :, 2], vmin = median - 2*std, vmax = median + 3*std)
    #ax3.imshow(tr[:, :, 2], origin='upper', cmap=plt.cm.bone)
    ax3.title.set_text('Difference')

    if show_fig:
        plt.show()
    return fig

def get_ra_dec_radius(objName, object_lists, wcsmode=1):
    if wcsmode == 1:
        targets = load_targets(object_lists)
        idx = np.where(objName == targets["objectID"])[0]
        if len(idx) == 0:
            ra = None
            dec = None
            radius = None
            print("%s missing from observation folder, please add."%objName)
        else:
            row = targets[idx][0]
            ra, dec = row["ra"], row["dec"]
            print ("%s, ra=%.5f, dec=%.5f"%(objName, ra, dec))
            radius = 0.5# 30 arcmin
    else:
        ra = None
        dec = None
        radius = None
    return ra, dec, radius


def makemovie(movieDir, fitsfiles, halfwidth=50, moviemode=1, aper_size=10, sky_in = 11, sky_out = 30):
    """
    Make a movie of a set of fits (each with multiple exposures)
    """
    from matplotlib.ticker import NullLocator
    fitsfiles = np.array(fitsfiles)
    nfiles = len(fitsfiles)
    mjds = np.zeros(nfiles)
    # sort the fits files by time
    for ii in range(nfiles):
        fitsfile = fitsfiles[ii]
        header = fits.open(fitsfile)[1].header
        if "GPS_TIME" in header:
            timeobsend = Time(header["GPS_TIME"])
        elif "DATE" in header:
            timeobsend = Time(header["DATE"])
        mjds[ii] = timeobsend.mjd
    fitsfiles = fitsfiles[np.argsort(mjds)]
    
    apthetas = np.linspace(0, 2*np.pi, 100)
    
    # loop over to plot
    cnt = 0
    for ii in range(nfiles):
    #for ii in range(1):
        fitsfile = fitsfiles[ii]
        fitsfileSplit = fitsfile.split("/")[-1].replace("_proc_regis","").replace(".fits","")
        print("%d/%d: %s"%(ii+1, nfiles, fitsfileSplit))
        hdulist = fits.open(fitsfile)
        nhdu = len(hdulist)
        for jj in range(nhdu):
        #for jj in range(20):
            if jj == 0: 
                continue
            if jj%50 == 0:
                print("  %d/%d"%(jj, nhdu))
            if (moviemode > 1) and (jj % moviemode !=1):
                continue
            header = hdulist[jj].header
            data = hdulist[jj].data
            xobj = header["X_OBJ"]
            yobj = header["Y_OBJ"]
            ystart = int(max(0, np.floor(yobj-halfwidth)))
            yend = int(min(data.shape[0], np.ceil(yobj+halfwidth)))
            xstart = int(max(0, np.floor(xobj-halfwidth)))
            xend = int(min(data.shape[0], np.ceil(xobj+halfwidth)))
            subdata = data[ystart:yend, xstart:xend]
            
            vmin = np.percentile(subdata,0.5)
            vmax = np.percentile(subdata,99.5)
            apxs = xobj +  np.cos(apthetas) * aper_size
            apys = yobj +  np.sin(apthetas) * aper_size
            skyapxs0 = xobj +  np.cos(apthetas) * sky_in
            skyapxs1 = xobj +  np.cos(apthetas) * sky_out
            skyapys0 = yobj +  np.sin(apthetas) * sky_in
            skyapys1 = yobj +  np.sin(apthetas) * sky_out
            
            plotName = os.path.join(movieDir,'image_%04d.png'%cnt)
            plt.figure(figsize=(5,5))
            plt.imshow(data, vmin=vmin,vmax=vmax, cmap='gray', origin = "lower")
            plt.xlim([xobj-halfwidth, xobj+halfwidth])
            plt.ylim([yobj-halfwidth, yobj+halfwidth])
            # plot aperture center
            plt.plot(xobj, yobj, 'rx', markersize=5)
            # plot the aperture
            plt.plot(apxs, apys, 'r', alpha = 0.5)
            # plot the sky region
            plt.plot(skyapxs0, skyapys0, 'b', alpha = 0.5)
            plt.plot(skyapxs1, skyapys1, 'b', alpha = 0.5)
            if jj==1:
                # write down the file name
                plt.text(xobj-halfwidth*0.9, yobj-halfwidth*0.9, fitsfileSplit[14:], color='r', fontsize= fs)
            plt.gca().set_axis_off()
            plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
                                hspace = 0, wspace = 0)
            plt.margins(0,0)
            plt.gca().xaxis.set_major_locator(NullLocator())
            plt.gca().yaxis.set_major_locator(NullLocator())
            plt.savefig(plotName, dpi=200)             
            plt.close()
            cnt +=1
            
    print ("Genrating the mpg file...")
    moviefiles = os.path.join(movieDir,"image_%04d.png")
    filename = os.path.join(movieDir,"movie.mpg")
    ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
    os.system(ffmpeg_command)

    rm_command = "rm %s/*.png" % (movieDir)
    os.system(rm_command)
    
    """
    print ("Genrating the gif file...")
    filename = os.path.join(movieDir,"movie.gif")
    ffmpeg_command = 'ffmpeg -an -y -r 20 -i %s -b:v %s %s'%(moviefiles,'5000k',filename)
    os.system(ffmpeg_command)
    """
    
