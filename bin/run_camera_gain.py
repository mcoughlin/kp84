#!/usr/bin/python

# Copyright (C) 2013 Michael Coughlin
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""KP84 camera script.

Script for running the KP84 camera.

Comments should be e-mailed to michael.coughlin@ligo.org.

"""

import os, sys, optparse, warnings
import time
import pyfits
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from andor2 import Andor

__author__ = "Michael Coughlin <michael.coughlin@ligo.org>"
__version__ = 1.0
__date__    = "9/22/2013"

# =============================================================================
#
#                               DEFINITIONS
#
# =============================================================================

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("-o","--outputDir",default="../output_gain")
    parser.add_option("-p","--plotDir",default="../plots_gain")
    parser.add_option("-m","--tempDir",default="../temp_gain")

    parser.add_option("--doPlots",  action="store_true", default=False)
    parser.add_option("--doTemperature",  action="store_true", default=False)
    parser.add_option("--doImages",  action="store_true", default=False)

    parser.add_option("-t", "--temperature", help="Temperature.", default=-10,type=int)
    parser.add_option("-N", "--nimages", help="Number of images.", default=60,type=int)
    parser.add_option("-d", "--duration", help="Duration.", default=1.0,type=float)

    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Run verbosely. (Default: False)")

    opts, args = parser.parse_args()

    # show parameters
    if opts.verbose:
        print >> sys.stderr, ""
        print >> sys.stderr, "running pylal_seismon_run..."
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

fitsdir = opts.outputDir
if not os.path.isdir(fitsdir):
    os.makedirs(fitsdir)
plotdir = opts.plotDir
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)
tempdir = opts.tempDir
if not os.path.isdir(tempdir):
    os.makedirs(tempdir)

T = opts.temperature
n_images = opts.nimages
exposure_time = opts.duration*1000.0

cam = Andor()
cam.Detector.OutputAmp(1)

if opts.doGain:
print(cam.EM.gain)
    print("Range of EM gains: %d-%d"%cam.EM.range)
    cam.EM.gain = opts.gain

print(cam.EM.mode)
print(cam.EM.is_on)
print(cam.EM.on)
print(cam.EM.is_on)

print(dir(cam.Temperature))
print(stop)

if opts.doTemperature:
    cam.Temperature.setpoint = T  # start cooling
    cam.Temperature.cooler = True

    status = cam.Temperature.read
    stabilized = status["status"]
    current_temp = int(status["temperature"])
    last_temp = current_temp
    print("Starting temp: %d"%current_temp)
    while not stabilized == "DRV_TEMP_STABILIZED":
        status = cam.Temperature.read
        stabilized = status["status"]
        current_temp = int(status["temperature"])
        print("Current temp: %d"%current_temp)
        time.sleep(10)
        last_temp = current_temp

if opts.doImages:
    cam.exposure = exposure_time

    filename = os.path.join(plotdir,'temp.dat')
    for ii in xrange(n_images):
        fitsfile = os.path.join(fitsdir,'fits_%03d.fits'%ii)
        tempfile = os.path.join(tempdir,'fits_%03d.dat'%ii)
        if os.path.isfile(fitsfile) and os.path.isfile(tempfile): continue

        status = cam.Temperature.read
        current_temp = int(status["temperature"])
        print("Image: %d, Current temp: %d"%(ii,current_temp))

        cam.Acquire.Single()
        data = cam.Acquire.snap(wait=True,type=32)
        hdulist = pyfits.PrimaryHDU(data)
        hdulist.writeto(fitsfile)

        fid = open(tempfile,'w')
        fid.write('%d\n'%current_temp)
        fid.close()
 
if opts.doPlots:
    adus, temps = [], []
    for ii in xrange(n_images):
        fitsfile = os.path.join(fitsdir,'fits_%03d.fits'%ii)
        plotfile = os.path.join(plotdir,'fits_%03d.png'%ii)
        tempfile = os.path.join(tempdir,'fits_%03d.dat'%ii)

        hdulist = pyfits.open(fitsfile)
        data = hdulist[0].data[-1,:,:]

        xshape, yshape = data.shape
        xlow, xhigh = np.min(xshape/2) - 100, np.min(xshape/2) + 100
        ylow, yhigh = np.min(yshape/2) - 100, np.min(yshape/2) + 100

        adu = np.median(data[xlow:xhigh,ylow:yhigh])
        adus.append(adu)
        temp = np.loadtxt(tempfile)
        temps.append(temp)

        if os.path.isfile(plotfile): continue

        vmin, vmax = np.percentile(data,10), np.percentile(data,90)
        #vmin, vmax = 300, 500
        plt.figure()
        plt.imshow(data,origin='lower',vmin=vmin,vmax=vmax)
        cbar = plt.colorbar()
        plt.show()
        plt.savefig(plotfile)
        plt.close()

    adus, temps = np.array(adus), np.array(temps)
    idx = np.where(np.arange(len(adus))>=len(adus)/2)[0]

    plotfile = os.path.join(plotdir,'adus.png')

    fig, ax1 = plt.subplots()
    ax1.plot(adus,'bx')
    ax1.set_xlabel('Exposure Number')
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('ADUs', color='b')
    ax1.tick_params('y', colors='b')

    ax2 = ax1.twinx()
    ax2.plot(temps, 'ro')
    ax2.set_ylabel('Temperature [C]', color='r')
    ax2.tick_params('y', colors='r')

    fig.tight_layout()
    plt.title('Relative RMS: %.2f'%(100*np.std(adus)/np.mean(adus)))
    plt.show()
    plt.savefig(plotfile)
    plt.close()

    plotfile = os.path.join(plotdir,'adus_zoom.png')

    fig, ax1 = plt.subplots()
    ax1.plot(adus,'bx')
    ax1.set_xlabel('Exposure Number')
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('ADUs', color='b')
    ax1.tick_params('y', colors='b')
    ax1.set_xlim([len(adus)/2,len(adus)])
    ax1.set_ylim([2750,2800])

    ax2 = ax1.twinx()
    ax2.plot(temps, 'ro')
    ax2.set_ylabel('Temperature [C]', color='r')
    ax2.tick_params('y', colors='r')
    ax2.set_xlim([len(adus)/2,len(adus)])
    ax2.set_ylim([-33,-27])

    fig.tight_layout()
    plt.title('Relative RMS: %.2f'%(100*np.std(adus[idx])/np.mean(adus[idx])))
    plt.show()
    plt.savefig(plotfile)
    plt.close()

