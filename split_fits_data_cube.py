import astropy.io.fits as pyfits
import numpy as np
import argparse

'''
split_fits_data_cube.py

Splits a KPED 3D data cube into separate 2D fits files the time direction

Usage:

python split_fits_data_cube.py -if ZTF05yY_9_g_20190826_061720_sg_o.fits

results in there being n 2D frames:

ZTF05yY_9_g_20190826_061720_sg_o_1.fits
ZTF05yY_9_g_20190826_061720_sg_o_2.fits
ZTF05yY_9_g_20190826_061720_sg_o_3.fits
...
ZTF05yY_9_g_20190826_061720_sg_o_n.fits

must funpack the data first

'''

parser = argparse.ArgumentParser()
parser.add_argument("-if", "--infile", help="input fits data cube, ZTF05yY_9_g_20190826_061720_sg_o.fits")
args = parser.parse_args()

infile = str(args.infile)

file_suffix = '.fits'
underscore = '_'

hdulist = pyfits.open(infile)

file_name_base = infile[:infile.find('.fits')]#create generic name for HDU elements 

for i in range(1,99):#have to skip the first one
    new_name = file_name_base + underscore + str(i) + file_suffix
    try:
        new_header = hdulist[i].header
        new_data = hdulist[i].data
        pyfits.writeto(new_name,new_data.astype(np.float32),overwrite=True,header=new_header)
    except IndexError:
        print "\n End of HDU file, it only has " + str(i-1) + " elements."
        break
