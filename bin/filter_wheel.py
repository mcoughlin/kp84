#!/usr/bin/env python

import optparse
import numpy as np
import FLI

def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser()

    parser.add_option("-m","--mask",default=0,type=int)
    parser.add_option("-f","--filter",default=0,type=int)
    parser.add_option("--doPosition", action="store_true",default=False)
    parser.add_option("--doGetPosition", action="store_true",default=False)
    parser.add_option("--doGetFilterList", action="store_true",default=False)

    opts, args = parser.parse_args()

    return opts


def initialize_connection():
    fws = FLI.filter_wheel.USBFilterWheel.find_devices()

    for fw in fws:
        if fw.model == "CenterLine Filter Wheel":
            fw0 = fw
    return fw0


# Parse command line
opts = parse_commandline()

filters = ["clear","r","g","I","dark"]
masks = ["clear","U","B","V","R"]

if opts.doGetFilterList:
    cnt = 0
    print("Number, filter, mask")
    for filt, mask in zip(filters,masks):
        print("%d, %s, %s"%(cnt,filt,mask))
        cnt = cnt + 1
    exit(0)

fw0 = initialize_connection()

if opts.mask > 4 or opts.mask < 0:
    raise Exception("Mask position must be integer 0-4")
elif opts.filter > 4 or opts.filter < 0:
    raise Exception("Filter position must be integer 0-4")   

if opts.doGetPosition:
    pos = fw0.get_filter_pos()

    mask = pos/5
    filt = np.mod(pos,5)

    print("Mask: %d"%mask)
    print("Filter: %d"%filt)

if opts.doPosition:
    position = 5 * opts.mask + opts.filter
    fw0.set_filter_pos(position)
    pos = fw0.get_filter_pos()

    mask = pos/5
    filt = np.mod(pos,5)

    print("Mask: %d"%mask)
    print("Filter: %d"%filt)


