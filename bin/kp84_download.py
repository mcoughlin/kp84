#!/usr/bin/env python
import os
import optparse
import subprocess
import numpy as np

def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser()

    parser.add_option("--downloadType", default="analysis")
    parser.add_option("--objName",default="14min")
    parser.add_option("--day",default="20190427")
    parser.add_option("-n",default=-1,type=int)

    opts, args = parser.parse_args()

    return opts

# Parse command line
opts = parse_commandline()
downloadType = opts.downloadType
objName = opts.objName
day = opts.day

outputDir = os.path.join("/Users/yuhanyao/Documents/GitHub/kp84/output", day, objName)
if not os.path.isdir(outputDir):
    os.makedirs(outputDir)

if downloadType == "data":
    ls_command = "ssh -p 22221 kped@140.252.53.120 ls /Data/%s/%s*.fits.fz" % (day, objName)
    result = subprocess.run(ls_command.split(" "), stdout=subprocess.PIPE)
    fitsfiles = list(filter(None,result.stdout.decode().split("\n")))

    nums = []
    for fitsfile in fitsfiles:
        fitsfileSplit = fitsfile.replace(".fits.fz","").replace(".fits","").split("_")
        try:
            num = int(fitsfileSplit[-1])
        except:
            num = -1
        nums.append(num)

    idx = np.where(np.array(nums) == opts.n)[0]
    for ii in idx:
        filename = fitsfiles[ii]

        scp_command = "scp -P 22221 kped@140.252.53.120:%s %s" % (filename, outputDir) 
        result = subprocess.run(scp_command.split(" "), stdout=subprocess.PIPE)
        

elif downloadType == "analysis":
    objpath = "/home/roboao/Michael/kp84/output/%s/%s/%s" % (day, objName, "product")

    scp_command = "scp -P 22220 -r roboao@140.252.53.120:%s/ %s" % (objpath, outputDir)
    result = subprocess.run(scp_command.split(" "), stdout=subprocess.PIPE)

else:
    print("downloadType must be data or analysis")
