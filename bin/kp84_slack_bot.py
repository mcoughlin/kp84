# Import the relevant packages
from socket import gethostname
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.time import Time
from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord
import astropy.units as u

from slack import RTMClient, WebClient
import numpy as np
import logging
import matplotlib.pyplot as plt
import io
import os
import sys
import glob
from astropy.time import Time
import traceback
import time

import kp84.visualize_utils

slack_token = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".slack_access_token.txt")
with open(slack_token, "r") as f:
    access_token = f.read()
web_client = WebClient(token=access_token)

def upload_fig(fig, user, filename, channel_id):
    imgdata = io.BytesIO()
    fig.savefig(imgdata, format='png', dpi=600, transparent=True)
    imgdata.seek(0)
    #wc = WebClient(token=bot_access_token)
    web_client.files_upload(
        file=imgdata.getvalue(),
        filename=filename,
        channels=channel_id,
        text="<@{0}>, here's the file {1} I've uploaded for you!".format(user, filename)
    )
    #fig.close()

def run_reductions(channel_id, setupDir="", outputDir="", bypass=False,
                  no_plots=False):

    thread_ts = time.time()

    #response = web_client.chat_postMessage(
    #    channel=channel_id,
    #    text='testing')

    if not bypass:
        try:
            converations = web_client.conversations_list(
                channel=channel_id
            )
            channel_slack_id = channel_id
            for converation in converations:
                for chan in converation.data["channels"]:
                    if chan["name"] == channel_id.replace("#",""):
                        channel_slack_id = chan["id"]
        
            delay_thresh = 60.0
        
            payload = web_client.conversations_history(
                channel=channel_slack_id,
                oldest=thread_ts-delay_thresh
            )
        except:
            return   
 
        if len(payload["messages"]) == 0:
            return
   
        doReduce, day = False, Time.now().isot.split("T")[0].replace("-","")
        objType, objName = 'variable', 'tmp'
        for mess in payload["messages"]:
            message_ts = float(mess["ts"])
            print(mess['text'])
            if np.abs(message_ts - thread_ts) > delay_thresh:
                continue
            txt = mess['text']
            txtsplit = list(filter(None,txt.split(" ")))
            if len(txtsplit) == 0: continue
            if txtsplit[0] == "kped":
                doReduce = True
                if len(txtsplit) == 2:
                    todo = txtsplit[1]
                elif len(txtsplit) == 3:
                    todo = txtsplit[1]
                    objName = txtsplit[2]
                elif len(txtsplit) == 4:
                    todo = txtsplit[1]
                    objName = txtsplit[2]
                    objType = txtsplit[3]
            user = mess['user']
        if not doReduce:
            return
    else:
        user, message_ts = 'test', thread_ts
        day = Time.now().isot.split("T")[0].replace("-","")
        #day = "20210224"
        #todo, objName = 'reduce', '063528_ZTF20acozryr-r'
        #objType = 'transient'
        day = "20210417"
        #todo, objName = 'movie', 'all'
        #todo, objName = 'stack', 'all'
        #todo, objName = 'setup_reduce', 'all'
        todo, objName = 'reduce', 'all'
        #todo, objName = 'reduce', '090531_ZTFJ15395027'
    

    message = []
    message.append("Hi <@{0}>! You are interested in KPED reductions right? Let me get right on that for you.".format(user))
    message.append('Received request %.1f seconds ago...' % (np.abs(message_ts - thread_ts)))
    message.append("We are looking at analyzing %s for you" % day)

    web_client.chat_postMessage(
        channel=channel_id,
        text="\n".join(message)
    )

    baseoutputDir = os.path.join(setupDir,day)
 
    if "setup" in todo:
        if objName == "redo":
            rm_command = "rm -rf %s" %baseoutputDir
            os.system(rm_command)

        if not os.path.isdir(baseoutputDir):
            message = []
            message.append("%s setup starting... please be patient." % baseoutputDir)
            web_client.chat_postMessage(
                channel=channel_id,
                text="\n".join(message)
            )           
             
            setup_command = "python kp84_setup_reduction --day %s" % day
            os.system(setup_command)
            message = []
            message.append("%s setup complete." % baseoutputDir)
        else:
            message = []
            message.append("%s already exists... please delete if you want it to be re-setup." % baseoutputDir)
        web_client.chat_postMessage(
            channel=channel_id,
            text="\n".join(message)
        )       

    directories = glob.glob(os.path.join(baseoutputDir,"*_*"))
    objs = []
    for directory in directories:
        obj = directory.split("/")[-1]
        objs.append(obj)

    message = []
    message.append(f"Found {len(objs)} candidates")
    message.append(str(objs))
    web_client.chat_postMessage(
        channel=channel_id,
        text="\n".join(message)
    )

    if "reduce" in todo:
        if not objName == "all":
            objsreduce = [objName]
        else:
            objsreduce = objs
                        

        for objNameTmp in objsreduce:
     
            if ("rgd" in objNameTmp) or ("-r" in objNameTmp) or ("-g" in objNameTmp):
                objType = "transient"
            else:
                objType = "variable"
            objType = "variable"

            baseoutputDir = os.path.join(outputDir, day, objNameTmp) # the output directory of this object
            outputProDir = os.path.join(baseoutputDir, "product")
  
            if not os.path.isdir(outputProDir):
                if objType == "variable":
                    setup_command = "python kp84_photometric_reduction --day %s --objName %s --doMakeMovie --doDynamicAperture" % (day, objNameTmp)
                elif objType == "transient":
                    setup_command = "python kp84_photometric_reduction --day %s --objName %s --doMakeMovie --doDynamicAperture --doTransient --doSubtraction" % (day, objNameTmp)
                else:
                    message = []
                    message.append("Sorry, %s transient type unknown. Please try variable or transient." % objType)
                    web_client.chat_postMessage(
                        channel=channel_id,
                        text="\n".join(message)
                    )
                    continue

                os.system(setup_command)
                message = []
                message.append("%s reduction complete." % objNameTmp)
            else:
                message = []
                message.append("%s already exists... please delete if you want it to be re-reduced." % baseoutputDir)
            web_client.chat_postMessage(
                channel=channel_id,
                text="\n".join(message)
            )
   
            finalforcefile = os.path.join(outputProDir,"lightcurve.forced")
            if not os.path.isfile(finalforcefile):
                message = []
                message.append("%s reduction failed... likely missed the target. Sorry!" % objNameTmp)        
                web_client.chat_postMessage(
                    channel=channel_id,
                    text="\n".join(message)
                )
                continue
   
            
            web_client.files_upload(
                file=finalforcefile,
                filename=finalforcefile.split("/")[-1],
                channels=channel_id,
                text="<@{0}>, here's the file {1} I've uploaded for you!".format(user, finalforcefile.split("/")[-1])
            )
    
            moviefile = os.path.join(outputProDir,"movie.mpg")
            if os.path.isfile(moviefile):
                web_client.files_upload(
                    file=moviefile,
                    filename=moviefile.split("/")[-1],
                    channels=channel_id,
                    text="<@{0}>, here's the file {1} I've uploaded for you!".format(user, moviefile.split("/")[-1])
                )
    
            tblforced = ascii.read(finalforcefile)
            mjd_forced = tblforced['MJD'].data
            mag_forced, magerr_forced = tblforced['mag'].data, tblforced['magerr'].data
            flux_forced, fluxerr_forced = tblforced['flux'].data, tblforced['fluxerr'].data
            idx = np.where(np.isfinite(mag_forced))[0]
            if not len(idx) == 0:    
                timetmp = (mjd_forced-mjd_forced[0])*24
                fig, ax1 = plt.subplots(1, 1, figsize=(9,6))
                ax1.errorbar(timetmp,mag_forced,magerr_forced,fmt='ko')
                ax1.set_xlabel('Time [hrs]')
                ax1.set_ylabel('Magnitude [ab]')
                idx = np.where(np.isfinite(mag_forced))[0]
                ymed = np.nanmedian(mag_forced)
                y10, y90 = np.nanpercentile(mag_forced[idx],10), np.nanpercentile(mag_forced[idx],90)
                ystd = np.nanmedian(magerr_forced[idx])
                ymin = y10 - 3*ystd
                ymax = y90 + 3*ystd
                ax1.set_ylim([ymin,ymax])
                ax1.set_xlim([min(timetmp), max(timetmp)])
                ax1.invert_yaxis()
                upload_fig(fig, user, "flux.png", channel_id)
                plt.close(fig)

            if objType == "transient":
                filename = '%s/triplet.png'%(outputProDir)
                web_client.files_upload(
                    file=filename,
                    filename=filename.split('/')[-1],
                    channels=channel_id,
                    text="<@{0}>, here's the file {1} I've uploaded for you!".format(user, filename.split('/')[-1])
                )

    if todo == "stack":
        if not objName == "all":
            objsreduce = [objName]
        else:
            objsreduce = objs

        for objName in objsreduce:
            
            dataDir = os.path.join(setupDir, day, objName)
            baseoutputDir = os.path.join(outputDir, day, objName) # the output directory of this object
            outputStackDir = os.path.join(baseoutputDir, "stack")

            fitsfiles = os.path.join(dataDir,'processing','*.fits')
            fitsstack = os.path.join(baseoutputDir,"stack.fits")
            pngstack = os.path.join(baseoutputDir,"stack.png")
            if not os.path.isdir(fitsstack) or not os.path.isdir(pngstack):
                setup_command = "python kp84_stack --inputfiles '%s' --outputfile %s --doOverwrite" % (fitsfiles, fitsstack)
                os.system(setup_command)
                message = []
                message.append("%s stack complete." % objName)
                web_client.chat_postMessage(
                    channel=channel_id,
                    text="\n".join(message)
                )

            web_client.files_upload(
                file=fitsstack,
                filename=fitsstack.split("/")[-1],
                channels=channel_id,
                text="<@{0}>, here's the file {1} I've uploaded for you!".format(user, fitsstack.split("/")[-1])
            )

            web_client.files_upload(
                file=pngstack,
                filename=pngstack.split("/")[-1],
                channels=channel_id,
                text="<@{0}>, here's the file {1} I've uploaded for you!".format(user, pngstack.split("/")[-1])
            )

    if todo == "movie":
        if not objName == "all":
            objsreduce = [objName]
        else:
            objsreduce = objs

        for objName in objsreduce:

            dataDir = os.path.join(setupDir, day, objName)
            baseoutputDir = os.path.join(outputDir, day, objName) # the output directory of this object
            outputStackDir = os.path.join(baseoutputDir, "stack")
            fitsfiles = os.path.join(dataDir,'processing','*.fits')

            moviefile = os.path.join(baseoutputDir,"movie.mpg")
            if not os.path.isdir(moviefile):
                setup_command = "python kp84_make_movie --inputfiles '%s' --outputfile %s" % (fitsfiles, moviefile)
                os.system(setup_command)
                message = []
                message.append("%s movie complete." % objName)
                web_client.chat_postMessage(
                    channel=channel_id,
                    text="\n".join(message)
                )

            web_client.files_upload(
                file=moviefile,
                filename=moviefile.split("/")[-1],
                channels=channel_id,
                text="<@{0}>, here's the file {1} I've uploaded for you!".format(user, moviefile.split("/")[-1])
            )

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--channel", type=str, default="reductions")
    parser.add_argument("-d", "--debug", action="store_true", default=False)
    parser.add_argument("-np", "--noplots", action="store_true", default=False)
    parser.add_argument("--setupDir", default = "/Backup/Data/archive_kped/data/reductions/")
    parser.add_argument("--outputDir", default = "/Backup/Data/archive_kped/data/photometry/")

    cfg = parser.parse_args()

    if cfg.channel == 'reductions':
        channel = 'C01EJMRKYKF'
    else:
        print('Sorry, I do not know that channel...')
        exit(0)

    if cfg.debug:
        run_reductions(channel, bypass=True,
                     no_plots=cfg.noplots,
                     setupDir=cfg.setupDir,
                     outputDir=cfg.outputDir)
        exit(0)

    while True:
        #try:
        print('Looking for some reducing to do!')
        run_reductions(channel, 
                           no_plots=cfg.noplots,
                           setupDir=cfg.setupDir,
                           outputDir=cfg.outputDir)
        #except:
        #    pass
        time.sleep(15)
