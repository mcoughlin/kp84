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
            user = mess['user']
        if not doReduce:
            return
    else:
        user, message_ts = 'test', thread_ts
        day = Time.now().isot.split("T")[0].replace("-","")
        todo, objName = 'reduce', '051317_ZTF-J012747.63+525813.0'

    message = []
    message.append("Hi <@{0}>! You are interested in KPED reductions right? Let me get right on that for you.".format(user))
    message.append('Received request %.1f seconds ago...' % (np.abs(message_ts - thread_ts)))
    message.append("We are looking at analyzing %s for you" % day)

    web_client.chat_postMessage(
        channel=channel_id,
        text="\n".join(message)
    )

    baseoutputDir = os.path.join(setupDir,day)
    if todo == "setup":
        if not os.path.join(baseoutputDir):
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

    if todo == "reduce":
        outputDir = os.path.join(outputDir, day, objName) # the output directory of this object
        outputProDir = os.path.join(outputDir, "product")

        if not os.path.join(outputProDir):
            setup_command = "python kp84_photometric_reduction --day %s --objName %s --doMakeMovie --doDynamicAperture" % (day, objName)
            os.system(setup_command)
            message = []
            message.append("%s reduction complete." % objName)
        else:
            message = []
            message.append("%s already exists... please delete if you want it to be re-reduced." % outputDir)
        web_client.chat_postMessage(
            channel=channel_id,
            text="\n".join(message)
        )

        finalforcefile = os.path.join(outputProDir,"lightcurve.forced")
        if not os.path.isfile(finalforcefile):
            message = []
            message.append("%s reduction failed... likely missed the target. Sorry!" % objName)        
            web_client.chat_postMessage(
                channel=channel_id,
                text="\n".join(message)
            )
            return

        web_client.files_upload(
            file=finalforcefile,
            filename=finalforcefile.split("/")[-1],
            channels=channel_id,
            text="<@{0}>, here's the file {1} I've uploaded for you!".format(user, finalforcefile.split("/")[-1])
        )

        moviefile = os.path.join(outputProDir,"movie.mpg")
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
        ax1.set_xlim(min(timetmp), max(timetmp))
        ax1.invert_yaxis()
        upload_fig(fig, user, "flux.png", channel_id)
        plt.close(fig)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--channel", type=str, default="reductions")
    parser.add_argument("-d", "--debug", action="store_true", default=False)
    parser.add_argument("-np", "--noplots", action="store_true", default=False)
    parser.add_argument("--setupDir", default = "/Data3/archive_kped/data/reductions/")
    parser.add_argument("--outputDir", default = "../output")

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
        try:
            print('Looking for some reducing to do!')
            run_reductions(channel, 
                           no_plots=cfg.noplots,
                           setupDir=cfg.setupDir,
                           outputDir=cfg.outputDir)
        except:
            pass
        time.sleep(15)
