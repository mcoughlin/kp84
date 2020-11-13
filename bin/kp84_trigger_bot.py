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

def run_triggering(channel_id, outfile="", object_lists="", requests="",
                   bypass=False,
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
   
        doTrigger = False
        for mess in payload["messages"]:
            message_ts = float(mess["ts"])
            print(mess['text'])
            if np.abs(message_ts - thread_ts) > delay_thresh:
                continue
            txt = mess['text']
            txtsplit = list(filter(None,txt.replace("\xa0"," ").split(" ")))
            if len(txtsplit) == 0: continue
            if txtsplit[0] == "kped":
                print(txtsplit)
                doTrigger = True
                if len(txtsplit) == 2:
                    todo = txtsplit[1]
                elif len(txtsplit) == 3:
                    todo = txtsplit[1]
                    objName = txtsplit[2]
                elif len(txtsplit) == 6:
                    todo = txtsplit[1]
                    objName = txtsplit[2]
                    ra = float(txtsplit[3])
                    dec = float(txtsplit[4])
                    priority = float(txtsplit[5])
            user = mess['user']
        if not doTrigger:
            return
    else:
        user, message_ts = 'test', thread_ts
        todo, objName = 'trigger', 'ZTF20acozryr'
        ra, dec, priority = 42.1846, 12.1372, 10

    message = []
    message.append("Hi <@{0}>! You are interested in KPED triggering right? Let me get right on that for you.".format(user))
    message.append('Received request %.1f seconds ago...' % (np.abs(message_ts - thread_ts)))
    message.append("We are looking at analyzing %s for you" % objName)

    web_client.chat_postMessage(
        channel=channel_id,
        text="\n".join(message)
    )

    if todo == "trigger":
        message = []
        message.append("We are looking at triggering %s for you" % objName)

        filenames = glob.glob(os.path.join(object_lists,'transients','%s_*' % objName))
        if len(filenames) > 0:
            message.append('Sorry... %s already has a pending observation.' % objName)
        else:
            filename = os.path.join(object_lists,'transients','%s_%d.dat' % (objName, Time.now().gps))
            fid = open(filename, 'w')
            fid.write('%s %.5f %.5f %.5f\n' % (objName, ra, dec, priority))
            fid.close()
    
            message.append("%s triggered..." % objName)
        web_client.chat_postMessage(
            channel=channel_id,
            text="\n".join(message)
        )

    elif todo == "delete":
        message = []
        filenames = glob.glob(os.path.join(object_lists,'transients','%s_*' % objName))
        if len(filenames) > 0:
            for filename in filenames:
                rm_command = "rm %s" % (filename)
                os.system(rm_command) 
                message.append('Deleting %s...' % filename)
        else:               
            message.append('Sorry... %s has no pending observations.' % objName)
        web_client.chat_postMessage(
            channel=channel_id,
            text="\n".join(message)
        )

    # check current observation
    scheduler_command = "python kp84_scheduler"
    os.system(scheduler_command)

    lines = [line.rstrip('\n') for line in open(outfile)]
    message = []
    message.append('Current observational priority:')
    for line in lines:
        message.append(line)
    web_client.chat_postMessage(
        channel=channel_id,
        text="\n".join(message)
    )

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--channel", type=str, default="trigger")
    parser.add_argument("-d", "--debug", action="store_true", default=False)
    parser.add_argument("-np", "--noplots", action="store_true", default=False)
    parser.add_argument("-f","--outfile",default="/home/kped/Software/Queue/Michael_queue/queue_target.dat")
    parser.add_argument("-l","--object_lists",default="/home/kped/KP84/object_lists")
    parser.add_argument("-r","--requests",default="/home/kped/Software/Queue/Michael_queue/requests/")

    cfg = parser.parse_args()

    if cfg.channel == 'trigger':
        channel = 'C01EN3TM18A'
    else:
        print('Sorry, I do not know that channel...')
        exit(0)

    if cfg.debug:
        run_triggering(channel, bypass=True,
                     no_plots=cfg.noplots,
                     outfile=cfg.outfile,
                     object_lists=cfg.object_lists,
                     requests=cfg.requests)
        exit(0)

    while True:
        #try:
        print('Looking for some triggering to do!')
        run_triggering(channel, 
                       no_plots=cfg.noplots,
                       outfile=cfg.outfile,
                       object_lists=cfg.object_lists,
                       requests=cfg.requests)
        #except:
        #    pass
        time.sleep(15)
