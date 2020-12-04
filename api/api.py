
import os
import glob
import time
import random
import string

import flask
from flask import request, jsonify, Response

from astropy.time import Time

import urllib
import requests
import jwt

app = flask.Flask(__name__)
app.config["DEBUG"] = True


secret_key_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".secret_key")
with open(secret_key_file, "r") as f:
    secret_key = f.read().rstrip()

object_lists="/home/kped/KP84/object_lists"
requests="/home/kped/Software/Queue/Michael_queue/requests"

@app.route('/api/queues', methods=['PUT','DELETE'])
def queues():

    if request.method == "PUT":
        payload = jwt.decode(request.data, secret_key, algorithms=['HS256'])
        cnt = 0
        for ii, target in enumerate(payload['targets']):
            request_id = target['id'] 
            objName = target['name']
            ra, dec = target['ra'], target['dec']
            priority = int(target['priority'])
            filt = target['filter']
            exp_time = target['exposure_time']
            exp_count = target['exposure_counts']
            program_pi = target['program_pi']

            for jj in range(exp_count):
                filename = os.path.join(object_lists,'fritz','%s-%s-%s-%d_%d.dat' % (objName, filt, request_id, cnt, Time.now().gps))
                fid = open(filename, 'w')
                fid.write('%s %.5f %.5f %s %.5f %.5f %s\n' % (objName, ra, dec, filt, exp_time, priority, program_pi))
                fid.close()
                cnt = cnt + 1

        res = Response('success',
                       status=201)
        return res
    elif request.method == "DELETE":
        payload = jwt.decode(request.data, secret_key, algorithms=['HS256'])

        cnt = 0
        for ii, target in enumerate(payload['targets']):
            request_id = target['id']
            objName = target['name']
            ra, dec = target['ra'], target['dec']
            priority = int(target['priority'])
            filt = target['filter']
            exp_time = target['exposure_time']
            exp_count = target['exposure_counts']
            program_pi = target['program_pi']

            for jj in range(exp_count):
                filenames = glob.glob(os.path.join(object_lists,'fritz','%s-%s-%s-%d_*.dat' % (objName, filt, request_id, cnt)))
                cnt = cnt + 1
                for filename in filenames:
                    rm_command = "rm %s" % filename
                    os.system(rm_command)
        res = Response('success',
                       status=201)

        return res

app.run()

