import datetime
import json
import os
import io
import urllib.parse
import math
import re
import requests
import shutil
import tempfile

import aplpy
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import time
import astropy.units as u
from astropy.table import Table
import pandas as pd
import matplotlib.style
import pkg_resources
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from flask import (
    abort, flash, jsonify, make_response, redirect, render_template, request,
    Response, url_for)
from flask_caching import Cache
from flask_login import (
    current_user, login_required, login_user, logout_user, LoginManager)
from wtforms import (
    BooleanField, FloatField, RadioField, TextField, IntegerField)
from wtforms_components.fields import (
    DateTimeField, DecimalSliderField, SelectField)
from wtforms import validators
from wtforms_alchemy.fields import PhoneNumberField
from passlib.apache import HtpasswdFile
from sqlalchemy.orm.exc import MultipleResultsFound, NoResultFound

from kp84.config import app
from kp84 import models

from PIL import Image
#
#
# From http://wtforms-alchemy.readthedocs.io/en/latest/advanced.html#using-wtforms-alchemy-with-flask-wtf  # noqa: E501
from flask_wtf import FlaskForm
from wtforms_alchemy import model_form_factory
# The variable db here is a SQLAlchemy object instance from
# Flask-SQLAlchemy package


BaseModelForm = model_form_factory(FlaskForm)



class ModelForm(BaseModelForm):
    @classmethod
    def get_session(cls):
        return models.db.session
#
#
#

# Server-side cache for rendered view functions.
cache = Cache(app, config={
    'CACHE_DEFAULT_TIMEOUT': 86400,
    'CACHE_REDIS_HOST': 'redis',
    'CACHE_TYPE': 'redis'})

def one_or_404(query):
    # FIXME: https://github.com/mitsuhiko/flask-sqlalchemy/pull/527
    rv = query.one_or_none()
    if rv is None:
        abort(404)
    else:
        return rv


def human_time(*args, **kwargs):
    secs = float(datetime.timedelta(*args, **kwargs).total_seconds())
    units = [("day", 86400), ("hour", 3600), ("minute", 60), ("second", 1)]
    parts = []
    for unit, mul in units:
        if secs / mul >= 1 or mul == 1:
            if mul > 1:
                n = int(math.floor(secs / mul))
                secs -= n * mul
            else:
                n = secs if secs != int(secs) else int(secs)
            parts.append("%s %s%s" % (n, unit, "" if n == 1 else "s"))
    return parts[0]


@app.route('/')
def index():

    objs = models.db.session.query(models.Object).all()
    objnames = []
    for obj in objs:
        objnames.append(obj.objname)
    objs = [x for _, x in sorted(zip(objnames, objs))]

    exposures = models.db.session.query(models.Exposure).all()

    days = []
    for exp in exposures:
        if exp.dateshort not in days:
            days.append(exp.dateshort)
    days = sorted(days)

    return render_template(
        'index.html',
        objs=objs,
        days=days)

@app.route('/obj/<objname>/')
def object(objname):

    query = models.db.session.query(models.Object).filter_by(objname=objname)

    try:
        obj = query.first()
    except NoResultFound:
        abort(404)
    exposures = models.db.session.query(models.Exposure).filter_by(objname=objname).all()

    return render_template(
            'obj.html',
            obj=obj,
            exposures=exposures)

@app.route('/obj/<objname>/image_id/<int:image_id>/exposure.png')
def exposure(objname, image_id):

    query = models.db.session.query(models.Cube).filter_by(objname=objname,
                                                           image_id=image_id)

    try:
        cubes = query.all()
    except NoResultFound:
        abort(404)

    fitsfiles = []
    for cube in cubes:
        fitsfiles.append(cube.filename)
    fitsfiles = sorted(fitsfiles)

    fig = Figure()
    f1 = aplpy.FITSFigure(fitsfiles[0],figure=fig)
    f1.show_grayscale(invert=False, stretch='power')
    fig.canvas.draw()

    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')


@app.route('/obj/<objname>/image_id/<int:image_id>/reduction.png')
def reduction(objname, image_id):

    query = models.db.session.query(models.Reduction).filter_by(objname=objname,
                                                                image_id=image_id)

    try:
        reduction = query.first()
    except NoResultFound:
        abort(404)
    mjd = reduction.mjd
    mag, magerr = reduction.mag, reduction.magerr
    flux, fluxerr = reduction.flux, reduction.fluxerr

    mjd = np.array(mjd)
    mag, magerr = np.array(mag), np.array(magerr)
    flux, fluxerr = np.array(flux), np.array(fluxerr)

    fig = Figure()
    ax = fig.add_subplot(1, 1, 1)
    timetmp = (mjd-mjd[0])*24
    ax.errorbar(timetmp,mag,magerr,fmt='ko')
    ax.set_xlabel('Time [hrs]')
    ax.set_ylabel('Magnitude [ab]')
    idx = np.where(np.isfinite(mag))[0]
    ymed = np.nanmedian(mag)
    y10, y90 = np.nanpercentile(mag[idx],10), np.nanpercentile(mag[idx],90)
    ystd = np.nanmedian(magerr[idx])
    ymin = y10 - 3*ystd
    ymax = y90 + 3*ystd
    ax.set_ylim([ymin,ymax])
    ax.set_xlim(min(timetmp), max(timetmp))
    ax.invert_yaxis()

    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

@app.route('/date/<day>/')
def date(day):
    
    query = models.db.session.query(models.Exposure).filter_by(dateshort=day)

    try:
        date = query.first()
    except NoResultFound:
        abort(404)

    exposures = models.db.session.query(models.Exposure).filter_by(dateshort=day).all()
    objects = []
    for exposure in exposures:
        objects.append(models.db.session.query(models.Object).filter_by(objname=exposure.objname).one())

    return render_template(
        'date.html',
        date=date,
        objects=objects,
        exposures=exposures)

@app.route('/calendar/')
def calendar():
    expall = models.db.session.query(models.Exposure)
    days=[]
    for exp in expall:
        if exp.dateshort not in days:
            days.append(exp.dateshort)
    days = sorted(days)
    days.reverse()
    objs = {}
    for day in days:
        objs[day] = []
        exposures = models.db.session.query(models.Exposure).filter_by(dateshort=day).all()
        for exp in exposures:
            objs[day].append(exp.objname)
        objs[day] = sorted(list(set(objs[day])))
    return render_template(
            'cal.html',
            days=days,
            objs=objs)
