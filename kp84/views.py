import datetime
import json
import os
import urllib.parse
import math
import re
import requests
import shutil
import tempfile

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import time
import astropy.units as u
from astropy.table import Table
import pandas as pd
import matplotlib.style
import pkg_resources

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

    ims = models.db.session.query(models.Image).all()

    for im in ims:
        print(im.filename)
    
    return render_template(
        'index.html',
        ims=ims)
