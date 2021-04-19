"""
Database schema.
"""

from datetime import datetime, date
import simplejson as json
import enum
import os
import glob
import time
import copy
import configparser

from astropy import table
from astropy import coordinates
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.io import fits
import astropy.io.ascii as asci
import pkg_resources
import numpy as np
import pandas as pd

import sqlalchemy as sa
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy.orm import sessionmaker, scoped_session, relationship
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy import func
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.schema import DropConstraint
from arrow.arrow import Arrow

from kp84.config import app

from flask_login.mixins import UserMixin
from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy(app)

DBSession = scoped_session(sessionmaker())
EXECUTEMANY_PAGESIZE = 50000
utcnow = func.timezone('UTC', func.current_timestamp())

data_types = {
    int: 'int',
    float: 'float',
    bool: 'bool',
    dict: 'dict',
    str: 'str',
    list: 'list'
    }


class Encoder(json.JSONEncoder):
    """Extends json.JSONEncoder with additional capabilities/configurations."""
    def default(self, o):
        if isinstance(o, (datetime, Arrow, date)):
            return o.isoformat()

        elif isinstance(o, bytes):
            return o.decode('utf-8')

        elif hasattr(o, '__table__'):  # SQLAlchemy model
            return o.to_dict()

        elif o is int:
            return 'int'

        elif o is float:
            return 'float'

        elif type(o).__name__ == 'ndarray': # avoid numpy import
            return o.tolist()

        elif type(o).__name__ == 'DataFrame':  # avoid pandas import
            o.columns = o.columns.droplevel('channel')  # flatten MultiIndex
            return o.to_dict(orient='index')

        elif type(o) is type and o in data_types:
            return data_types[o]

        return json.JSONEncoder.default(self, o)

def to_json(obj):
    return json.dumps(obj, cls=Encoder, indent=2, ignore_nan=True)

class BaseMixin(object):
    query = DBSession.query_property()
    id = db.Column(db.Integer, primary_key=True)
    created_at = db.Column(db.DateTime, nullable=False, default=utcnow)
    modified = db.Column(db.DateTime, default=utcnow, onupdate=utcnow,
                         nullable=False)

    @declared_attr
    def __tablename__(cls):
        return cls.__name__.lower() + 's'

    __mapper_args__ = {'confirm_deleted_rows': False}

    def __str__(self):
        return to_json(self)

    def __repr__(self):
        attr_list = [f"{c.name}={getattr(self, c.name)}"
                     for c in self.__table__.columns]
        return f"<{type(self).__name__}({', '.join(attr_list)})>"

    def to_dict(self):
        if db.inspection.inspect(self).expired:
            DBSession().refresh(self)
        return {k: v for k, v in self.__dict__.items() if not k.startswith('_')}

    @classmethod
    def get_if_owned_by(cls, ident, user, options=[]):
        obj = cls.query.options(options).get(ident)

        if obj is not None and not obj.is_owned_by(user):
            raise AccessError('Insufficient permissions.')

        return obj

    def is_owned_by(self, user):
        raise NotImplementedError("Ownership logic is application-specific")

    @classmethod
    def create_or_get(cls, id):
        obj = cls.query.get(id)
        if obj is not None:
            return obj
        else:
            return cls(id=id)

Base = declarative_base(cls=BaseMixin)

# The db has to be initialized later; this is done by the app itself
# See `app_server.py`
def init_db(user, database, password=None, host=None, port=None):
    url = 'postgresql://{}:{}@{}:{}/{}'
    url = url.format(user, password or '', host or '', port or '', database)

    conn = sa.create_engine(url, client_encoding='utf8')
#                            executemany_mode='values',
#                            executemany_values_page_size=EXECUTEMANY_PAGESIZE)

    DBSession.configure(bind=conn)
    Base.metadata.bind = conn

    return conn


class Object(Base):
    """Object information"""

    objname = db.Column(
        db.String,
        nullable=False,
        unique=True,
        comment='Object Name')

    RA = db.Column(
        db.Float,
        nullable=False,
        comment='Right Ascension of the object')

    Dec = db.Column(
        db.Float,
        nullable=False,
        comment='Declination of the object')

    exposures = relationship(lambda: Exposure)

class Exposure(Base):
    """Exposure information"""

    objname = sa.Column(
        sa.String,
        sa.ForeignKey(Object.objname),
        nullable=False,
        comment='Object Name')

    date = db.Column(
        db.DateTime,
        nullable=False,
        unique=True,
        comment='UTC event timestamp',
        index=True)

    dateshort = db.Column(
        db.String,
        nullable=False,
        comment='Year-Month-Day time')

    image_id = db.Column(
        db.BigInteger,
        nullable=False,
        unique=True,
        comment='Image ID')

    cubes = relationship(lambda: Cube)
    
    reductions = relationship(lambda: Reduction)

class Cube(Base):
    """Cube information"""

    objname = sa.Column(
        sa.String,
        sa.ForeignKey(Object.objname),
        nullable=False,
        comment='Object Name')

    date = db.Column(
        db.DateTime,
        nullable=False,
        comment='UTC event timestamp',
        index=True)

    image_id = db.Column(
        sa.ForeignKey(Exposure.image_id),
        nullable=False,
        comment='Image ID')

    filename = db.Column(
        db.String,
        nullable=False,
        comment='Filename')

    exposure_time = db.Column(
        db.Float,
        nullable=False,
        comment='Exposure Time')


class Reduction(Base):
    """Reduction information."""

    objname = sa.Column(
        sa.String,
        sa.ForeignKey(Object.objname),
        nullable=False,
        comment='Object Name')

    image_id = db.Column(
        sa.ForeignKey(Exposure.image_id),
        nullable=False,
        unique=True,
        comment='Image ID')

    mjd = db.Column(
        db.ARRAY(db.Float),
        nullable=False,
        comment='MJD values')

    mag = db.Column(
        db.ARRAY(db.Float),
        nullable=False,
        comment='magnitude values')

    magerr = db.Column(
        db.ARRAY(db.Float),
        nullable=False,
        comment='magnitude error values')

    flux = db.Column(
        db.ARRAY(db.Float),
        nullable=False,
        comment='flux values')

    fluxerr = db.Column(
        db.ARRAY(db.Float),
        nullable=False,
        comment='flux error values')

def ingest_images(config, lookback, repeat=False):

    lookbackTD = TimeDelta(lookback,format='jd')

    folders = glob.glob(os.path.join(config["images"]["directory"],"20*"))
    for folder in folders:
        folderSplit = folder.split("/")
        yyyymmdd = folderSplit[-1]

        date = Time("%s-%s-%sT00:00:00"%(yyyymmdd[:4],
                                         yyyymmdd[4:6],
                                         yyyymmdd[6:8]),
                    format='isot', scale='utc')

        if Time.now() - date > lookbackTD: continue

        filenamesall = glob.glob(os.path.join(folder,"*.fits")) + glob.glob(os.path.join(folder,"*.fits.fz"))
        exposures = []
        for filename in filenamesall:
            filenameSplit = filename.split('/')
            filenameSplit = filenameSplit[-1].split('_')
            if filenameSplit[0] == "kped":
                exposure = filenameSplit[3]
            else:
                continue
            if exposure not in exposures:
                exposures.append(exposure)

        for exposure in exposures:
            fitsfiles = sorted(glob.glob('%s/*%s*.fit*'%(folder,exposure)))
            for filename in fitsfiles:
                hdul = fits.open(filename)
                filenameSplit = filename.split("/")[-1].split("_")
                imgid = int(filenameSplit[1] + filenameSplit[2])
                objname = filenameSplit[3]

                gpstime_start = Time(hdul[1].header['GPS_TIME'],
                                     format='isot', scale='utc')
                gpstime_end = Time(hdul[-1].header['GPS_TIME'],
                                     format='isot', scale='utc') 
                RA = hdul[0].header['RAD']
                Dec = hdul[0].header['DecD']
                exposure_time = (gpstime_end - gpstime_start).sec
                dateshort = str(gpstime_start.datetime)[:10]

                image_data = fits.getdata(filename, ext=1)

                exists = db.session.query(Object.objname).filter_by(objname=objname).first() is not None
                if not exists:
                    db.session().merge(Object(RA=RA,
                                              Dec=Dec,
                                              objname=objname))

                
                exists = Exposure.query.filter_by(image_id=imgid,
                                                  objname=objname).first() is not None
                if not exists:
                    print(imgid, objname)
                    db.session().merge(Exposure(objname=objname,
                                                image_id=imgid,
                                                date=gpstime_start.datetime,
                                                dateshort=dateshort))
                db.session().merge(Cube(filename=filename,
                                         image_id=imgid,
                                         objname=objname,
                                         date=gpstime_start.datetime,
                                         exposure_time=exposure_time))

                print('Ingested filename: %s' % filename)
                db.session().commit()


def ingest_reductions(config, lookback, repeat=False):

    lookbackTD = TimeDelta(lookback,format='jd')

    folders = glob.glob(os.path.join(config["reductions"]["directory"],"20*"))
    for folder in folders:
        folderSplit = folder.split("/")
        yyyymmdd = folderSplit[-1]

        date = Time("%s-%s-%sT00:00:00"%(yyyymmdd[:4],
                                         yyyymmdd[4:6],
                                         yyyymmdd[6:8]),
                    format='isot', scale='utc')

        if Time.now() - date > lookbackTD: continue

        objectsall = glob.glob(os.path.join(folder,"*-*"))
        for objdir in objectsall:
            productdir = os.path.join(objdir, 'product')
            lcurve = os.path.join(productdir, 'lightcurve.forced')
            if not os.path.isfile(lcurve): continue

            objdirSplit = objdir.split("/")
            imgid = int(objdirSplit[-2] + objdirSplit[-1].split("_")[0])
            objname = objdirSplit[-1].split("_")[1]

            reduction = Reduction.query.filter_by(image_id=imgid,
                                                  objname=objname).first()
            if reduction is not None:
                continue

            tblforced = asci.read(lcurve)
            mjd_forced = tblforced['MJD'].data
            mag_forced, magerr_forced = tblforced['mag'].data, tblforced['magerr'].data
            flux_forced, fluxerr_forced = tblforced['flux'].data, tblforced['fluxerr'].data

            db.session().merge(Reduction(objname=objname,
                                         image_id=imgid,
                                         mjd=mjd_forced,
                                         mag=mag_forced,
                                         magerr=magerr_forced,
                                         flux=flux_forced,
                                         fluxerr=fluxerr_forced))
            print('Ingested filename: %s' % lcurve)
            db.session().commit()            

def run_kped(init_db=False):

    if init_db:
        ingest_images(config, args.lookback, repeat=True)
        ingest_reductions(config, args.lookback, repeat=True)
    else:
        ingest_images(config, args.lookback, repeat=True)
        ingest_reductions(config, args.lookback)

    exposures = Exposure.query.all()
    for exposure in exposures:
        print(exposure.objname, exposure.date)

    reductions = Reduction.query.all()
    for reduction in reductions:
        print(reduction.objname, reduction.image_id)

if __name__ == "__main__":

    from argparse import ArgumentParser
    parser = ArgumentParser()

    parser.add_argument('-i', '--init_db', action='store_true', default=False)
    parser.add_argument('-C', '--config', default='input/config.yaml')
    parser.add_argument('-l', '--lookback', default=7, help='lookback in days')
    parser.add_argument("-d", "--debug", action="store_true", default=False)

    args = parser.parse_args()

    if not os.path.isfile(args.config):
        print('Missing config file: %s' % args.config)
        exit(1)

    config = configparser.ConfigParser()
    config.read(args.config)

    conn = init_db(config['database']['user'],
                   config['database']['database'],
                   password=config['database']['password'],
                   host=config['database']['host'],
                   port=config['database']['port'])

    if args.init_db:
        print(f'Creating tables on database {conn.url.database}')
        Base.metadata.reflect()
        Base.metadata.drop_all()
        Base.metadata.create_all()

        print('Refreshed tables:')
        for m in Base.metadata.tables:
            print(f' - {m}')

    if args.debug:
        run_kped(init_db=args.init_db)
        exit(0)

    while True:
        #try:
        print('Looking for some images to analyze!')
        run_kped(init_db=args.init_db)
        #except:
        #    pass
        time.sleep(15)
