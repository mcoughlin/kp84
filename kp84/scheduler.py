#!/usr/bin/env python

import os, sys, optparse, shutil
import glob
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from datetime import datetime
from time import strptime

import astropy.table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
from astropy.coordinates import Angle
from astropy.table import Table, vstack
from astropy.time import Time, TimeDelta
from astropy.io import fits
from astropy.io import ascii

from astroplan import Observer
from astroplan import FixedTarget
from astroplan import ObservingBlock
from astroplan.constraints import TimeConstraint
from astroplan.constraints import AtNightConstraint, AirmassConstraint
from astroplan.scheduling import Transitioner
from astroplan.scheduling import SequentialScheduler
from astroplan.scheduling import PriorityScheduler
from astroplan.scheduling import Schedule
from astroplan.plots import plot_schedule_airmass
from astroplan import observability_table

def convert_to_hex(val, delimiter=':', force_sign=False):
    """
    Converts a numerical value into a hexidecimal string

    Parameters:
    ===========
    - val:           float
                     The decimal number to convert to hex.

    - delimiter:     string
                     The delimiter between hours, minutes, and seconds
                     in the output hex string.

    - force_sign:    boolean
                     Include the sign of the string on the output,
                     even if positive? Usually, you will set this to
                     False for RA values and True for DEC

    Returns:
    ========
    A hexadecimal representation of the input value.
    """
    s = np.sign(val)
    s_factor = 1 if s > 0 else -1
    val = np.abs(val)
    degree = int(val)
    minute = int((val  - degree)*60)
    second = (val - degree - minute/60.0)*3600.
    if degree == 0 and s_factor < 0:
        return '-00{2:s}{0:02d}{2:s}{1:.2f}'.format(minute, second, delimiter)
    elif force_sign or s_factor < 0:
        deg_str = '{:+03d}'.format(degree * s_factor)
    else:
        deg_str = '{:02d}'.format(degree * s_factor)
    return '{0:s}{3:s}{1:02d}{3:s}{2:.2f}'.format(deg_str, minute, second, delimiter)

def load_ilaria_objects(filename):

    columns = ['name', 'ra', 'dec', 'duration']

    observations = astropy.io.ascii.read(filename,format='csv',
                                         data_start=0,
                                         names = columns,
                                         delimiter = " ")

    names = ["requestID", "programID", "objectID", "ra_hex", "dec_hex", "epoch", "ra_rate", "dec_rate", "mag", "exposure_time", "filter", "mode", "pi", "comment"]
    targets = []
    cnt = 0

    filt = 1
    program_id = 8
    program_pi = "caiazzo"

    for ii, observation in enumerate(observations.filled()):
        name = observation["name"].replace(" ","-")
        ra, dec = observation["ra"], observation["dec"]
        duration = observation["duration"]

        if duration > 7200:
            exposure_time = 7200
        else:
            exposure_time = duration

        ra_hex, dec_hex = convert_to_hex(ra*24/360.0,delimiter=':'), convert_to_hex(dec,delimiter=':')

        if filt == 1:
            filt = "FILTER_SLOAN_G"
        elif filt == 2:
            filt = "FILTER_SLOAN_R"
        elif filt == "3":
            filt = "FILTER_JOHNSON_I"

        comment = "%.0f_%.5f" % (10, exposure_time)

        targets.append([name,
                        program_id,
                        name,
                        ra_hex, dec_hex, 2000, 0, 0,
                        0, exposure_time,
                        filt, 13, program_pi,
                        comment])

    targets = Table(rows=targets, names=names)

    sigs, periods = [], []
    coords, target = [], []
    ras, decs = [], []
    for row in targets:
        comment = row["comment"]
        commentSplit = comment.split("_")
        sig, period = float(commentSplit[0]), float(commentSplit[1])
        sigs.append(sig)
        periods.append(period)

        ra_hex, dec_hex = row["ra_hex"], row["dec_hex"]

        ra  = Angle(ra_hex, unit=u.hour).deg
        dec = Angle(dec_hex, unit=u.deg).deg

        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        tar = FixedTarget(coord=coord, name=row["objectID"])
        coords.append(coord)
        target.append(tar)
        ras.append(ra)
        decs.append(dec)

    targets["sig"] = sigs
    targets["periods"] = periods
    targets["coords"] = coords
    targets["target"] = target
    targets["ra"] = ras
    targets["dec"] = decs
    targets["programID"] = 1
    targets["priority"] = sigs
    targets["redo"] = 1
    targets["delta_redo"] = 100

    targets.sort("sig")
    targets.reverse()
    targets = astropy.table.unique(targets, keys=["objectID"])

    targets['ra_rate'].dtype= np.float64
    targets['dec_rate'].dtype= np.float64

    print(targets)

    return targets

def load_kevin_objects(filename):

    columns = ['type', 'name', 'ra', 'dec',
               'period_days', 'period_minutes', 'new']

    observations = astropy.io.ascii.read(filename,format='csv',
                                         data_start=0,
                                         names = columns)
    observations['new'].fill_value = 0

    names = ["requestID", "programID", "objectID", "ra_hex", "dec_hex", "epoch", "ra_rate", "dec_rate", "mag", "exposure_time", "filter", "mode", "pi", "comment"]
    targets = []
    cnt = 0

    filt = 1
    program_id = 8
    program_pi = "burdge"

    for ii, observation in enumerate(observations.filled()):
        name = observation["name"].replace(" ","-")
        ra, dec = observation["ra"], observation["dec"]
        period = observation["period_days"]
        new = observation["new"]

        exposure_time = int(1.1*period*86400)

        if exposure_time < 1800:
            exposure_time = 1800

        ra_hex, dec_hex = convert_to_hex(ra*24/360.0,delimiter=':'), convert_to_hex(dec,delimiter=':')

        if filt == 1:
            filt = "FILTER_SLOAN_G"
        elif filt == 2:
            filt = "FILTER_SLOAN_R"
        elif filt == "3":
            filt = "FILTER_JOHNSON_I"

        if new == '0':
            comment = "%.0f_%.5f" % (0, period)
        else:
            comment = "%.0f_%.5f" % (1e23, period)

        targets.append([name,
                        program_id,
                        name,
                        ra_hex, dec_hex, 2000, 0, 0,
                        0, exposure_time,
                        filt, 13, program_pi,
                        comment])

    targets = Table(rows=targets, names=names)

    sigs, periods = [], []
    coords, target = [], []
    ras, decs = [], []
    for row in targets:
        comment = row["comment"]
        commentSplit = comment.split("_")
        sig, period = float(commentSplit[0]), float(commentSplit[1])
        sigs.append(sig)
        periods.append(period)

        ra_hex, dec_hex = row["ra_hex"], row["dec_hex"]

        ra  = Angle(ra_hex, unit=u.hour).deg
        dec = Angle(dec_hex, unit=u.deg).deg

        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        tar = FixedTarget(coord=coord, name=row["objectID"])
        coords.append(coord)
        target.append(tar)
        ras.append(ra)
        decs.append(dec)

    targets["sig"] = sigs
    targets["periods"] = periods
    targets["coords"] = coords
    targets["target"] = target
    targets["ra"] = ras
    targets["dec"] = decs
    targets["programID"] = 1
    targets["priority"] = sigs
    targets["redo"] = 1
    targets["delta_redo"] = 0.5

    targets.sort("sig")
    targets.reverse()
    targets = astropy.table.unique(targets, keys=["objectID"])

    targets['ra_rate'].dtype= np.float64
    targets['dec_rate'].dtype= np.float64

    return targets

def load_jan_objects(filename):

    observations = astropy.io.ascii.read(filename,format='csv')

    names = ["requestID", "programID", "objectID", "ra_hex", "dec_hex", "epoch", "ra_rate", "dec_rate", "mag", "exposure_time", "filter", "mode", "pi", "comment"]
    targets = []
    cnt = 0

    filt = 1
    program_id = 7
    program_pi = "vanRoestel"

    for ii, observation in enumerate(observations):

        ra, dec = observation["Ra"], observation["Dec"]
        period = observation["Period (d)"]
        exposure_time = int(1.1*period*86400)
        Gmag = observation["G"]

        # only want the bright stuff
        if Gmag > 18: continue

        ra_hex, dec_hex = convert_to_hex(ra*24/360.0,delimiter=':'), convert_to_hex(dec,delimiter=':')

        if filt == 1:
            filt = "FILTER_SLOAN_G"
        elif filt == 2:
            filt = "FILTER_SLOAN_R"
        elif filt == "3":
            filt = "FILTER_JOHNSON_I"

        comment = "%.0f_0" % (1e15 - Gmag)
        targets.append([observation['Shortname'],
                        program_id,
                        observation['Shortname'],
                        ra_hex, dec_hex, 2000, 0, 0,
                        1e5, exposure_time,
                        filt, 9, program_pi,
                        comment])

    targets = Table(rows=targets, names=names)

    sigs, periods = [], []
    coords, target = [], []
    ras, decs = [], []
    for row in targets:
        comment = row["comment"]
        commentSplit = comment.split("_")
        sig, period = float(commentSplit[0]), float(commentSplit[1])
        sigs.append(sig)
        periods.append(period)

        ra_hex, dec_hex = row["ra_hex"], row["dec_hex"]

        ra  = Angle(ra_hex, unit=u.hour).deg
        dec = Angle(dec_hex, unit=u.deg).deg

        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        tar = FixedTarget(coord=coord, name=row["objectID"])
        coords.append(coord)
        target.append(tar)
        ras.append(ra)
        decs.append(dec)

    targets["sig"] = sigs
    targets["periods"] = periods
    targets["coords"] = coords
    targets["target"] = target
    targets["ra"] = ras
    targets["dec"] = decs
    targets["programID"] = 1
    targets["priority"] = sigs

    targets.sort("sig")
    targets.reverse()
    targets = astropy.table.unique(targets, keys=["objectID"])

    targets['ra_rate'].dtype= np.float64
    targets['dec_rate'].dtype= np.float64

    return targets

def load_observations(filename):

    lines = [line.rstrip('\n') for line in open(filename)]
    data = {}
    for line in lines:
        lineSplit = line.split("=")
        data[lineSplit[0]] = lineSplit[1]

    return data

def load_variables(filename):

    names = ["objectID", "ra", "dec", "exposure_time", "filter", "priority"]
    observations = astropy.io.ascii.read(filename,names=names)

    names = ["requestID", "programID", "objectID", "ra_hex", "dec_hex", "epoch", "ra_rate", "dec_rate", "mag", "exposure_time", "filter", "mode", "pi", "comment"]
    targets = []
    cnt = 0
    for ii, observation in enumerate(observations):

        ra, dec = observation["ra"], observation["dec"]
        ra_hex, dec_hex = convert_to_hex(ra*24/360.0,delimiter=':'), convert_to_hex(dec,delimiter=':')

        if observation["filter"] == "g":
            filts = ["FILTER_SLOAN_G"] 
        elif observation["filter"] == "r":
            filts = ["FILTER_SLOAN_R"]
        elif observation["filter"] == "V":
            filts = ["FILTER_JOHNSON_V"]
        for filt in filts:
            comment = "%d_0" % (observation["priority"])
            targets.append(["%s-%s"%(observation["objectID"],filt[-1].lower()),
                            3,
                            "%s-%s"%(observation["objectID"],filt[-1].lower()),
                            ra_hex, dec_hex, 2000, 0, 0,
                            1e5, observation["exposure_time"],
                            filt, 9, "Coughlin",
                            comment])
            cnt = cnt + 1

    targets = Table(rows=targets, names=names)

    sigs, periods = [], []
    coords, target = [], []
    ras, decs = [], []
    for row in targets:
        comment = row["comment"]
        commentSplit = comment.split("_")
        sig, period = float(commentSplit[0]), float(commentSplit[1])
        sigs.append(sig)
        periods.append(period)

        ra_hex, dec_hex = row["ra_hex"], row["dec_hex"]

        ra  = Angle(ra_hex, unit=u.hour).deg
        dec = Angle(dec_hex, unit=u.deg).deg

        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        tar = FixedTarget(coord=coord, name=row["objectID"])
        coords.append(coord)
        target.append(tar)
        ras.append(ra)
        decs.append(dec)

    targets["sig"] = sigs
    targets["periods"] = periods
    targets["coords"] = coords
    targets["target"] = target
    targets["ra"] = ras
    targets["dec"] = decs
    targets["programID"] = 1
    targets["priority"] = sigs

    targets.sort("sig")
    targets.reverse()
    targets = astropy.table.unique(targets, keys=["objectID"])

    targets['ra_rate'].dtype= np.float64
    targets['dec_rate'].dtype= np.float64

    return targets

def load_transients(filename):

    names = ["objectID", "ra", "dec", "priority"]
    observations = astropy.io.ascii.read(filename,names=names)

    names = ["requestID", "programID", "objectID", "ra_hex", "dec_hex", "epoch", "ra_rate", "dec_rate", "mag", "exposure_time", "filter", "mode", "pi", "comment"]
    targets = []
    cnt = 0
    for ii, observation in enumerate(observations):

        ra, dec = observation["ra"], observation["dec"]
        ra_hex, dec_hex = convert_to_hex(ra*24/360.0,delimiter=':'), convert_to_hex(dec,delimiter=':')

        filts = ["FILTER_SLOAN_G", "FILTER_SLOAN_R"]
        for filt in filts:
            comment = "10000000000000000000_0"
            targets.append(["%s-%s"%(observation["objectID"],filt[-1].lower()),
                            3,
                            "%s-%s"%(observation["objectID"],filt[-1].lower()),
                            ra_hex, dec_hex, 2000, 0, 0,
                            1e5, 180.0,
                            filt, 9, "Coughlin",
                            comment])
            cnt = cnt + 1

    targets = Table(rows=targets, names=names)

    sigs, periods = [], []
    coords, target = [], []
    ras, decs = [], []
    for row in targets:
        comment = row["comment"]
        commentSplit = comment.split("_")
        sig, period = float(commentSplit[0]), float(commentSplit[1])
        sigs.append(sig)
        periods.append(period)

        ra_hex, dec_hex = row["ra_hex"], row["dec_hex"]

        ra  = Angle(ra_hex, unit=u.hour).deg
        dec = Angle(dec_hex, unit=u.deg).deg

        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        tar = FixedTarget(coord=coord, name=row["objectID"])
        coords.append(coord)
        target.append(tar)
        ras.append(ra)
        decs.append(dec)

    targets["sig"] = sigs
    targets["periods"] = periods
    targets["coords"] = coords
    targets["target"] = target
    targets["ra"] = ras
    targets["dec"] = decs
    targets["programID"] = 1
    targets["priority"] = sigs

    targets.sort("sig")
    targets.reverse()
    targets = astropy.table.unique(targets, keys=["objectID"])

    targets['ra_rate'].dtype= np.float64
    targets['dec_rate'].dtype= np.float64

    return targets

def load_NEO(filename):

    filenameSplit = filename.split("_")
    names = ["requestID", "programID", "objectID", "ra_hex", "dec_hex", "epoch", "ra_rate", "dec_rate", "mag", "exposure_time", "filter", "mode", "pi", "comment"]
    targets = astropy.io.ascii.read(filename,names=names,format='csv')
    targets['priority'] = 1e15
    targets['ra_rate'].dtype= np.float64
    targets['dec_rate'].dtype= np.float64
    targets = targets.filled()

    coords, target = [], []
    ras, decs = [], []
    for row in targets:
        ra_hex, dec_hex = row["ra_hex"], row["dec_hex"]

        ra  = Angle(ra_hex, unit=u.hour).deg
        dec = Angle(dec_hex, unit=u.deg).deg

        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        tar = FixedTarget(coord=coord, name=row["objectID"])
        coords.append(coord)
        target.append(tar)
        ras.append(ra)
        decs.append(dec)

    targets["coords"] = coords
    targets["target"] = target
    targets["ra"] = ras
    targets["dec"] = decs


    mjd_now = Time.now().mjd
    mjds = []
    for row in targets:
        if not row["comment"]: continue
        tt = datetime(*strptime(row["comment"], '%Y-%b-%d %H:%M:%S.%f')[:6])
        mjds.append(Time(tt, format='datetime').mjd)
    mjds = np.array(mjds)
    idx = np.argmin(np.abs(mjd_now-mjds))
    targets = targets[idx]

    return targets

def load_time_sensitive(filename,buffertime=15./60/24):

    filenameSplit = filename.split("_")
    priority = float(filenameSplit[-1].replace(".dat",""))

    names = ["requestID", "programID", "objectID", "ra_hex", "dec_hex", "epoch", "ra_rate", "dec_rate", "mag", "exposure_time", "filter", "mode", "pi", "comment"]
    targets = astropy.io.ascii.read(filename,names=names)

    mjd_now = Time.now().mjd

    sigs, periods = [], []
    coords, target = [], []
    ras, decs = [], []
    for row in targets:
        comment = row["comment"]
        commentSplit = comment.split("_")
        mjd_start, mjd_end = float(commentSplit[0]), float(commentSplit[1])

        sig, period = priority, 0
 
        dt = mjd_start - mjd_now
        #dt_hrs = dt*24.0
        #if (dt_hrs < 0) or (dt_hrs > 300.0):
        if dt < 0 or dt > buffertime:        
            continue

        sigs.append(sig)
        periods.append(period)

        ra_hex, dec_hex = row["ra_hex"], row["dec_hex"]

        ra  = Angle(ra_hex, unit=u.hour).deg
        dec = Angle(dec_hex, unit=u.deg).deg

        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        tar = FixedTarget(coord=coord, name=row["objectID"])
        coords.append(coord)
        target.append(tar)
        ras.append(ra)
        decs.append(dec)

    if len(sigs) == 0:
        return []

    targets["sig"] = sigs
    targets["periods"] = periods
    targets["coords"] = coords
    targets["target"] = target
    targets["ra"] = ras
    targets["dec"] = decs
    targets["programID"] = 1
    targets["priority"] = np.array(sigs) + priority

    targets.sort("sig")
    targets.reverse()
    targets = astropy.table.unique(targets, keys=["objectID"])

    targets['ra_rate'].dtype= np.float64
    targets['dec_rate'].dtype= np.float64

    return targets

def load_periodic_time_sensitive(filename,buffertime=15./60/24):

    filenameSplit = filename.split("_")
    priority = float(filenameSplit[-1].replace(".dat",""))

    names = ["requestID", "programID", "objectID", "ra_hex", "dec_hex", "epoch", "ra_rate", "dec_rate", "mag", "exposure_time", "filter", "mode", "pi", "comment"]
    targets = astropy.io.ascii.read(filename,names=names)

    jd_now = Time.now().jd

    sigs, periods = [], []
    coords, target = [], []
    ras, decs = [], []
    exposures = []
    for row in targets:
        comment = row["comment"]
        commentSplit = comment.split("_")
        # units are in days, format is HJD
        p,t0,window = float(commentSplit[0]), float(commentSplit[1]), float(commentSplit[2])

        sig, period = priority, 0
 
        # convert current time to heliocentre
        ut_helio = JD2HJD(jd_now,ra,dec) 

        # calculate the number of eclipses since t0
        N = (ut_helio-t0).value//p

        tn = t0+p*(N+1) - ut_helio # time until next window midpoint
        dt = tn - 0.5*window # time to start of next window
        if dt > buffertime or dt < 0:
            # if next window is too far away or already ongoing, skip
            continue

        sigs.append(sig)
        periods.append(period)

        ra_hex, dec_hex = row["ra_hex"], row["dec_hex"]

        ra  = Angle(ra_hex, unit=u.hour).deg
        dec = Angle(dec_hex, unit=u.deg).deg

        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        tar = FixedTarget(coord=coord, name=row["objectID"])
        coords.append(coord)
        target.append(tar)
        ras.append(ra)
        decs.append(dec)
        exposures.append((window+dt)*3600*24) # set the exposuretime to windowsize + time untill window start

    if len(sigs) == 0:
        return []

    targets["sig"] = sigs
    targets["periods"] = periods
    targets["coords"] = coords
    targets["target"] = target
    targets["ra"] = ras
    targets["dec"] = decs
    targets["exposure_time"] = exposures
    targets["programID"] = 1
    targets["priority"] = np.array(sigs) + priority

    targets.sort("sig")
    targets.reverse()
    targets = astropy.table.unique(targets, keys=["objectID"])

    targets['ra_rate'].dtype= np.float64
    targets['dec_rate'].dtype= np.float64

    return targets

def JD2HJD(JD,ra,dec):
    # JD2HJD for Kittpeak!
    location = EarthLocation.from_geodetic(-111.5967*u.deg, 31.9583*u.deg,
                                           2096*u.m)
    kp = Observer(location=location, name="Kitt Peak",timezone="US/Arizona")


    # convert JD to HJD
    target = coord.SkyCoord(ra*u.deg,dec*u.deg, frame='icrs')
    #tsite = coord.EarthLocation.of_site(site)
    times = Time(JD, format='jd',
                      scale='utc', location=location)
    ltt_helio = times.light_travel_time(target, 'heliocentric')

    HJD = JD+ltt_helio

    return HJD

def load_backup(filename):

    filenameSplit = filename.split("_")
    priority = float(filenameSplit[-1].replace(".dat",""))
    targets = astropy.io.ascii.read(filename)
    ncolumns = len(targets.columns)

    if ncolumns == 16:
        names = ["requestID", "programID", "objectID", "ra_hex", "dec_hex", "epoch", "ra_rate", "dec_rate", "mag", "exposure_time", "filter", "mode", "pi", "comment","redo","delta_redo"]
    elif ncolumns == 15:
        names = ["requestID", "programID", "objectID", "ra_hex", "dec_hex", "epoch", "ra_rate", "dec_rate", "mag", "exposure_time", "filter", "mode", "pi", "comment","redo"]
    elif ncolumns == 14:
        names = ["requestID", "programID", "objectID", "ra_hex", "dec_hex", "epoch", "ra_rate", "dec_rate", "mag", "exposure_time", "filter", "mode", "pi", "comment"]
    targets = astropy.io.ascii.read(filename,names=names)

    sigs, periods = [], []
    coords, target = [], []
    ras, decs = [], []
    for row in targets:
        comment = row["comment"]
        commentSplit = comment.split("_")
        sig, period = float(commentSplit[0]), float(commentSplit[1])
        sigs.append(sig)
        periods.append(period)
    
        ra_hex, dec_hex = row["ra_hex"], row["dec_hex"]
    
        ra  = Angle(ra_hex, unit=u.hour).deg
        dec = Angle(dec_hex, unit=u.deg).deg
    
        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        tar = FixedTarget(coord=coord, name=row["objectID"])
        coords.append(coord)
        target.append(tar)
        ras.append(ra)
        decs.append(dec)
    
    targets["sig"] = sigs
    targets["periods"] = periods
    targets["coords"] = coords
    targets["target"] = target
    targets["ra"] = ras
    targets["dec"] = decs
    targets["programID"] = 1
    targets["priority"] = np.array(sigs) + priority    
    targets["redo"] = 1

    targets.sort("sig")
    targets.reverse()
    targets = astropy.table.unique(targets, keys=["objectID"])

    targets['ra_rate'].dtype= np.float64
    targets['dec_rate'].dtype= np.float64

    return targets

def load_rgd_objects(filename):

    names = ["requestID", "programID", "objectID", "ra_hex", "dec_hex", "ra_rate", "dec_rate", "exposure_time", "filter", "mode", "pi", "comment"]
    targets = astropy.io.ascii.read(filename,format='no_header', delimiter=',',names=names)

    sigs, periods = [], []
    coords, target = [], []
    ras, decs = [], []
    for row in targets:
        sigs.append(0)
        periods.append(0)

        ra_hex, dec_hex = row["ra_hex"], row["dec_hex"]

        ra  = Angle(ra_hex, unit=u.hour).deg
        dec = Angle(dec_hex, unit=u.deg).deg

        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        tar = FixedTarget(coord=coord, name=row["objectID"])
        coords.append(coord)
        target.append(tar)
        ras.append(ra)
        decs.append(dec)

    targets["sig"] = sigs
    targets["periods"] = periods
    targets["coords"] = coords
    targets["target"] = target
    targets["ra"] = ras
    targets["dec"] = decs
    targets["programID"] = 1
    targets["epoch"] = 2000
    targets["priority"] = np.array(sigs)

    targets.sort("sig")
    targets.reverse()
    targets = astropy.table.unique(targets, keys=["objectID"])

    targets['ra_rate'].dtype= np.float64
    targets['dec_rate'].dtype= np.float64
    #targets['requestID'].dtype = np.float64

    return targets

def load_GW(filename):
    lines = [line.rstrip('\n') for line in open(filename)]
    data = eval(lines[0])

    queue_name = data['queue_name']
    observations = data['targets']

    names = ["requestID", "programID", "objectID", "ra_hex", "dec_hex", "epoch", "ra_rate", "dec_rate", "mag", "exposure_time", "filter", "mode", "pi", "comment"]
    targets = []
    for ii, observation in enumerate(observations):

        ra, dec = observation["ra"], observation["dec"]
        ra_hex, dec_hex = convert_to_hex(ra*24/360.0,delimiter=':'), convert_to_hex(dec,delimiter=':')

        filt = observation["filter_id"]
        if filt == 1:
            filt = "FILTER_SLOAN_G"
        elif filt == 2:
            filt = "FILTER_SLOAN_R"
        elif filt == "3":
            filt = "FILTER_JOHNSON_I"

        comment = "10000000000000000000_0"
        targets.append([str(observation['request_id']),
                        observation['program_id'],
                        str(observation['request_id']),
                        ra_hex, dec_hex, 2000, 0, 0,
                        1e5, observation['exposure_time'],
                        filt, 9, observation['program_pi'],
                        comment])

    targets = Table(rows=targets, names=names)

    sigs, periods = [], []
    coords, target = [], []
    ras, decs = [], []
    for row in targets:
        comment = row["comment"]
        commentSplit = comment.split("_")
        sig, period = float(commentSplit[0]), float(commentSplit[1])
        sigs.append(sig)
        periods.append(period)

        ra_hex, dec_hex = row["ra_hex"], row["dec_hex"]

        ra  = Angle(ra_hex, unit=u.hour).deg
        dec = Angle(dec_hex, unit=u.deg).deg

        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        tar = FixedTarget(coord=coord, name=row["objectID"])
        coords.append(coord)
        target.append(tar)
        ras.append(ra)
        decs.append(dec)

    targets["sig"] = sigs
    targets["periods"] = periods
    targets["coords"] = coords
    targets["target"] = target
    targets["ra"] = ras
    targets["dec"] = decs
    targets["programID"] = 1
    targets["priority"] = sigs

    targets.sort("sig")
    targets.reverse()
    targets = astropy.table.unique(targets, keys=["objectID"])

    targets['ra_rate'].dtype= np.float64
    targets['dec_rate'].dtype= np.float64

    return targets

def load_targets(object_lists):
    targets_all = []
    programs = {"variables": 1e18,
                "backup": 1,
                "GW": -100000000,
                "transients": 1e21,
                "NEO": 10,
                "time_sensitive": 10000,
                "jan_objects": 20,
                "rgd_objects": 1e20,
                "periodic_time_sensitive": 30,
                "kevin_objects": 1e23,
                "ilaria_objects": 1e24}

    for program in programs:
        priority = programs[program]
        object_list_dir = os.path.join(object_lists, program)
        filenames = glob.glob(os.path.join(object_list_dir, '*'))
        for filename in filenames:
            if program == "variables":
                targets = load_variables(filename)
            elif program == "backup":
                targets = load_backup(filename)
            elif program == "GW":
                targets = load_GW(filename) 
            elif program == "transients":
                targets = load_transients(filename)
            elif program == "NEO":
                targets = load_NEO(filename)
            elif program == "time_sensitive":
                targets = load_time_sensitive(filename)
            elif program == "periodic_time_sensitive":
                targets = load_periodic_time_sensitive(filename)
            elif program == "jan_objects":
                targets = load_jan_objects(filename)
            elif program == "kevin_objects":
                targets = load_kevin_objects(filename)
            elif program == "ilaria_objects":
                targets = load_ilaria_objects(filename)
            elif program == "rgd_objects":
                targets = load_rgd_objects(filename)
            else:
                print("How do I load objects from the %s program?" % program)
                exit(0)
            if len(targets) == 0: continue
            targets['priority'] = targets['priority'] + priority 
   
            if not 'redo' in targets.columns:
                targets['redo'] = np.zeros(targets["ra"].shape)
                #targets['redo'] = np.ones(targets["ra"].shape)
 
            if not 'delta_redo' in targets.columns:
                targets['delta_redo'] = np.zeros(targets["ra"].shape)
            
            targets.remove_column("coords")
            targets_all.append(targets)
    
    targets = vstack(targets_all)

    return targets

