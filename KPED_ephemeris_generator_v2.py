import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy import wcs
import astropy.units as u
import argparse
from astroquery.jplhorizons import Horizons

'''
KPED_ephemeris_generator_v2.py

Generates KPED-friendly commands for tracking and observing Solar System objects with the KPED system.

Python modules required:

argparse
numpy
astropy
astroquery

directions for use:
ipython -i -- finder_chart_bolin_edits.py @obj "2020 AV2" @ta 2458878.5 @pt "0,20,1" 


req. ID, prog. ID, obj. ID, RA, Dec, epoch, dRA/dt ("/hr), ddec/dt ("/hr), V, exp. (s), filter, camera mode, PI name, time of obs'
08XXXXXX,1,ZTF05yY,20:45:46.000,+02:52:07.00,2000.0,817.20,-648.00,18.8,30.000,FILTER_SLOAN_R,9,Bryce Bolin, 2019-Aug-26 06:08:16.793

#sample executions

python KPED_ephemeris_generator_v2.py @ri 08XXXXXX @pid 8 @on "2020 AV2" @fi FILTER_SLOAN_R @exp 10.000 @pi bbolin @cm 9 @sd 2017-10-01 @ed 2017-10-02 @ts 10m


'''
parser = argparse.ArgumentParser(prefix_chars='@')
parser.add_argument("@ri", "@@request_id", help="request id, e.g., 08XXXXXX")
parser.add_argument("@pid", "@@program_id", help="program id, e.g., 8")
parser.add_argument("@on", "@@object_name", help="object name, e.g., 2020 AV2")
parser.add_argument("@fid", "@@filter", help="filter, e.g., FILTER_SLOAN_R")
parser.add_argument("@exp", "@@exposure_time_s", help="exposure time in s, e.g., 10")
parser.add_argument("@pi", "@@pi_name", help="pi name, e.g., Bryce_Bolin")
parser.add_argument("@cm", "@@camera_mode", help="camera mode, e.g., 9")
parser.add_argument("@sd", "@@start_date", help="date to start the ephemeris, e.g., '2017-10-01'")
parser.add_argument("@ed", "@@end_date", help="date to end the ephemeris, e.g., '2017-10-02'")
parser.add_argument("@ts", "@@time_step", help="time step to produce the ephemeris in minutes, e.g., '10m'")
args = parser.parse_args()


request_id = str(args.request_id)
program_id = str(args.program_id)
object_name = str(args.object_name)
filter = str(args.filter)
exposure_time_s = str(args.exposure_time_s)
pi_name = str(args.pi_name)
camera_mode = str(args.camera_mode)
start_date = str(args.start_date)
end_date = str(args.end_date)
time_step = str(args.time_step)

hours_to_deg = 15.0
years_to_months  = 12.
years_to_days = 365.
days_to_hours = 24.
hours_to_minutes = 60.
minutes_to_seconds = 60.
deg_to_arcsec = 3600.
djcal_precision = 5
astrores_seconds_rounding_format_RA = 2
astrores_seconds_rounding_format_Dec = 1
equinox= "2000.0"
space = ' '
dash = '-'

#functions

def dictionary_month_string_two_digit(numerical_month_two_digit_string):
    month = dict([['01','Jan'], ['02','Feb'], ['03','Mar'], ['04','Apr'], ['05','May'], ['06','Jun'], ['07','Jul'], ['08','Aug'], ['09','Sep'], ['10','Oct'], ['11','Nov'], ['12','Dec']])
    return month[numerical_month_two_digit_string]

def get_current_pos(object_name, start_date, end_date, time_step):
    if 'COMET' in object_name:
        obj = Horizons(id=object_name[object_name.find(' ')+1:], epochs={'start':start_time, 'stop':end_time,'step':step}, id_type='id')
    if not 'COMET' in object_name:
        obj = Horizons(id=object_name, epochs={'start':start_date, 'stop':end_date,'step':time_step})
    eph = obj.ephemerides()
    pos = SkyCoord(eph['RA'], eph['DEC'], unit='deg')
    vec = obj.vectors()
    times = Time(np.asarray(vec['datetime_jd']), format='jd', scale='utc')
    return np.asarray(pos.ra),np.asarray(pos.dec),times, np.asarray(eph['RA_rate']), np.asarray(eph['DEC_rate']), np.asarray(eph['V'])

def convert_deg_to_hms_RA(deg):
    decimal_m, h = np.modf(deg/hours_to_deg )
    decimal_s, m = np.modf(np.abs(decimal_m) * hours_to_minutes)
    s = np.round(decimal_s*minutes_to_seconds,astrores_seconds_rounding_format_RA)
    return h,m,s

def convert_deg_to_dms_Dec(deg):
    decimal_m, d = np.modf(deg)
    decimal_s, m = np.modf(np.abs(decimal_m) * hours_to_minutes)
    s = np.round(decimal_s*minutes_to_seconds,astrores_seconds_rounding_format_Dec)
    return d,m,s

ra_deg, dec_deg, times,ra_rate_arcsec_hr, dec_rate_arcsec_hr, V_mag = get_current_pos(object_name, start_date, end_date, time_step)
t_isot = times.isot

entire_date_string = np.copy(times.mjd.astype('string'))
for i in range(0,len(entire_date_string)):
    t_date_storage = t_isot[i][:t_isot[i].find('T')]
    year_storage = t_date_storage[:t_date_storage.find('-')]
    month_storage = dictionary_month_string_two_digit(t_date_storage[t_date_storage.find('-')+1: t_date_storage.find('-')+3])
    day_storage = t_date_storage[t_date_storage.find('-')+4:]
    time_storage = t_isot[i][t_isot[i].find('T')+1:]
    entire_date_string[i] = year_storage + dash + month_storage + dash + day_storage + space + time_storage

RA_H,RA_M,RA_S = convert_deg_to_hms_RA(ra_deg)
Dec_d,Dec_m,Dec_s = convert_deg_to_dms_Dec(dec_deg)

RA_H_string,RA_M_string,RA_S_string = RA_H.astype('int').astype('string'),RA_M.astype('int').astype('string'),np.round(RA_S,2).astype('string')
Dec_d_string,Dec_m_string,Dec_s_string =Dec_d.astype('int').astype('string'),Dec_m.astype('int').astype('string'),np.round(Dec_s,1).astype('string')

Dec_sign = np.copy(Dec_d_string)
Dec_sign[np.where(Dec_d>=0.0)] = '+'
Dec_sign[np.where(Dec_d<0.0)] = '-'
Dec_d_string = np.abs(Dec_d).astype('int').astype('string')


#times_string_isot_no_t = np.array([s.replace('T','-') for s in t.isot])
print 'req. ID, prog. ID, obj. ID, RA, Dec, epoch, dRA/dt ("/hr), ddec/dt ("/hr), V, exp. (s), filter, camera mode, PI name, time of obs'


for i in range(0, len(entire_date_string)):
    print '%s,%s,%s,%s,%s,%s,%.2f,%.2f,%.1f,%s,%s,%s,%s,%s'%(request_id, program_id, object_name, RA_H_string[i].zfill(2)+":"+RA_M_string[i].zfill(2)+":"+RA_S_string[i].zfill(5), Dec_sign[i].zfill(1)+Dec_d_string[i].zfill(2)+":"+Dec_m_string[i]+":"+Dec_s_string[i].zfill(4), equinox, ra_rate_arcsec_hr[i], dec_rate_arcsec_hr[i], V_mag[i], exposure_time_s, filter, camera_mode, pi_name, entire_date_string[i])
