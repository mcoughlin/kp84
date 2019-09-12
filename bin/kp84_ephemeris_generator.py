import sys
import pyslalib.slalib as sla
import os
import argparse
import numpy as np
from astropy.time import Time


'''

KPED_ephemeris_generator.py

Generates KPED-friendly commands for tracking and observing Solar System objects with the KPED system.

Python modules required:

sys
pyslalib
os
argparse
numpy
time
astropy

directions for use:

directly copy output from https://www.projectpluto.com/ephem.htm for obvservables using:

Number of steps : 2160
Step size: 0.00138888 (2 minutes)

pbpaste > project_pluto_test; python KPED_ephemeris_generator.py -if project_pluto_test -ri 08XXXXXX -pi 8 -on p 10.000 -pn bbolin -cm 9 > ZTF064m_ephemerides_2019_08_27

input into project_plut_test from project pluto is like:

Date (UTC).dddddddd   RA              Dec         delta   r     elong  mag  RA '/hr dec  " sig PA
---- --------------  --------------  -----------  -----  ------ -----  --- ------------- --------
2019 08 16.07225723  02 36 23.275   +04 00 07.47  .78944 1.4306 104.4 20.5   1.32   0.01 1.64 125

Arguments are:

-if: name of file containing project pluto output from https://www.projectpluto.com/ephem.htm
-ri: request id, e.g., 08XXXXXX we will get a better formal one
-pi: id number for PI, e.g. 8, we will get a formal one later
-on: object name, e.g., ZTF05yY
-fi: filter id, FILTER_SLOAN_R is good for Solar System objects.
-exp: total integration time in data cube
-pn: PI name, e.g., bbolin
-cm: camera mode, e.g., 9 from:

1: science mode, bin 1, fastest frames (a bit over 1s)
2: science mode, bin 1, 2s frames
3: science mode, bin 1, 3s frames
4: science mode, bin 1, 5s frames
5: science mode, bin 1, 10s frames
6: science mode, bin 2, 1s frames
7: science mode, bin 2, 2s frames
8: science mode, bin 2, 3s frames
9: science mode, bin 2, 10s frames
10: rapid readout mode, bin 1, fastest frames (~8.6Hz) (Gain = 300)
11: rapid readout mode, bin 1, 1s frames  (Gain = 300)

output format is in the following format:


req. ID, prog. ID, obj. ID, RA, Dec, epoch, dRA/dt ("/hr), ddec/dt ("/hr), V, exp. (s), filter, camera mode, PI name, time of obs'
08XXXXXX,1,ZTF05yY,20:45:46.000,+02:52:07.00,2000.0,817.20,-648.00,18.8,30.000,FILTER_SLOAN_R,9,Bryce Bolin, 2019-Aug-26 06:08:16.793

'''

#functions

def dictionary_month_string_two_digit(numerical_month_two_digit_string):
    month = dict([['01','Jan'], ['02','Feb'], ['03','Mar'], ['04','Apr'], ['05','May'], ['06','Jun'], ['07','Jul'], ['08','Aug'], ['09','Sep'], ['10','Oct'], ['11','Nov'], ['12','Dec']])
    return month[numerical_month_two_digit_string]

parser = argparse.ArgumentParser()
parser.add_argument("-if", "--infile", help="location of input ephemeris from Project Pluto e.g., /Users/bolin/NEO/Follow_up/APO_observing/project_pluto_test")
parser.add_argument("-ri", "--request_id", help="request id, e.g., 08XXXXXX")
parser.add_argument("-pi", "--program_id", help="program id, e.g., 8")
parser.add_argument("-on", "--object_name", help="object name, e.g., ZTF05iT")
parser.add_argument("-fi", "--filter", help="filter, e.g., FILTER_SLOAN_R")
parser.add_argument("-exp", "--exposure_time_s", help="exposure time in s, e.g., 30")
parser.add_argument("-pn", "--pi_name", help="pi name, e.g., Bryce_Bolin")
parser.add_argument("-cm", "--camera_mode", help="object name, e.g., ZTF05iT")
args = parser.parse_args()

infile = str(args.infile)
request_id = str(args.request_id)
program_id = str(args.program_id)
object_name = str(args.object_name)
filter = str(args.filter)
exposure_time_s = str(args.exposure_time_s)
pi_name = str(args.pi_name)
camera_mode = str(args.camera_mode)

equinox= "2000.0"
space = ' '
dash = '-'

os.system('sed \'/----/,$!d\' ' + infile + ' | tail -n+2 > ' + object_name +'_tmp')
year,month,date,RA_H,RA_M,RA_S,Dec_d,Dec_m,Dec_s,V_mag,ra_rate_arcmin_hr,dec_rate_arcmin_hr = np.loadtxt(object_name +'_tmp',usecols =(0,1,2,3,4,5,6,7,8,12,13,14)).T
os.system('rm ' + filename_out)



date_mjd_UTC = np.array([map(sla.sla_caldj,year,month,date)]).T[0].T[0] + np.modf(date)[0]
t = Time(date_mjd_UTC, format='mjd', scale='utc')
t_isot = t.isot
entire_date_string = np.copy(date_mjd_UTC.astype('string'))
for i in range(0,len(entire_date_string)):
    t_date_storage = t.isot[i][:t.isot[i].find('T')]
    year_storage = t_date_storage[:t_date_storage.find('-')]
    month_storage = dictionary_month_string_two_digit(t_date_storage[t_date_storage.find('-')+1: t_date_storage.find('-')+3])
    day_storage = t_date_storage[t_date_storage.find('-')+4:]
    time_storage = t_isot[i][t_isot[i].find('T')+1:]
    entire_date_string[i] = year_storage + dash + month_storage + dash + day_storage + space + time_storage


RA_H_string,RA_M_string,RA_S_string = RA_H.astype('int').astype('string'),RA_M.astype('int').astype('string'),RA_S.astype('string')
Dec_d_string,Dec_m_string,Dec_s_string =Dec_d.astype('int').astype('string'),Dec_m.astype('int').astype('string'),Dec_s.astype('string')
ra_rate_arcsec_hr, dec_rate_arcsec_hr =ra_rate_arcmin_hr*60., dec_rate_arcmin_hr*60.#arcmin_hr to arcsec_hr

Dec_sign = np.copy(Dec_d_string)
Dec_sign[np.where(Dec_d>=0.0)] = '+'
Dec_sign[np.where(Dec_d<0.0)] = '-'
Dec_d_string = np.abs(Dec_d).astype('int').astype('string')


#times_string_isot_no_t = np.array([s.replace('T','-') for s in t.isot])
print('req. ID, prog. ID, obj. ID, RA, Dec, epoch, dRA/dt ("/hr), ddec/dt ("/hr), V, exp. (s), filter, camera mode, PI name, time of obs')


for i in range(0, len(entire_date_string)):
    print('%s,%s,%s,%s,%s,%s,%.2f,%.2f,%.1f,%s,%s,%s,%s,%s'%(request_id, program_id, object_name, RA_H_string[i]+":"+RA_M_string[i]+":"+RA_S_string[i], Dec_sign[i]+Dec_d_string[i]+":"+Dec_m_string[i]+":"+Dec_s_string[i], equinox, ra_rate_arcsec_hr[i], dec_rate_arcsec_hr[i], V_mag[i], exposure_time_s, filter, camera_mode, pi_name, entire_date_string[i]))
