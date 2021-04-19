
import os
from astropy.time import Time, TimeDelta
import astropy.units as u

start_time = Time('2020-10-29T00:00:00', format='isot')
ndays = 170

for ii in range(ndays):
    time = start_time + TimeDelta(ii*u.day)
    system_command = "python kp84_slack_bot.py -d --day %s" % time.isot.split("T")[0].replace("-","")
    os.system(system_command)
