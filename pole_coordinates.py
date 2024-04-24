#
# Time-dependent South Pole coordinate conversion library
# This is an approximation with errors of order ~1 ft.
#

import datetime
import numpy as np
import pygeodesy as geo
from pygeodesy.ellipsoidalKarney import LatLon

M_PER_FT = 0.3048

# Fit to pole marker movement, 1989–2020
V_M_PER_YR = 10.04
THETA_FLOW_DEG = 129.72

# T0 is fit to pole marker 2000
# Pole northing, easting in ft
POLE2000_N = 50809.76
POLE2000_E = 49491.12

# Convert pole (northing, easting) in feet to
# UPS in meters, given a datetime.date
def pole2ups(northing, easting, date, falsed=False):

    dt_years = (date-datetime.date(2000, 1, 1)).days/365.24
    ups_n_ft = northing + np.sin(np.deg2rad(THETA_FLOW_DEG))*V_M_PER_YR*FT_PER_M * dt_years - POLE2000_N
    ups_e_ft = easting + np.cos(np.deg2rad(THETA_FLOW_DEG))*V_M_PER_YR*FT_PER_M * dt_years - POLE2000_E

    ups_n_m = ups_n_ft*M_PER_FT
    ups_e_m = ups_e_ft*M_PER_FT
    if (falsed):
        ups_n_m += 2000000
        ups_e_m += 2000000

    return (ups_n_m, ups_e_m)

# Convert UPS in meters to pole (northing, easting) in feet
# given a datetime.date
def ups2pole(ups_n_m, ups_e_m, date, falsed=False):
    if (falsed):
        ups_n_m -= 2000000
        ups_e_m -= 2000000

    ups_n_ft = ups_n_m/M_PER_FT
    ups_e_ft = ups_e_m/M_PER_FT

    dt_years = (date-datetime.date(2000, 1, 1)).days/365.24
    northing = ups_n_ft - np.sin(np.deg2rad(THETA_FLOW_DEG))*V_M_PER_YR*FT_PER_M * dt_years + POLE2000_N
    easting = ups_e_ft - np.cos(np.deg2rad(THETA_FLOW_DEG))*V_M_PER_YR*FT_PER_M * dt_years + POLE2000_E

    return (northing, easting)

#
# Convert pole (northing, easting) in feet to
# Lat / Lon, given a datetime.date.  Requires
# pygeodesy
def pole2latlon(northing, easting, date):
    (ups_n, ups_e) = pole2ups(northing, easting, date, falsed=False)
    u = geo.Ups(northing = ups_n, easting = ups_e, pole='S', falsed=False)
    ll = u.toLatLon()
    return(ll.lat, ll.lon)

# Convert lat/lon in fractional degrees (SP at -90.0, 0)
# pole northing, easting in feet, given a
# datetime.date. South pole is at -90, 0.
def latlon2pole(lat, lon, date):
    ups = geo.ups.toUps8(LatLon(lat,lon), pole='S', falsed=False)
    return ups2pole(ups.northing, ups.easting, date, falsed=False)
