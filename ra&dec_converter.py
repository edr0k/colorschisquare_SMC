# -*- coding: utf-8 -*-
"""
Created on Sat Aug 17 21:05:15 2019

@author: souza
"""

from astropy import units as u
from astropy.coordinates import SkyCoord

# c = SkyCoord('21 33 27.02 -00 49 23.7', unit=(u.hourangle, u.deg))
# print c
ra_dec = input('Type Ra (hour) and Dec (deg) between quotation marks: ')
print
ra_dec
c = SkyCoord(ra_dec, unit=(u.hourangle, u.deg))
print
c
