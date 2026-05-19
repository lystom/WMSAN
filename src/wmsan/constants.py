import numpy as np
from math import log, radians

"""
Constants for the wmsan package.
"""

# Package version
__version__ = "2026.1.0"

## General constants
R_E = 6.371*1e6 # radius of the earth in meters
LG10 = log(10) # log of 10
g = 9.81


## Constants for the primary microseisms

## Constants for secondary microseisms
RES_MOD = radians(0.5) # resolution of WW3 P2L model files