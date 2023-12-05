from .utilsdev import *  
from .starmask import *
from .cftools import * 
from .fiberassign import *
from .mockimg import * 
from .mcompl import *
from .random import * 
from .pipe import * 
# import numpy as np
# import logging 
# import matplotlib.pyplot as plt
# from astropy.cosmology import FlatLambdaCDM
# from scipy.interpolate import interp1d

def xyz2ruv(x, y, z, unpack = True):
    '''
    球坐标-> 笛卡尔坐标
    u0 = 0, ==> u = [-180, 180]; u0 = 180, ==> u = [0, 360]
    v0 = 0, ==> v = [-90, 90];   v0 = 90, ==> v = [0, 180]
    '''
    r = x**2 + y**2 + z**2
    u = np.arctan2(y,x)
    v = np.arctan2(z,np.sqrt(x**2+y**2))
    u = u*180/np.pi;
    v = v*180/np.pi;
    if unpack: 
        return r, u, v
    else: 
        return np.vstack([r, u, v]).T

def ruv2xyz(r, u, v, unpack = True):
    x = r*np.cos(v/180.0*np.pi)*np.cos(u/180.0*np.pi)
    y = r*np.cos(v/180.0*np.pi)*np.sin(u/180.0*np.pi)
    z = r*np.sin(v/180.0*np.pi)
    if unpack: 
        return x, y, z
    else: 
        return np.vstack([x, y, z]).T
