import numpy as np
import pandas as pd
import multiprocessing as mp
from itertools import repeat
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from scipy.spatial import cKDTree, KDTree
from sklearn.neighbors import BallTree 

def masscomplete(z, mag, lmass, mag_limit = 25, completeness = 95, zbin = None, right= True, include_lowest = True, silence = True): 
    if not silence: print( '#################################################################' ) 
    if not silence: print( '######################### masscomplete --- mag_limit = %0.3f '% mag_limit )
    if not silence: print( '######################### completeness = %0.3f '% completeness + '%' )
    tab = np.array([z, mag, lmass]).T
    nsize   = np.shape(tab)[0]
    max_mag = np.max(mag)
    min_mag = np.min(mag)
    if not silence: print( '### The size of input sample: %i'%nsize)
    df   = pd.DataFrame(tab, columns = ['z', 'mag', 'lmass'])
    if not silence: print( '### After appling magnitude limit and redshift range')
    df   = df[(df.mag < mag_limit)&(df.z <= max(zbin))&(df.z >= min(zbin))]
    if not silence: print( df.describe() )
    df['zslices'], bins  = pd.cut(df.z, zbin, retbins = True, right= right, include_lowest = include_lowest)
    zslices           = df.zslices.value_counts().to_frame().sort_index()
    zslices           = zslices.rename(columns={'zslices': 'ncount'})
    zslices['lmcomp'] = np.nan
    zslices['z_mean']       =  [0.5*(zbin[ii]+zbin[ii+1]) for ii in range(len(zbin)-1)]
    zslices['z_up']         =  [zbin[ii+1] for ii in range(len(zbin)-1)]
    zslices['z_lo']         =  [zbin[ii]   for ii in range(len(zbin)-1)]

    groupby = df.groupby('zslices') 
    for zslice_index in zslices.index: 
        df_ = groupby.get_group(zslice_index)
        mag_threhold = np.percentile(df_.mag, 80, interpolation = 'midpoint')
        df_ = df_[df_.mag >= mag_threhold]
        df_['lmass_scaled'] = df_.lmass + 0.4*(df_.mag - mag_limit)
        zslices.loc[zslice_index, 'lmcomp'] =   np.percentile(df_.lmass_scaled, completeness,  interpolation = 'midpoint')
    if not silence: print( zslices )
    return zslices