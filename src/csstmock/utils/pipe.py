import time
import numpy  as np
import healpy as hp 
import Corrfunc
from matplotlib import pyplot as plt
from scipy.spatial import KDTree
import os 
from astropy.table import Table 
class pipe(object):  
    def pipe_asfunc(ra, dec, iszspec, targetid = None, zspec = None, zphot = None, outputdir = './output/', outputfmt = 'txt'):
        '''
        Calculate angular selection function. 
		'''
        ra      = np.atleast_1d(ra) 
        dec     = np.atleast_1d(dec)  
        iszspec = np.atleast_1d(iszspec) 

        nsides = [2**ii for ii in range(6,13)]
        #
		#---- healpy weight
		#
        from .cftools import upweight_healpy
        for nside in nsides: 
            ipix,  new_w, new_iszspec = upweight_healpy(ra, dec, iszspec, nside)
            subdir   = os.path.join(outputdir, 'upweight_healpix%s'%nside)
            if not os.path.isdir(subdir): os.makedirs(subdir) 
            pipe.pipe_save([ipix,  new_w, new_iszspec.astype('int') ], ['ipix', 'wht', 'select'], 'upwht', dir = subdir, outputfmt = outputfmt)

        if zspec is not None: 
            # nearest neigbour 
            pass 

        pass

    def pipe_save(vals, colnames, name, dir = './', outputfmt = 'txt'): 
        '''
        Save function for the pipelines. 
        '''
        outputfmt = outputfmt.replace('.', '') 
        validfmts = ['txt', 'fits']
        if not outputfmt in validfmts: 
            msg = '"%s" is not in the valid list. Valid list is %s'%(outputfmt, validfmts)
            raise ValueError(msg) 
        else: 
            savename = os.path.join(dir, name + '.' +  outputfmt)
            

            t = Table() 
            for val, col in zip(vals, colnames): t[col] = val 
            if outputfmt == 'fits': t.write(savename, overwrite = True)  
            if outputfmt == 'txt':  t.write(savename, format = 'ascii.fixed_width', delimiter = None,  overwrite=True)
            del t
            return 0  

        