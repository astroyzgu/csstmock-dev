from .skycov_healpy import *
from .mask  import *
from .fiber import * 
from .visual_old import * 
import numpy as np 
import healpy as hp 

def pix2border(nside, pixid, nest = False, outlier = True ): 
    '''
    Roughly estimate the boundary of a set of pixels. 
    
    paramter: 
    ---------
    nside: int, nside of healpy, 1,2,4,8,16, ..., 2^
    pixid: array_like, the set of pixels
    nest: bool
    return: 
    ---------
    pixid_edge: array_like, the set of pixels on the envelope (edge)
    '''
    # 扩大一圈外部轮廓并填充空心区域
    pixid_neighs  = hp.get_all_neighbours(nside, pixid, nest = nest)
    pixid_neighs  = pixid_neighs[pixid_neighs >= 0]
    pixid_enlarge = np.union1d(pixid_neighs, pixid);
    # 再次扩大一圈外部的轮廓， 
    # 1. 没有交集的区域即是外围轮廓
    pixid_neighs  = hp.get_all_neighbours(nside, pixid_enlarge, nest = nest); #print(np.shape(pixid_neighs) )
    pixid_neighs  = pixid_neighs[pixid_neighs >= 0]
    pixid_edge    = np.setdiff1d(pixid_neighs, pixid_enlarge, assume_unique = False )
    # 2. 有交集的区域是轮廓
    if outlier is False: 
        pixid_edge = hp.get_all_neighbours(nside, pixid_edge, nest = nest)
        pixid_edge = pixid_edge[pixid_edge >= 0]
        pixid_edge = hp.get_all_neighbours(nside, pixid_edge, nest = nest)
        pixid_edge = pixid_edge[pixid_edge >= 0]
        pixid_edge = np.intersect1d(pixid_edge, pixid, assume_unique = False )
    return pixid_edge 

def fracmap(nside, a, d, mask, zbins = None, z = None): 
    '''
    parameter:
    ----------
    nside: int 
        nside of healpix, all sky is divided into 12*n^2 regions with equal area  (n=1,2,3, ...). 
    a: array-like
        Right Ascension
    d: array-like
        Declination
    bins: ndarray
        The serial number of bins, start from 0 to nbin - 1
    mask: array-like, bool 
        if true, tract as numerator.
    z: array-like
        thrid quantity (e.g., magnitude), if wanting a map of fraction as a function z.
    zbins: array-like
        The bin of thrid quantity (e.g., magnitude).
    return: 
    ----------
    wmap: 
        the map of fraction (1D or 2D if thrid quantity is available)
    pix: array-like
        the region id of healpix 
    w: array-like 
        the weighted value of input galaxies 
    '''

    if zbins is not None: 
        zbins = np.atleast_1d(zbins)
        dbins = zbins[1:] - zbins[:-1]
        dbins, nbins = np.unique(dbins, return_counts = True)
        if len(dbins)!=1: raise('The interval of input zbins is not same.')
        iseq    = np.empty( len(z), dtype = 'int64' )
        iseq[:] = (z-zbins[0])/dbins; nbins = nbins[0]
        wmap = np.zeros( shape = (12*nside*nside, nbins) ) + np.nan 
    else: 
        nbins = 1 
        wmap = np.zeros(12*nside*nside) + np.nan 
#----------------------------------------
    pix = hp.ang2pix(nside, a, d, lonlat = True)
    for ii in range(nbins):  
        if nbins != 1: pix_ = pix[ii  == iseq]
        if nbins != 1: mask_= mask[ii == iseq]
        if nbins == 1: pix_ = pix 
        if nbins == 1: mask_= mask  
        w1 = np.zeros(12*nside*nside) + np.nan # 如果分母没有星系，wmap  = np.nan
        w2 = np.zeros(12*nside*nside) # 如果分子没有找到，分母找到了，wmap = 0
        pix1, count1 = np.unique(pix_, return_counts = True) 
        pix2, count2 = np.unique(pix_[mask_], return_counts = True) 
        w1[pix1] = count1 
        w2[pix2] = count2
        if nbins != 1: wmap[:,ii]  = 1.0*w2/w1
        if nbins == 1: wmap[:]     = 1.0*w2/w1
        
    if nbins != 1: 
        iseq[iseq < 0      ] = 0 
        iseq[iseq > nbins-1] = nbins-1
    else: 
        iseq  = 0
    indx  = pix*nbins + iseq; 
    return wmap, pix, wmap.flatten()[indx] 
