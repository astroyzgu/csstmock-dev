from .skycov_healpy import *
from .mask  import *
from .fiber import *
from .LSScats import * 

# from .visual_old import * 
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

def fracmap1d_plot(X, indx_spec, bins, ipix = None, labels = None): 
    nlabel = X.shape[1]; nsample = X.shape[0]
    if ipix is None: ipix = np.zeros(nsample)
    print(ipix)
    indx_spec = indx_spec.astype('float')
    if labels is None: labels = ['Feature %s'%ii for ii in range(nlabel)]
    fig, axs = plt.subplots(1, nlabel, figsize = (5*nlabel, 5)) 
    for ii in range(X.shape[1]): 
        ax = axs[ii]
        if isinstance(bins, int): bins_in = [bins]
        if isinstance(bins,list): bins_in = [bins[ii]]  
        frac1d, edges = fracmapdd(X[:,ii], indx_spec, bins = bins_in, ipix = ipix) 
        x = (edges[0][1:] + edges[0][:-1])*0.5 
        h,  _ = np.histogram(X[:,ii], bins = edges[0])
        h1, _ = np.histogram(X[:,ii][indx_spec!=0], bins = edges[0]) 
        ax.plot(x,  h/np.sum(h), drawstyle = 'steps-mid', color = 'b', label = 'Data')
        ax.plot(x, h1/np.sum(h), drawstyle = 'steps-mid', color = 'r', label = 'Data')
        ax.plot(x, frac1d,color = 'gray', alpha = 0.5) 
        ax.plot(x, h1/h, drawstyle = 'steps-mid', marker = '*', color = 'k', label = 'Data')
        ax.set_xlabel(labels[ii])
        ax.set_ylabel('Fraction of Spec') 
    plt.show() 

def cartview(vmap):
    from healpy.newvisufunc import projview, newprojplot
    nside = hp.npix2nside(len(vmap)) 
    projview(
        vmap,
        coord=["G"],
        graticule=True,
        graticule_labels=True,
        unit="surface number density [deg^-2]",
        xlabel="longitude",
        ylabel="latitude",
        title='nside=%d'%nside, 
        norm ='log', 
        cb_orientation="vertical",
        projection_type="cart",
    ); 
    plt.tight_layout()
    plt.show()