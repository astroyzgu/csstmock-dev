from .skycov_healpy import *
from .mask  import *
from .fiber import *
from .LSScats import * 

# from .visual_old import * 
import numpy as np 
import healpy as hp 

class BinMap:
    def __init__(self, Xbins=None, Cbins=None, interp=None):
        '''
        Class for creating a bin map based on input data.

        Parameters:
        -----------
        Xbins: list, optional
            List of bin edges for the X variables.
        Cbins: list, optional
            List of bin edges for the C variables.
        '''
        self.Xbins = Xbins
        self.Cbins = Cbins
        self.interp = interp

    def interpdd(self, data, method='cubic'):
        from scipy import interpolate
        z = data.copy()
        valid_index  = np.where(~np.isnan(data))
        valid_value  = data[valid_index]
        indices      = np.where( np.isnan(data))
        z[indices]   = interpolate.griddata(valid_index, valid_value, indices, method=method)
        return z

    def assure2d(self, val):
        '''
        Helper function to ensure that the input array is 2-dimensional.

        Parameters:
        -----------
        val: array-like
            Input array.

        Returns:
        --------
        val_2d: array-like
            2-dimensional version of the input array.
        '''
        val = np.array(val)
        ndim = np.ndim(val)
        if ndim == 0:
            return np.empty((0, 2))
        if ndim == 1:
            return val[:, np.newaxis]
        if ndim > 2:
            return val.reshape(-1, val.shape[-1])
        else:
            return val

    def autobins(self, X=None, C=None):
        '''
        Automatically determine the bin edges based on the input data.

        Parameters:
        -----------
        X: array-like, optional
            Input X variables.
        C: array-like, optional
            Input C variables.

        Returns:
        --------
        bins: list
            List of bin edges for X and C variables.
        '''
        X = self.assure2d(X)
        if self.Cbins is None: 
            C  = self.assure2d(C)
            if np.prod(np.shape(C)) != 0:
                Cmaxs = np.max(self.assure2d(C), axis=0)
                self.Cbins = [np.arange(Cmax + 2) for Cmax in Cmaxs]
            else: 
                self.Cbins = []

        if self.Xbins is None: 
            X = self.assure2d(X) 
            if np.prod(np.shape(X)) != 0:
                self.Xbins = [np.histogram_bin_edges(X[:, ii], bins=10, range=None, weights=None) for ii in range(X.shape[1])]
            else:
                self.Xbins = [] 
        self.bins = self.Xbins + self.Cbins
        return self.bins

    def fit(self, X, y, C=None):
        '''
        Fit the bin map using the input data.

        Parameters:
        -----------
        X: array-like
            Input X variables.
        y: array-like
            Target variable.
        C: array-like, optional
            Input C variables.
        '''
        if C is None:
            data = self.assure2d(X)
        elif X is None:
            data = self.assure2d(C)
        else:
            data = np.hstack([self.assure2d(X), self.assure2d(C)])
        bins = self.autobins(X, C)
        histdd1, _ = np.histogramdd(data, bins=bins, weights=y)
        histdd2, _ = np.histogramdd(data, bins=bins, weights=None)
        self.histdd1 = histdd1
        self.histdd2 = histdd2 
        self.valmap = np.zeros_like(histdd2) + np.nan
        self.shape = histdd2.shape
        self.valmap[histdd2 != 0] = histdd1[histdd2 != 0] / histdd2[histdd2 != 0]
        if self.interp is not None: 
            self.valmap = self.interpdd(self.valmap, method=self.interp) 
            self.valmap_pad = np.pad(self.valmap, 1, mode='edge')   
        else: 
            self.valmap_pad = np.pad(self.valmap, 1, mode='constant', constant_values=np.nan)

    def predict_proba(self, X, C=None):
        '''
        Predict the probabilities based on the input data.

        Parameters:
        -----------
        X: array-like
            Input X variables.
        C: array-like, optional
            Input C variables.

        Returns:
        --------
        proba: array-like
            Predicted probabilities.
        '''
        if C is None:
            data = self.assure2d(X)
        elif X is None:
            data = self.assure2d(C)
        else:
            data = np.hstack([self.assure2d(X), self.assure2d(C)])
        queryslice = tuple(np.digitize(data[:, ii], bins=self.bins[ii], right=False) for ii in range(data.shape[1]))
        proba      = self.valmap_pad[queryslice]
        self.index_outer  = np.isnan(proba) 
        num_of_outer = np.sum(self.index_outer)  
        if num_of_outer > 0: 
            if self.interp is None: 
                print('Warning: %s values are out of range, returning to nan.'%num_of_outer ) 
            else: 
                print('Warning: %s values are out of range, returning to use interped values.'%num_of_outer)
        return proba
    
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
    from astropy.coordinates import SkyCoord 
    theta = np.linspace(0, 360, 101)
    gp  = SkyCoord(theta, theta*0, frame='galactic', unit='deg')
    newprojplot(gp.icrs.ra.degree , gp.icrs.dec.degree, color = 'k', lonlat = True)
    # hp.graticule(dpar=30, dmer=30) 
    plt.tight_layout()
    plt.show()