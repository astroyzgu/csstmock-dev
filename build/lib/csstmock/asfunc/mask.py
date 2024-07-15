import healpy as hp 
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt 
import time

def digitize(X, edges, ipix = None):
    '''
    digitize the data into the given bins.

    Parameters:
    - X (array): The input data with shape (n_samples, n_features).
    - edges (list): The bin edges for each feature.
    - ipix (array, optional): The pixel indices for the data. Default is None.

    Returns:
    - slice (tuple): The tuple of indices for each bin. 
    '''
    if ipix is None: 
        ndim = len(edges)
    else: 
        ndim = len(edges) - 1
    index = []; 
    # print(len(edges), ndim) 
    for ii in range(len(edges)):
        bin = edges[ii]
        if ii <= ndim - 1: xxx = X[:,ii]; 
        if ii == ndim: xxx = ipix; 
        ind = np.digitize(xxx, bin, right=False) 
        ind[ind==0] = 1
        ind[ind==len(bin)] = len(bin) - 1
        # print(ii, len(bin), len(ind) )
        index.append(ind-1)
    index  = np.array(index).T
    irange = np.arange(np.shape(index)[1]); print(irange)
    slice  = tuple(index[:, ii].astype(int) for ii in irange )
    return slice 

def fracmapdd(X, y, bins = 10, ipix = None): 
    '''
    Compute the fraction map of the given data set. 

    Parameters:
    - X (array): The feature matrix with the shape (n_samples, n_features).
    - y (array): The target values with the shape (n_samples,). 
    - bins (int or list): The number of bins or the bin edges for the histogram. 
    - ipix (array, optional): The pixel indices for the data points. Default is None.

    Returns:
    - fracdd (array): The nD array of the fraction map.
    - frac1d (list): The list of 1D arrays of the fraction map.
    - edges (list): The list of bin edges for each dimension of the histogram.
    '''
    if np.ndim(X) == 1: X = X.reshape(-1, 1)
    # size = np.shape(X)
    hists1, edges = np.histogramdd(X, bins= bins) 
    if not ipix is None: 
        ipix_unique, ipix_counts = np.unique(ipix, return_counts=True)
    else:
        ipix = np.array( [0]*np.sum(hists1).astype(int) ) 
        ipix_unique = [0] 
        ipix_counts = np.sum(hists1)  
    edges.append( np.hstack([ipix_unique,  ipix_unique[-1]+1]) )
    # print(np.shape(X), np.shape(ipix.reshape(-1,1)) )
    X  = np.concatenate([X, ipix.reshape(-1,1)], axis = 1)
   #  X  = np.hstack([X, ipix.reshape(-1,1)] ).T 
    hists1, edges = np.histogramdd(X, bins= edges) 
    hists2, edges = np.histogramdd(X, bins= edges, weights=y) 
    fracdd = hists1*0.0+np.nan
    fracdd[hists1!=0] = 1.0*hists2[hists1!=0]/hists1[hists1!=0]
    # frac1d = []
    # idim   = range(len(edges))
    # for ii in idim[::-1]:  
    #     axes = tuple( np.delete(np.arange(len(edges)), ii)  )
    #     y2 = np.sum(hists2.T, axis=axes)
    #     y1 = np.sum(hists1.T, axis=axes)
    #     yy = y1*0.0+np.nan
    #     yy[y1!=0] = y2[y1!=0]/y1[y1!=0]; # yerr = np.sqrt(yy*(1-yy)/y1)
    #     # print(ii, axes, np.shape(yy) )
    #     frac1d.append(yy) 
    return fracdd, edges 

def __circle_counting(ra,  dec, u, v, a, w = None): 
    '''
    Count the number of points within a circle of radius a centered at (u, v) in the (ra, dec) coordinates.

    Parameters:
    - ra (array-like): The right ascension coordinates of the points.
    - dec (array-like): The declination coordinates of the points.
    - u (float or array-like): The right ascension of the center of the circle.
    - v (float or array-like): The declination of the center of the circle.
    - a (float): The radius of the circle in degrees.
    - w (array-like, optional): The weights of the points. Default is None.

    Returns:
    - indx (array): The number of points within the circle for each center (u, v). 
    '''
    a = np.atleast_1d(a)
    xyzstar = hp.ang2vec(u,  v,   lonlat = True).reshape(-1, 3)
    xyzgala = hp.ang2vec(ra, dec, lonlat = True).reshape(-1, 3)
    indx =  np.zeros( len(ra) )*0.0
    from scipy.spatial import KDTree
    radius  = 2*np.sin(a/2*np.pi/180)
    kd_tree = KDTree(xyzgala); 
    indice  = kd_tree.query_ball_point(xyzstar, r = radius)
    if w is None: w = np.ones(len(u))
    counts  = np.array([len(ind) for ind in indice]) 
    wht     = np.repeat(w, counts)
    indice  = np.hstack(indice) 
    unq, inv = np.unique( indice,return_inverse=True)
    wcount   = np.bincount(inv, wht.reshape(-1))
    if len(unq) != 0:  indx[unq] = wcount
    return indx

def __ellipse_counting(ra,  dec, u, v, a, pa, ba, w = None): 
    '''
    Count the number of points within an ellipse centered at (u, v) in the (ra, dec) coordinates.

    Parameters:
    - ra (array-like): The right ascension coordinates of the points.
    - dec (array-like): The declination coordinates of the points.
    - u (float or array-like): The right ascension of the center of the ellipse.
    - v (float or array-like): The declination of the center of the ellipse.
    - a (float): The semi-major axis of the ellipse in degrees.
    - pa (float): The position angle of the ellipse in degrees.
    - ba (float): The axis ratio of the ellipse.
    - w (array-like, optional): The weights of the points. Default is None.

    Returns:
    - indx (array): The number of points within the ellipse for each center (u, v). 
    '''

    xyzstar = hp.ang2vec(u,  v,   lonlat = True).reshape(-1, 3)
    xyzgala = hp.ang2vec(ra, dec, lonlat = True).reshape(-1, 3) 
    a = np.atleast_1d(a)
    pa= np.atleast_1d(pa)*np.pi/180.0
    ba= np.atleast_1d(ba)
    from scipy.spatial import KDTree
    radius  = 2*np.sin(a/2*np.pi/180) 
    kd_tree = KDTree(xyzgala); 
    indice  = kd_tree.query_ball_point(xyzstar, r = radius)
    # indice  = np.unique(np.hstack(indice) ).astype(int)
    # indx[indice] = True 
    istar =  np.arange(len(u) ) 
    indx  =  np.zeros( len(ra) )*0.0
    if w is None: w = np.ones(len(u))
    for istar_, x_, y_, a_, pa_, ba_, w_ in zip(istar, u, v, a, pa, ba, w): 
        igala_ = indice[istar_]; 
        if len(igala_) == 0: continue 
        igala_ = np.array(igala_)
        # Transform the point to the ellipse's coordinates
        ra_  = ra[igala_]
        dec_ = dec[igala_]
        ra_prime  =  (ra_  - x_) * np.cos(pa_) * np.cos(dec_*np.pi/180) - (dec_ - y_) * np.sin(pa_)
        dec_prime =  (ra_  - x_) * np.sin(pa_) * np.cos(dec_*np.pi/180) + (dec_ - y_) * np.cos(pa_)
        a_prime  =  a_*ba_ # *np.cos(dec*np.pi/180)
        b_prime  =  a_
        ang      = ra_prime**2/a_prime**2 + dec_prime**2/b_prime**2
        igala_   = igala_[ang < 1]
        if len(igala_) == 0: continue
        indx[igala_] = indx[igala_] + w_
    return indx

def ellipse_counting(ra, dec, u, v, a, pa, ba, w = None):  
    ''' 
    Count the number of points within an ellipse centered at (u, v) in the (ra, dec) coordinates.

    parameters
    ----------
    ra, dec:  float, scalars or array-like
	Angular coordinates of input targets on the sphere
    u, v: float, scalars or array-like
	Angular coordinates of central point of ellipse on the sphere
    a: float, scalars or array-like
	The length of Semi-major axis of ellipse [degree] 
    pa:  float, scalars or array-like
	Position angle [degree]. PA == 0 is North (ΔDEC>0), PA = 90 is WEST (ΔRA > 0). 
    ba:  float, scalars or array-like
	b/a [0,1]. if ba == 1, the shape is circle, namely the cap on the sphere. 
    w: array-like, optional
    Weights of the points. Default is None. 

    Returns
    -------
    indx: array
        Number of points within the ellipse for each center (u, v).
    ''' 
    u   = np.atleast_1d(u) 
    v   = np.atleast_1d(v) 
    if w is None: w = 1
    if isinstance(a,  (int, float)): a   = u*0.0 + a 
    if isinstance(pa, (int, float)): pa  = u*0.0 + pa 
    if isinstance(ba, (int, float)): ba  = u*0.0 + ba
    if isinstance(w,  (int, float)): w   = u*0.0 + w 
    a   = np.atleast_1d(a) 
    pa  = np.atleast_1d(pa) 
    ba  = np.atleast_1d(ba) 

    indx     = np.zeros(len(ra))#.astype(bool) 
    iscircle = (ba==1)
    if np.sum( iscircle) != 0: 
        # print('Number of circle masking:', np.sum(iscircle))
        indx1     = __circle_counting(ra, dec, 
                                   u[iscircle], v[iscircle], a[iscircle], w[iscircle])
        indx      = indx + indx1
    if np.sum(~iscircle) != 0: 
        # print('Number of ellipse masking:', np.sum(~iscircle))
        indx2    = __ellipse_counting(ra, dec, 
                                   u[~iscircle],  v[~iscircle], a[~iscircle], 
                                   pa[~iscircle], ba[~iscircle],w[~iscircle]) 
        indx      = indx + indx2
    return indx

def __circle_masking(ra,  dec, u, v, a): 
    a = np.atleast_1d(a)
    xyzstar = hp.ang2vec(u,  v,   lonlat = True).reshape(-1, 3)
    xyzgala = hp.ang2vec(ra, dec, lonlat = True).reshape(-1, 3)
    indx =  np.zeros( len(ra) ).astype(bool)
    from scipy.spatial import KDTree
    radius  = 2*np.sin(a/2*np.pi/180)
    kd_tree = KDTree(xyzgala); 
    indice  = kd_tree.query_ball_point(xyzstar, r = radius)
    indice  = np.unique(np.hstack(indice) ).astype(int)
    indx[indice] = True
    return indx

def __ellipse_masking(ra,  dec, u, v, a, pa, ba): 
    xyzstar = hp.ang2vec(u,  v,   lonlat = True).reshape(-1, 3)
    xyzgala = hp.ang2vec(ra, dec, lonlat = True).reshape(-1, 3) 
    a = np.atleast_1d(a)
    pa= np.atleast_1d(pa)*np.pi/180.0
    ba= np.atleast_1d(ba)
    from scipy.spatial import KDTree
    radius  = 2*np.sin(a/2*np.pi/180) 
    kd_tree = KDTree(xyzgala); 
    indice  = kd_tree.query_ball_point(xyzstar, r = radius)
    # indice  = np.unique(np.hstack(indice) ).astype(int)
    # indx[indice] = True
    istar = np.arange(len(u) ) 
    indx =  np.zeros( len(ra) ).astype(bool)
    for istar_, x_, y_, a_, pa_, ba_ in zip(istar, u, v, a, pa, ba): 
        igala_ = indice[istar_]; 
        # print(istar_, igala_)
        if len(igala_) == 0: continue 
        igala_ = np.array(igala_)
        # Transform the point to the ellipse's coordinates
        ra_  = ra[igala_]
        dec_ = dec[igala_]
        ra_prime  =  (ra_  - x_) * np.cos(pa_) * np.cos(dec_*np.pi/180) - (dec_ - y_) * np.sin(pa_)
        dec_prime =  (ra_  - x_) * np.sin(pa_) * np.cos(dec_*np.pi/180) + (dec_ - y_) * np.cos(pa_)
        a_prime  =  a_*ba_ # *np.cos(dec*np.pi/180)
        b_prime  =  a_
        ang      = ra_prime**2/a_prime**2 + dec_prime**2/b_prime**2
        igala_   = igala_[ang < 1]
        if len(igala_) == 0: continue
        indx[igala_] = True
    return indx 

def ellipse_plot(x, y, a, pa, ba, lonlat = True): 
    theta = np.arange(0, 2*np.pi, 0.01)
    x_prime= a*ba * np.cos(theta)#*np.cos(y*np.pi/180)
    y_prime= a    * np.sin(theta)
    pa     = -pa*np.pi/180.0
    rotation_matrix = np.array([[np.cos(pa), -np.sin(pa)],
                                [np.sin(pa),  np.cos(pa)]])
    rotated_ellipse = np.dot(rotation_matrix, np.array([x_prime, y_prime])) 
    dec  = rotated_ellipse[1, :] + y 
    if lonlat: 
        ra   = rotated_ellipse[0, :]/np.cos(dec*np.pi/180) + x
    else: 
        ra   = rotated_ellipse[0, :] + x
    return ra, dec

def ellipse_masking(ra, dec, u, v, a, pa, ba):  
    ''' 
    Returns a boolean array of the same shape as ra that is True where the position 
(ra, dec) is in the ellipse shape. 

    parameters
    ----------
    ra, dec:  float, scalars or array-like
	Angular coordinates of input targets on the sphere
    u, v: float, scalars or array-like
	Angular coordinates of central point of ellipse on the sphere
    a: float, scalars or array-like
	The length of Semi-major axis of ellipse [degree] 
    pa:  float, scalars or array-like
	Position angle [degree]. PA == 0 is North (ΔDEC>0), PA = 90 is WEST (ΔRA > 0). 
    ba:  float, scalars or array-like
	b/a [0,1]. if ba == 1, the shape is circle, namely the cap on the sphere. 
    Returns
    -------
    vetomap: ndarray
	 An array with size = 12*nside*nside. 1.0 means w/i survey; 0.0 means w/o survey.  
    nside: int 
         nside of the healpy pixelization
    ''' 

    u   = np.atleast_1d(u) 
    v   = np.atleast_1d(v) 
    if isinstance(a,  (int, float)): a   = u*0.0 + a 
    if isinstance(pa, (int, float)): pa  = u*0.0 + pa 
    if isinstance(ba, (int, float)): ba  = u*0.0 + ba 
    a   = np.atleast_1d(a) 
    pa  = np.atleast_1d(pa) 
    ba  = np.atleast_1d(ba) 

    indx     = np.zeros(len(ra)).astype(bool) 
    iscircle = (ba==1)
    if np.sum( iscircle) != 0: 
        # print('Number of circle masking:', np.sum(iscircle))
        indx1     = __circle_masking(ra, dec, 
                                   u[iscircle], v[iscircle], a[iscircle])
        indx[indx1] = True
    if np.sum(~iscircle) != 0: 
        # print('Number of ellipse masking:', np.sum(~iscircle))
        indx2    = __ellipse_masking(ra, dec, 
                                   u[~iscircle],  v[~iscircle], a[~iscircle], 
                                   pa[~iscircle], ba[~iscircle]) 
        indx[indx2] = True
    return indx

def sphrand_uniform(nrand, ramin, ramax, decmin, decmax): 
    '''
	Draw a random sample with uniform distribution on a sphere

    Parameters
    ----------
    nrand : int
        the number of the random points to return
    ramin, ramax: float
        Right Ascension between ramin and ramax [degrees] 
    decmin, decmax: float 
        Declination between decmin and decmax [degrees]
    Returns
    -------
    ra, dec : ndarray
        the random sample on the sphere within the given limits.
        arrays have shape equal to nrand.
    skycov: float 
        sky coverage [deg^2].
    '''
    zmax = np.sin( np.asarray(decmax) * np.pi / 180.)
    zmin = np.sin( np.asarray(decmin) * np.pi / 180.)

    z   = np.random.uniform(zmin, zmax,  nrand)
    DEC = (180. / np.pi) * np.arcsin(z) 
    RA  = np.random.uniform(ramin, ramax, nrand)

    skycov = (zmax - zmin)*180/np.pi *(ramax - ramin)
    return RA, DEC, skycov

def sphrand_healpy(nrand, nside, pix):  #, ramin = None, ramax = None, decmin = None, decmax = None): 
    '''
	Draw a random sample with uniform distribution on the given region of a sphere defined by healpix. 

    Parameters
    ----------
    nrand : int
        the number of the random points to return
    nside: int 
        nside of the healpy pixelization 
    pix: ndarray, int 
        pixel number(s)
    Returns
    -------
    ra, dec : ndarray
        the random sample on the sphere within the given region defined by healpix.
        arrays have shape equal to nrand.
    skycov: float 
        sky coverage [deg^2].
    ''' 
    pix      = np.asarray(pix).astype(int)
    pixarea  = hp.nside2pixarea(nside, degrees = True)
    skycov_healpy = pixarea*len(pix) 
    
    lon, lat = hp.pix2ang(nside, pix, lonlat = True)
    indx_box = np.hstack([np.argmax(lon), np.argmin(lon), np.argmax(lat), np.argmin(lat)])
    vec      = hp.boundaries(nside, pix[indx_box], step = 1)
    lon, lat = hp.vec2ang(vec, lonlat = True) 
    ramax  = np.max(lon); ramin  = np.min(lon) 
    decmax = np.max(lat); decmin = np.min(lat)
    nrand_ = 0
    RA = []; DEC = []
    while nrand_ < nrand:
        arand, drand, __skycov__ = sphrand_uniform( int(nrand*1.2), ramin, ramax, decmin, decmax)
        pix_rand = hp.ang2pix(nside, arand, drand, lonlat = True) 
        indx     = np.isin(pix_rand, pix)
        nrand_ = nrand_ + np.sum(indx)
        RA.append( arand[indx]) 
        DEC.append(drand[indx]) 
        # print('Generating %s random points. Targeting --> %s.'%(nrand_, nrand) )
    RA   = np.hstack(RA)
    DEC  = np.hstack(DEC)
    indx = np.arange(nrand).astype('int') # , replace = False)
    return RA[indx], DEC[indx], skycov_healpy
    