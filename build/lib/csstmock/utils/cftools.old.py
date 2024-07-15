import time
import numpy  as np
import healpy as hp 
import Corrfunc
from matplotlib import pyplot as plt
from scipy.spatial import KDTree
import os 

def angular_density(a, d, max_sep = None, k = None):
    from scipy.spatial import KDTree
    x = 1.0*np.cos(d/180.0*np.pi)*np.cos(a/180.0*np.pi)
    y = 1.0*np.cos(d/180.0*np.pi)*np.sin(a/180.0*np.pi)
    z = 1.0*np.sin(d/180.0*np.pi)

    xyz = np.vstack([x,y,z]).T
    kd  = KDTree(xyz) 
    if max_sep is not None: 
        area = 2*np.pi*(1-np.cos(max_sep*np.pi/180))/(np.pi/180)**2
        max_sep = 2*np.sin(0.5*max_sep*np.pi/180); # 角度转为弧度转为对应笛卡尔坐标下的直线距离
        number = kd.query_ball_point(xyz, r = max_sep, return_length = True) 
        a_den = 1.0*number/area
    else: 
        d, i= kd.query(xyz, k + 1)
        d = d[:, -1]
        d = np.arcsin(0.5*d)*2*180/np.pi
        print( np.min(d),  np.max(d) )
        area  = 2*np.pi*(1-np.cos(d*np.pi/180))/(np.pi/180)**2
        a_den = (k+1)/area
    return a_den

def upweight_nearest(gala, indx, max_sep = 60.0/3600): 
    '''
    gala: array-like 
        a, d, z, w 
    indx: array-like, bool
        galaxies are assigned, if True
    '''
    u = gala[:,0]; v = gala[:,1]; r = gala[:,2]; w = gala[:,3]
    x = 1.0*np.cos(v/180.0*np.pi)*np.cos(u/180.0*np.pi)
    y = 1.0*np.cos(v/180.0*np.pi)*np.sin(u/180.0*np.pi)
    z = 1.0*np.sin(v/180.0*np.pi)
    xyz  = np.vstack([x, y, z]).T 
    max_sep = 2*np.sin(0.5*max_sep*np.pi/180); # 角度转为弧度转为对应笛卡尔坐标下的直线距离
    #------------------------------------------------------------
    ngal  = np.shape(gala)[0]
    index = np.zeros(ngal).astype(bool) 
    igal = np.arange(ngal)
    index[indx]  = True
    igal1 = igal[ index]; #  有光谱的样本编号（igal1）
    igal2 = igal[~index]; # 没有光谱的样本编号(igal2)
    
    kd   = KDTree(  xyz[ index])
    d, i = kd.query(xyz[~index], k = 1, workers = -1)
    select = d < max_sep;
    igal2_ = igal2[   select  ] # 符合条件的，有光谱的样本编号
    igal1_ = igal1[ i[select] ] # 符合条件的，没有有光谱的样本编号
    #------------------------------------------------------------
    # upweighting by nearby untargeted objects 
    uniq1, count1 = np.unique(igal1_, return_counts = True) 
    w[uniq1] = w[uniq1] + count1 # <-- upweighting有光谱的样本
    gala[:,3]      = w
    gala[~index,3] = 0
    # if np.sum(w) != ngal: print('Warning: the sum of w (%.2f) != input number (%.2f)'%(np.sum(w), ngal) )  
    return gala

    
def assignz_nearest(gala, indx, max_sep = 60.0/3600): 
    u = gala[:,0]; v = gala[:,1]; r = gala[:,2]; w = gala[:,3]
    
    x = 1.0*np.cos(v/180.0*np.pi)*np.cos(u/180.0*np.pi)
    y = 1.0*np.cos(v/180.0*np.pi)*np.sin(u/180.0*np.pi)
    z = 1.0*np.sin(v/180.0*np.pi)
    xyz  = np.vstack([x, y, z]).T 
    max_sep = 2*np.sin(0.5*max_sep*np.pi/180);  # 角度转为弧度转为对应笛卡尔坐标下的直线距离
    #------------------------------------------------------------
    ngal  = np.shape(gala)[0]
    index = np.zeros(ngal).astype(bool) 
    igal = np.arange(ngal)
    index[indx] = True
    igal1 = igal[ index]; #  有光谱的样本编号
    igal2 = igal[~index]; # 没有光谱的样本编号
    
    kd   = KDTree(  xyz[ index])
    d, i = kd.query(xyz[~index], k = 1, workers = -1)
    select = d < max_sep
    igal2_ = igal2[   select  ] 
    igal1_ = igal1[ i[select] ]
    #------------------------------------------------------------
    r[igal2_]     = r[ igal1_ ] # assign 距离（也就红移）
    index[igal2_] = True    
    w[igal2_]     = 1 # <--- assign 权重为1 
    
    gala[:,2]= r
    gala[:,3]= w
    gala[~index,3] = 0 # 剩下的权重赋为1（既不是输入的光谱样本， 也不没有assign的）
    #if np.sum(~select) != 0: print('Warning: %.2f targets are not assgined'%np.sum(~select))
    return gala 


import numpy as np 
from astropy.table import Table 
def DD_angularupweight(gala, sbin, angularupweight):
    '''
    the pair-counts of autocorrelation with angular upweighting 
    
    '''
    #-------------------------------------------------------------------------------------
    ngal  = np.shape(gala)[0]
    x = gala[:,2]*np.cos(gala[:,1]/180.0*np.pi)*np.cos(gala[:,0]/180.0*np.pi)
    y = gala[:,2]*np.cos(gala[:,1]/180.0*np.pi)*np.sin(gala[:,0]/180.0*np.pi)
    z = gala[:,2]*np.sin(gala[:,1]/180.0*np.pi)
    w = gala[:,3]
    gala3d = np.vstack([x,y,z]).T 
    x = 1.0*np.cos(gala[:,1]/180.0*np.pi)*np.cos(gala[:,0]/180.0*np.pi)
    y = 1.0*np.cos(gala[:,1]/180.0*np.pi)*np.sin(gala[:,0]/180.0*np.pi)
    z = 1.0*np.sin(gala[:,1]/180.0*np.pi)
    gala2d = np.vstack([x,y,z]).T 
    
    #-------------------------------------------------------------------------------------
    from scipy.interpolate import interp1d
    import timeit
    upwht   = interp1d( np.log10(angularupweight[0]), angularupweight[1], kind = 'nearest', bounds_error=False, fill_value= (angularupweight[1][0], angularupweight[1][-1]))  
    

    #---- too slow 
    from scipy.spatial import cKDTree as KDTree 

    smax   = np.max(sbin);
    smin   = np.min(sbin); 
    
    size = 60
    igal = np.array_split(np.arange(ngal), size)
    
    nloop   = len(sbin) - 1; 
    npairs = np.zeros( (nloop, size, size) )*0.0
    anglew = np.zeros( (nloop, size, size) )*0.0
    pairsw = np.zeros( (nloop, size, size) )*0.0
    
    import tqdm
    for kk  in tqdm.tqdm(range(size*size)): 
        
        rank1 = int(kk%size)
        rank2 = int(kk/size)
        igal1     = igal[rank1] 
        igal2     = igal[rank2] 
        product =  w[igal1]*w[igal2, np.newaxis]
        sep3d = np.sum((gala3d[igal1] - gala3d[igal2,np.newaxis,:])**2,axis = -1) 
        sep2d = np.sum((gala2d[igal1] - gala2d[igal2,np.newaxis,:])**2,axis = -1) 
        sep2d = 2*np.arcsin(0.5*sep2d)*180/np.pi # 换成角度制 
        pairwht = sep2d.copy()*0; nonzero = sep2d != 0 
        pairwht[nonzero] = upwht( np.log10(sep2d[nonzero]) );             
        print(len(igal1), len(igal2), np.shape(product), np.shape(pairwht), np.shape(sep3d) )   

        for ii in range(nloop): 
            indx   = (sep3d>sbin[ii])&(sep3d<=sbin[ii+1])
            npairs[ii, rank1, rank2] = np.sum( indx ) # pair count 
            anglew[ii, rank1, rank2] = np.sum( pairwht[indx] ) # pair count with angular upweight
            pairsw[ii, rank1, rank2] = np.sum( product[indx] ) # product 
        
    npairs = np.sum(npairs, axis = 2).sum(axis = 1) 
    anglew = np.sum(anglew, axis = 2).sum(axis = 1) 
    pairsw = np.sum(pairsw, axis = 2).sum(axis = 1) 
    nonzero = npairs != 0 
    anglew[nonzero] = anglew[nonzero]/npairs[nonzero]
    pairsw[nonzero] = pairsw[nonzero]/npairs[nonzero]
    anglew[~nonzero] = 0 
    pairsw[~nonzero] = 0 
    return npairs, anglew, pairsw


    #-------------------------------------------------------------------------------------
    #---- out of memeroy scipy.cKDTree query_pairs  
#     from scipy.spatial import cKDTree as KDTree 
    
#     method1 = 'scipy.cKDTree query_pairs'
#     method2 = 'scipy.cKDTree query_pairs'

#     smax   = np.max(sbin); 
#     smin   = np.min(sbin); 
#     kd3d   = KDTree(gala3d);
#     indexes = kd3d.query_pairs(smax, output_type = 'ndarray') 
#     print( np.shape(indexes) ) 
#     sep3d   =  np.sqrt( np.sum( (gala3d[indexes[:,0]] - gala3d[indexes[:,1]])**2, axis = 1) )
#     sep2d   =  np.sqrt( np.sum( (gala2d[indexes[:,0]] - gala2d[indexes[:,1]])**2, axis = 1) )
#     product =  w[indexes[:,0]]*w[indexes[:,1]]
    
#     selectpairs = sep3d > smin
#     sep3d       = sep3d[selectpairs]
#     sep2d       = sep2d[selectpairs]
#     product     = product[selectpairs]

#     sep2d = 2*np.arcsin(0.5*sep2d)*180/np.pi # 换成角度制 
#     pairwht = sep2d.copy()*0 + 1; 
#     nonzero = sep2d != 0 
#     pairwht[nonzero] = upwht( np.log10(sep2d[nonzero])); 
    
#     nloop = len(sbin) - 1; 
#     DD    = np.zeros(nloop)*0.0 
#     WW    = np.zeros(nloop)*0.0 
#     for ii in range(nloop): 
#         indx   = (sep3d>sbin[ii])&(sep3d<=sbin[ii+1])
#         DD[ii] = np.sum(  pairwht[indx] )
#         WW[ii] = np.mean( product[indx] )
#     DD    = DD*2 # query_pairs的数目 N*(N-1)/2, 所以计数需要补上
    
    return DD, WW


class cftools(object):
    def __init__(): 
        pass 
    
    def fastinput(data_uvz, data_idx, randfactor, Om0): 
        '''
        
        rand_uvz, rand_idx, wdata, wrand = cftools.astinput(data_uvz, data_idx, randfactor = 10, Om0 = 0.315)
        '''

        ndata1 = np.sum( data_idx)
        ndata2 = np.sum(~data_idx)
        ndata  = ndata1 + ndata2
        if np.shape(data_uvz)[0] != ndata: raise ValueError('len(data_idx)!=np.shape(data_uvz)[0]')

        wdata        = np.zeros( (ndata, 5) )*0.0
        wdata[:,0]   = 1
        wdata[ data_idx,1] = 1
        wdata[ data_idx,2] = 1
        wdata[~data_idx,3] = 1
        wdata[~data_idx,4] = 1

        #---- generate rand sample (adz) and convert to xyz (assuming ΛCDM Om0 = 0.315)
        #---- randfactor   = 10

        nrand1 = randfactor*ndata1
        nrand2 = randfactor*ndata2
        nrand  = nrand1 + nrand2
        rand_idx = np.zeros(nrand).astype(bool) 

        rand_uvz1, Nz1 = cftools.generand(data_uvz[ data_idx], int(1.2*nrand1), dzwindow = 0.1, percentage = True)
        rand_uvz2, Nz2 = cftools.generand(data_uvz[~data_idx], int(1.2*nrand2), dzwindow = 0.1, percentage = True)
        rand_uvz1    = rand_uvz1[np.arange(nrand1)]
        rand_uvz2    = rand_uvz2[np.arange(nrand2)]

        rand_uvz     = np.vstack([rand_uvz1, rand_uvz2])
        rand_idx[np.arange(nrand1)] = True

        wrand        = np.zeros( (randfactor*ndata, 5) )
        wrand[:,0]     = 1
        wrand[ rand_idx,1] = 1
        wrand[~rand_idx,2] = 1
        wrand[ rand_idx,3] = 1
        wrand[~rand_idx,4] = 1

        return rand_uvz, rand_idx, wdata, wrand

    def cf_jk(DD, DR, RR, Ndata_jk, Nrand_jk, axis = 2):

        ndim   = np.ndim(DD); shape = np.array( np.shape(DD) ) 
        if (ndim <= 2)|(ndim!=np.ndim(DR))|(ndim!=np.ndim(RR)):
            raise('Error: dim of DD, DR, RR is not same or <= 2.')

        axis = np.atleast_1d(axis)
        axes = np.arange(ndim)
        axes = tuple( axes[~np.isin(axes, axis)])

        Ndata  = np.sum(Ndata_jk)
        Nrand  = np.sum(Nrand_jk)

        njk           = np.shape(DD)[0]
        output_shape  = shape[axis]
        DD_ = 1.0*np.sum(DD, axis = axes)/Ndata/(Ndata-1) 
        DR_ = 1.0*np.sum(DR, axis = axes)/Ndata/Nrand
        RR_ = 1.0*np.sum(RR, axis = axes)/Nrand/(Nrand-1)
        #print(np.shape(DD_), np.shape(DR_), np.shape(RR_))

        xis = np.zeros(output_shape)*0.0
        nonzeros      =  RR_ != 0
        xis[nonzeros] = (DD_[nonzeros] - 2*DR_[nonzeros] + RR_[nonzeros])/RR_[nonzeros]

        xis_jk = np.zeros((njk, *output_shape))*0.0

        for ijk in np.arange(njk): 
            Ndata_ = Ndata - Ndata_jk[ijk]
            Nrand_ = Nrand - Nrand_jk[ijk]

            DD_ = DD.copy(); DD_[ijk,:] = 0; DD_[:,ijk] = 0;
            DR_ = DR.copy(); DR_[ijk,:] = 0; DR_[:,ijk] = 0;
            RR_ = RR.copy(); RR_[ijk,:] = 0; RR_[:,ijk] = 0;

            DD_ = 1.0*np.sum(DD_, axis = axes)/Ndata_/(Ndata_-1) 
            DR_ = 1.0*np.sum(DR_, axis = axes)/Ndata_/Nrand_ 
            RR_ = 1.0*np.sum(RR_, axis = axes)/Nrand_/(Nrand_-1)  

            # print(DD_[:3], DR_[:3], RR_[:3])
            nonzeros       = RR_ != 0
            xis_jk[ijk,nonzeros]  = (DD_[nonzeros] - 2*DR_[nonzeros] + RR_[nonzeros])/RR_[nonzeros] 

        njk         = np.shape(xis_jk)[0]; 
        mean_xis_jk = np.mean( xis_jk, axis = 0)

        nbin        = np.prod(output_shape) 
        del_xis_jk  = xis_jk.reshape(njk, nbin) - mean_xis_jk.reshape(nbin) 
        # Cij         =  1.0*(njk - 1)/njk*np.dot(del_xis_jk.T, del_xis_jk) 
        Cij         = 1.0/njk*np.dot(del_xis_jk.T, del_xis_jk) 

        # or var_xis_jk  = (njk-1)*np.var(xis_jk, axis = 0)
        Cij         = Cij.reshape((nbin, nbin) )
        iii = np.arange(nbin) 
        var_xis_jk  = Cij[iii, iii]
        var_xis_jk  = var_xis_jk.reshape(output_shape)  
        sig_xi_jk   = np.sqrt(var_xis_jk) 
        return xis, mean_xis_jk, sig_xi_jk, xis_jk, Cij

    
    def jkang(data, rand, njk = 100, nside = 128): 
        from sklearn.cluster import KMeans
        import healpy as hp
        pixdata    = hp.vec2pix(nside, data[:,0], data[:,1], data[:,2])

        uniq_pixdata = np.unique(pixdata)
        uniq_vecdata = np.vstack( hp.pix2vec(nside, uniq_pixdata) ).T

        kmeans  = KMeans(n_clusters = njk, n_init = 'auto', init='k-means++',  random_state=0).fit( uniq_vecdata)
        jk_num       = np.zeros(12*nside*nside) - 1
        jk_num[uniq_pixdata] = kmeans.labels_
        ijk_data   = jk_num[pixdata]
    
        pixrand    = hp.vec2pix(nside, rand[:,0], rand[:,1], rand[:,2])
        ijk_rand   = jk_num[pixrand]

        return ijk_data, ijk_rand


    def uvz2xyz(uvz, Om0 = 0.315): 
        if (np.ndim(uvz)!=2)|(np.shape(uvz)[-1]!=3): raise('Error: uvz format is not (N,3)')    
        from astropy.cosmology import FlatLambdaCDM
        cosmo   = FlatLambdaCDM(100, Om0=0.315);
        '''convert a, d, DC to x, y, z'''
        r = cosmo.comoving_distance(uvz[:,2]).value 
        u = uvz[:,0]
        v = uvz[:,1]
        x = r*np.cos(v/180.0*np.pi)*np.cos(u/180.0*np.pi)
        y = r*np.cos(v/180.0*np.pi)*np.sin(u/180.0*np.pi)
        z = r*np.sin(v/180.0*np.pi)
        return np.vstack([x,y,z]).T

    def generand(data, Nrand, ranges = None, dzwindow = 0.1, percentage = False):
    
        if (np.ndim(data)!=2)|(np.shape(data)[-1]!=3): raise('Error: data format is not (N,3)')    
        amax, dmax, zmax = np.max(data, axis = 0)
        amin, dmin, zmin = np.min(data, axis = 0)

        if ranges is None: 
            ranges  = np.array([amin, amax, dmin, dmax, zmin, zmax])
        if len(ranges) != 6: raise('Error: Parameter of generand, len(ranges)!=6')
        #
        #----- unifrom distribution on the sphere 
        #
        # Lambert azimuthal equal-area projection
        # 从球面上的点（0,0,-1） ==> 投影 ==> z = 0的平面
        # α = arctan(y/x); r = sqrt(x^2+y^2)
        # δ = 2arcsin(r);  
        # r<1的区域生成n个随机数
        ranges  = np.array(ranges)
        a_range = ranges[0:2]*np.pi/180
        r_range = np.sin( 0.5*(90 + ranges[2:4])*np.pi/180 )**2  # (0,1) 
        alpha = np.random.uniform(*a_range, Nrand)
        r     = np.sqrt( np.random.uniform(*r_range, Nrand) )
        delta = 2*np.arcsin(r) - np.pi/2# arcsin==>(0,90)

        #
        #----- sampling redshift distribution 
        #
        zlo = ranges[4]; zup = ranges[5]

        nzbin= 1001
        zbin = np.linspace(zlo, zup, nzbin)
        hist, edge = np.histogram(data[:,2], bins = zbin)

        if percentage: dzwindow = (zup - zlo)*dzwindow

        if dzwindow > (zup - zlo):
            print('warning: dzwindow (%s) larger than the redshift range of data.'%dzwindow)
            dzwindow = (zup - zlo)*0.1
            print('set dzwindow == 0.05*(zmax - zmin) = %s'%dzwindow)

        zref = 0.5*(zbin[1:]+zbin[:-1])
        prob = hist/np.sum(hist)/(edge[1]-edge[0])  
        window_len  = int( np.abs( (nzbin-1)*dzwindow/(zup - zlo)) )
        #smomth_prob = np.convolve(prob, np.ones((window_len,))/window_len, mode='same')
        from scipy.signal import savgol_filter
        smomth_prob = savgol_filter(prob, window_len, 3, mode ='nearest') 

        # print(len(smomth_prob), len(zref), len(hist), window_len) 

        from scipy.interpolate import interp1d
        Nz     = interp1d(zref, smomth_prob, kind = 'linear', bounds_error=False, fill_value=(0,0) )
        smomth_probmax  = np.max( Nz(zref) )

        resamples = []
        while len(resamples) < Nrand:  
            redshift = np.random.uniform(zlo, zup, Nrand)
            prob     = np.random.uniform(0, smomth_probmax*1.05, Nrand) 
            redshift = redshift[prob <=  Nz(redshift)]
            resamples  = np.hstack([resamples, redshift])
        redshift_new = np.random.choice(resamples, Nrand, replace = False) 
        rand  = np.vstack([alpha/np.pi*180, delta/np.pi*180, redshift_new]).T
        return rand, Nz

    def cf(gala, rand, bins, nthreads = 12, angular = False, 
           DD = None, DR = None, RR = None, 
           verbose = 0, return_counts = False, angularupweight = None): 
        '''
        if angular is True, 
            gala/rand ==> ra, dec, 1, wht 
            return: ω(θ),  angular correlation function 

        if angular is False
            gala/rand ==> ra, dec, co-distance, wht 
            return: ξ(s) correlation function
        '''

        gala = gala.copy()
        rand = rand.copy() 
        if verbose != 0: print(np.shape(gala), np.shape(rand)) 
        mu_max   = 1.0
        nmu_bins = 1
        if angular == True: 
            gala[:,2] = 1; 
            rand[:,2] = 1;
            bins = 2*np.sin(0.5*bins*np.pi/180.) 
        if np.shape(gala)[-1] == 3: gala = np.hstack([gala, gala[:,0].reshape(-1,1)*0 +1])
        if np.shape(rand)[-1] == 3: rand = np.hstack([rand, rand[:,0].reshape(-1,1)*0 +1])
        
        t2 = time.time()
        
        Ngala = np.sum(gala[:,3]) 
        Nrand = np.sum(rand[:,3]) 
        # DD_counts 
        if DD is None: 
            DD_counts = Corrfunc.mocks.DDsmu_mocks(1, 1, nthreads, mu_max, nmu_bins, bins,
                             gala[:,0], gala[:,1],gala[:,2], weights1 = gala[:,3], 
                             is_comoving_dist=True, weight_type = "pair_product")
            DD = np.array( [ DD[4]*DD[5] for DD in DD_counts] )
            
        if DR is None: 
            DR_counts = Corrfunc.mocks.DDsmu_mocks(0, 1, nthreads, mu_max, nmu_bins, bins,
                                 RA1=gala[:,0], DEC1=gala[:,1], CZ1=gala[:,2], weights1=gala[:,3], 
                                 RA2=rand[:,0], DEC2=rand[:,1], CZ2=rand[:,2], weights2=rand[:,3],
                                 is_comoving_dist=True, weight_type = "pair_product")
            DR = np.array( [ DR[4]*DR[5] for DR in DR_counts] )
        if RR is None: 
            RR_counts = Corrfunc.mocks.DDsmu_mocks(1, 1, nthreads,  mu_max, nmu_bins, bins,
                                 rand[:,0], rand[:,1], rand[:,2], weights1 = rand[:,3], 
                                 is_comoving_dist=True, weight_type = "pair_product")
            RR = np.array( [ RR[4]*RR[5] for RR in RR_counts] )
            
        #print('Ngala Nrand', Ngala, Nrand)
        #print('DD', DD, DD.dtype, np.shape(DD) )
        #print('DR', DR, DR.dtype, np.shape(DR) )
        #print('RR', RR, RR.dtype, np.shape(RR) )
        if return_counts: 
            return DD, DR, RR; 
        
        DD = 1.0*DD/(Ngala*(Ngala-1)); 
        DR = 1.0*DR/(Ngala*Nrand) 
        RR = 1.0*RR/(Nrand*(Nrand-1))
        corr = 0.0*RR; nonzero = RR != 0.0
        corr[nonzero] = 1.0*(DD[nonzero]-2*DR[nonzero]+RR[nonzero])/RR[nonzero]
        #corr = 1.0*DD/RR 
        t1 = time.time()
        #if verbose != 0: print('It takes %.2f s'% (t1-t2) )
        if not return_counts: 
            return corr
        
    def wp(data, rand, bins, pimax = 20, nthreads = 12):
        '''
        return: χ(rp, pi)
        '''
        data = data.copy()
        rand = rand.copy() 
        DD_counts = DDrppi_mocks(1, 1, nthreads, pimax, bins, 
                                 data[:,0], data[:,1], data[:,2], weights1 = data[:,3],
                                 is_comoving_dist=True, weight_type = "pair_product")
        DR_counts = DDrppi_mocks(0, 1, nthreads, pimax, bins,
                                 RA1=data[:,0], DEC1=data[:,1], CZ1=data[:,2], weights1 = data[:,3],
                                 RA2=rand[:,0], DEC2=rand[:,1], CZ2=rand[:,2], weights2 = rand[:,3],
                                 is_comoving_dist=True, weight_type = "pair_product")
        RR_counts = DDrppi_mocks(1, 1, nthreads, pimax, bins, 
                                 rand[:,0], rand[:,1], rand[:,2], weights1 = rand[:,3],
                                 is_comoving_dist=True, weight_type = "pair_product") 
        Ndata   = np.shape(data)[0]; 
        Nrand   = np.shape(rand)[0]; 
        nrpbins = len(bins) - 1; 
        npibin  = int( len(DD_counts)/nrpbins )

        DD = np.array([DD[4] for DD in DD_counts])
        DR = np.array([DR[4] for DR in DR_counts])
        RR = np.array([RR[4] for RR in RR_counts])

        DD = 1.0*DD/(Ndata*(Ndata-1)); 
        DR = 1.0*DR/(Ndata* Nrand)
        RR = 1.0*RR/(Nrand*(Nrand-1)); 
        corr = 0.0*RR; nonzero = RR != 0.0

        corr[nonzero] = 1.0*(DD[nonzero]-2*DR[nonzero]+RR[nonzero])/RR[nonzero]
        corr = corr.reshape(nrpbins, npibin)
        # wp
        corr = 2*np.sum( corr.reshape(nrpbins, npibin), axis = 1 )

        # convert_rp_pi_counts_to_wp
        #     wp = convert_rp_pi_counts_to_wp(Ndata, Ndata, Nrand, Nrand,
        #                                  DD_counts, DR_counts,
        #                                  DR_counts, RR_counts, nrpbins, pimax)
        # cf = convert_3d_counts_to_cf(Ndata, Ndata, Nrand, Nrand,
        #                             DD_counts, DR_counts, DR_counts, RR_counts) 
        #
        # wp = np.empty(nrpbins)
        # for i in range(nrpbins):
        #     wp[i] = 2.0 * 1.0 * np.sum( cf[i * npibin:(i + 1) * npibin])
        return corr 
        

    def kde_sampling(z, N, method = 'kde', zlo = None, zup = None): 
        '''
        z: input catalogue of redshift 
        N: number need sampling
        '''
        import numpy as np 
        from scipy.interpolate import interp1d
        from scipy.stats import gaussian_kde
        import emcee

        if zlo is None: zlo = np.min(z)
        if zup is None: zup = np.max(z)
        kde  = gaussian_kde(z)
        zref = np.linspace(zlo, zup, 1001) 
        #--- gaussian_kde is too slow, using interp1d to speed up
        prob = kde(zref) 
        Nz   = interp1d( zref, prob, kind = 'linear', bounds_error=False, fill_value=(0,0) )
        logNz= interp1d( zref, np.log(prob), kind = 'linear', bounds_error=False, fill_value=(-np.inf,-np.inf) )

        if method == 'pdf':
            samples = []
            while len(samples) < N:  
                redshift = np.random.uniform(zlo, zup, 10*N)
                probmax  = np.max( Nz(redshift) )
                prob     = np.random.uniform(0, probmax*1.05, 10*N)
                redshift = redshift[prob <=  Nz(redshift)]
                samples  = np.hstack([samples, redshift])

        if method == 'kde':  
            ndim    = 1
            nwalkers= 32;   nrun = int(N/31-500) 
            if nrun < 1000: nrun = 1000
            p0      = np.random.uniform(0, 1, size = (nwalkers, ndim) ) 
            sampler = emcee.EnsembleSampler(nwalkers, ndim, logNz )
            state   = sampler.run_mcmc(p0, nrun)
            samples = sampler.get_chain()
            samples = samples[:500,:,:].flatten()

        if method == 'resample': 
            samples = z

        redshift = samples[(samples >= zlo)&(samples <= zup)] 
        replace = False
        if N > len(redshift): replace = True    
        redshift = np.random.choice(redshift, N, replace = replace) 
        return redshift, Nz
    
    def DRang_IW_jk(data, rand, s_bin, a_bin, 
                    W_data = None, W_rand   = None, 
                  ijk_data = None, ijk_rand = None, 
                  isDD = 0,
                  ijk_data_task = None, 
                  ijk_rand_task = None, 
                  savedir = None): 

        '''return DRang_jk = (ijk_data, ijk_rand, isbin, iabin, iw)'''

        mu_max    = 1;
        nmu_bins  = 1;
        ns_bins   = len(s_bin) + 1;
        na_bins   = len(a_bin) + 1;

        njk = len(np.unique(ijk_data));

        if np.ndim(W_data) == 0: W_data = data[:,0]*0 + 1.0
        if np.ndim(W_data) == 1: W_data = W_data.reshape(-1, 1)
        if np.ndim(W_data) == 2: nw = np.shape(W_data)[1]
        if np.ndim(W_data) >  2: raise ValueError('Error shape of W_data, ndim > 2')

        if np.ndim(W_rand) == 0: W_rand = rand[:,0]*0 + 1.0
        if np.ndim(W_rand) == 1: W_rand = W_rand.reshape(-1, 1)
        if np.ndim(W_rand) == 2: nw = np.shape(W_rand)[1]
        if np.ndim(W_rand) >  2: raise ValueError('Error shape of W_rand, ndim > 2')

        if isDD: 
            if np.shape(data) != np.shape(rand): 
                raise ValueError('isDD is True, but data and rand is not same shape.')
            else: 
                if not np.allclose(data, rand): 
                    raise ValueError('isDD is True, values of data and rand is not same.')
        else: 
            pass 

        if ijk_data is None: ijk_data = data[:,0]*0
        if ijk_rand is None: ijk_rand = rand[:,0]*0

        from csstmock.ext import DRang_IW
        DRang_IW_jkmat  = np.zeros( (njk, njk, ns_bins, na_bins, nw) )

        if ijk_data_task is None: ijk_data_task   = np.arange(njk)
        if ijk_rand_task is None: ijk_rand_task   = np.arange(njk)
        for ijk_data_ in ijk_data_task: 
            for ijk_rand_ in ijk_rand_task: 
                
                if isDD:
                    if(ijk_data_< ijk_rand_): continue
                    if(ijk_data_>=ijk_rand_): pass
                else: 
                    pass

                data_ = data[ijk_data_ == ijk_data,:]
                rand_ = rand[ijk_rand_ == ijk_rand,:]
                W_data_ = W_data[ijk_data_ == ijk_data,:]
                W_rand_ = W_rand[ijk_rand_ == ijk_rand,:]
                mat     = DRang_IW(data_, rand_, s_bin, a_bin, W_data_, W_rand_, isDD = 0)
                    
                if not savedir is None: 
                    if not os.path.isdir(savedir): os.makedirs(savedir) 
                    
                    filename = os.path.join(savedir, 'paircounts-jk-%04i-%04i.npy'%(ijk_data_, ijk_rand_) ) 
                    np.save(filename, mat)

                    filename = os.path.join(savedir, 'paircounts-jk-%04i-%04i.npy'%(ijk_rand_, ijk_data_) ) 
                    if isDD&(ijk_data_>ijk_rand_): np.save(filename, mat) 
                    
                DRang_IW_jkmat[ijk_data_, ijk_rand_] = mat
                if isDD&(ijk_data_>ijk_rand_): 
                    DRang_IW_jkmat[ijk_rand_, ijk_data_] = mat

        return DRang_IW_jkmat

class upwht(object):
    def __init__(): 
        pass 
    def upweight_healpy(data, indx, nside): 
        '''
        
        ''' 
        if (np.ndim(data)==2): 
            if np.shape(data)[1] != 3: raise ValueError('The shape of data is not (N,3)') 
        else: 
            raise ValueError('The shape of data is not (N,3)') 

        hpwht1 = np.zeros( 12*nside*nside)*0.0
        hpwht2 = np.zeros( 12*nside*nside)*0.0
        hpwht  = np.zeros( 12*nside*nside)*0.0
        
        ipix1  = hp.vec2pix(nside, data[indx,0], data[indx,1], data[indx,2]) 
        ipix2  = hp.vec2pix(nside, data[:,0], data[:,1], data[:,2]) 
        
        uniq_ipix1, count1 = np.unique(ipix1, return_counts = True )
        uniq_ipix2, count2 = np.unique(ipix2, return_counts = True )
        
        hpwht1[uniq_ipix1] = count1
        hpwht2[uniq_ipix2] = count2

        nonzeros = hpwht1! = 0.0
        hpwht[nonzeros] = 1.0*hpwht2[nonzeros]/hpwht1[nonzeros]
        neww  = hpwht[ipix2]

        return neww, indx
    
    def upweight_nearest(data, indx, w = None, max_sep = 180): 
        '''
        data: xyz 
        indx: array-like, bool
            galaxies are assigned, if True
        w: weight, if None, all galaxies are weighted by 1. 
        '''
        #u = gala[:,0]; v = gala[:,1]; r = gala[:,2]; w = gala[:,3]
        #x = 1.0*np.cos(v/180.0*np.pi)*np.cos(u/180.0*np.pi)
        #y = 1.0*np.cos(v/180.0*np.pi)*np.sin(u/180.0*np.pi)
        #z = 1.0*np.sin(v/180.0*np.pi)
        #xyz  = np.vstack([x, y, z]).T
        if (np.ndim(data)==2): 
            if np.shape(data)[1] != 3: raise ValueError('The shape of data is not (N,3)') 
        else: 
            raise ValueError('The shape of data is not (N,3)') 
            
        if w is None: w = data[:,0]*0+1 
        r   = np.sqrt( np.sum(data**2, axis = 1) ) 
        xyz = data/r.reshape(-1,1)
        # print(np.shape(xyz), np.shape(r))
        max_sep = 2*np.sin(0.5*max_sep*np.pi/180); # 角度转为弧度转为对应笛卡尔坐标下的直线距离
        #------------------------------------------------------------
        ngal  = np.shape(data)[0]
        
        index = np.zeros(ngal).astype(bool) 
        igal = np.arange(ngal)
        index[indx]  = True
        igal1 = igal[ index]; #  有光谱的样本编号（igal1）
        igal2 = igal[~index]; # 没有光谱的样本编号(igal2)

        kd   = KDTree(   xyz[ index])
        d, i = kd.query( xyz[~index], k = 1, workers = -1)
        select = d < max_sep;
        igal2_ = igal2[   select  ] # 符合条件的，没有有光谱的样本编号, indx[igal2_] === False
        igal1_ = igal1[ i[select] ] # 符合条件的，有光谱的样本编号,    indx[igal1_] === True 

        # upweighting by nearby untargetwed objects (igal1_, igal2_)
        #w      = np.array([1,1,0.5,1])*1.0
        #igal1_ = np.array([0,0,1,3])#*1.0
        #igal2_ = np.array([1,1,2,1])#*1.0
        
        newwht    = w.copy()
        for i1, i2 in zip(igal1_,igal2_):
            newwht[i1] = newwht[i1] + w[i2]
        return newwht, index

    def assignz_nearest(data, indx, z = None, max_sep = 180): 
        if (np.ndim(data)==2): 
            if np.shape(data)[1] != 3: raise ValueError('The shape of data is not (N,3)') 
        else: 
            raise ValueError('The shape of data is not (N,3)')

        r   = np.sqrt( np.sum(data**2, axis = 1) ) 
        xyz = data/r.reshape(-1,1)
        # print(np.shape(r), np.shape(xyz) )
        max_sep = 2*np.sin(0.5*max_sep*np.pi/180);  # 角度转为弧度转为对应笛卡尔坐标下的直线距离
        #------------------------------------------------------------
        ngal  = np.shape(data)[0]
        index = np.zeros(ngal).astype(bool) 
        igal  = np.arange(ngal)
        index[indx] = True
        igal1 = igal[ index]; #  有光谱的样本编号
        igal2 = igal[~index]; # 没有光谱的样本编号

        kd   = KDTree(  xyz[ index])
        d, i = kd.query(xyz[~index], k = 1, workers = -1)
        select = d < max_sep
        igal2_ = igal2[   select  ] 
        igal1_ = igal1[ i[select] ]
        #------------------------------------------------------------
        if z is None: 
            newr = r.copy()
        else: 
            if np.shape(r) == np.shape(z): 
                newr = z
            else: 
                raise ValueError('np.shape(data)[0] != len(z) ') 
        
        newr[igal2_]  = newr[igal1_] # assign 距离（也就红移）
        index[igal2_] = True    
        return newr, index 