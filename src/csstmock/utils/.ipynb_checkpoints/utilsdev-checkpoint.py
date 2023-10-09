import numpy as np
import logging 
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d

def fracndim(data, bin_list = None, mask = None):
    return fracndim2(data, bin_list = bin_list, mask = mask)

def fracndim2(data, bin_list = None, mask = None):
    '''
    outershape: 循环层数和每层的次数
    '''
    if bin_list is None: 
        return np.sum(mask)/np.prod(np.shape(mask))
    
    if isinstance(data, list): data = np.array(data).T
    
    outershape = [len(bins)-1 for bins in bin_list] 
    for iloop, d in enumerate(outershape): 
        if d == -1: 
            bin_new = np.sort( np.unique(data[:,iloop])  )
            bin_list[iloop]  = np.append(bin_new, bin_new[-1]+1) 
            
    hist2, edges = np.histogramdd(data, bins=bin_list)
    if mask is not None: hist1, edges = np.histogramdd(data[mask,:], bins=edges)
    if mask is not None: 
        fracndim = 1.0*hist1/hist2
    else: 
        fracndim = 1.0*hist2
    return fracndim, edges

def fracndim1(data, bin_list = None, mask = None):
    '''
    outershape: 循环层数和每层的次数
    '''
    if bin_list is None: 
        return np.sum(mask)/np.prod(np.shape(mask))
    
    if isinstance(data, list): data = np.array(data).T
    
    outershape = [len(bins)-1 for bins in bin_list] 
    for iloop, d in enumerate(outershape): 
        if d == -1: 
            bin_new = np.sort( np.unique(data[:,iloop])  )
            bin_list[iloop]  = np.append(bin_new, bin_new[-1]+1) 
            outershape[iloop]= len(bin_list[iloop]) - 1
    outer = np.zeros( outershape ) 
    loops = [0 for _ in outershape]
    count = 0 
    while loops[-1]!=outershape[-1]: 
        for iloop in range( len(outershape) ): 
            lo = bin_list[iloop][loops[iloop]]
            up = bin_list[iloop][loops[iloop] + 1] 
            indx_ = (data[:,iloop] >= lo)&(data[:,iloop] < up)
            if iloop == 0: index = indx_
            if iloop != 0: index = index&indx_
        print( loops )
    #--------------------------------------------
        if mask is not None: 
            outer[tuple(loops)] = 1.0*np.sum(index&mask)/np.sum(index)
        else: 
            outer[tuple(loops)] = 1.0*np.sum(index)
            
    #---------------------------------------------
    # control the loops
        loops[0]+=1
        count+=1
        for iloop in range( len(outershape) - 1 ): 
            if loops[iloop] == outershape[iloop] : 
                loops[iloop+1]+= 1
                loops[iloop]   = 0
        # print(count, loops, outershape)
        iloop=0
    #---------------------------------------------
    return outer, bin_list


class lftools(object):
    """
    Update:  calzmax2 - zmax[np.inan(zmax)] = 0 

    """
    def __init__(self, H0, Om0, precision=6): 
        self.cosmo   = FlatLambdaCDM(H0, Om0=Om0); 
        self.zgrid   = self.build_zgrid(zup = 15, zlo = 10**(-precision), precision = precision, Omega_M = Om0, H0 = H0)

    def __init__0(self, zgal, mag_app, kcorr, kcorr_func = None, 
                 skycov = None, weight = None, igrp = None, zgrp = None, rank = None, mag_cut = None,
                 Msun = 4.83, Omega_M = 0.315, H0 = 67.4, zmaxiter = 3, verbose = 1):
        '''
        
        '''
        self.verbose = verbose
        self.zgal    = zgal
        self.mag_app = mag_app
        self.kcorr   = kcorr
        self.kcorr_func = kcorr_func
        self.ngal    = len(zgal)
        self.cosmo   = FlatLambdaCDM(H0, Om0=Omega_M); 
        self.zgrid   = self.build_zgrid(precision = 6, Omega_M = Omega_M, H0 = Omega_M)
        self.zmaxgrid= self.build_zmaxgrid(kcorr_func=kcorr_func, zgrid = self.zgrid) 
        self.Msun    = Msun
        self.mode    = self.initialize(skycov = skycov, weight = weight, igrp = igrp, zgrp = zgrp, rank = rank, mag_cut = mag_cut, verbose = self.verbose) 
        self.loglum  = self.app2abs(zgal, mag_app, kcorr = kcorr, Msun = Msun, use_h = True, zgrid = self.zgrid, verbose = self.verbose)
        self.zmax    = self.calzmax(zgal, mag_app, mag_cut = mag_cut, zmaxgrid = self.zmaxgrid, zgrid = self.zgrid)
        # print(np.shape(self.zmax) )

    def app2abs(self, zgal, mag_app, kcorr, Msun = None, use_h = True):  
        '''
        Convet apparent magnitude to absolute magnitude to a given redshift. 
        mag_abs = mag_app - 5*np.log10(DL[Mpc]) - 25 - kcorr
        mag_abs - 5*log10(h) = mag_app - 5*np.log10(DL[Mpc/h]) - 25 - kcorr - 5*log10(h)

        mag_abs - Msun = -2.5*log10( L/Lsun )
        mag_abs - 5*log10(h) - Msun = -2.5*log10( L/Lsun ) - 5*log10(h)
                                    = -2.5*log10( L*h*h/Lsun )
        log10( L/(Lsun/h^2) ) = -0.4(mag_abs - 5*log10(h) - Msun) 

        Parameter: 
        ----------
        mag_app: array_like
            apparent magnitude of galaxies
        z: array_like
            redshift of galaxies
        kcorr: array_like, or constant
        Msun: array_like, or None
            if Msun is None (defaults), return magnitude
            if Msun is settled, e.g., 4.83, return luminosity
        use_h: bool
            if use_h is True (defults), the results of absolute magnitude (or luminosity) is Mabs - 5*log10(h) (or L[Lsun/h^2])
            if use_h is False, the results of absolute magnitude (or luminosity) is Mabs (or L[Lsun])
        Results: 
        ----------- 
        magnitude/log10(luminosity): array_like 
            use Msun and use_h to control the output
        '''
        zgrid   = self.zgrid # lftools.build_zgrid(precision = 4, Omega_M = Omega_M, H0 = H0)
        zgal    = np.atleast_1d(zgal)
        mag_app = np.atleast_1d(mag_app)
        kcorr = np.atleast_1d(kcorr)
        ngal    = len(zgal)
        DL      = zgrid.z2DL(zgal)   # luminosity distance, Unit: X/h [Mpc]， e.g., 100 [Mpc/h]
        if use_h == False: DL = DL/zgrid.h # 100 [Mpc/h] --> 142.86 Mpc if h = 0.7 
        DM   = 5*np.log10(DL) + 25 # DISTANCE MODULE
        mag_abs = mag_app - DM - kcorr
        if Msun is not None: 
            loglum = -0.4*(mag_abs - Msun)  
            return np.array(loglum)  # if use_h == True, Unit: Lsun/h^2 
        else: 
            return np.array(mag_abs) # if use_h == True, Mabs - 5*log10(h) 

    def calzmax(self, zgal, mag_app, mag_cut, kcorr_func = None):
        '''
            mag_abs = mag_app - 5*np.log10(DL(z)[Mpc]) - 25 - kcorr(z)
            mag_abs = mag_cut - 5*np.log10(DL(zmax)[Mpc]) - 25 - kcorr(zmax)
    
            ==> mag_cut - mag_app =  func(zmax) - func(z),   func(z) = kcorr(z) + 5*np.log10(DL(z)[Mpc])
        '''
        zgal    = np.atleast_1d(zgal)
        mag_app = np.atleast_1d(mag_app)
        zgrid    = self.zgrid # lftools.build_zgrid(precision = 4, Omega_M = Omega_M, H0 = H0)
        zmaxgrid = lftools.build_zmaxgrid(kcorr_func, zgrid = zgrid) 
        dmag = mag_cut - mag_app + zmaxgrid.func1(zgal)
        zmax = zmaxgrid.func2(dmag)
        return zmax

    def callf(self, zlo, zup, zgal, zmax, loglum, wht, loglum_bin = np.arange(6, 13+0.1, 0.1), bootstrap = 200, skycov = 41252.96): 
        indx_z   = (zlo < zgal)&(zgal <= zup)&(zmax>=zgal)
        #if wht is None: wht = loglum*0+1.0
        loglum = loglum[indx_z]
        zgal = zgal[indx_z]
        zmax = zmax[indx_z]; zmax[zmax > zup] = zup 
        wht  = wht[indx_z]
        vmax = self.zgrid.z2VC(zmax) - self.zgrid.z2VC(zlo)
        weight = 41252.96/skycov/vmax; # 1/vmax
        weight[np.isinf(weight)|np.isnan(weight)] = 0
        
        loglum_x = 0.5*(loglum_bin[1:]+loglum_bin[:-1])  
        numden_y = self.bootstrap(loglum, loglum_bin, weight = weight*wht, bootstrap = bootstrap )
        lf0 = numden_y[:,0] 
        ave = np.mean(numden_y, axis = 1 )
        std = np.std( numden_y, ddof = 1, axis = 1)
        #-- print(np.shape(lf0), np.shape(ave), np.shape(std) ) 
        logelo = np.abs(   np.log10(ave) - np.log10( np.abs(ave-std) ) )
        logeup = np.abs( - np.log10(ave) + np.log10( np.abs(ave+std) ) )
        hist,edge = np.histogram( loglum, loglum_bin)
        res = np.vstack([loglum_x, np.log10(ave), logelo, logeup, hist, lf0]).T
        return res 

    def bootstrap(self, loglum, loglum_bin, weight, bootstrap = 1):
        '''
        bootstrap1 is used to calculate the luminosity function and the uncertainty, 
        which do the boostrap by random sampling of galaxies. 
        '''
        loglumc    = 0.5*( loglum_bin[:-1] + loglum_bin[1:] )
        numden     = np.zeros( (len(loglumc), bootstrap)) 
        weight[np.isnan(weight)|np.isinf(weight)] = 0
        ngal = np.shape(loglum)[0]; 
        nwht = np.zeros( (ngal, bootstrap) ) # number weight 
        rand = np.random.randint(0, ngal, size = (ngal, bootstrap))
        #
        #------ bootstrap  
        #
        for ii in range(bootstrap): 
            uniq, counts = np.unique(rand[:, ii], return_counts = True) 
            nwht[uniq, ii]  = counts
            if ii == 0: nwht[:, ii] = 1
            
        for ii in range(len(numden)): 
            indx = (loglum_bin[ii] < loglum)&(loglum < loglum_bin[ii+1]) 
            if np.sum(indx) == 0: 
                numden[ii, :] = np.nan
            else: 
                numden[ii, :] = np.sum(nwht[indx, :]*weight[indx].reshape(-1, 1), axis = 0) 
        return numden/np.abs(loglum_bin[1]-loglum_bin[0]) 


    def initialize(self, skycov = None, weight = None, igrp = None, zgrp = None, rank = None, mag_cut = None, verbose = 1): 
        '''
        
        '''
        self.weight = weight
        self.skycov = skycov
        self.zgrp = zgrp
        self.igrp = igrp
        self.rank = rank
        self.mag_cut = mag_cut
        if self.mag_cut is None: self.mag_cut = np.max(self.mag_app) 
        if self.weight  is None: 
            self.weight = 1
            #print('Weight of galaxies is not specific, set it to %s. '% self.weight)
        if (self.rank is None): 
            if verbose != 0: print('lftools: Calculate luminosity function') 
            # print('Cannot calculate conditional luminosity function. (id, redshift, c/s) of groups is not provided.')
            self.mode = 'lf' 
        else: 
            # self.skycov = 4*np.pi*(180/np.pi)**2
            #print('sky coverage of galaxies is not specific, set it to %s deg^2. '%self.skycov) 
            if verbose != 0: print('lftools: Calculate conditional luminosity function') 
            self.mode = 'clf'
        return self.mode 
#    @classmethod
#    def data(cls, data, mag_cut, kcorr_func, skycov, Msun):        
#        data[]
#        return cls
    def check(self, index = None): 
        ngal = len(self.zgal); 
        if index is None: index = np.arange(ngal)
        nsel = len(self.zgal[index])
        zmax = np.max(self.zgal[index]); zmin = np.min(self.zgal[index]) 
        print('#--- select %s/%s galaxies at %s <= z <= %s'%(nsel, ngal, zmax, zmin) ) 
        print('#--- luminosity: %s <= log(L) <= %s'% (np.min(self.loglum), np.max(self.loglum)))
        
    @classmethod 
    def build_zgrid(cls, zup = 15, zlo = 0, precision = 3, Omega_M = 0.315, H0 = 67.4, scale = 'log'):
        '''
        Build the gird of redshift according to the precision. 
        Generate the corrsponding distance and volume as a reference. 
        '''
        from scipy.interpolate import interp1d
        cosmo   = FlatLambdaCDM(H0, Om0=Omega_M); 
        #---------- redshift grid --> DC, DL (Mpc/h)
        cls.Omega_M = Omega_M
        cls.H0    = H0
        cls.cosmo = cosmo 
        cls.zup   = zup
        cls.zlo   = zlo
        cls.h     = cosmo.H0.value/100.0 # 67.4 km / (Mpc s)
        cls.precision = precision
        if scale == 'log': 
            cls.log1pz= np.arange(np.log10(1+zlo), np.log10(1+zup) + 10**(-precision), 10**(-precision))
            cls.zref  = 10**cls.log1pz - 1
        elif scale == 'linear': 
            cls.zref  = np.arange(zlo, zup+ 10**(-precision), 10**(-precision)) 
            cls.log1pz= np.log10(cls.zref + 1)
        #--------
        #log1pz = np.arange(np.log10(1+zlo), np.log10(1+0.1) + 10**(-precision), 10**(-precision))
        #zref1  = 10**log1pz - 1
        #zref2  = np.arange(0.1+10**(-precision), zup+10**(-precision), 10**(-precision) )        
        #cls.zref  = np.hstack(zref1, zref2)
        #--------
        
        # np.arange(zlo, zup+10**(-precision), 10**(-precision)); # zref = [0., 0.001, 0.002, ..]
        # 142.86 [Mpc] * 0.7 = 100 [Mpc/h]
        cls.DCref = np.array(cosmo.comoving_distance(cls.zref).value)*cls.h; # comoving distance, Unit: Mpc/h
        cls.DLref =  (1+cls.zref)*cls.DCref    # luminosity distance, Unit: Mpc/h
        cls.VCref = 4.0/3.0*np.pi*cls.DCref**3 # comoving   volume, Mpc^3/h^3
        cls.z2DC  = interp1d(cls.zref, cls.DCref, kind = 'linear', bounds_error=False, fill_value="extrapolate")
        cls.z2DL  = interp1d(cls.zref, cls.DLref, kind = 'linear', bounds_error=False, fill_value="extrapolate")
        cls.z2VC  = interp1d(cls.zref, cls.VCref, kind = 'linear', bounds_error=False, fill_value="extrapolate")
        cls.DC2z  = interp1d(cls.DCref, cls.zref, kind = 'linear', bounds_error=False, fill_value="extrapolate")
        cls.DL2z  = interp1d(cls.DLref, cls.zref, kind = 'linear', bounds_error=False, fill_value="extrapolate")
        cls.VC2z  = interp1d(cls.VCref, cls.zref, kind = 'linear', bounds_error=False, fill_value="extrapolate")
        return cls

    @classmethod 
    def build_ngrid1(cls, z, i = None, precision = 4, zgrid = None, Omega_M = 0.315, H0 = 67.4, **kwargs): 
        '''
        Input group redshift and group id, return the number as a function of redshift. 
        Use zref as redshift bin to count. 
        '''
        if zgrid is None: zgrid = lftools.build_zgrid(precision = precision, Omega_M = Omega_M, H0 = H0) 
        if i is not None: 
            uniq, index = np.unique(i, return_index=True)
            z  = z[index]; i = i[index]
        hist, edge  = np.histogram(z, bins = zgrid.zref)
        cumhist = np.cumsum(hist)
        cumhist = np.append(cumhist, cumhist[-1])
        z2ngrid = interp1d(edge, cumhist, kind = 'previous', bounds_error=False, fill_value="extrapolate")
        return z2ngrid
    @classmethod
    def bootstrap_ngal(cls, ngal, bootstrap = 1, igrp = None): 
        '''
        #---- 对ngal个星系进行随机抽样， 或者对星系团以及其中的星系进行随机抽样
        #---- 返回：numgal_rand, 星系在每一次抽样中出现的次数
        '''
        numgal_rand = np.zeros((ngal, bootstrap)) # 星系在每一次抽样中出现的次数
        if igrp is None: 
            idxgal     = np.random.randint(0, ngal, (ngal, bootstrap))
            for ii in range(bootstrap): 
                uniq, counts = np.unique(idxgal[:, ii], return_counts = True) 
                numgal_rand[uniq,ii] = counts
                if ii == 0: numgal_rand[:,0] = 1 
        else: 
            if ngal != np.shape(igrp)[0]: logging.error( 'input igrp did not have same size with number of galaxies') 
            uniq, inverse, rich = np.unique(igrp, return_inverse = True, return_counts = True)
            ngrp = len(uniq) # 星系群的总数
            #print('星系群编号', uniq, '， 共%s个星系群'%ngrp)
            #print('成员星系数', rich) 
            numgrp_rand = np.zeros((ngrp, bootstrap)) # 星系团在每一次抽样中出现的次数
            #--- 开始随机抽样
            idxgrp     = np.random.randint(0, ngrp, (ngrp, bootstrap))
            for ii in range(bootstrap): 
                idxgrp_uniq, counts = np.unique(idxgrp[:,ii], return_counts = True)
                numgrp_rand[idxgrp_uniq,ii] = counts
                numgal_rand[:,ii] = numgrp_rand[inverse,ii]
                # print('星系团在抽样中出现的次数:', numgrp_rand[:,ii], '星系在抽样中出现的次数:', numgal_rand[:,ii] )
                if ii == 0: numgal_rand[:,0] = 1 
        return numgal_rand 
    @classmethod 
    def build_ngrid2(cls, z, i = None, precision = 10):
        '''
        Input group redshift and group id, return the number as a function of redshift. 
        Use input redshift to calculate the cumulative counts. 
        '''
        if i is not None: 
            i, index = np.unique(i, return_index=True)
            z = z[index]
        z   = np.sort(z); z = np.append(z, z[-1]+10**(-precision) ) 
        num = np.arange( len(z) ) + 1
        z2ngrid= interp1d(z-10**(-precision), num, kind = 'next', bounds_error=False, fill_value=(num[0], num[-1]) )
        return z2ngrid 

    @classmethod
    def build_ngrid(cls, *args, **kwargs):
        '''
        mag_app = -2.5log10(L) + 5log10(DL1)
        mag_cut = -2.5log10(L) + 5log10(DL2) 

        '''
        return lftools.build_ngrid2(*args, **kwargs)


    @classmethod 
    def build_zmaxgrid(cls, kcorr_func = None, zgrid = None, Omega_M = 0.315, H0 = 67.4,): 
        if kcorr_func is None: 
            dmag   = 5*np.log10(zgrid.DLref) 
        else:
            dmag   = kcorr_func(zgrid.zref) + 5*np.log10(zgrid.DLref) 
        cls.dmag    = dmag
        cls.func1  = interp1d(zgrid.zref, dmag, kind = 'linear', bounds_error=False, fill_value="extrapolate")
        cls.func2  = interp1d(dmag, zgrid.zref, kind = 'linear', bounds_error=False, fill_value="extrapolate")        
        return cls 

    
    @classmethod 
    def calzmax2(cls, zgal, mag_app, mag_cut, kcorr_func = None, zgrid = None, zmaxgrid = None, Omega_M = 0.315, H0 = 67.4, verbose = 1, **kwargs):
        '''
            mag_abs = mag_app - 5*np.log10(DL(z)[Mpc]) - 25 - kcorr(z)
            mag_abs = mag_cut - 5*np.log10(DL(zmax)[Mpc]) - 25 - kcorr(zmax)
    
            ==> mag_cut - mag_app =  func(zmax) - func(z),   func(z) = kcorr(z) + 5*np.log10(DL(z)[Mpc])
        '''
        if verbose != 0: print('lftools: Calculate zmax of galaxies (v2)')
        zgal    = np.atleast_1d(zgal)
        mag_app = np.atleast_1d(mag_app)
        if zgrid    is None: zgrid    = lftools.build_zgrid(precision = 6, Omega_M = Omega_M, H0 = H0) 
        if zmaxgrid is None: zmaxgrid = lftools.build_zmaxgrid(kcorr_func, zgrid = zgrid) 
        dmag = mag_cut - mag_app + zmaxgrid.func1(zgal)
        zmax = zmaxgrid.func2(dmag)
        zmax[np.isnan(zmax)] = 0 
        return zmax

    @classmethod
    def calzmax1(cls, zgal, mag_app, mag_cut, kcorr_func = None, iter = 0, zgrid = None, Omega_M = 0.315, H0 = 67.4, verbose = 1):
        '''
        mag_app = -2.5log10(L) + 5log10(DL1)
        mag_cut = -2.5log10(L) + 5log10(DL2) 

        '''
        if verbose != 0: print('lftools: Calculate zmax of galaxies (v1)')

        if zgrid is None: zgrid   = lftools.build_zgrid(precision = 6, Omega_M = Omega_M, H0 = H0) 
        zgal    = np.atleast_1d(zgal)
        mag_app = np.atleast_1d(mag_app)
        DL1  = zgrid.z2DL(zgal)  # luminosity distance, Unit: Mpc/h
        DL2  = mag_cut - mag_app + 5*np.log10(DL1)
        DL2  = 10**(DL2/5.0)
        z1   = zgrid.DL2z(DL2)
        for  ii in range(iter):
            DL2  = mag_cut - mag_app + 5*np.log10(DL1) + kcorr_func(zgal) - kcorr_func(z1)  
            DL2  = 10**(DL2/5.0) 
            z1   = zgrid.DL2z(DL2) 
        return z1
    
    @classmethod
    def calzmax0(cls, *args, **kwargs):
        '''
        mag_app = -2.5log10(L) + 5log10(DL1)
        mag_cut = -2.5log10(L) + 5log10(DL2) 

        '''
        return lftools.calzmax2(*args, **kwargs)

    
    @classmethod 
    def app2abs1(cls, zgal, mag_app, kcorr=0, Msun = None, use_h = True, zgrid = None, Omega_M = 0.315, H0 = 67.4, verbose = 1):  
        '''
        Convet apparent magnitude to absolute magnitude to a given redshift. 
        mag_abs = mag_app - 5*np.log10(DL[Mpc]) - 25 - kcorr
        mag_abs - 5*log10(h) = mag_app - 5*np.log10(DL[Mpc/h]) - 25 - kcorr - 5*log10(h)

        mag_abs - Msun = -2.5*log10( L/Lsun )
        mag_abs - 5*log10(h) - Msun = -2.5*log10( L/Lsun ) - 5*log10(h)
                                    = -2.5*log10( L*h*h/Lsun )
        log10( L/(Lsun/h^2) ) = -0.4(mag_abs - 5*log10(h) - Msun) 

        Parameter: 
        ----------
        mag_app: array_like
            apparent magnitude of galaxies
        z: array_like
            redshift of galaxies
        kcorr: array_like, or constant
        Msun: array_like, or None
            if Msun is None (defaults), return magnitude
            if Msun is settled, e.g., 4.83, return luminosity
        use_h: bool
            if use_h is True (defults), the results of absolute magnitude (or luminosity) is Mabs - 5*log10(h) (or L[Lsun/h^2])
            if use_h is False, the results of absolute magnitude (or luminosity) is Mabs (or L[Lsun])
        Results: 
        ----------- 
        magnitude/log10(luminosity): array_like 
            use Msun and use_h to control the output
        '''
        if verbose != 0: print('lftools: Convet apparent magnitude to absolute magnitude (or luminosity)')
        if zgrid is None: zgrid = lftools.build_zgrid(precision = 4, Omega_M = Omega_M, H0 = H0)
        zgal    = np.at_least1d(zgal )
        mag_app = np.at_least1d(mag_app) 
        kcorr   = np.at_least1d(kcorr) 
        ngal    = len(zgal)
        DL      = zgrid.z2DL(zgal)   # luminosity distance, Unit: X/h [Mpc]， e.g., 100 [Mpc/h]
        if use_h == False: DL = DL/zgrid.h # 100 [Mpc/h] --> 142.86 Mpc if h = 0.7 
        DM   = 5*np.log10(DL) + 25 # DISTANCE MODULE
        mag_abs = mag_app - DM.value - kcorr
        if Msun is not None: 
            loglum = -0.4*(mag_abs - Msun)  
            return np.array(loglum.value)  # if use_h == True, Unit: Lsun/h^2 
        else: 
            return np.array(mag_abs.value) # if use_h == True, Mabs - 5*log10(h) 

    @classmethod
    def bootstrap1(cls, loglum, loglum_bin, weight, bootstrap = 1, verbose = 1):
        '''
        bootstrap1 is used to calculate the luminosity function and the uncertainty, 
        which do the boostrap by random sampling of galaxies. 
        '''
        if verbose != 0: print('lftools: bootstrap to calculate the luminosity function and the uncertainty') 
        loglumc    = 0.5*( loglum_bin[:-1] + loglum_bin[1:] )
        numden     = np.zeros( (len(loglumc), bootstrap)) 
        weight[np.isnan(weight)|np.isinf(weight)] = 0
        ngal = np.shape(loglum)[0]; 
        nwht = np.zeros( (ngal, bootstrap) ) # number weight 
        rand = np.random.randint(0, ngal, size = (ngal, bootstrap))
        #
        #------ bootstrap  
        #
        for ii in range(bootstrap): 
            uniq, counts = np.unique(rand[:, ii], return_counts = True) 
            nwht[uniq, ii]  = counts
            if ii == 0: nwht[:, ii] = 1
            
        for ii in range(len(numden)): 
            indx = (loglum_bin[ii] < loglum)&(loglum < loglum_bin[ii+1]) 
            if np.sum(indx) == 0: 
                numden[ii, :] = np.nan
            else: 
                numden[ii, :] = np.sum(nwht[indx, :]*weight[indx].reshape(-1, 1), axis = 0) 
        return loglumc, numden/np.abs(loglum_bin[1]-loglum_bin[0])  
    
    @classmethod 
    def bootstrap3(cls, loglum, bins, weight, rank, igrp = None, bootstrap = 1, verbose = 1):
        '''
        bootstrap2 is used to calculate the conditional luminosity function and the uncertainty, 
        which do the boostrap by random sampling of groups. 
        '''
        if verbose != 0: print('lftools: bootstrap to calculate the uncertainty of (conditional) luminosity function') 
        xc     = 0.5*( bins[:-1] + bins[1:] )
        ngal   = len(loglum) # 输入星系的总数
        weight[np.isnan(weight)|np.isinf(weight)] = 0 # 预处理，如权重为异常值， 则不考虑后续的记数
        #------------------------------------------------------------------------------------------------------
        numgal_rand = lftools.bootstrap_ngal(ngal, bootstrap = bootstrap, igrp = igrp) 
        
        #------------------------------------------------------------------------------------------------------
        #
        #--- bootstrap 
        #
        yboots0 = np.zeros( (len(xc), bootstrap)) # 用于存放central结果
        yboots1 = np.zeros( (len(xc), bootstrap)) # 用于存放satellite结果

        for ii in range(len(xc)): 
            # central CLF
            indx   = (bins[ii] < x)&(x < bins[ii+1])&(rank == 0)
            if np.sum(indx) == 0: 
                yboots0[ii, :] = 0
            else: 
                yboots0[ii, :] = np.sum(numgal_rand[indx,:]*weight[indx].reshape(-1, 1), axis = 0) 
            # satellite CLF
            indx   = (bins[ii] < x)&(x < bins[ii+1])&(rank != 0)
            if np.sum(indx) == 0: 
                yboots1[ii, :] = 0
            else: 
                yboots1[ii, :] = np.sum(numgal_rand[indx,:]*weight[indx].reshape(-1, 1), axis = 0) 
                
        return xc, yboots0/np.abs(bins[1]-bins[0]), yboots1/np.abs(bins[1]-bins[0]) 
    
    @classmethod 
    def bootstrap2(cls, x, bins, weight, igrp = None, bootstrap = 1, verbose = 1):
        '''
        bootstrap2 is used to calculate the conditional luminosity function and the uncertainty, 
        which do the boostrap by random sampling of groups. 
        '''
        if verbose != 0: print('lftools: bootstrap to calculate the uncertainty of (conditional) luminosity function') 
        xc     = 0.5*( bins[:-1] + bins[1:] )
        nx     = np.shape(x)[0]; 
        ngal = len(x) # 输入星系的总数
        weight[np.isnan(weight)|np.isinf(weight)] = 0 # 预处理，如权重为异常值， 则不考虑后续的记数
        #------------------------------------------------------------------------------------------------------
        numgal_rand = lftools.bootstrap_ngal(ngal, bootstrap = bootstrap, igrp = igrp) 
        #------------------------------------------------------------------------------------------------------
        #
        #--- bootstrap 
        #
        yboots     = np.zeros( (len(xc), bootstrap)) # 用于存放结果
        for ii in range(len(xc)): 
            indx   = (bins[ii] < x)&(x < bins[ii+1])
            if np.sum(indx) == 0: 
                yboots[ii, :] = 0
            else: 
                yboots[ii, :] = np.sum(numgal_rand[indx,:]*weight[indx].reshape(-1, 1), axis = 0) 
                # print(('[%.2f, %.2f]')%(bins[ii], bins[ii+1]), ':', 1/weight[indx][:3]) 
        return xc, yboots/np.abs(bins[1]-bins[0]) 
        
        
    def run(self, zlo, zup, loglum_bin = np.arange(8, 13+0.1, 0.1), index = None, bootstrap = 1, mode = None): 
        if index is None: index = np.ones(self.ngal).astype(bool)
        if self.mode == 'lf':  
            indx_z   = (zlo < self.zgal)&(self.zgal < zup)&(self.zmax>=self.zgal)
            
        if self.mode == 'clf': 
            # rank: 必要，区分中心星系和卫星星系
            # zgrp: 星系群的红移，计算星系团的累积数目和红移的关系numgrp(<z)。有，直接使用；无，利用中心星系的红移。 
            # igrp: 星系群的编号，1.用于找到在红移范围内的星系团对应的成员星系; 2. 用于boostrap是对星系团随机抽样
            #          有: 选择在红移范围内的星系群， 若成员星系超出了红移的范围的，也纳入其中
            #          无: 只能选择在红移范围内的星系，只能直接对星系随机抽样
            # 
            #--- 构建红移z与小于z的星系团数目之间的关系 numgrp_lt_zmax = z2ngrid(zmax)
            #
            if self.zgrp is not None: 
                # 如果group的红移已经提供了， 直接使用
                indx_z  = (zlo < self.zgrp)&(self.zgrp < zup)&(self.zmax>=self.zgal)
                zgrp    = self.zgrp[indx_z&(self.rank ==0)]
            else: 
                # 如果group的红移没有提供了， 则使用中心星系的红移作为group的红移 
                indx_z  = (self.zmax>=self.zgal)
                zgrp    =  self.zgal[indx_z&(self.rank ==0)]
            # 星系群的红移，计算星系团的累积数目和红移的关系numgrp(<z)
            z2ngrid  = lftools.build_ngrid( zgrp[(zlo < zgrp)&(zgrp < zup)])
            #
            #--- 选择出用于计算CLF的星系样本 
            #
            if self.igrp is not None:  
                # 
                # 如果提供了groupid, 则选择满足红移范围的的group。 若成员星系超出了红移的范围的，也纳入其中
                # 
                if self.zgrp is not None: 
                    igrp     = self.igrp[(self.rank ==0)&(zlo < self.zgrp)&(self.zgrp < zup)&(self.zmax>=self.zgal)]
                if self.zgrp is None:     
                    igrp     = self.igrp[(self.rank ==0)&(zlo < self.zgal)&(self.zgal < zup)&(self.zmax>=self.zgal)]
                uniq_igrp= np.unique( igrp ) 
                indx_z   = np.isin(self.igrp, uniq_igrp) 
                indx_z   = indx_z & (self.zmax>=self.zgal)
            else: 
                if self.zgrp is not None: indx_z   = (zlo < self.zgrp)&(self.zgrp < zup)&(self.zmax>=self.zgal) #
                if self.zgrp is None:     indx_z   = (zlo < self.zgal)&(self.zgal < zup)&(self.zmax>=self.zgal) #

        index   = index & indx_z
        zmax    = self.zmax[index]
        loglum  = self.loglum[index]
        
        if self.weight is not None: weight = self.weight[index]
        
        if self.mode == 'lf': 
            #print('Running to calculate the luminosity function')
            allsky= 4*np.pi*(180/np.pi)**2; 
            omega = self.skycov/allsky
            zmax[zmax > zup]  = zup; 
            vmax  = omega*(self.z2VC(zmax) - self.z2VC(zlo) ) 
            loglumc, numden = lftools.bootstrap1(loglum, loglum_bin, weight = 1.0/weight/vmax, bootstrap = bootstrap, verbose = self.verbose)

        if self.mode == 'clf': 
            if self.igrp is not None: 
                igrp = self.igrp[index]
                #zgal = self.zgal[index]; # print( np.max(zgal), np.min(zgal) )
            else: 
                igrp = None 
            #print('Running to calculate the conditional luminosity function')
            zmax[zmax > zup]  = zup; 
            numgrp_lt_zmax  = z2ngrid(zmax) - z2ngrid(zlo); 
            
            loglumc, numden = lftools.bootstrap2(loglum, loglum_bin, weight = 1.0/weight/numgrp_lt_zmax, igrp = igrp, bootstrap=bootstrap, verbose = self.verbose) 
            #loglumc, numden = lftools.bootstrap1(loglum, loglum_bin, weight = 1.0/weight/numgrp_lt_zmax, bootstrap=bootstrap, verbose = self.verbose) 

        #-------- summary bootstrap results
        lfdata   = np.zeros( (len(loglumc), 6) ) + np.nan 
        for ii in range(len(loglumc)):
            numden_ = numden[ii, :] 
            numden_[ np.isnan(numden_)|np.isinf(numden_) ] = 0
            numden_ = numden_[numden_ != 0] 
            ave = np.mean(numden_ ) 
            err = np.std( numden_, ddof = 1, axis = 1) 
            lfdata[ii,0] = loglumc[ii]
            lfdata[ii,1] = np.log10(ave)
            lfdata[ii,2] = np.log10(err)
            lfdata[ii,3] = np.abs(   np.log10(ave) - np.log10( np.abs(ave-err) ) )
            lfdata[ii,4] =  - np.log10(ave) + np.log10( np.abs(ave+err) ) 
            
        hist,edge= np.histogram(  loglum, loglum_bin)
        lfdata[:,5] = hist

        # numden = np.ma.masked_array(numden, mask = np.isnan(numden) ) 
        # ave = np.mean(numden, axis = 1) 
        # err_lo = np.log10(ave) - np.log10(ave-err)
        # err_hi = np.log10(ave+err) - np.log10(ave) 
        # lfdata = np.vstack([loglumc, np.log10(ave), np.log10(err), err_lo, err_hi]).T 
        return lfdata, numden

    


import os
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
from astropy.io  import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord 

class rgbimg(object): 
    '''
    The requirements for the images are to 
        1. present an homogeneous background, 
        2. distinguish bright to faint objects
    '''
    def __init__(self, imgs): 
        self.imgs = imgs
        
    @classmethod 
    def from_cube(cls, cube): 
        pass
    
    @classmethod 
    def from_fits(cls, filenames): 
        if len(filenames) != 3: print('inputs != 3') 
        imgs, hdrs = rgbimg.readimgs(filenames, aslist = True)
        return cls(imgs)
    #--- end __init__ end ---# 
    
    @classmethod
    def readimg_worker(cls, filename, ihdu = 0): 
        '''
        write for mutliprocessing
        read image array from a file (jpg or fits). 
        '''
        split = filename.split('.') 
        if split[-1] not in ['jpg', 'jpeg', 'png', 'fits']: 
            logging.error(u" Unsupported filename: %s"% filename)
        if split[-1] in ['jpg', 'jpeg', 'png']: 
            # import imageio
            import imageio.v2 as imageio
            try:
                img = imageio.imread(filename) 
                img = img/255.0
                img = img[::-1]
                #print(img.shape)
                #img = np.transpose(img,(1,2,0))
                hdr = None 
            except: 
                print('rm -rf %s'%filename)
                img = None; hdr = None
        if split[-1]=='fits': 
            from astropy.io import fits
            try:
                hdu = fits.open(filename)
                img = np.array(hdu[ihdu].data)
                hdr = hdu[ihdu].header
                hdu.close()
            except: 
                print('rm -rf %s'%filename)
                img = None; hdr = None
        return [img, hdr]
    
    @classmethod
    def readimgs(cls, inputs, ihdu = 0, aslist = False): 
        '''
        load all the images and convert them into a tensor (indx, xpix, ypix, 1).  
        '''
        inputs = np.array(inputs)
        imgs = []; hdrs = []
        for inputfile in inputs: 
            img_, hdr_ = rgbimg.readimg_worker(inputfile, ihdu = ihdu) 
            imgs.append(img_)
            hdrs.append(hdr_)
            
        # print('### Finish reading %i files'%(len(inputs)))
        if aslist == False: 
            if(len(inputs)!=1): imgs = np.stack(imgs, axis = 0)
            if(len(inputs)==1): imgs = imgs[0][np.newaxis,:, :]
            print('###-----------------------------------------------------------------')
            print('### The shape of image array = ', imgs.shape)
            return imgs, hdrs
        else:
            print('###-----------------------------------------------------------------')
            print('### return as list')
            return imgs, hdrs
    @classmethod
    def forward(cls, x, Q = 10.0, minimum = 0.0, stretch = 0.5 ): 
        return 1.0/Q*np.arcsinh( Q*(x-minimum)/stretch )
    @classmethod
    def inverse(cls, x, Q = 10.0, minimum = 0.0, stretch = 0.5): 
        return np.sinh(x*Q)*stretch/Q + minimum
    @classmethod 
    def easy_bkg(cls, data): 
        from astropy.stats import sigma_clipped_stats, SigmaClip
        from photutils.segmentation import detect_threshold, detect_sources
        from photutils.utils import circular_footprint
        sigma_clip  = SigmaClip(sigma=3.0, maxiters=10)
        threshold   = detect_threshold(data, nsigma=2.0, sigma_clip=sigma_clip)
        segment_img = detect_sources(data, threshold, npixels=10)
        footprint = circular_footprint(radius=10)
        mask = segment_img.make_source_mask(footprint=footprint)
        mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
        return median, std
    
    def run(self):
        from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog
        image_r = self.imgs[0]; bkg_r, std_r = rgbimg.easy_bkg(image_r); 
        image_g = self.imgs[1]; bkg_g, std_g = rgbimg.easy_bkg(image_g); 
        image_b = self.imgs[2]; bkg_b, std_b = rgbimg.easy_bkg(image_b); 
        r = image_r - bkg_r#/std_r
        g = image_g - bkg_g#/std_g
        b = image_b - bkg_b#/std_b
        # 4.9,5.7,7.8 
        from astropy.visualization import SqrtStretch
        from astropy.visualization import ZScaleInterval
        from astropy.visualization import make_lupton_rgb
        lo_val, up_val = np.percentile(np.hstack((r.flatten(), g.flatten(), b.flatten())), (0.05, 99.5)) 
        image = make_lupton_rgb(0.7*r, 1.1*g, 1.8*b, stretch=0.05, Q=8, minimum = -0.005, filename = None )        
        image = image/255
        return image

class sex(object): 
    def __init__(self, img, hdr = None, verbose = 0): 
        import photutils
        # print(photutils.__version___) # == 1.3.0 ) 
        npixels = 10
        self.img       = img
        self.hdr       = hdr
        self.verbose   = verbose
        # self.mask      = detect.make_source_mask(self.img, nsigma=1, npixels=npixels, dilate_size=10)
        from astropy.convolution import convolve
        from astropy.convolution import Gaussian2DKernel
        from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog
        
        self.threshold = detect_threshold(self.img, nsigma=0.5) #, mask=self.mask)
        sigma          = 1.5; npixels = 10
        kernel         = Gaussian2DKernel(sigma, x_size=3, y_size=3); 
        convolved_data = convolve(self.img, kernel) # kernel.normalize()
        self.seg       = detect_sources(self.img, self.threshold, npixels=npixels) # , kernel=kernel)
        if self.seg is not None: 
            self.seg       = deblend_sources(self.img, self.seg, npixels=npixels, nlevels=32, contrast=0.01, progress_bar=False )
            self.tab       = SourceCatalog(  self.img, self.seg).to_table()
            self.tab       = self.tab[np.argsort(-self.tab['segment_flux']) ]
            if self.hdr is None:  
                from astropy import wcs
                wcs2d           = wcs.WCS(self.hdr)
                ra, dec = wcs2d.wcs_pix2world( self.tab['xcentroid'], self.tab['ycentroid'], 0) 
                self.tab['ra']  = ra 
                self.tab['dec'] = dec
            if(self.verbose!=0): print('Eesy Source Extract get %s potential sources'%len(self.tab) )
        else: 
            self.tab = None

    def return_mainobj(self): 
        if self.tab is None: 
            segm = self.img.copy().astype(bool) 
            segm = 1
            bkg  = ~segm
        elif  len(self.tab)  == 1: 
            segm     = self.seg.copy()
            segm = segm.data.astype(bool)
            bkg  = ~self.seg.data.astype(bool) 
        elif len(self.tab) > 1: 
            x0  = 0.5*np.shape(self.img)[0]; y0 = 0.5*np.shape(self.img)[1]
            x   = np.atleast_1d(self.tab['xcentroid'])
            y   = np.atleast_1d(self.tab['ycentroid'])
            pix_dist = (x - x0)**2 + (y - y0)**2
            ic = np.argmin(pix_dist.value) 
            t = self.tab.copy(); t.remove_row(ic)
            segm = self.seg.copy()  
            segm.remove_label(label = t['label'])
            segm = segm.data.astype(bool)
            bkg  = ~self.seg.data.astype(bool) 
        return  segm, bkg    
 
    def return_stamps(self): 
        from astropy.nddata import Cutout2D
        from astropy.io import fits 
        # positions = [(ix, iy) for ix, iy in zip(self.tab['ra'],self.tab['dec'] )]
        positions = [(ix, iy) for ix, iy in zip(self.tab['xcentroid'],self.tab['ycentroid'] )]
        bbox_ysize = self.tab['bbox_xmax'] - self.tab['bbox_xmin']
        bbox_xsize = self.tab['bbox_ymax'] - self.tab['bbox_ymin'] 
        sizes = [(xsize, ysize) for xsize, ysize in zip(bbox_xsize, bbox_ysize )]
        iths       = range(len(positions))
        self.new_imgs = []; self.new_hdrs = [] 
        for ith, position, size in zip(iths, positions, sizes): 
            hdr_wcs2d  = wcs.WCS(self.hdr.copy()) 
            cutout     = Cutout2D(self.img, position, size, hdr_wcs2d)
            new_img    = cutout.data
            hdu        = fits.ImageHDU(new_img)
            hdu.header = self.hdr.copy()
            hdu.header.update( cutout.wcs.to_header())
            new_hdr    = hdu.header
            self.new_imgs.append( new_img )
            self.new_hdrs.append( new_hdr )
        return self.new_imgs, self.new_hdrs

    def return_Rectangles(self):  
        self.stamp_mask = np.array( self.img.copy()*0.0, dtype = bool )
        self.patches = []
        for ii in range(len(self.tab)): 
            w_ = wcs.WCS( self.new_hdrs[ii] )
            w  = wcs.WCS( self.hdr )
            ra_min, dec_min = w_.wcs_pix2world(0, 0, 0)
            ix_min, iy_min  = w.wcs_world2pix(ra_min, dec_min, 0)
            ix_min, iy_min  = np.around([ix_min, iy_min],0).astype(int) 
            xy     = (ix_min, iy_min); 
            width = self.new_imgs[ii].shape[1]; height = self.new_imgs[ii].shape[0]
            ix_max = ix_min + width
            iy_max = iy_min + height
            self.patches.append(Rectangle(xy, width, height, angle=0.0)) 
            self.stamp_mask[iy_min:iy_max, ix_min:ix_max] = True
	# p = PatchCollection(patches, edgecolors='r', facecolors='None', linewidth = 3 ) 
	# ax.add_collection(p) 
        return self.stamp_mask, self.patches

class mockimg(object):
    '''
    build mock image
    
    parameter 
    ---------
    a1:
    d1:
    stps: 
    
    '''
    def __init__(self, a = None, d = None, stps = None, mask = None, wcs = None, survey = 'desi', cache = './__cache__/', enlarge = 1.0, verbose = 0): 
        self.ngal = np.shape(a1)[0]
        self.cache= cache
        self.a = np.atleast_1d(a)
        self.d = np.atleast_1d(d)
        self.stps = stps
        self.mask = mask 
        self.wcs    = wcs 
        self.survey = survey
        self.enlarge= enlarge 
        if not os.path.exists(cache): os.makedirs(cache)
        
            
    @classmethod
    def readimg_worker(cls, filename, ihdu = 0): 
        '''
        write for mutliprocessing
        read image array from a file (jpg or fits). 
        '''
        split = filename.split('.') 
        if split[-1] not in ['jpg', 'jpeg', 'png', 'fits']: 
            logging.error(u" Unsupported filename: %s"% filename)
        if split[-1] in ['jpg', 'jpeg', 'png']: 
            # import imageio
            import imageio.v2 as imageio
            try:
                img = imageio.imread(filename) 
                img = img/255.0
                hdr = None 
            except: 
                print('rm -rf %s'%filename)
                img = None; hdr = None
        if split[-1]=='fits': 
            from astropy.io import fits
            try:
                hdu = fits.open(filename)
                img = np.array(hdu[ihdu].data)
                hdr = hdu[ihdu].header
                hdu.close()
            except: 
                print('rm -rf %s'%filename)
                img = None; hdr = None
        return [img, hdr]

    @classmethod
    def readimgs(cls, inputs, ihdu = 0, aslist = False ): 
        '''
        load all the images and convert them into a tensor (indx, xpix, ypix, 1).  
        '''
        inputs = np.array(inputs)
        imgs = []; hdrs = []
        for inputfile in inputs: 
            img_, hdr_ = mockimg.readimg_worker(inputfile, ihdu = ihdu) 
            imgs.append(img_)
            hdrs.append(hdr_)
            
        # print('### Finish reading %i files'%(len(inputs)))
        if aslist == False: 
            if(len(inputs)!=1): imgs = np.stack(imgs, axis = 0)
            if(len(inputs)==1): imgs = imgs[0][np.newaxis,:, :]
            print('###-----------------------------------------------------------------')
            print('### The shape of image array = ', imgs.shape)
            return imgs, hdrs
        else:
            print('###-----------------------------------------------------------------')
            print('### return as list')
            return imgs, hdrs


    @classmethod
    def ruv2xyz(cls, r, u, v):
        x = r*np.cos(v/180.0*np.pi)*np.cos(u/180.0*np.pi)
        y = r*np.cos(v/180.0*np.pi)*np.sin(u/180.0*np.pi)
        z = r*np.sin(v/180.0*np.pi)
        return x, y, z

    @classmethod
    def xyz2ruv(cls, x, y, z):  
        r = x**2 + y**2 + z**2
        u = np.arctan2(y,x)
        v = np.arctan2(z, np.sqrt(x**2+y**2))
        u = u*180/np.pi;   
        v = v*180/np.pi
        return r, u, v

    @classmethod
    def autowcs(cls, a, d, pixscale, rho, enlarge = 1.01): 
        a_ = np.array([np.min(a), np.max(a)])
        d_ = np.array([np.min(d), np.max(d)])
        x, y, z = mockimg.ruv2xyz(1, a_, d_)
        x0 = np.mean(x)
        y0 = np.mean(y)
        z0 = np.mean(z)
        _,a0,d0 = mockimg.xyz2ruv(x0,y0,z0)
        obj0 = SkyCoord(a0, d0, unit='deg')
        obj  = SkyCoord(a,  d, unit='deg')
        sep  = obj0.separation(obj).degree
        size = 2*np.max(sep)*enlarge
        from astropy.wcs import WCS
        alpha_center = a0; NAXIS1 = int(size/pixscale); xpix_center  = 0.5 + (NAXIS1)*0.5;  
        delta_center = d0; NAXIS2 = int(size/pixscale); ypix_center  = 0.5 + (NAXIS2)*0.5;
        w = WCS(naxis=2);  
        w.wcs.cd    = [[-pixscale*np.cos(rho), pixscale*np.sin(rho)], 
                       [ pixscale*np.sin(rho), pixscale*np.cos(rho)]]
        w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        w.wcs.crval = [alpha_center, delta_center] 
        w.wcs.crpix = [ xpix_center,  ypix_center]
        w.pixel_shape = [NAXIS1, NAXIS2]
        print('#--- autowcs by input coordinates') 
        print()
        # print('#--- pixscale = %f degree/pixel'%pixscale )
        print(w) 
        return w

    @classmethod
    def isinBOX(cls, box, old_rec):
        '''
        parameters
        ----------
        box: array_like
             box_xmin, box_ymin, box_xmax, box_ymax
        old_rec: array_like, shape = [4, :]
             old rectangle, xmin, ymin, xmax, ymax
        return: 
        indx: array_like, bool
              x == -1, no corner in the box; 
              x ==  0, all corners fully in the box 
              x ==  1, left  bottom corner in the box
              x ==  2, right bottom corner in the box
              x ==  3, right upper  corner in the box
              x ==  4, left  upper  corner in the box
        new_rec: array_like, shape = [4, :] 
             new rectangle of joint region, xmin, ymin, xmax, ymax
        #-------------------------------------------------------
        #                          +------------(jxmax,jymax)
        #                          |                  |   
        #                +---------|--(ixmax,iymax)   |   xmin = max( ixmin, jxmin )
        #                |         |        |         |   ymin = max( iymin, jymin )
        #            height  (jxmin,jymin)------------+  
        #                |                  |             xmax = min( ixmax, jxmax )
        #      (ixmin,iymin)---- width -----+             ymax = min( iymax, jymax )
        #       
        #       if (xmax > xmin)&(ymax > ymin):  is  in the box 
        #       if (xmax < xmin)|(ymax < ymin):  not in the box 
        #-------------------------------------------------------

        '''
        rec = old_rec.copy()
        ixmin, iymin, ixmax, iymax = box 
        rec[0,:][ rec[0,:] < ixmin ] = ixmin 
        rec[1,:][ rec[1,:] < iymin ] = iymin
        rec[2,:][ rec[2,:] > ixmax ] = ixmax 
        rec[3,:][ rec[3,:] > iymax ] = iymax
        indx = np.zeros( np.shape(old_rec)[-1] )
        indx = (rec[2,:] <= rec[0,:])|(rec[3,:] <= rec[1,:])
        return ~indx, rec

    @classmethod
    def fillin(cls, img, stps, xc, yc, return_fill = False ): 
        ngal =  len(xc)
        xs   =  np.array( [stp.shape[1] for stp in stps] ) 
        ys   =  np.array( [stp.shape[0] for stp in stps] ) 
        ori_rec = np.array([xc-0.5*xs, yc-0.5*ys, xc+0.5*xs, yc+0.5*ys] ) 
        img_rec = np.array([-0.5, -0.5, img.shape[1], img.shape[0] ]) 
        indx, stp_rec = mockimg.isinBOX(img_rec, ori_rec) 
        #----------------------------------------------------------------
        ori_rec  = np.array(ori_rec +0.5).astype(int); 
        stp_rec  = np.array(stp_rec +0.5).astype(int); 
        img_rec  = np.array(img_rec +0.5).astype(int); 
        fill = np.zeros( (8, ngal) )  -1
        for ii in range(ngal): 
            x0, y0, x1, y1   = stp_rec[:,ii]
            xref, yref, _, _ = ori_rec[:,ii]
            x0_ = x0 - xref
            x1_ = x1 - xref
            y0_ = y0 - yref
            y1_ = y1 - yref
            if indx[ii]: fill[:,ii] = [x0, y0, x1, y1, x0_, y0_, x1_, y1_ ]
            if indx[ii]: img[y0:y1, x0:x1] = stps[ii][y0_:y1_, x0_:x1_] 
        return img, fill

    @classmethod
    def download_from_desidr9():
        pass
    
    @classmethod
    def wcs2hdu(cls, w):
        NAXIS1, NAXIS2 = w.pixel_shape
        img = np.zeros( (NAXIS2, NAXIS1) ).astype(float)
        hdu = fits.ImageHDU(img); 
        hdu.header.update( w.to_header()) 
        hdr = hdu.header
        return img, hdr

    @classmethod
    def save_as_fits(cls, outputname, img, hdr): 
        hdu = fits.ImageHDU( img ); 
        hdu.header = hdr
        hdul= fits.HDUList([fits.PrimaryHDU(), hdu])
        hdul.writeto( outputname, overwrite =  True)
        hdul.close()
        
    @classmethod
    def save_as_jpeg(cls, filename, img):
        split = filename.split('.')
        if split[-1]=='jpeg': 
            import imageio
            imageio.imsave(filename, img) 
    
    @classmethod 
    def wcs2reg(cls, w, npt = 1, return_world = True ): 
        '''
        Input: w, wcs object of astropy.wcs
        npt = 1 # Number of elements for each side of the pixel.
        Return: 
           The profile in the celestial coordinate system (ra, dec). 
        '''
        npt   = npt + 1 
        nxpix = w.pixel_shape[0]
        nypix = w.pixel_shape[1]
        nxrange = np.linspace(0, nxpix, npt) - 0.5
        nyrange = np.linspace(0, nypix, npt) - 0.5
        xbox = np.hstack([ nxrange[0:-1:1],       [nxrange[-1]]*(npt-1),   nxrange[-1:0:-1],     [nxrange[0]]*(npt-1) ] )#.astype( int )
        ybox = np.hstack([[nyrange[0]]*(npt-1),   nyrange[0:-1:1],         [nyrange[-1]]*(npt-1), nyrange[-1:0:-1] ] )#.astype( int )
        if return_world: 
            ra_box, dec_box = w.wcs_pix2world(xbox, ybox, 0); 
            return ra_box, dec_box
        else: 
            return xbox, ybox
    @classmethod  
    def downloads(cls, a, d, cache = './__cache__/', as_shell = None, survey = 'desi', pixscale=0.262, size = None, suffix = 'fits', bands = 'grz'):
        if survey == 'desi': 
            if pixscale is None: pixscale = 0.262
            if size is None: size = [30, 30];
        
        outputnames = [] 
        if as_shell is not None: f = open(as_shell, 'w+')
        else:
            if not os.path.exists(cache): os.makedirs(cache)
        for ii in range( len(a) ):
            a_ = a[ii]
            d_ = d[ii]
            outputname =  os.path.join(cache, '%s.'%ii + suffix)
            if  (os.path.exists(outputname)) : 
                pass
                #print('"%s" exists'%outputname) 
            elif (suffix!='jpeg')&(suffix!='fits'): 
                print('Suffix of outputname "%s" is not correct. '%outputname)
                print('Suffix should be "jpeg" or "fits"')
            else:  
                # os.system('rm -rf %s'%outputname)
                cmd = 'wget -O %s --no-check-certificate '% outputname
                web = ' "https://www.legacysurvey.org/viewer/%s-cutout?ra=%f&dec=%f&height=%s&width=%s&layer=ls-dr9&pixscale=%s&bands=%s" '%(suffix, a_, d_, size[0], size[1], pixscale, bands)
                if as_shell is not None: 
                    f.writelines([cmd + web + ' \r\n']) 
                else: 
                    os.system( cmd + web)
            outputnames.append(outputname)
        if as_shell is not None: f.close()  
        return outputnames
    
    @classmethod 
    def easy_bkg(cls, data): 
        from astropy.stats import sigma_clipped_stats, SigmaClip
        from photutils.segmentation import detect_threshold, detect_sources
        from photutils.utils import circular_footprint
        sigma_clip  = SigmaClip(sigma=3.0, maxiters=10)
        threshold   = detect_threshold(data, nsigma=2.0, sigma_clip=sigma_clip)
        segment_img = detect_sources(data, threshold, npixels=10)
        footprint = circular_footprint(radius=10)
        mask = segment_img.make_source_mask(footprint=footprint)
        mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
        return median, std
    @classmethod 
    def rgb(self, imgs, survey = 'desi'):
        from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog
        image_r = imgs[0]; bkg_r, std_r = mockimg.easy_bkg(image_r); 
        image_g = imgs[1]; bkg_g, std_g = mockimg.easy_bkg(image_g); 
        image_b = imgs[2]; bkg_b, std_b = mockimg.easy_bkg(image_b); 
        r = image_r - bkg_r#/std_r
        g = image_g - bkg_g#/std_g
        b = image_b - bkg_b#/std_b
        # 4.9,5.7,7.8 
        from astropy.visualization import SqrtStretch
        from astropy.visualization import ZScaleInterval
        from astropy.visualization import make_lupton_rgb
        # lo_val, up_val = np.percentile(np.hstack((r.flatten(), g.flatten(), b.flatten())), (0.05, 99.5)) 
        if survey == 'desi': 
            image = make_lupton_rgb(0.7*r, 1.1*g, 1.8*b, stretch=0.05, Q=8, minimum = -0.005, filename = None ) 
        else: 
            image = make_lupton_rgb(r, g, b, stretch=0.05, Q=10, minimum = 0, filename = None )
        image = image/255
        return image
    
    def run(self):
        if self.mask is None: 
            self.mask = np.ones(ngal).astype(bool)  
        if self.wcs is None:  
            if self.survey == 'desi': pixscale = 0.262/3600; rho = 0
            self.wcs = mockimg.autowcs(self.a[self.mask], 
                                       self.d[self.mask], 
                                        pixscale = pixscale, rho = 0, enlarge = self.enlarge )
            
        #--- 预处理 
        img, hdr = mockimg.wcs2hdu(self.wcs) 
        stprgbs    = self.stps
        print() 
        nchannel = len(stprgbs)
        if nchannel == 1: print('#--- Run as single channel')  
        if nchannel != 1: print('#--- Run as %s channel'%nchannel ) 
        imgrgb = np.array([img for ii in range(nchannel)]) 
        mskrgb = imgrgb.copy()*0.0
        
        self.xc,  self.yc = self.wcs.wcs_world2pix(self.a, self.d, 0); 
        for nc in range(nchannel): 
        #--- 确定图像的区域（Imaging strategy）
            img      = imgrgb[nc]
            img_msk  = mskrgb[nc]
        #--- 准备需要的子图（prepare stamps & Locations） 
            stps    = stprgbs[nc]            
            stp_msks= [stp*0+1 for stp in stps]
        #--- 抽取目标和背景(extract object mask & background noise)
            noises = [];
            for ii, stp in enumerate(stps):
                stp_msk, bkg_msk = sex(stp).return_mainobj() 
                stp_msks[ii] = stp_msk
                noises.append(stp[bkg_msk]) 
                stps[ii][~stp_msk] = 0 

            noises = np.hstack(noises) 
            
        #--- 叠加准备的子图（stamp superposition）
            img,     fill= mockimg.fillin(img,     stps,     self.xc,  self.yc)
            img_msk, fill= mockimg.fillin(img_msk, stp_msks, self.xc,  self.yc) 
        #--- 添加背景噪声（background noise)
            n_fill = np.count_nonzero(fill[0,:] != -1) 
            n_noise = np.count_nonzero(img_msk==0)
            if n_fill != 0: 
                img[img_msk==0] = np.random.choice( noises, n_noise, replace = True)
            imgrgb[nc] = img
            mskrgb[nc] = img_msk
        if nchannel == 1: 
            imgrgb = imgrgb[0]
            mskrgb = mskrgb[0]
        if nchannel != 1: 
            imgrgb = np.transpose(imgrgb, (1,2,0)) 
            mskrgb = np.transpose(mskrgb, (1,2,0)) 
        print('#--- return shape:', np.shape(imgrgb) ) 
        return imgrgb, mskrgb, hdr

