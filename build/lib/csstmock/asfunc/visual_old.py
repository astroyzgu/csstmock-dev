#------------------------------------------------------------------------------
#  This block provides a way to show the skymap based on the matplotlib.basemap
#
#  -- 0 < ra < 360; -90 < dec < 90
#
#  pix2border 
#  visual(ax = ax).moll() 
#  visual(ax = ax).cyl() 
#------------------------------------------------------------------------------
# Written by Gu Yizhou
# Refer:  https://github.com/kadrlica/skymap
#         https://matplotlib.org/basemap/api/basemap_api.html
#------------------------------------------------------------------------------

try:
    import matplotlib 
    from mpl_toolkits.basemap import Basemap
except:
    print () 
    print('Current matplotlib version: ', matplotlib.__version__) 
    print('module asfunc_old requires basemap=1.3.6. ')
    print('Degrade matplotlib to 3.4.3: conda install matplotlib=3.4.3')
    print('Then, using "pip install basemap" to install')
    print () 
    print('---------------------------------------------------------')
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp 
import tqdm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import os
from scipy.spatial import Delaunay


def add_colorbar(ax, label, ticks, vmax, vmin, cmap, position = "bottom" ): 
    '''
    parameter: 
    position: {"left", "right", "bottom", "top"}
        Where the new axes is positioned relative to the main axes.
    ''' 
    from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
    ax1 = make_axes_locatable(ax);
    cax = ax1.append_axes(position, size="3%", pad="6%")
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sm   = plt.cm.ScalarMappable(cmap=plt.get_cmap('plasma'), norm=norm) 

    if position in ["bottom", "top"]: orientation = 'horizontal' 
    if position in ["left", "right"]: orientation = 'vertical' 

    cb1 = plt.colorbar(sm, cax=cax, orientation = orientation ) 
    cb1.set_label(label) # , fontsize = 14)
    cb1.set_ticks(ticks) 



class Skymap(Basemap): 
    '''
    Draw skymap by matplotlib.Basemap
    '''
    def __init__(self, *arg, **kwargs): 
        super(Skymap, self).__init__(*arg, **kwargs)
        self.resolution = kwargs.pop('resolution', 'l') 
        self.projection = kwargs.pop('projection', None) 
        self.lon_0      = kwargs.pop('lon_0', 0) 
        if self.projection == 'moll':  
           self.proj_xmin = self.lon_0 - 180.0; 
           self.proj_xmax = self.lon_0 + 180.0; 
           self.proj_ymin = -90;  self.proj_ymax  = 90
        if self.projection == 'cyl':  
           axis            = kwargs.pop('axis', None) 
           if axis == None: axis = [0.0, 360.0, -90, 90]
           self.proj_xmin = axis[0]; self.proj_xmax = axis[1]
           self.proj_ymin = axis[2]; self.proj_ymax = axis[3]

    def fill_boundary(self, boundary, wht = None, **kwargs): 
        #print(wht)
        #print(len(wht))
        wrapleft = self.lon_0 - 180.00000; wrapright = self.lon_0 + 180
        patches = []; new_wht = [] 
        for ith, ipix, val in tqdm.tqdm(zip(range(len(boundary)), boundary.keys(),boundary.values()), total = len(boundary))  :
            ra, dec = split_into_wraps(val[::2], val[1::2], wrapleft = wrapleft, wrapright = wrapright)
            for ra_, dec_ in zip(ra, dec): 
                ra_[ra_ <= wrapleft + 1E-08]  = wrapleft   + 1E-8
                ra_[ra_ >= wrapright- 1E-08]  = wrapright  - 1E-8
                x_, y_  = self(ra_, dec_)
                if((len(x_)>2)): 
                    polygon = Polygon(np.vstack((x_, y_)).T, True)
                    patches.append(polygon) 
                    #print(ith, wht[ith-1]) 
                if((len(x_)>2)&(not wht is None)): new_wht.append(wht[ith])
        if wht is None: 
           print('# plot_boundary')
           edgecolor   = kwargs.pop('edgecolor', 'k') 
           facecolor   = kwargs.pop('facecolor', 'w')
           p = PatchCollection(patches, edgecolor = edgecolor, facecolor = facecolor, **kwargs) 
        else: 
           print('# fill_boundary')
           cm   = kwargs.pop('cmap', plt.get_cmap('rainbow')) 
           vmin = kwargs.pop('vmin', 1.0*np.min(wht))
           vmax = kwargs.pop('vmax', 1.0*np.max(wht))
           p = PatchCollection(patches, cmap = cm, **kwargs) 
           p.set_array(np.array(new_wht))
           p.set_clim([vmin, vmax])
        self.ax.add_collection(p)


    def hpviewer(self, vmap, nest = False, npt = 1, **kwargs): 
        nside = np.sqrt( len(vmap)/12.0)
        if not nside in 2**np.arange(0, 30): 
            ValueError("Length %s of vmap must equal to 12*Nside^2, Nside = 1,2,3,4,...,29."% len(vmap) )
        nside = int(nside); # print(nside)
        ipix  = np.arange( 12*nside*nside )
        pixid = ipix[~np.isnan(vmap)] 
        vmap  = vmap[~np.isnan(vmap)] 
        ra,  dec = boundaries_ang(nside, pixid, nest, npt = npt)
        wrapleft = self.lon_0 - 180.00000; wrapright = self.lon_0 + 180
        patches = []; new_wht = [] 
        for ith, ipix in enumerate(pixid):
            ra__, dec__ = split_into_wraps(ra[ith], dec[ith], wrapleft = wrapleft, wrapright = wrapright)
            for ra_, dec_ in zip(ra__, dec__): 
                ra_[ra_ <= wrapleft + 1E-08]  = wrapleft   + 1E-8
                ra_[ra_ >= wrapright- 1E-08]  = wrapright  - 1E-8
                if len(ra_) <= 2: continue 
                x_, y_  = self(ra_, dec_) 
                polygon = Polygon( np.vstack((x_, y_)).T, True )
                patches.append(polygon) 
                new_wht.append(vmap[ith])
        if ('edgecolor' in kwargs.keys())|('facecolor' in kwargs.keys()): 
            edgecolor   = kwargs.pop('edgecolor', None) 
            facecolor   = kwargs.pop('facecolor', None)
            p = PatchCollection(patches, edgecolor = edgecolor, facecolor = facecolor, **kwargs) 
            self.ax.add_collection(p)
            return p 
 
        if ('cmap' in kwargs.keys()): 
            cm   = kwargs.pop('cmap', plt.get_cmap('rainbow')) 
            vmin = 1.0*np.min(vmap)
            vmax = 1.0*np.max(vmap)
            vmin = kwargs.pop('vmin', vmin)
            vmax = kwargs.pop('vmax', vmax)
            norm = plt.Normalize(vmin=vmin, vmax=vmax)
            norm = kwargs.pop('vmax', norm)
            p = PatchCollection(patches, cmap = cm, **kwargs) 
            p.set_array(np.array(new_wht))
            p.set_clim([norm.vmin, norm.vmax])
            p.set_norm( norm )
            self.ax.add_collection(p)
            return p


    def hpviewer_fast(self, vmap, nest = False, npt = 2, **kwarg): 
        nside = np.sqrt( len(vmap)/12.0)
        if not nside in 2**np.arange(0, 30): 
           raise ValueError("Length %s of vmap must equal to 12*Nside^2, Nside = 1,2,3,4,...,29."% len(vmap) )
        nside = int(nside); # print(nside)
        ipix  = np.arange( 12*nside*nside )
        vmap[np.isnan(vmap)] = -np.inf  
        from scipy.spatial import Delaunay
        ra,  dec = boundaries_ang(nside, ipix, nest, npt = npt)
        ra  = np.hstack(ra); #print(np.shape(ra)) 
        dec = np.hstack(dec);
        ra[ra <= self.proj_xmin]  = ra[ra  <= self.proj_xmin] + 360
        ra[ra >= self.proj_xmax]  = ra[ra  >= self.proj_xmax] - 360
        ra[(ra <= self.proj_xmin + 1E-8)]   = self.proj_xmin + 1E-8 
        ra[(ra >= self.proj_xmax - 1E-8)]   = self.proj_xmax - 1E-8 
        ra_app  = [self.proj_xmin+1E-8]*18 + [self.proj_xmax-1E-8]*18 
        dec_app = np.append( np.linspace(self.proj_ymin, self.proj_ymax, 18), np.linspace(self.proj_ymin, self.proj_ymax, 18))
        ra  = np.append( ra,  ra_app ) 
        dec = np.append(dec, dec_app ) 
        #--------------------------------------
        points_ang, index = np.unique( np.vstack([ra, dec]).T , return_index = True, axis = 0 ) 
        #print( ' points grid :', np.shape(points_ang) )
        #points_xy, index = np.unique( np.vstack([x, y]).T , return_index = True, axis = 0 ) 
        tri    =  Delaunay( points_ang )
        amid   =  tri.points[:,0][tri.simplices].mean(axis=1)
        dmid   =  tri.points[:,1][tri.simplices].mean(axis=1)
        ipix   =  hp.ang2pix(nside, amid, dmid, nest = nest, lonlat = True)
        zfaces = vmap[ipix] 
        x, y   = self( tri.points[:,0],  tri.points[:,1])
        mask   = np.zeros(np.shape(zfaces)).astype(bool)
        mask[np.isinf(zfaces)] = True
        #cm = self.ax.triplot(  x, y , tri.simplices, color = 'k', **kwarg)
        #cm = self.ax.tripcolor( x, y , tri.simplices, facecolors= zfaces, **kwarg)
        p = self.ax.tripcolor( x, y , tri.simplices, mask = mask, facecolors= zfaces, **kwarg)
        p.set_edgecolor('none')

        #--------------------------------------
        #for xmid_, ymid_, ipix_ in zip(amid, dmid, ipix): 
        #    xmid_, ymid_ = self(xmid_, ymid_)
        #    self.ax.scatter( xmid_, ymid_, color = 'k' )
        #    self.ax.text( xmid_, ymid_,  ipix_ )
        #ra = 20; dec = 30; xm, ym = self(ra, dec)
        #elf.ax.scatter( xm, ym, color = 'r' )
        #ra, dec   = self(xm, ym, inverse = True )
        #print(ra, dec, self.proj_xmax, self.lon_0, self.proj_xmin) 
        #phase = int( (ra - ra%(self.proj_xmax-self.proj_xmin))/(self.proj_xmax-self.proj_xmin) )
        #ra    = ra - phase*(self.proj_xmax-self.proj_xmin); 
        #xm, ym = self(ra, dec)
        #self.ax.scatter( xm, ym, color = 'g' )
        #print(ra, dec, self.proj_xmax, self.lon_0, self.proj_xmin, phase) 
        return p 
 
    def plot(self, ra, dec, **kwargs): 
        wrapleft = self.lon_0 - 180.00000; wrapright = self.lon_0 + 180
        ra, dec  = split_into_wraps(ra, dec, wrapleft = wrapleft, wrapright = wrapright)
        for ra_, dec_ in zip(ra, dec): 
            ra_[ra_ <= wrapleft + 1E-08]  = wrapleft   + 1E-8
            ra_[ra_ >= wrapright- 1E-08]  = wrapright  - 1E-8
            xy = self(ra_, dec_)
            self.ax.plot(*xy, **kwargs)
              
    def plot_poly(self, a, d, npt = 11, **kwargs): 
        for ii in range(-1, len(a)-1):
             a1 = a[ii]; a2 = a[ii+1]
             d1 = d[ii]; d2 = d[ii+1]
             xx = np.linspace(a1, a2, npt)
             yy = np.linspace(d1, d2, npt)
             self.plot(xx, yy, **kwargs)
      
    def text(self, x, y, text, **kwargs): 
        xpt, ypt = self(x, y)
        self.ax.text(xpt, ypt, text, **kwargs)

    def scatter(self, ra, dec, **kwargs):
        xy = self(ra, dec)
        self.ax.scatter( *xy, **kwargs) 

    def circle(self,a,d,r, npt = 101, **kwargs): 
        x,  y,  z  = ruv2xyz(1.0, a, d)
        x0, y0, z0 = ruv2xyz(1.0, a, d-r)
        t          = np.linspace(0, 360, npt) 
        x_, y_, z_ = rot_by_spinaxis(x0, y0, z0, x, y, z, t) 
        r_, a_, d_ = xyz2ruv(x_, y_, z_)
        a_[a_ < 0]  =  a_[a_ < 0]  + 360 
        xy = self(a_,d_)
        self.ax.plot(*xy, **kwargs)

    def set_parallels(self, ticks = [-60, -30, 0, 30, 60],  shift = 0,  color='gray', labels = [1,0,0,0], **kwargs): 
        if(self.projection != 'moll'): 
            self.drawparallels(ticks,color=color,dashes=[1,3], labels= labels, labelstyle = '+/-', **kwargs)
        if(self.projection == 'moll'): 
            for tick in ticks: 
                wrapleft = self.lon_0 - 179.9; wrapright = self.lon_0 + 180
                x, y = self(wrapleft + shift, tick)
                if tick >=0: va = 'center' 
                if tick < 0: va = 'top' 
                self.ax.text(x+10000, y, r'$ %i^\circ$'%tick, ha = 'right', va= va)
                self.drawparallels(ticks,color=color,dashes=[1,3], labels= [0,0,0,0], labelstyle = '+/-', **kwargs)

         
    def set_meridians(self, ticks = np.arange(0, 360, 30), shift = 0,  color='gray', labels = [0,0,0,1],  **kwargs): 
        if(self.projection == 'moll'):
            for xlabel in ticks: 
                mx, my = self(xlabel, shift) 
                self.ax.text( mx, my, r'$ %i^\circ$'%xlabel,  **kwargs)
                self.drawmeridians(ticks,color=color,dashes=[1,3],labels = [0,0,0,0], **kwargs)
        if(self.projection != 'moll'):
            self.drawmeridians(ticks,color=color,dashes=[1,3], labels=labels, labelstyle = '+/-',  **kwargs)

    def set_galactic(self, b0 = 0, **kwargs):
        lon_0= self.lon_0  
        l0   = np.linspace(-180, 179.99, 360); b0   = b0 + 0.0*l0
        obj0 = SkyCoord(l0, b0, frame='galactic', unit='deg')
        ra0  = obj0.icrs.ra.degree  
        dec0 = obj0.icrs.dec.degree
        ra0[ra0 < lon_0 -180] = ra0[ra0 < lon_0 -180] + 360
        ra0[ra0 > lon_0 +180] = ra0[ra0 > lon_0 +180] - 360
        idx = np.argmin(ra0) 
        ra0  = np.roll( ra0, -idx)
        dec0 = np.roll(dec0, -idx)
        return self.plot(ra0, dec0, **kwargs)
        
    def set_ecliptic(self, b0 = 0, **kwargs):
        lon_0= self.lon_0  
        l0   = np.linspace(-180, 179.99, 360); b0   = b0 + 0.0*l0
        obj0 = SkyCoord(l0, b0, frame='geocentricmeanecliptic', unit='deg')
        ra0  = obj0.icrs.ra.degree  
        dec0 = obj0.icrs.dec.degree
        ra0[ra0 < lon_0 -180] = ra0[ra0 < lon_0 -180] + 360
        ra0[ra0 > lon_0 +180] = ra0[ra0 > lon_0 +180] - 360
        idx = np.argmin(ra0) 
        ra0  = np.roll( ra0, -idx)
        dec0 = np.roll(dec0, -idx)
        return self.plot(ra0, dec0, **kwargs)   
        #isort = np.argsort(ra0)
        #ra0 = ra0[isort] ; dec0 = dec0[isort]
        #xpt, ypt = self(ra0, dec0 ) # degree
        #self.ax.plot( xpt, ypt, color = 'r', linewidth = 2, label = 'galactic plane') 

    def set_gc(self, **kwargs): 
        obj0 = SkyCoord([0, 0], [-90, 90], frame='galactic', unit='deg')
        ra0  = obj0.icrs.ra.degree
        dec0 = obj0.icrs.dec.degree
        xpt, ypt = self(ra0, dec0 ) # degree
        self.ax.scatter( xpt, ypt, **kwargs ) 


    def fill_poly(self, ra, dec, **kwargs): 
        edgecolor   = kwargs.pop('edgecolor', 'k') 
        facecolor   = kwargs.pop('facecolor', 'w') 
        wrapleft = self.lon_0 - 180.00000; wrapright = self.lon_0 + 180
        ra, dec = split_into_wraps(ra, dec, wrapleft = wrapleft, wrapright = wrapright)
        patches = [] 
        for ra_, dec_ in zip(ra, dec): 
            ra_[ra_ <= wrapleft + 1E-08]  = wrapleft   + 1E-8
            ra_[ra_ >= wrapright- 1E-08]  = wrapright  - 1E-8
            x_, y_ = self(ra_, dec_)
            #self.ax.fill(*xy, **kwargs)
            polygon = Polygon(np.array([x_, y_]).T, True)
            patches.append(polygon) 
        p = PatchCollection(patches, **kwargs) # cmap=plt.get_cmap('rainbow'), edgecolor='k') 
        self.ax.add_collection(p)

    def draw_poly(self, ra, dec, **kwargs): 
        edgecolor   = kwargs.pop('edgecolor', 'k') 
        facecolor   = kwargs.pop('facecolor', 'w') 
        patches = [] 
        x_, y_  = self(ra, dec)
        polygon = Polygon(np.array([x_, y_]).T, True)
        patches.append(polygon) 
        p = PatchCollection(patches, edgecolor = edgecolor, facecolor = facecolor, **kwargs) # cmap=plt.get_cmap('rainbow'), edgecolor='k') 
        self.ax.add_collection(p)
        #p.set_array(whts)
        #p.set_clim([min(whts),max(whts)])

    def draw_box(self, axis, **kwargs): 
        llcrnrlon = axis[0] 
        urcrnrlon = axis[1] 
        llcrnrlat = axis[2]
        urcrnrlat = axis[3]
        ra  = [llcrnrlon, urcrnrlon, urcrnrlon, llcrnrlon, llcrnrlon]
        dec = [llcrnrlat, llcrnrlat, urcrnrlat, urcrnrlat, llcrnrlat]
        self.draw_path(ra, dec, **kwargs)


    def fill_boundary(self, boundary, wht = None, **kwargs): 
        #print(wht)
        #print(len(wht))
        parallel  = kwargs.pop('parallel', 1) 
        wrapleft = self.lon_0 - 180.00000; wrapright = self.lon_0 + 180
        patches = []; new_wht = [] 
        if( parallel != 1):  
            print('Starting Parallels with %s processes'%parallel)
            pool = mp.Pool(parallel)
            p    = tqdm.trange(len(ipixs))
            results = []
            print('aaa')
        else: 
            for ith, ipix, val in tqdm.tqdm(zip(range(len(boundary)), boundary.keys(),boundary.values()), total = len(boundary))  :
                ra, dec = split_into_wraps(val[::2], val[1::2], wrapleft = wrapleft, wrapright = wrapright)
                for ra_, dec_ in zip(ra, dec): 
                    ra_[ra_ <= wrapleft + 1E-08]  = wrapleft   + 1E-8
                    ra_[ra_ >= wrapright- 1E-08]  = wrapright  - 1E-8
                    x_, y_  = self(ra_, dec_)
                    if((len(x_)>2)): 
                        polygon = Polygon(np.vstack((x_, y_)).T, True)
                        patches.append(polygon) 
                        #print(ith, wht[ith-1]) 
                    if((len(x_)>2)&(not wht is None)): new_wht.append(wht[ith])
        if wht is None: 
           print('# plot_boundary')
           edgecolor   = kwargs.pop('edgecolor', 'k') 
           facecolor   = kwargs.pop('facecolor', 'w')
           p = PatchCollection(patches, edgecolor = edgecolor, facecolor = facecolor, **kwargs) 
        else: 
           print('# fill_boundary')
           cm   = kwargs.pop('cmap', plt.get_cmap('rainbow')) 
           vmin = kwargs.pop('vmin', 1.0*np.min(wht))
           vmax = kwargs.pop('vmax', 1.0*np.max(wht))
           p = PatchCollection(patches, cmap = cm, **kwargs) 
           p.set_array(np.array(new_wht))
           p.set_clim([vmin, vmax])
        self.ax.add_collection(p)






class visual(): 
    def __init__(self, *arg, **kwargs):
        self.ax         = kwargs.pop('ax', None )

    def moll(self, lon_0 = 180, xticks = np.arange(0, 360, 30), yticks = np.arange(-60, 61, 30)): 
        m = Skymap(projection = 'moll', lon_0 = lon_0, ax = self.ax, resolution='l', celestial = True)
        self.lon_0 = lon_0
        #m.set_parallels(yticks)
        #m.set_meridians(xticks)
        # m.set_galactic()
        # m.set_gc()
        # self.ax.set_title("moll")
        return m 

    def cyl(self, axis = None, xticks = None, yticks = None ):
        if axis == None: axis = [0, 360, -90, 90]
        if axis[1] - axis[0] <= 60: d_ticks = 5
        if axis[1] - axis[0] >  60: d_ticks = 10
        if axis[1] - axis[0] > 120: d_ticks = 30
        llcrnrlon = axis[0] 
        urcrnrlon = axis[1] 
        llcrnrlat = axis[2]
        urcrnrlat = axis[3]
        self.lon_0 = 0.5*(llcrnrlon + urcrnrlon)
        m = Skymap(projection = 'cyl', ax = self.ax, resolution='l', lon_0 = self.lon_0,   
        llcrnrlon = llcrnrlon,  llcrnrlat = llcrnrlat,  urcrnrlon = urcrnrlon,  urcrnrlat = urcrnrlat)
        #xticks = np.arange(llcrnrlat, urcrnrlat+0.01, 30)
        #yticks = np.arange(llcrnrlon, urcrnrlon+0.01, 30)
        if xticks == None: xticks  = np.arange(0, 360, d_ticks); 
        if yticks == None: yticks = np.arange(-90, 91, d_ticks) 
        m.set_parallels(ticks = yticks, labels = [0,0,0,0])
        m.set_meridians(ticks = xticks, labels = [0,0,0,0])
        self.ax.set_yticks( yticks )
        self.ax.set_xticks( xticks )
        self.ax.axis(axis)
        self.ax.set_xlabel('RA')
        self.ax.set_ylabel('DEC')
        #self.ax.set_title("cyl")
        return m 

#--- start boundaries_ang
#-----------------------------------------------------------------------------------------------
def boundaries_ang(nside, ipixs, nest = False, npt = 1): 
    ra_collect = []; dec_collect = []
    ipixs = np.atleast_1d(ipixs)
    for ipix in ipixs: 
        x, y, z = hp.boundaries(nside, ipix, step = npt, nest = nest)
        a0, d0  = hp.pix2ang(nside, ipix, nest =  nest, lonlat= True)
        a0 = a0%360
        _, a, d = xyz2ruv(x, y, z) 
        a[a - a0 > 180] = a[a - a0 > 180] - 360
        a[a - a0 < -180] = a[a - a0 < -180] + 360
        if d[0] == 90: 
            a[0] = a[1]; 
            a = np.append(a, a[-1]); d = np.append(d, 90)
        if d[npt*2] == -90: 
            a = np.roll(a, 2*npt)
            d = np.roll(d, 2*npt)
            a[0] = a[1]; 
            a = np.append(a, a[-1]); d = np.append(d, -90)
        ra_collect.append(a)
        dec_collect.append(d)
    return ra_collect, dec_collect
def boundaries_ang_(nside, ipixs, nest = False, npt = 1): 
    ra_collect = []; dec_collect = []
    ipixs = np.atleast_1d(ipixs)
    for ipix in ipixs: 
        x, y, z    = hp.boundaries(nside, ipix, step = npt, nest = nest)  
        ra0, dec0  = hp.pix2ang(nside, ipix, nest =  nest, lonlat= True) 
        _, ra, dec = xyz2ruv(x, y, z)
        wrapleft = ra0 - 180; wrapright = ra0 + 180
        ra[ra <= wrapleft ]  = ra[ra <= wrapleft]  + 360
        ra[ra >= wrapright]  = ra[ra >= wrapright] - 360
        ra[(ra <= 1E-8)&(ra >= -1E-8)]  = 0 
        ra[(ra <= 360+1E-8)&(ra >= 360-1E-8)] = 360
        if(((ipix+1)%(nside*nside)==0)&( ipix< 4*nside*nside)): 
           ra[0]=ra[1]
           ra  = np.append( ra[:2*npt+1:1], 2*ra[2*npt]-ra[2*npt::-1])
           dec = np.append(dec[:2*npt+1:1], dec[2*npt::-1])
        if( (ipix   %(nside*nside)==0)&((ipix>=8*nside*nside))): 
           ra   = np.roll(ra, 2*npt)
           dec  = np.roll(dec,2*npt)
           ra[0]=ra[1]
           ra  = np.append( ra[:2*npt+1:1], 2*ra[2*npt]-ra[2*npt::-1])
           dec = np.append(dec[:2*npt+1:1], dec[2*npt::-1])
        ra_collect.append(ra)
        dec_collect.append(dec)
    return ra_collect, dec_collect
#--- end boundaries_ang
#-----------------------------------------------------------------------------------------------

### --- need module
def supply_at_wrap(x_, y_, wrap = 360): 
    diff = x_ - np.roll(x_, -1)
    supply_indx = np.where((diff == 0.0)&(x_ == wrap))[0]
    interp_index = []; interp_value = []; n =  len(x_)
    for ii in supply_indx: 
        for y_insert in np.linspace(y_[ii], y_[ii-n+1],  11):
            interp_index.append(ii-n+1)
            interp_value.append(y_insert)
    #print('interp:', interp_index)
    x0 = np.insert(x_, interp_index, wrap)
    y0 = np.insert(y_, interp_index, interp_value)
    return x0, y0
def interp_at_wrap(x0, y0, wrap = 360): 
    x0 = np.atleast_1d(x0)
    y0 = np.atleast_1d(y0)
    diff   = x0 - wrap
    diff   = diff*np.roll(diff, -1)
    interp_index = []; interp_value = []
    n = len(diff)
    for i, d in enumerate(diff): 
        if(d>=0): continue
        if(d< 0): 
            interp_index.append(i+1)
            k = (y0[i] -  y0[i-n+1])/(x0[i] -  x0[i-n+1])
            y = k*(wrap - x0[i]) + y0[i]
            interp_value.append( y ) 
    x0 = np.insert(x0, interp_index, wrap)
    y0 = np.insert(y0, interp_index, interp_value)
    return x0, y0
def split_at_wrap(x0, y0, wrap = 360, wrapvalue = -360): 
    if(wrapvalue < 0): left_value = 0.0; right_value = wrapvalue;
    if(wrapvalue > 0): left_value = wrapvalue; right_value = 0.0;
    #---------------------------------------------------
    split_index = np.where(x0 <= wrap)[0]
    x_left = x0[split_index] 
    y_left = y0[split_index]
    x_left, y_left  = supply_at_wrap(x_left, y_left, wrap = wrap)
    x_left = x_left  + left_value
    #---------------------------------------------------
    split_index = np.where(x0 >= wrap)[0]
    x_right = x0[split_index] 
    y_right = y0[split_index]
    x_right, y_right = supply_at_wrap(x_right, y_right, wrap = wrap)
    x_right = x_right + right_value
    return [x_left, x_right], [y_left, y_right]
def split_into_wraps(x0, y0, wrapleft = 0.0, wrapright = 360.0): 
    x0 = np.atleast_1d(x0)
    y0 = np.atleast_1d(y0)
    wrap_value = wrapright - wrapleft
    if (x0 >=wrapright).any()&(x0<=wrapleft).any(): 
        #print(x0, (x0 >= wrapright).any(), wrapleft)
        raise Exception('Input data out of the wrapleft and wrapright')
    elif (x0 <= wrapleft).any(): 
        x0, y0 = interp_at_wrap(x0, y0, wrap = wrapleft)
        x0, y0 = split_at_wrap( x0, y0, wrap = wrapleft,  wrapvalue =    wrap_value) 
    elif (x0 >= wrapright).any(): 
        x0, y0 = interp_at_wrap(x0, y0, wrap = wrapright)
        x0, y0 = split_at_wrap( x0, y0, wrap = wrapright, wrapvalue = -1*wrap_value)  
    else: 
        x0 = [x0]
        y0 = [y0]
    return x0, y0
#--- end visual module 
#----------------------------------------------------------------------------------------------- 

##### 笛卡尔三轴坐标操作 #####
# 点乘 numpy.matmul(x1, x2)
# 叉乘 numpy.cross(x1, x2)
def xyz2ruv(x, y, z): # , u0 = 0): 
    '''
    球坐标-> 笛卡尔坐标
    u0 = 0, ==> u = [-180, 180]; u0 = 180, ==> u = [0, 360]
    v0 = 0, ==> v = [-90, 90];   v0 = 90, ==> v = [0, 180]
    '''
    r = x**2 + y**2 + z**2
    u = np.arctan2(y,x)
    v = np.arctan2(z,np.sqrt(x**2+y**2))
    #if (u0 != 180)&(u0 != 0): raise('u0 must be 0 or 180')
    #if (v0 != 90)&(u0 != 0): raise('v0 must be 0 or 90')
    u = u*180/np.pi;  #u[u > u0 -180]=u[u > u0 -180]+360 
    v = v*180/np.pi
    return r, u, v

def ruv2xyz(r, u, v):
    '''
    笛卡尔坐标-> 球坐标
    u = [0, 360]
    v = [-90, 90]
    '''
    x = r*np.cos(v/180.0*np.pi)*np.cos(u/180.0*np.pi)
    y = r*np.cos(v/180.0*np.pi)*np.sin(u/180.0*np.pi)
    z = r*np.sin(v/180.0*np.pi)
    return x, y, z

def rot_by_spinaxis(x0, y0, z0, x, y, z, t): 
    '''
    向量(x0, y0, z0)绕轴(x,y,z)旋转t度
    '''
    t = t/180.0*np.pi
    r = np.sqrt(x**2 + y**2 + z**2)
    x = 1.0*x/r; y = 1.0*y/r; z = 1.0*z/r
    cs  = np.cos(t); sn = np.sin(t)
    #rot = np.array( [[ cs+(1-cs)*x**2, (1-cs)*x*y-sn*z,(1-cs)*x*z+sn*y], 
    #                 [(1-cs)*x*y+sn*z,  cs+(1-cs)*y**2,(1-cs)*y*z-sn*x],
    #                 [(1-cs)*x*z-sn*y, (1-cs)*y*z+sn*x, cs+(1-cs)*z**2]])
    #arr = np.matmul(rot, arr0)
    # x_ = arr[0]; y_ = arr[1]; z_ = arr[2]
    x_ = ( cs+(1-cs)*x**2)*x0 + ( (1-cs)*x*y-sn*z )*y0 + ((1-cs)*x*z+sn*y)*z0
    y_ = ((1-cs)*x*y+sn*z)*x0 + (  cs+(1-cs)*y**2 )*y0 + ((1-cs)*y*z-sn*x)*z0
    z_ = ((1-cs)*x*z-sn*y)*x0 + ( (1-cs)*y*z+sn*x )*y0 + ( cs+(1-cs)*z**2)*z0
    return x_, y_, z_
