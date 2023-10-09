import os
import subprocess
import numpy as np
from astropy.io  import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord 
import imageio

class sex(object): 
    def __init__(self, img, hdr = None,  verbose = 0 ): 
        import photutils
        # print(photutils.__version___) # == 1.3.0 ) 
        self.img       = img
        self.hdr       = hdr
        self.verbose   = verbose
        # self.mask      = detect.make_source_mask(self.img, nsigma=1, npixels=npixels, dilate_size=10)
        from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog
        from photutils.background import Background2D, MedianBackground
        bkg_estimator = MedianBackground()
        bkg       = Background2D(self.img, (20, 20), filter_size=(3, 3), bkg_estimator=bkg_estimator)
        self.img -= bkg.background  # subtract the background
        self.threshold = bkg.background_rms        


    def return_mainobj(self, npixels = 10, sigma = 3.0, kernel_size = 5, threshold = 0.5, contrast = 0.3):
        for kk in [0,]: 
            from astropy.convolution import convolve
            from astropy.convolution import Gaussian2DKernel
            from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog, SourceFinder 
            # cold mode to detect large objects
            if kk == 0: 
                sigma          = 3.0; 
                kernel         = Gaussian2DKernel(sigma, x_size=kernel_size, y_size=kernel_size); 
                convolved_data = convolve(self.img, kernel) # kernel.normalize()
                self.seg  = detect_sources(convolved_data, threshold*self.threshold, npixels=npixels) # , kernel=kernel)
                self.seg  = deblend_sources(convolved_data, self.seg, npixels=npixels, nlevels = 32, contrast = contrast, progress_bar=False )
            # hot  mode to detect more small objects
            if kk == 1: 
                sigma          = 1.0; npixels = 10;
                kernel         = Gaussian2DKernel(sigma, x_size=5, y_size=5); 
                convolved_data = convolve(self.img, kernel) # kernel.normalize()
                self.seg  = detect_sources(convolved_data, 0.5*self.threshold, npixels=npixels) # , kernel=kernel)             
                self.seg  = deblend_sources(convolved_data, self.seg, npixels=npixels, nlevels = 32, contrast = 0.01, progress_bar=False, mode = 'linear' )
            self.tab   = SourceCatalog(  self.img, self.seg).to_table()
            x0  = 0.5*np.shape(self.img)[0]; 
            y0  = 0.5*np.shape(self.img)[1]
            x   = np.atleast_1d(self.tab['xcentroid'])
            y   = np.atleast_1d(self.tab['ycentroid'])
            pix_dist = (x - x0)**2 + (y - y0)**2
            imain    = np.argmin(pix_dist.value) 
            signal   = self.tab['max_value'][imain]
            noise    = np.mean(self.threshold)
            if (signal > 0.2)|(signal/noise > 50): break 
            #signal   = self.tab['segment_flux'][imain]/self.tab['area'][imain]
            #print(signal.value, noise, signal.value/noise) 
            #if (signal.value/noise > 3)&(kk==0): break 

        #     if self.hdr is None:  
        #         from astropy import wcs
        #         wcs2d           = wcs.WCS(self.hdr)
        #         ra, dec = wcs2d.wcs_pix2world( self.tab['xcentroid'], self.tab['ycentroid'], 0) 
        #         self.tab['ra']  = ra 
        #         self.tab['dec'] = dec
        #     if(self.verbose!=0): print('Eesy Source Extract get %s potential sources'%len(self.tab) )
        # else: 
        #     self.tab = None
            
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
    
    
class mockimg(): 
    
    def __init__(self, outputdir = './output/', survey = 'desidr9', band = 'z'): 
        if not os.path.isdir(outputdir): os.makedirs(outputdir)
        self.outputdir = outputdir 
        self.survey    = survey
        self.band      = band

        
    def downloads(self, a, d, as_shell = None, survey = None, pixscale=0.262, size = None, bands = None, suffix = 'fits'):
        cache  = os.path.join(self.outputdir, 'downloads')
        
        if survey is None: survey = self.survey
        if bands  is None: bands  = self.band

        if  self.survey  == 'desidr9':
            if pixscale is None: pixscale = 0.262
            if size is None: size = [30, 30];
        else: 
            ValueError('Invalid survey name: %s'%survey ) 
        
        outputnames = [] 
        if not os.path.exists(cache): os.makedirs(cache)
            
        for ii in range( len(a) ):
            a_ = a[ii]
            d_ = d[ii]
            outputname =  os.path.join(cache, '%s-%s_%.6f_%.6f.'%(self.survey, bands, a_, d_) + suffix)
            # print(outputname) 
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
                    print(ii, cmd + web )
                    subprocess.run(cmd + web, shell=True, capture_output=True) 
                    # os.system( cmd + web + ' > /dev/null' )
            outputnames.append(outputname)
        if as_shell is not None: f.close()  
        return outputnames

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
    def readimgs(cls, inputs, ihdu = 0, aslist = False, return_headers = False ): 
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
            if(len(inputs)==1): hdrs = hdrs[0]
            #print('###-----------------------------------------------------------------')
            #print('### The shape of image array = ', imgs.shape)
        
        if return_headers: 
            return imgs, hdrs
        else:
            return imgs
    @classmethod
    def wcs2hdu(cls, w):
        NAXIS1, NAXIS2 = w.pixel_shape
        img = np.zeros( (NAXIS2, NAXIS1) ).astype(float)
        hdu = fits.ImageHDU(img); 
        hdu.header.update( w.to_header()) 
        hdr = hdu.header
        return img, hdr
    
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

    def autowcs(self, a, d, enlarge = 1.01):
        '''
        autowcs by input coordinates
        '''
        # pixscale = 0.1, rho = 0, 
        if  self.survey  == 'desidr9': pixscale = 0.262/3600; rho = 0; 
            
        a_ = np.array([np.min(a), np.max(a)])
        d_ = np.array([np.min(d), np.max(d)])
        x, y, z = mockimg.ruv2xyz(1, a_, d_)
        x0 = np.mean(x)
        y0 = np.mean(y)
        z0 = np.mean(z)
        _,a0,d0 = mockimg.xyz2ruv(x0,y0,z0)
        obj0 = SkyCoord(a0, d0, unit='deg')
        obj  = SkyCoord(a,  d,  unit='deg')
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
        # print('#--- pixscale = %f degree/pixel'%pixscale )
        # print(w) 
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
            if indx[ii]: img[y0:y1, x0:x1] = img[y0:y1, x0:x1] + stps[ii][y0_:y1_, x0_:x1_] 
        return img, fill
    
    
    @classmethod 
    def run(cls, a, d, stps, wcs, stp_msks = None, bkg_msks = None): 
            
        #--- 预处理 
        if np.ndim(stps[0]) !=2: ValueError('Error format of stamps')
        
        #--- 确定图像的区域（Imaging strategy）
        img,hdr = mockimg.wcs2hdu(wcs)
        img_msk = img.copy()*0.0
        
        #--- 准备需要的子图（prepare stamps & Locations）         
        xc,  yc = wcs.wcs_world2pix(a, d, 0); 

        #--- 抽取目标和背景(extract object mask & background noise)
        if stp_msks is None: stp_msks= [stp*0+1 for stp in stps]
        
        noises = [];        
        for ii, stp in enumerate(stps):
            if (stp_msks is None)|(bkg_msks is None): 
                stp_msk, bkg_msk = sex(stp).return_mainobj()
                stp_msks[ii] = stp_msk
            else: 
                stp_msk = stp_msks[ii]
                bkg_msk = bkg_msks[ii]
            noises.append( stp[bkg_msk] ) 
            stps[ii][~stp_msk] = 0
        noises = np.hstack(noises) 
            
        #--- 叠加准备的子图（stamp superposition）
        img,     fill= mockimg.fillin(img,     stps,     xc,  yc)
        img_msk, fill= mockimg.fillin(img_msk, stp_msks, xc,  yc) 
        #--- 添加背景噪声（background noise)
        n_fill = np.count_nonzero(fill[0,:] != -1) 
        n_noise = np.count_nonzero(img_msk==0)
        if n_fill != 0: 
            img[img_msk==0] = np.random.choice( noises, n_noise, replace = True)
        return img, img_msk, hdr
    
    def save(self, img, hdr, msk = None): 
        outputdir = self.outputdir 
        band      = self.band 
        #--- save as fits 
        hdu = fits.ImageHDU( img ); 
        hdu.header = hdr
        hdul= fits.HDUList([fits.PrimaryHDU(), hdu])
        hdul.writeto( os.path.join(outputdir, 'mockimg-%s.fits'%band), overwrite =  True)
        hdul.close()
        #--- save as png 
        import imageio
        imageio.imsave( os.path.join(outputdir, 'mockimg-%s.jpeg'%band), img) 
        if not msk is None: 
            hdu = fits.ImageHDU( msk ); 
            hdu.header = hdr 
            hdul= fits.HDUList([fits.PrimaryHDU(), hdu])
            hdul.writeto( os.path.join(outputdir, 'mockmsk-%s.fits'%band), overwrite =  True)
            hdul.close()
            imageio.imsave( os.path.join(outputdir, 'mockmsk-%s.jpeg'%band), msk) 
        return 1

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
    
    def rgb(self, rgbbands):
        from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog
        fitsnames = [os.path.join(self.outputdir, 'mockimg-%s.fits'%band) for band in rgbbands ]
        imgs    = mockimg.readimgs(fitsnames, ihdu = 1)
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
        if self.survey == 'desidr9': 
            image = make_lupton_rgb(0.7*r, 1.1*g, 1.8*b, stretch=0.05, Q=8, minimum = -0.005, filename = None )
        else: 
            image = make_lupton_rgb(r, g, b, stretch=0.05, Q=10, minimum = 0, filename = None )
        image = image/255
        imageio.imsave( os.path.join(self.outputdir, 'mockrgb.jpeg'), image)
        return image

    
def download_from_sdssdr16_worker(ra, dec, outputname, size = [30, 30], pixscale=0.396127, verbose=0, overwrite = False):
    suffix = outputname.split('.')[-1] 
    if  (os.path.exists(outputname))&(overwrite == False) : 
        if(verbose!=0): print('"%s" exists'%outputname) 
        return 1 
    elif (suffix!='jpeg'): 
        print('Suffix of outputname "%s" is not correct. '%outputname)
        print('Suffix should be "jpeg"') 
        return 1
    else:         
        if(verbose!=0): print(ra, dec, outputname, size, pixscale)
        os.system('rm -rf %s'%outputname)
        cmd = 'wget -O %s'% outputname
        web = ' "http://skyserver.sdss.org/dr16/SkyServerWS/ImgCutout/getjpeg?ra=%f&dec=%f&scale=%f&width=%s&height=%s" '%(ra, dec, pixscale, size[0], size[1])
        if(verbose!=0): print( cmd + web )
        return os.system( cmd + web)

def download_from_sdssdr16(ra, dec, outputname, size = [30, 30], pixscale=0.262, verbose=0, overwrite = False, parallel = 1): 
    ra = np.atleast_1d(ra)
    dec= np.atleast_1d(dec)
    outputname=np.atleast_1d(outputname) 
    size = np.atleast_1d(size)
    total = len(ra); 
    if (size.ndim == 1 )&(size.shape[-1] == 2): size = np.repeat(size[np.newaxis,:], total, axis = 0)
    if (size.ndim != 2 )|(size.shape[-1] != 2): print('Error in size; the shape of size is %s'% np.shape(size))
    with parallel_backend("loky") :
        parallel = Parallel(n_jobs=parallel) 
        res      = parallel(delayed(download_from_sdssdr16_worker)(ra[ii], dec[ii], outputname[ii], size = size[ii], pixscale = pixscale, verbose = verbose, overwrite = overwrite) for ii in range(total) ) 
    return res

def download_from_desidr9_worker(ra, dec, outputname, size = [30, 30], pixscale=0.262, verbose=0, overwrite = False):
    suffix = outputname.split('.')[-1] 
    if  (os.path.exists(outputname))&(overwrite == False) : 
        if(verbose!=0): print('"%s" exists'%outputname) 
        return 1 
    elif (suffix!='jpeg')&(suffix!='fits'): 
        print('Suffix of outputname "%s" is not correct. '%outputname)
        print('Suffix should be "jpeg" or "fits"')
        return 1
    else:         
        if(verbose!=0): print(ra, dec, outputname, size, pixscale)
        os.system('rm -rf %s'%outputname)
        cmd = 'wget -O %s'% outputname
        web = ' "https://www.legacysurvey.org/viewer/%s-cutout?ra=%f&dec=%f&height=%s&width=%s&layer=ls-dr9&pixscale=%s&bands=grz" '%(suffix, ra, dec, size[0], size[1], pixscale)
        if(verbose!=0): print( cmd + web )
        return os.system( cmd + web)

def download_from_desidr9(ra, dec, outputname, size = [30, 30], pixscale=0.262, verbose=0, overwrite = False, parallel = 1): 
    ra = np.atleast_1d(ra)
    dec= np.atleast_1d(dec)
    outputname=np.atleast_1d(outputname) 
    size = np.atleast_1d(size)
    total = len(ra); 
    if (size.ndim == 1 )&(size.shape[-1] == 2): size = np.repeat(size[np.newaxis,:], total, axis = 0)
    if (size.ndim != 2 )|(size.shape[-1] != 2): print('Error in size; the shape of size is %s'% np.shape(size))
    # print(np.shape(size)) 
    #for ii in range(total):  
    #    print( ra[ii], dec[ii], outputname[ii], size[ii], pixscale ) 
    with parallel_backend("loky") :
        parallel = Parallel(n_jobs=parallel) 
        res      = parallel(delayed(download_from_desidr9_worker)(ra[ii], dec[ii], outputname[ii], size = size[ii], pixscale = pixscale, verbose = verbose, overwrite = overwrite) for ii in range(total) ) 
    return res


def wcs2reg(w, npt = 1, return_world = True ): 
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
        
