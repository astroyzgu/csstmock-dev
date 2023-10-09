import os 
import numpy as np 
import healpy as hp
import glob
from astropy.coordinates import SkyCoord

def skycov_avail():
    LOCATION    = os.path.dirname(os.path.abspath(__file__)) 
    print(LOCATION) 
    search_path = os.path.join(LOCATION, './skycov*.npy')
    path_available = glob.glob(search_path) 
    survey_available = []
    nside_available = []
    for ii in range(len(path_available)):
        filename = path_available[ii]; # print(filename)
        filename = filename.split('/')[-1]
        nside    = filename.split('-')[0].replace('skycov', '')
        survey   = filename.replace('skycov%s-'%nside, '').replace('.npy', '')
        survey_available.append(survey)
        nside_available.append(int(nside))
    survey_available = np.array(survey_available)
    nside_available  = np.array(nside_available)
    return survey_available, nside_available

def skycov_healpy(surveyname):
    ''' 
    Return the vetomap and nside of preset healpy pixelization of some given surveys.

    parameters
    ----------
    surveyname: str 
	name of preset survey, 'desidr9', 'hscdr3', 'csstv0'
	- surveyname == 'desidr9':   18350.26 deg^2, nside = 256, (pixarea = 5.246e-02 deg^2/pixel) 
	- surveyname == 'lsdr9-ngc': , nside = 256, (pixarea = 5.246e-02 deg^2/pixel) 
	- surveyname ==  'hscdr3':   967.39 deg^2, nside = 512, (pixarea = 1.311e-02 deg^2/pixel)
	- surveyname ==  'csstv0':  16787.94 deg^2, nside = 512  (pixarea = 1.311e-02 deg^2/pixel)  
    Returns
    -------
    vetomap: ndarray
	 An array with size = 12*nside*nside. 1.0 means w/i survey; 0.0 means w/o survey.  
    nside: int 
         nside of the healpy pixelization
    ''' 
    survey_available, nside_available = skycov_avail() 
    # print(survey_available)
    # print(nside_available)
    # survey_available =   ['desidr9', 'hscdr3', 'csstv0', 'lsdr9ngc']
    if not surveyname in survey_available: 
        raise ValueError('%s is not available. Only the following is avaible:', survey_available)
    else:
        nside = nside_available[surveyname == survey_available][0]
        LOCATION    = os.path.dirname(os.path.abspath(__file__)) 
        wmap_healpy = os.path.join(LOCATION, './skycov%s-%s.npy'%(nside, surveyname))

    pix      = np.load(wmap_healpy).astype('int64')
#    nside    = wmap_healpy.split('/')[-1].split('-')[0].replace('skycov', '')
#    nside    = int(nside)
    if not nside in 2**np.arange(0, 30):
        raise ValueError( "Running as HEALPY; %s pixels are found."% indx.shape[0] + \
            "However, corresponding nside (%s) must equal to 2^N, N = 1,2,3,4,...,29. \n"% nside + \
            "Thus, corresponding nside must larger than nside >= 2^%s (%s)"%( int( np.log2(nside) ), 2**int( np.log2(nside) ) ) ) 
            
    w = np.arange(12*nside*nside).astype('float64')*0.0
    w[pix] = 1.0
    
    npix = len(pix)
    pixarea = hp.nside2pixarea(nside, degrees= True) 
    print( 'sky coverage is %.2f deg^2 using nside = %s (pixarea = %9.3e deg^2/pixel)'%(pixarea*npix, nside, pixarea) )
    ra, dec = hp.pix2ang(nside, pix, lonlat = True)
    coord   = SkyCoord(ra, dec, unit = 'deg', frame='icrs')  # using degrees directly
    npix = int( np.sum( coord.galactic.b < 0) )
    print( 'South (b<0) coverage = %.2f deg^2'%(pixarea*npix ) )
    npix = int( np.sum( coord.galactic.b > 0) )
    print( 'North (b>0) coverage = %.2f deg^2'%(pixarea*npix ) )
    return w, nside

def assignwht_healpy(x, y, w): 
    w        = np.array(w)
    nside    = np.sqrt( np.shape(w)[0]/12.0 )
    if not nside in 2**np.arange(0, 30): 
            raise ValueError( "Running as HEALPY; %s pixels are found."% indx.shape[0] + \
            "However, corresponding nside (%s) must equal to 2^N, N = 1,2,3,4,...,29. \n"% nside + \
            "Thus, corresponding nside must larger than nside >= 2^%s (%s)"%( int( np.log2(nside) ), 2**int( np.log2(nside) ) ) ) 
    nside = int(nside) 
    ipix = hp.ang2pix(nside, x, y, lonlat = True)
    return w[ipix]

# import matplotlib.pyplot as plt
# if '__main__' == __name__: 
#     x = np.array([1,40])
#     y = np.array([1,40]) 
#     #wht, nside = skycov_healpy('desidr9') 
#     #wht, nside = skycov_healpy('hscdr3') 
#     wht, nside = skycov_healpy('csstv0') 
#     ipix = np.arange(12*nside*nside)
#     ipix = ipix[wht==1]
#     w          = assignwht_healpy(x,y, wht) 
#     print(w)
#     hp.mollview(wht, rot = [118, 0, 0])
#     hp.projscatter(x, y, lonlat = True, )
#     plt.show() 




