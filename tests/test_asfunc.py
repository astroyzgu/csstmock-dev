import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import csstmock.asfunc as asfunc 
def test_asfunc(): 
    surveyavail, _  = asfunc.skycov_avail() 
    print(surveyavail)
    assert len(surveyavail) != 0
    vetomap, nside = asfunc.skycov_healpy('desidr9')

    # generate random sample inside the pix 
    pix     = np.arange(12*nside*nside)
    pix_ngc = pix[vetomap>0]
    lon, lat, skycov = asfunc.sphrand_healpy(100, nside, pix_ngc)

    # plot 
    # vetomap[vetomap <= 0] = -np.inf # to make other region blank when plot. 
    # hp.mollview(vetomap, rot = [200, 0, 0], cbar = False)
    # hp.projscatter(lon, lat, s = 1, color = 'r', lonlat = True, label = 'random')
    # for lon in range(0, 360, 30): hp.projtext(lon+6, 0+2, lon, lonlat = True) 

    # hp.graticule()
    # plt.legend()
    # plt.savefig('1.png') 
    # plt.show() 
