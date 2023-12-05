import numpy as np 
class sphrand():
    def generand_range(n, range = [-180, 180, -90, 90], seed = None):
        #
        #----- unifrom distribution on the sphere
        # Lambert azimuthal equal-area projection
        # 从球面上的点（0,0,-1） ==> 投影 ==> z = 0的平面
        # α = arctan(y/x); r = sqrt(x^2+y^2)
        # δ = 2arcsin(r);
        # r<1的区域生成n个随机数
        if seed is not None:  np.random.seed(seed)
        range     = np.atleast_1d(range)
        a_range   = range[0:2]*np.pi/180
        r_range   = np.sin( 0.5*(90 + range[2:4])*np.pi/180 )**2  # (0,1)
        alpha     = np.random.uniform(*a_range, n)
        r         = np.sqrt( np.random.uniform(*r_range, n) )
        delta     = 2*np.arcsin(r) - np.pi/2# arcsin==>(0,90)
        #data      = np.vstack([alpha/np.pi*180, delta/np.pi*180]).T
        #if xyz == True:
        #    x = np.cos(alpha) * np.cos(delta)
        #    y = np.sin(alpha) * np.cos(delta)
        #    z = np.sin(delta)
        #    data  = np.vstack([x,y,z]).T
        return alpha/np.pi*180, delta/np.pi*180
