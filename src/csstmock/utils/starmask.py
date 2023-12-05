import matplotlib.pyplot as plt 
import numpy as np 
class starmask(): 
    #def __init__(self, star, gala): 
    #    self.star = star
    #    self.gala = gala
    def xyz2ruv(x, y, z):  
        '''
        球坐标-> 笛卡尔坐标
        u0 = 0, ==> u = [-180, 180]; u0 = 180, ==> u = [0, 360]
        v0 = 0, ==> v = [-90, 90];   v0 = 90, ==> v = [0, 180]
        '''
        r = x**2 + y**2 + z**2
        u = np.arctan2(y,x)
        v = np.arctan2(z,np.sqrt(x**2+y**2))
        u = u*180/np.pi; 
        v = v*180/np.pi;
        return r, u, v

    def ruv2xyz(r, u, v):
        x = r*np.cos(v/180.0*np.pi)*np.cos(u/180.0*np.pi)
        y = r*np.cos(v/180.0*np.pi)*np.sin(u/180.0*np.pi)
        z = r*np.sin(v/180.0*np.pi)
        return x, y, z
    
    def uniform_sph(RAlim, DEClim, size):
        """Draw a uniform sample on a sphere

        Parameters
        ----------
        RAlim : tuple
            select Right Ascension between RAlim[0] and RAlim[1]
            units are degrees
        DEClim : tuple
            select Declination between DEClim[0] and DEClim[1]
        size : int (optional)
            the size of the random arrays to return (default = 1)

        Returns
        -------
        RA, DEC : ndarray
            the random sample on the sphere within the given limits.
            arrays have shape equal to size.
        """
        zlim = np.sin(np.pi * np.asarray(DEClim) / 180.)

        z = zlim[0] + (zlim[1] - zlim[0]) * np.random.random(size)
        DEC = (180. / np.pi) * np.arcsin(z)
        RA = RAlim[0] + (RAlim[1] - RAlim[0]) * np.random.random(size)
        return RA, DEC

    def sphrand(n, box = [-180, 180, -90, 90], xyz = False, seed = None): 
        #
        #----- unifrom distribution on the sphere 
        # Lambert azimuthal equal-area projection
        # 从球面上的点（0,0,-1） ==> 投影 ==> z = 0的平面
        # α = arctan(y/x); r = sqrt(x^2+y^2)
        # δ = 2arcsin(r);  
        # r<1的区域生成n个随机数
        if seed is not None:  np.random.seed(seed)
        box     = np.atleast_1d(box)
        a_range = box[0:2]*np.pi/180
        r_range = np.sin( 0.5*(90 + box[2:4])*np.pi/180 )**2  # (0,1) 
        alpha = np.random.uniform(*a_range, n)
        r     = np.sqrt( np.random.uniform(*r_range, n) )
        #x = r*np.cos(alpha) 
        #y = r*np.sin(alpha) 
        #data  = np.vstack([x, y]).T
        delta = 2*np.arcsin(r) - np.pi/2# arcsin==>(0,90)
        data  = np.vstack([alpha/np.pi*180, delta/np.pi*180]).T
        if xyz == True: 
            x = np.cos(alpha) * np.cos(delta)  
            y = np.sin(alpha) * np.cos(delta)  
            z = np.sin(delta) 
            data  = np.vstack([x,y,z]).T
        return data
    @classmethod 
    def crosscorr(cls, star, gala, r = 1): 
        from scipy.spatial import KDTree
        KDs = KDTree(star)
        KDg = KDTree(gala)

        leng = KDg.query_ball_point(star, r = 1/180*np.pi, return_length = True)
        idxs = np.hstack( [[ii]*leng[ii] for ii in range(len(leng)) ] )

        indx = KDg.query_ball_point(star, r = 1/180*np.pi, return_sorted = True) 
        idxg = np.hstack(indx)
        return idxs, idxg 
    
    @classmethod
    def mv2dxyz(cls, vec, from_vec, to_vec):   
        '''
        move points following the rotation from 1 point (from_vec) to another (to_vec)
        '''
        from numpy.linalg import norm
        from_vec = np.atleast_2d(from_vec) 
        to_vec   = np.atleast_2d(to_vec) 
        vec      = np.atleast_2d(vec) 
        from_vec  = from_vec/norm(from_vec, axis = 1)[:, np.newaxis]
        to_vec    = to_vec/norm(to_vec, axis = 1)[:, np.newaxis]
        #norm1  = norm(  from_vec, axis = 1)     
        #norm2  = norm(  to_vec, axis = 1) 
        #print(norm1, norm2) 
        print(from_vec[:2], to_vec[:2])  
        cs        = np.sum(from_vec*to_vec, axis = 1); # print(cos.shape, norm1.shape, norm2.shape)
        vec_cross = np.cross(from_vec,to_vec, axis = 1)
        sn        = norm(vec_cross, axis = 1); 
        #angle     = np.arctan2(sn, cs)
        #print(cs, np.cos(angle) )
        #print(sn, np.sin(angle) )

        vec_cross  = vec_cross/norm(vec_cross, axis = 1)[:, np.newaxis]

        x = vec_cross[:,0]; x0 = vec[:,0]
        y = vec_cross[:,1]; y0 = vec[:,1]
        z = vec_cross[:,2]; z0 = vec[:,2]
        x_ = ( cs+(1-cs)*x**2)*x0 + ( (1-cs)*x*y-sn*z )*y0 + ((1-cs)*x*z+sn*y)*z0
        y_ = ((1-cs)*x*y+sn*z)*x0 + (  cs+(1-cs)*y**2 )*y0 + ((1-cs)*y*z-sn*x)*z0
        z_ = ((1-cs)*x*z-sn*y)*x0 + ( (1-cs)*y*z+sn*x )*y0 + ( cs+(1-cs)*z**2)*z0
        return np.vstack([x_, y_, z_]).T
    @classmethod
    def run(self, gala, star, r = 1):
        idxs, idxg = starmask.crosscorr(star, gala, r = r)
        dxyz       = starmask.mv2dxyz(gala[idxg], from_vec = star[idxs], to_vec = [1,0,0] )
        _, a_, d_  = starmask.xyz2ruv(dxyz[:,0],dxyz[:,1],dxyz[:,2])
        return a_, d_
    
