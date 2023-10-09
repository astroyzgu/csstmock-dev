from astropy.table import Table, join, join_skycoord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import KDTree

def quickfa(xyzgal, xyzfib, r = 1.48/60, xyz = True):
    '''
    quick fiberassign 
    parameter: 
    xyzgal: array-like 
        coordinates of galaxy 
    xyzfib: array-like 
        coordinates of  fiber
    r: search radius in degee 
    xyz: if xyz is True:  coordinates means x,y,z.
         if xyz is False: coordinates means r,u,v.
    return: iround, ifiber 
    iround >= 1, round id be assigned. nround == 2 means that the galaxy is assigned in 2nd round. 
    iround =  0, not assigned. 
    iround = -1, galaxies which fiber is impossible to point to.
    
    ifiber >= 0, fiber id at that round be assigned. 
    ifiber = -1, galaxies which are not assigned. 
    '''
    if xyz: 
        pass
    else: 
        x = 1.0*np.cos(xyzgal[:,2]/180.0*np.pi)*np.cos(xyzgal[:,1]/180.0*np.pi)
        y = 1.0*np.cos(xyzgal[:,2]/180.0*np.pi)*np.sin(xyzgal[:,1]/180.0*np.pi)
        z = 1.0*np.sin(xyzgal[:,2]/180.0*np.pi)
        xyzgal = np.vstack([x, y, z]).T
        for iround_, xyzfib_ in enumerate(xyzfib):
            x = 1.0*np.cos(xyzfib_[:,2]/180.0*np.pi)*np.cos(xyzfib_[:,1]/180.0*np.pi)
            y = 1.0*np.cos(xyzfib_[:,2]/180.0*np.pi)*np.sin(xyzfib_[:,1]/180.0*np.pi)
            z = 1.0*np.sin(xyzfib_[:,2]/180.0*np.pi)
            xyzfib[iround_] = np.vstack([x, y, z]).T 

    rgal = np.sqrt( np.sum(xyzgal**2, axis = 1) )
    xyzgal = xyzgal/rgal.reshape(-1,1)
    for iround_, xyzfib_ in enumerate(xyzfib):
        rfib_ = np.sqrt( np.sum(xyzfib_**2, axis = 1) )
        xyzfib[iround_] = xyzfib_/rfib_.reshape(-1,1)
        # print('round:', iround_+1, 'fiber:', np.shape(xyzfib[iround_])) 

    
    ngal   = np.shape(xyzgal)[0]; # 总的星系数量
    nround = len(xyzfib); # 总的轮次
    t1 = Table() 
    t1['igal']   = np.arange(ngal) # 给输入星系编号
    t1['iround'] = - 1             # 在哪一轮观测到？
    t1['ifiber'] = - 1             # 安排了哪一个光纤? 
    t1['ngal_of_fib'] = np.inf # 安排了光纤，它搜索范围内星系的数目
    
    rsearch = 2.0*np.sin(0.5*r/180*np.pi) # >>> 换为直线距离
    
    '''
    找到不在光纤搜索范围内的星系。 
    '''
    
    index_avail  = np.zeros(len(t1)).astype(bool)
    for iround_, xyzfib_ in enumerate(xyzfib):
        kd_fib_        = KDTree(xyzfib_);
        nfib_per_gal_  = kd_fib_.query_ball_point(xyzgal, r = rsearch,  return_length = True)
        index_avail[ nfib_per_gal_ > 0 ] = True 
    t1['iround'][index_avail]  = 0
    if np.sum(~index_avail) == ngal: 
        iround = np.array(t1['iround']).astype(int) 
        ifiber = np.array(t1['ifiber']).astype(int) 
        return iround, ifiber
    
    '''
    按轮次分配光纤给能够指向的星系 
    轮次：iround_+1；该轮光纤的搜索中心：xyzfib_；
    '''
    #print('#--- start ') 
    for iround__, xyzfib__ in enumerate(xyzfib): 
        nfib__   = np.shape(xyzfib__)[0]; # 在这轮中光纤的数目            
        ifib__   = np.arange(nfib__);    # 给这轮中的光纤编号
        #print('#--- New round') 
        for istep in [0, 1]:
            # 选出还没有使用的光纤
            index_fib_   = np.ones(nfib__).astype(bool)
            if istep !=0: 
                fibid_used_  = np.unique( t1['ifiber'][(t1['ifiber']!=-1)&(t1['iround']==iround__+1)].value) 
                index_fib_[fibid_used_]  =  False
            ifib_        = ifib__[index_fib_] 
            xyzfib_      = xyzfib__[index_fib_] 

            # 选出还没有安排光纤的星系
            indx_   = t1['iround']==0; 
            igal_   = np.array( t1['igal'][ t1['iround']==0] ).astype(int) ;
            xyzgal_ = xyzgal[indx_]
            
            #print('round:', iround__+1, 'fiber:', np.shape(xyzfib_)) 
            #print('  igal:', igal_[:10])
            #print('  ifib:', ifib_[:10])
            #print(t1)

            # 查看每一个光纤搜索半径内星系的数目ngal_per_fib
            kd_gal        = KDTree(xyzgal_);
            ngal_per_fib_ = kd_gal.query_ball_point(xyzfib_, r = rsearch, return_length = True )
            igal_per_fib_ = kd_gal.query_ball_point(xyzfib_, r = rsearch, return_sorted = False)
            #--- 完全随机给每一个光纤确定一个星系
            #--- bug: 如果星系同时在两个光纤的搜索半价内， 在同一轮观测中有可能两个光纤都指向了这个星系。
            #--- fix: 加入第二轮 istep == 1 （调整环节） 
            for fibid, galid in enumerate(igal_per_fib_): 
                if len(galid) == 0: continue
                galid = np.array(galid)
                galid_final =  igal_[ np.random.choice(galid, 1)] 
                if ngal_per_fib_[fibid] < t1['ngal_of_fib'][galid_final]: 
                    t1['ifiber'][galid_final] = ifib_[fibid]
                    t1['ngal_of_fib'][galid_final] = ngal_per_fib_[fibid]
                t1['iround'][galid_final] = iround__ + 1
    iround = np.array(t1['iround']).astype(int)
    ifiber = np.array(t1['ifiber']).astype(int)
    return iround, ifiber            
