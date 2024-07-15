import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import logging
from astropy.table import Table, join, join_skycoord
from scipy.spatial import KDTree

import csstmock.utils as utils 
def quickfa1pass(ra_tgt, dec_tgt, ra_fib, dec_fib, rsearch, mask_tgt = None, subpriority = None, finetune = False): 
    '''
    quick fiberassign 
    parameter
    ===========
    ra_tgt, dec_tgt: array-like 
        coordinates of targets  
    ra_fib, dec_fib: list of array-like 
        coordinates of pointing centers of robotic fibers. 
    rsearch: float 
        search radius in degee 
    subpriority: array-like
        subpriority of each targets
    mask_tgt: array-like 
        False (or 0), waiting to be assigned; True (or 1), masked; 
    
    return
    =========== 
    ifiber: array-like
        ifiber >= 0, fiber id at that round be assigned. 
        ifiber = -1, galaxies which are not assigned. 
        ifiber = -2, galaxies which are not reachable. 
    '''
    xyztgt = utils.ruv2xyz(1.0, ra_tgt, dec_tgt, unpack=False) 
    xyzfib = utils.ruv2xyz(1.0, ra_fib, dec_fib, unpack=False) 
    sep3d  = 2.0*np.sin(0.5*rsearch/180*np.pi) # >>> 换为直线距离

    ntgt   = np.shape(xyztgt)[0]; # 目标的总数量
    nfib   = np.shape(xyzfib)[0]; # 光纤的总数量
    if mask_tgt is None: 
        mask_tgt = np.zeros(ntgt).astype(bool)
    else: 
        mask_tgt = np.array(mask_tgt).astype(bool)
    itgt   = np.zeros(ntgt) - 2; 
    ifib   = np.zeros(ntgt) - 2; 

    # 查看每一个光纤搜索半径内星系的数目ngal_per_fib和编号itgt_per_fib
    kd            = KDTree(xyztgt)

    #--- ntgt_per_fib = kd.query_ball_point(xyzfib, r = sep3d, return_length = True)
    itgt_per_fib = kd.query_ball_point(xyzfib, r = sep3d, return_sorted = False)
    itgt_avail   = np.unique(np.hstack(itgt_per_fib)).astype(int) 
    itgt[itgt_avail] = -1
    ifib[itgt_avail] = -1

    #--- 为每一个光纤完全随机的指派一个星系
    if subpriority is None: subpriority = np.ones(ntgt, dtype='float32')
    for fibid_, tgtid_ in enumerate(itgt_per_fib):
        tgtid_        = np.array(tgtid_)
        if len(tgtid_)== 0:   continue 
        
        mask_tgt_     = mask_tgt[tgtid_]
        if np.sum(~mask_tgt_)==0:  continue

        tgtid_        = tgtid_[~mask_tgt_] 
        subpriority_  = subpriority[tgtid_]
        tgtid_        = tgtid_[np.max(subpriority_) == subpriority_]
        tgtid_        = np.random.choice(tgtid_, 1)
        itgt[tgtid_]  = tgtid_
        ifib[tgtid_]  = fibid_

    return itgt, ifib
    #--- finetune 查看有没有其他光纤有能力指向这个星系 
    # itgt_assign, time_assign = np.unique(itgt[itgt>=0], return_counts=True) 
    # itgt_finetune = itgt_assign[time_assign >= 2] 
    # indx_finetune = np.isin(itgt, itgt_finetune) 
    # ifib_finetune = ifib[indx_finetune] 
    # ntgt_per_fib[ifib_finetune]
                # #--- 查看有没有其他光纤有能力指向这个星系 
                # if (t1['ngal_of_fib'][galid_final] == -1)|(ngal_per_fib_[fibid] < t1['ngal_of_fib'][galid_final]): 
                #     # 1. 如果这个星系还没有安排上光纤， 那先安排上； 后面在遍历所有光纤时，如果有其他光纤也指向了这个星系，则转入2。 
                #     # 2. 如果当前光纤搜索半径内星系的数目较少， 这优先用这根光纤。
                #     t1['ifiber'][galid_final]      = ifib_[fibid]
                #     t1['ngal_of_fib'][galid_final] = ngal_per_fib_[fibid]
                # t1['iround'][galid_final] = iround__





def quickfa(xyztgt, xyzfib, r = 1.48/60, xyz = True):
    '''
    quick fiberassign 
    parameter: 
    xyztgt: array-like 
        coordinates of galaxy 
    xyzfib: array-like 
        coordinates of  fiber
    r: search radius in degee 
    xyz: if xyz is True:  coordinates means x,y,z.
         if xyz is False: coordinates means r,u,v.
    return: iround, ifiber 

    iround >= 0, round id be assigned. nround == 1 means that the galaxy is assigned in 2nd round. 
    iround = -1, not assigned. 
    iround = -2, galaxies which fiber is impossible to point to.
    
    ifiber >= 0, fiber id at that round be assigned. 
    ifiber = -1, galaxies which are not assigned. 
    ifiber = -2, galaxies which are not assigned. 
    '''
    formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO) 
    ch.setFormatter(formatter)
    # 定义日志记录器，并将以上两者添加到记录器中
    logger = logging.getLogger() # 'mylogger')
    logger.addHandler(ch)

    if xyz: 
        pass
    else: 
        x = 1.0*np.cos(xyztgt[:,2]/180.0*np.pi)*np.cos(xyztgt[:,1]/180.0*np.pi)
        y = 1.0*np.cos(xyztgt[:,2]/180.0*np.pi)*np.sin(xyztgt[:,1]/180.0*np.pi)
        z = 1.0*np.sin(xyztgt[:,2]/180.0*np.pi)
        xyztgt = np.vstack([x, y, z]).T
        for iround_, xyzfib_ in enumerate(xyzfib):
            x = 1.0*np.cos(xyzfib_[:,2]/180.0*np.pi)*np.cos(xyzfib_[:,1]/180.0*np.pi)
            y = 1.0*np.cos(xyzfib_[:,2]/180.0*np.pi)*np.sin(xyzfib_[:,1]/180.0*np.pi)
            z = 1.0*np.sin(xyzfib_[:,2]/180.0*np.pi)
            xyzfib[iround_] = np.vstack([x, y, z]).T



    logger.info('') 
    logger.info( 'Search radius of each fiber positioner is %s deg.'% r ) 
    rsearch = 2.0*np.sin(0.5*r/180*np.pi) # >>> 换为直线距离

    logger.info('') 
    ngal   = np.shape(xyztgt)[0]; # 总的星系数量
    logger.info( 'Normalizing input positions of %s targets.'% ngal ) 
    rgal   = np.sqrt( np.sum(xyztgt**2, axis = 1) )
    logger.info( 'Normalizing positions of fiber positioner.') 
    xyztgt = xyztgt/rgal.reshape(-1,1)
    nround = 0 
    nfiber = []
    for iround_, xyzfib_ in enumerate(xyzfib):
        nfib_ = np.shape(xyzfib[iround_])[0] 
        logger.info( ' - SET %s has %s fiber positioners'%(iround_, nfib_) ) 
        rfib_ = np.sqrt( np.sum(xyzfib_**2, axis = 1) )
        xyzfib[iround_] = xyzfib_/rfib_.reshape(-1,1)
        nround = nround + 1 
        nfiber.append(nfib_) 



    t1 = Table() 
    t1['igal']   = np.arange(ngal) # 给输入星系编号
    t1['iround'] = -2             # 在哪一轮观测到？
    t1['ifiber'] = -2             # 安排了哪一个光纤? 
    t1['ngal_of_fib'] = -1  # 安排了光纤，它搜索范围内星系的数目

    '''
    找到不在光纤搜索范围内的星系。 
    '''
    
    index_avail  = np.zeros(len(t1)).astype(bool)
    for iround_, xyzfib_ in enumerate(xyzfib):
        kd_fib_        = KDTree(xyzfib_);
        nfib_per_gal_  = kd_fib_.query_ball_point(xyztgt, r = rsearch,  return_length = True)
        index_avail[ nfib_per_gal_ > 0 ] = True 
    t1['iround'][index_avail]  = -1
    t1['ifiber'][index_avail]  = -1

    if np.sum(~index_avail) == ngal: 
        iround = np.array(t1['iround']).astype(int) 
        ifiber = np.array(t1['ifiber']).astype(int) 
        return iround, ifiber
    
    '''
    按轮次分配光纤给能够指向的星系 
    轮次：iround___；该轮光纤的搜索中心：xyzfib__；
    '''
    logger.info('') 
    logger.info('') 
    logger.info('#--- starting the processes of fiber assignment')
    #--- bug: 如果星系同时在两个光纤的搜索半价内， 在同一轮观测中有可能两个光纤都指向了这个星系。
    #--- fix: 加入第二轮 istep == 1 （调整环节）, 如果两光纤都指向这个星系，优先安排光纤搜索范围内星系少的那根纤；剩下的光纤重新安排。  
    for iround__, xyzfib__ in enumerate(xyzfib): 
        nfib__   = np.shape(xyzfib__)[0]; # 在这轮中光纤的数目            
        ifib__   = np.arange(nfib__);     # 给这轮中的光纤编号
        logger.info('') 
        logger.info('  %sth set of fibers containing %s'%(iround__, nfib__) ) 
        logger.info('') 
        logger.info('   - running...') 
        #print('#--- New round') 
        for istep in [0, 1]: 
            # 选出还没有使用的光纤
            if istep == 0: index_fib_   = np.ones(nfib__).astype(bool)
            if istep == 1: logger.info('   - fine-tune...') #  %s '% istep ) 

            ifib_        =   ifib__[index_fib_] 
            xyzfib_      = xyzfib__[index_fib_] 

            # 选出还没有安排光纤的星系
            indx_   = t1['iround'] < 0; 
            igal_   = np.array( t1['igal'][indx_] ).astype(int) 
            xyztgt_ = xyztgt[indx_]

            # 查看每一个光纤搜索半径内星系的数目ngal_per_fib
            kd            = KDTree(xyztgt_)
            ngal_per_fib_ = kd.query_ball_point(xyzfib_, r = rsearch, return_length = True )
            igal_per_fib_ = kd.query_ball_point(xyzfib_, r = rsearch, return_sorted = False)

            #--- 为每一个光纤完全随机的指派一个星系
            for fibid, galid in enumerate(igal_per_fib_): 
                if len(galid) == 0: continue
                galid       =  np.array(galid)
                galid_final =  igal_[ np.random.choice(galid, 1)]

                #--- 查看有没有其他光纤有能力指向这个星系 
                if (t1['ngal_of_fib'][galid_final] == -1)|(ngal_per_fib_[fibid] < t1['ngal_of_fib'][galid_final]): 
                    # 1. 如果这个星系还没有安排上光纤， 那先安排上； 后面在遍历所有光纤时，如果有其他光纤也指向了这个星系，则转入2。 
                    # 2. 如果当前光纤搜索半径内星系的数目较少， 这优先用这根光纤。
                    t1['ifiber'][galid_final]      = ifib_[fibid]
                    t1['ngal_of_fib'][galid_final] = ngal_per_fib_[fibid]
                t1['iround'][galid_final] = iround__
            
            # 选出还没有使用的光纤 
            fibid_used_, counts = np.unique( t1['ifiber'][(t1['ifiber'] >= 0)&(t1['iround'] == iround__)].value, return_counts = True)
            if np.sum(counts > 1) != 0:  logger.info('    - %sth set: Repeating fibers appear.' % iround__ )
            index_fib_[fibid_used_]   =  False
            ngal_unassigned = np.sum( t1['ifiber'] < 0 ) 
            logger.info('    - %sth set: Used fibers  %s'%(iround__,  np.sum( ~index_fib_ )  ) )
            logger.info('    - %sth set: Not assigned targets  %s, %.2f precent'%(iround__, ngal_unassigned,  100.0*ngal_unassigned/ngal) )
    iround = np.array(t1['iround']).astype(int)
    ifiber = np.array(t1['ifiber']).astype(int)
    return iround, ifiber            
