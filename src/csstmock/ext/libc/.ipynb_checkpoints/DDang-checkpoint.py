import numpy as np 

def DDang(data, sbin, abin, xyz = False):  # , sbin, angularupweight): 
    '''
    data:  xyz or adr
    lsbin: log of s bins
    labin: log of angular bins 

    return: 
    DD: counts within s and a bins

    '''
    #---- prepare input/output data 
    
    sbin = np.array(sbin.copy() ).astype("float64")
    abin = np.array(abin.copy() ).astype("float64")
    data =np.ascontiguousarray(data.copy() )
    lsbin = np.log10(sbin)
    labin = np.log10(abin)

    if xyz is False:
        #print('input is a, d, r')
        #import Corrfunc
        #DD_counts = Corrfunc.mocks.DDsmu_mocks(1, 1, 12, 1, 1, 10**lsbin,
        #             data[:,0], data[:,1],data[:,2], data[:,2]*0+1, 
        #             is_comoving_dist=True, weight_type = "pair_product")
        #DD_counts = np.array([ DD[4] for DD in DD_counts])
        #print(DD_counts, np.shape(DD_counts), np.sum(DD_counts) )

        x = data[:,2]*np.cos(data[:,1]/180.0*np.pi)*np.cos(data[:,0]/180.0*np.pi)
        y = data[:,2]*np.cos(data[:,1]/180.0*np.pi)*np.sin(data[:,0]/180.0*np.pi)
        z = data[:,2]*np.sin(data[:,1]/180.0*np.pi)
        data = np.vstack([x,y,z]).T 
    data = np.array(data).astype("float64")

    #from sklearn.neighbors import KDTree
    #KDT_D     = KDTree(data)
    #counts_DD = KDT_D.two_point_correlation(data, 10**lsbin)
    #counts_DD = np.diff( counts_DD)
    #counts_DD = np.array(counts_DD).astype('int64')
    #print(counts_DD, np.shape(counts_DD), np.sum(counts_DD) )
    
    r    = np.sqrt( np.sum(data**2, axis = 1) )
    data1= data/r.reshape(-1, 1)  # a, d, 1 in xyz space 
    ngal = np.shape(data)[0]
    dim2 = np.shape(data)[1]

    #---- link to function of cpp  
    import os
    from ctypes import CDLL, POINTER
    from ctypes import c_double, c_int, c_long
    from numpy.ctypeslib import ndpointer
    
    dirname, filename = os.path.split(os.path.abspath(__file__)) 
    libc = CDLL( os.path.join(dirname,'libc/DRang_IW_c.so')) 
    libc.DDang_C.argtypes = ( ndpointer(dtype=np.float64), 
                              ndpointer(dtype=np.float64),
                              c_long, c_long,
                              ndpointer(dtype=np.float64), c_long,
                              ndpointer(dtype=np.float64), c_long,)

    nsbin = len(lsbin) + 1; # -1 <==[0,1,2]==> 3
    nabin = len(labin) + 1; # -1 <==[0,1,2]==> 3
    # libc.DDang_C.restype  =   POINTER(c_double*nsbin*nabin )
    # mat      = np.frombuffer(mat_ptr.contents)
    # mat      = mat.reshape(nsbin, nabin)
    libc.DDang_C.restype  =   ndpointer(dtype=np.float64, shape = (nsbin, nabin) ) 
    #print('c-----')
    mat  = libc.DDang_C( data, data1, ngal, dim2,
                             lsbin, len(lsbin), labin, len(labin) )
    # DD_counts = np.sum(mat[1:-1], axis = 1)
    # print(DD_counts, np.shape(DD_counts), np.sum(DD_counts) )
        
    return mat #.astype("int64")

