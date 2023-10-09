from ctypes import * 
from numpy.ctypeslib import ndpointer
from time import time 
import numpy as np 

        

def DRang_IW(data, rand, sbin, abin, W_data, W_rand, isDD = True): 
    '''
    data:  xyz or adr
    lsbin: log of s bins
    labin: log of angular bins 

    return: 
    DD: counts within s and a bins

    '''
    #---- prepare input/output data 
    
    sbin = np.array( sbin.copy() ).astype("float64")
    abin = np.array( abin.copy() ).astype("float64")
    lsbin = np.log10(sbin)
    labin = np.log10(abin)
    
    
    if not np.allclose( np.diff(lsbin), (lsbin[1]-lsbin[0])*np.ones(len(lsbin)-1)):
        print('lsbin:', lsbin)
        raise('Error: Inconsistent step spacing of lsbin:')
    if not np.allclose( np.diff(labin), (labin[1]-labin[0])*np.ones(len(labin)-1)): 
        print('labin:', labin)
        raise('Error: Inconsistent step spacing of labin:', labin )
    

    data   =np.ascontiguousarray(data.copy() ).astype("float64")
    rand   =np.ascontiguousarray(rand.copy() ).astype("float64")
    W_data =np.ascontiguousarray(W_data.copy() ).astype("float64")
    W_rand =np.ascontiguousarray(W_rand.copy() ).astype("float64")
    

    dim1 = Ndata = np.shape(data)[0]
    dim2 = Nrand = np.shape(rand)[0]
    
    if isDD: 
        if dim1 != dim2: raise ValueError('isDD is True, but Ndata != Nrand.')
    else: 
        pass 
    
    if np.ndim(W_data) == 0: nw = 1; W_data = data[:,0]*0 + 1.0
    if np.ndim(W_data) == 1: nw = 1; W_data = W_data.reshape(-1, 1)
    if np.ndim(W_data) == 2: nw = np.shape(W_data)[1]
    if np.ndim(W_data) >  2: raise ValueError('Error shape of W_data, ndim > 2')

    if np.ndim(W_rand) == 0: nw = 1; W_rand = rand[:,0]*0 + 1.0
    if np.ndim(W_rand) == 1: nw = 1; W_rand = W_rand.reshape(-1, 1)
    if np.ndim(W_rand) == 2: nw = np.shape(W_rand)[1]
    if np.ndim(W_rand) >  2: raise ValueError('Error shape of W_rand, ndim > 2')

    if isDD: 
        if np.shape(data) != np.shape(rand): 
            raise ValueError('isDD is True, but data and rand is not same shape.')
        else: 
            if not np.allclose(data, rand): 
                raise ValueError('isDD is True, values of data and rand is not same.')
    else: 
        pass 
        

    #---- link to function of cpp
    
    import os
    from ctypes import CDLL, POINTER
    from ctypes import c_double, c_int, c_long
    from numpy.ctypeslib import ndpointer
    
    dirname, filename = os.path.split(os.path.abspath(__file__)) 
    libc = CDLL( os.path.join(dirname,'libc/DRang_IW_c.so')) 
    libc.DRang_IW_c.argtypes = ( 
                              ndpointer(dtype=np.float64),     # // -- xyz1, data # (Ndata, 3)
                              ndpointer(dtype=np.float64),     # // -- xyz2, rand # (Nrand, 3)
                              ndpointer(dtype=np.float64),     # // -- w1, weight of data
                              ndpointer(dtype=np.float64),     # // -- w2, weight of rand
                              c_long,                          # // -- dim1, Ndata 
                              c_long,                          # // -- dim2, Nrand  
                              c_long,                          # // -- dim3, Nw, type number of individual weights 
                              ndpointer(dtype=np.float64),     # // -- lsbin, 
                              c_long,                          # // -- nsbin, 
                              ndpointer(dtype=np.float64),     # // -- labin, 
                              c_long,                          # // -- nabin 
                              c_long
                            ) 
    
    nsbin = len(lsbin) + 1; # -1 <==[0,1,2]==> 3
    nabin = len(labin) + 1; # -1 <==[0,1,2]==> 3
    
    if isDD: 
        isDD = 1
    else: 
        isDD = 0
    libc.DRang_IW_c.restype  =   ndpointer(dtype=np.float64, shape = (nsbin, nabin, nw) ) 
    print('data:', np.shape(data),  'W_data:', np.shape(W_data)) 
    print('rand:', np.shape(rand),  'W_rand:', np.shape(W_rand)) 
    print(dim1, dim2, nw, isDD) 
    mat  = libc.DRang_IW_c(data, rand, W_data, W_rand, 
                         dim1, dim2, nw,
                         lsbin, len(lsbin), labin, len(labin), isDD)
    return mat #.astype("int64")


