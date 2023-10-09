import os 
import matplotlib.pyplot as plt 
import numpy as np 

def mpi4py_dummy(snap, basedir):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    if rank == 0: print('#-------------------')
    if rank == 0: print('# run python function mpi4py_dummy')
    if rank == 0: print('# ')
    if rank == 0: print('snapshot: %s'% snap)
    if rank == 0: print('basedir: %s'% basedir)
    ngal = 100 
    larr = np.arange(ngal, dtype = 'i8') 
    indx = np.array_split( np.arange(ngal), size )[rank]
    
    darr = np.random.uniform( size = (2, len(indx))).astype('f8')
    
    darr   = comm.gather( darr, root = 0) 
    larr   = comm.gather( larr[indx], root = 0) 

    if  rank == 0:
        darr = np.hstack(darr) 
        larr = np.hstack(larr) 
        if 1 == len(np.shape(larr)): larr = larr[:, None].T
        if 1 == len(np.shape(darr)): darr = darr[:, None].T
        
        print( 'shape of darr:', np.shape(darr), darr.dtype )
        print( 'shape of larr:', np.shape(larr), larr.dtype )
        ngal, nd1 = np.shape(darr) 
        ngal, nl1 = np.shape(larr) 
        print( '1st  line of darr:',  darr[ :, 0]  ) 
        print( '2nd  line of darr:',  darr[ :, 1]  ) 
        print( 'last line of darr:',  darr[ :,-1]  ) 
        print('# ')
        print('# end mpi4py dummy')
        print('#-------------------')
        return darr, larr, ngal
    else:
        # print(rank)
        os._exit(0)


if __name__ == '__main__': 
    mpi4py_dummy(snap = 127, basedir = '/home/cossim/Jiutian/M1000/')
