import os 
import time
import numpy as np
from csstmock.dataio import io_jiutian

def test_io_jiutian():
    snapnum  = 127
    filenum  = io_jiutian.fullsky_z2_filenum(snapnum)
    arr      = io_jiutian.fullsky_z2(snapnum, np.arange(filenum)) 
    filenum  = io_jiutian.subhalo_filenum(snapnum) 
    arr      = io_jiutian.subhalo(snapnum, [0, 1]) 
    filenum, filelst  = io_jiutian.halo_filenum(snapnum, return_filelist = True) 
    hdr,grp,sub       = io_jiutian.halo(snapnum, [0, 1]) 
    subfind = io_jiutian.subfind_catalog(filelst[1]) 
    #print(subfind.sub_idbm[-5:])
    #print(sub['sub_idbm'][-5:] ) 
    assert filenum != 0
    assert len(arr)!= 0
if __name__ == '__main__': 
    test_io_jiutian() 
    
