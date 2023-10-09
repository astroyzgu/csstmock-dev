import os 
import time
import numpy as np
from csstmock.dataio import io_jiutian

def test_io_jiutian():
    snapnum  = 127
    filenum  = io_jiutian.fullsky_z2_filenum(snapnum)
    arr      = io_jiutian.fullsky_z2(snapnum, np.arange(filenum)) 
    assert filenum != 0
    assert len(arr)!= 0
