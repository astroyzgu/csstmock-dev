from my_plugin import ffi
import numpy as np
import myfunc 
import time 
from numpy import dtype 

ctype2dtype = {'uint8_t': dtype('uint8'), 'int': dtype('int32'), 'long': dtype('int64'), 'double': dtype('float64'), 'int32_t': dtype('int32'), 'float': dtype('float32'), 'int8_t': dtype('int8'), 'int16_t': dtype('int16'), 'uint32_t': dtype('uint32'), 'int64_t': dtype('int64'), 'uint64_t': dtype('uint64'), 'uint16_t': dtype('uint16')}

def asstring(ffi, ptr, length = 256, Nc = 1):
    string_bytes = ffi.string( ptr[0:length*length])
    string_ = str( string_bytes, encoding='UTF-8') 
    string_ = np.array( [string_[ii*length:(ii+1)*length].rstrip() for ii in range(Nc) ] )
    return string_ 
    

def asarray(ffi, ptr, shape, **kwargs):
    length = np.prod(shape)
    # Get the canonical C type of the elements of ptr as a string.
    T = ffi.getctype(ffi.typeof(ptr).item) 
    #print(T)
    if T not in ctype2dtype:
        raise RuntimeError("Cannot create an array for element type: %s" % T) 
    a = np.frombuffer(ffi.buffer(ptr, length * ffi.sizeof(T)), ctype2dtype[T])\
          .reshape(shape, **kwargs) 
    return a
#
# @ffi.def_extern() 
# def box_filenum_jiutian_c(snapnum, filenum): 
#     snapnum = asarray( ffi, snapnum, 1)   
#     filenum = asarray( ffi, filenum, 1)   
#     from csstmock.dataio import io_jiutian
#     filenum[0] = io_jiutian.box_filenum(snapnum[0])
# @ffi.def_extern() 
# def box_jiutian_c(snapnum, ifile, larr, darr, nmax): 
#     snapnum = asarray( ffi, snapnum, 1)   
#     ifile   = asarray( ffi, ifile,   1)   
#     nmax    = asarray( ffi, nmax,    1)   
#     larr    = asarray( ffi, larr, (nmax[0], 5) )
#     darr    = asarray( ffi, darr, (nmax[0],20) )

#------------------ start isin_survey_c
@ffi.def_extern()
def isin_survey_c(ra, dec, n, survey, veto): 
    n         = asarray( ffi, n, 1)   
    ra        = asarray( ffi, ra,  n[0])   
    dec       = asarray( ffi, dec, n[0])   
    veto      = asarray( ffi, veto,n[0]) 
    strsurvey = asstring(ffi, survey, length = 256, Nc = 1 )[0] 
    import csstmock.asfunc as asfunc
    import healpy as hp 
    print(strsurvey, '---') 
    survey_available, nside_available = asfunc.skycov_avail() 
    if strsurvey in survey_available: 
        wht, nside = asfunc.skycov_healpy(strsurvey) 
        ipix       = hp.ang2pix(nside, ra, dec, lonlat = True)
        veto[:] = wht[ipix]
    else: 
        logging.error('%s is not available. /n Available survey is %s'%(strsurvey, survey_available))

#------------------ end isin_survey_c   

#------------------ start fullsky_z2_jiutian_c  
@ffi.def_extern() 
def fullsky_z2_filenum_jiutian_c(snapnum, filenum): 
    snapnum = asarray( ffi, snapnum, 1)   
    filenum = asarray( ffi, filenum, 1)   
    from csstmock.dataio import io_jiutian
    filenum[0] = io_jiutian.fullsky_z2_filenum(snapnum[0])
@ffi.def_extern()
def fullsky_z2_jiutian_c(snapnum, ifile, larr, darr, nmax): 
    snapnum = asarray( ffi, snapnum, 1)   
    ifile   = asarray( ffi, ifile,   1)   
    nmax    = asarray( ffi, nmax,    1)   
    larr    = asarray( ffi, larr, (nmax[0], 5) )
    darr    = asarray( ffi, darr, (nmax[0],20) )
    # >>> start to read data 
    t1  = time.time() 
    from csstmock.dataio import io_jiutian
    arr = io_jiutian.fullsky_z2(snapnum, ifile[0])
    colnames = np.array(list(arr.dtype.names))
    nsub = arr.shape[0]
    col_i8type = ['trackID', 'hosthaloID', 'rank', 'snapNum', 'group_nr']
    idx_i8type = np.isin(colnames, col_i8type) 
    col_f8type = colnames[~idx_i8type]
    # <<< end 
    # >>> assign to the sharing memary 
    nmax[0] = nsub 
    from numpy.lib.recfunctions import structured_to_unstructured
    larr[:nsub,:] = structured_to_unstructured(arr[col_i8type], dtype='i8')
    darr[:nsub,:] = structured_to_unstructured(arr[col_f8type], dtype='f4') 
    t2 = time.time() 
    # print ('fullsky_z2_jiutian_c(snapnum = %s) takes %.2f s'%(snapnum[0], t2 - t1))
#------------------ end fullsky_z2_jiutian_c 


@ffi.def_extern() 
def dummydata_c(larr, darr, carr, Nl1, Nd1, Nc1, N2): 
    Nl1    = asarray( ffi, Nl1, 1)   
    Nd1    = asarray( ffi, Nd1, 1) 
    Nc1    = asarray( ffi, Nc1, 1)
    N2     = asarray( ffi, N2,  1)
    larr   = asarray( ffi, larr, (N2[0],Nl1[0]) )   
    darr   = asarray( ffi, darr, (N2[0],Nd1[0]) ) 
    carr   = asstring( ffi, carr, length = 256, Nc = Nc1[0]) 
    larr[0,:] = np.arange(Nl1[0]) 
    print ('#-------------------------')
    print ('# recive data from fortran')
    print( Nl1, Nd1, Nc1, N2)
    print( 'larr:', larr, larr.dtype )
    print( 'darr:', darr, darr.dtype )
    print( 'carr:', carr )
    print ('# end')
    print ('#-------------------------')
    # input dummy data
    # 
    # N2[:]  = 4;  
    #larr_dummy = np.arange(N2[0]*Nl1[0]).astype('i8').reshape( N2[0], Nl1[0])
    #darr_dummy = np.arange(N2[0]*Nd1[0]).astype('f8').reshape( N2[0], Nd1[0])
    #carr_dummy = np.random.choice(['this', 'is', 'a', 'test'], N2[0], replace=True).reshape(N2[0], Nc1[0]) 
    #larr[:N2[0],:] = larr_dummy[:N2[0], :] 
    #darr[:N2[0],:] = darr_dummy[:N2[0], :] 
    #carr[:N2[0],:] = carr_dummy[:N2[0], :] 

@ffi.def_extern() 
def dummydata_mpi_c(snap, basedir, darr, larr, nd1, nl1, ngal):
    print ('#-------------------------')
    snap    = asarray( ffi, snap, 1)   
    basedir = asstring(ffi, basedir)   
    print ('# recive the snapshot number: %s '% snap[0] )
    print ('# running mpi to load dummy data')
    darr_, larr_, ngal_ = myfunc.mpi4py_dummy( 
                          snap = snap[0], 
                          basedir = basedir[0] ) 
    nd1     = asarray( ffi, nd1, 1)   
    nl1     = asarray( ffi, nl1, 1)   
    ngal    = asarray( ffi,ngal, 1)   
    larr    = asarray(ffi, larr, (nl1[0], ngal[0]) )
    darr    = asarray(ffi, darr, (nd1[0], ngal[0]) )
    nl1[:]  = np.shape(larr_)[0]
    nd1[:]  = np.shape(darr_)[0]
    ngal[:] = ngal_
    darr[:nd1[0], :ngal[0]] = darr_[:, :]
    larr[:nl1[0], :ngal[0]] = larr_[:, :]

@ffi.def_extern() 
def readsnapshot_jiutian_c(isnap, sub_in,idsub_in, Nd1,Nd2max, Nin): 
    print ('#-------------------------')
    isnap    = asarray( ffi, isnap, 1)
    basedir = '/home/cossim/Jiutian/M1000/'
    basedir_hbt    = basedir + '/hbt/' 
    basedir_groups = basedir + '/groups/' 
    
    print ('# recive the snapshot number: %s '% isnap[0] )
    from io_jiutian import readsnapshot
    ngal_, larr_, darr_ = readsnapshot(
                          snapnum = int(isnap[0]), 
                          basedir_hbt = basedir_hbt, 
                          basedir_groups = basedir_groups)
    print(ngal_, darr_.shape, darr_.dtype, larr_.shape, larr_.dtype) 
    Nl1     = [3]
    Nd1     = asarray( ffi, Nd1, 1) # 10 
    Nd2max  = asarray( ffi, Nd2max, 1)
    Nin     = asarray( ffi, Nin, 1) 

    Nd1[0]  = darr_.shape[-1]
    idsub_in  = asarray(ffi, idsub_in, (Nd2max[0], Nl1[0],  ) )
    sub_in    = asarray(ffi,   sub_in, (Nd2max[0], Nd1[0],) )
    Nin[:]  = ngal_ 

    sub_in[:ngal_ ,  :Nd1[0],]   = darr_[:,:]
    idsub_in[:ngal_, 0]    = larr_[:, 0] # host id 
    idsub_in[:ngal_, 1]    = larr_[:, 1] # rank  
    idsub_in[:ngal_, 2]    = larr_[:, 2] # trackID
    #idsub_in[:ngal_, 3]    = larr_[:, 3] # next host id

@ffi.def_extern() 
def readlightcone_jiutian_c( sub_in,idsub_in, Nd1,Nd2max, Nin): 
    print ('#-------------------------')
    from io_jiutian import lightcone4csstmock
    ngal_, larr_, darr_  = lightcone4csstmock()
    Nd1     = asarray( ffi, Nd1, 1) # 10 
    Nd2max  = asarray( ffi, Nd2max, 1)
    Nin     = asarray( ffi, Nin, 1) 

    Nd1[0]  = darr_.shape[-1]
    Nl1     = [2]
    idsub_in  = asarray(ffi, idsub_in, (Nd2max[0], Nl1[0],  ) )
    sub_in    = asarray(ffi,   sub_in, (Nd2max[0], Nd1[0],) )
    Nin[:]  = ngal_ 

    sub_in[:ngal_ ,  :Nd1[0],]   = darr_[:,:]
    idsub_in[:ngal_, 0]    = larr_[:, 0] # haloid
    idsub_in[:ngal_, 1]    = larr_[:, 1] # rank 
    idsub_in[:ngal_, 2]    = larr_[:, 2] # TrackID 

@ffi.def_extern() 
def skycov_c( x, y, w, survey, Ngalmax, Ngal): 
    import csstmock.asfunc as asfunc
    survey   = asstring(ffi, survey)[0]
    Ngalmax  = asarray( ffi, Ngalmax, 1)[0]
    Ngal     = asarray( ffi, Ngal   , 1)[0]
    print ('#-----start python ---')
    print(  survey, "Ngalmax = ", Ngalmax, "Ngal = ", Ngal) 
    x = asarray( ffi, x, (Ngalmax,) )[:]
    y = asarray( ffi, y, (Ngalmax,) )[:]
    w = asarray( ffi, w, (Ngalmax,) )[:]
    print( 'First 6 of x:', x[:6] )
    print( 'First 6 of y:', y[:6] )
    print( 'Only assign the first Ngal (= %s) galaxies. The rest of them are invalid.'% Ngal )
    wht, nside = asfunc.skycov_healpy('desidr9') 
    w[:Ngal]   = asfunc.assignwht_healpy(x[:Ngal],y[:Ngal], wht) 
    print( 'First 6 of w:', w[:6] )
    print ('#-----end  python ---')
