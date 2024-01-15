import numpy as np
import os
import sys

class subfind_catalog:
  '''
  code for reading Subfind's subhalo_tab files
  '''
  def __init__(self, curfile): 
 
    #self.filebase = basedir + "/groups_" + str(snapnum).zfill(3) + "/subhalo_tab_" + str(snapnum).zfill(3) + "."
 
    #print()
    #print("reading subfind catalog for snapshot",snapnum,"of",basedir)
 
    have_veldisp = True
 
    #curfile = self.filebase + str(filenum)
    
    if (not os.path.exists(curfile)):
      print("file not found:", curfile)
      sys.exit()
    
    f = open(curfile,'rb')
    
    ngroups = np.fromfile(f, dtype=np.uint32, count=1)[0]    # Number of groups within this file
    totngroups = np.fromfile(f, dtype=np.uint64, count=1)[0] # Total number of groups for this snapshot.
    nids = np.fromfile(f, dtype=np.uint32, count=1)[0]       
    totnids = np.fromfile(f, dtype=np.uint64, count=1)[0]
    ntask   = np.fromfile(f, dtype=np.uint32, count=1)[0]
    nsubs   = np.fromfile(f, dtype=np.uint32, count=1)[0]
    totnsubs = np.fromfile(f, dtype=np.uint64, count=1)[0]
    
    self.ngroups= ngroups # Number of    groups within this file chunk.
    self.nids   = nids  # ???? print(self.nids)  
    self.nfiles = ntask # Total number of file chunks the group catalog is split between.
    self.nsubs  = nsubs # Number of subgroups within this file chunk.

    self.totngrp = totngroups # Total number of    groups for this snapshot.
    self.totnsub = totnsubs   # Total number of subgroups for this snapshot.


    self.group_len = np.empty(ngroups, dtype=np.uint32)
    self.group_offset = np.empty(ngroups, dtype=np.uint32)
    self.group_nr = np.empty(ngroups, dtype=np.uint64)
    self.group_cm = np.empty(ngroups, dtype=np.dtype((np.float32,3)))
    self.group_vel = np.empty(ngroups, dtype=np.dtype((np.float32,3)))
    self.group_pos = np.empty(ngroups, dtype=np.dtype((np.float32,3)))
    self.group_m_mean200 = np.empty(ngroups, dtype=np.float32)
    self.group_m_crit200 = np.empty(ngroups, dtype=np.float32)
    self.group_m_tophat200 = np.empty(ngroups, dtype=np.float32)
    self.group_veldisp = np.empty(ngroups, dtype=np.float32)
    if have_veldisp:
      self.group_veldisp_mean200 = np.empty(ngroups, dtype=np.float32)
      self.group_veldisp_crit200 = np.empty(ngroups, dtype=np.float32)
      self.group_veldisp_tophat200 = np.empty(ngroups, dtype=np.float32)
    self.group_nsubs = np.empty(ngroups, dtype=np.uint32)
    self.group_firstsub = np.empty(ngroups, dtype=np.uint32)
    
    if nsubs > 0: 
      self.sub_len = np.empty(nsubs, dtype=np.uint32)
      self.sub_offset = np.empty(nsubs, dtype=np.uint32)
      self.sub_grnr = np.empty(nsubs, dtype=np.uint64)
      self.sub_nr = np.empty(nsubs, dtype=np.uint64)
      self.sub_pos = np.empty(nsubs, dtype=np.dtype((np.float32,3)))
      self.sub_vel = np.empty(nsubs, dtype=np.dtype((np.float32,3)))
      self.sub_cm = np.empty(nsubs, dtype=np.dtype((np.float32,3)))
      self.sub_spin = np.empty(nsubs, dtype=np.dtype((np.float32,3)))
      self.sub_veldisp = np.empty(nsubs, dtype=np.float32)
      self.sub_vmax = np.empty(nsubs, dtype=np.float32)
      self.sub_vmaxrad = np.empty(nsubs, dtype=np.float32)
      self.sub_halfmassrad = np.empty(nsubs, dtype=np.float32)
      #self.sub_shape = np.empty(nsubs, dtype=np.dtype((np.float32,6)))
      self.sub_ebind = np.empty(nsubs, dtype=np.float32)
      self.sub_pot = np.empty(nsubs, dtype=np.float32)
      #self.sub_profile = np.empty(nsubs, dtype=np.dtype((np.float32,9)))
      self.sub_parent = np.empty(nsubs, dtype=np.uint32)
      self.sub_idbm = np.empty(nsubs, dtype=np.uint64)
    #--------------------------------------------------------------------
    self.group_len = np.fromfile(f, dtype=np.uint32, count=ngroups)
    # group_len - Integer counter of the total number of particles/cells of all types in this group.
    self.group_offset = np.fromfile(f, dtype=np.uint32, count=ngroups)
    self.group_nr = np.fromfile(f, dtype=np.uint64, count=ngroups)

    self.group_cm = np.fromfile(f, dtype=np.dtype((np.float32,3)), count=ngroups)
    self.group_vel = np.fromfile(f, dtype=np.dtype((np.float32,3)), count=ngroups)
    self.group_pos = np.fromfile(f, dtype=np.dtype((np.float32,3)), count=ngroups)
    # group_cm  (N,3) - Center of mass of the group (mass weighted relative coordinates of all particles/cells in the group)
    # group_vel (N,3) - Velocity of the group, The peculiar velocity is obtained by multiplying this value by 1/a. 
    # group_pos (N,3) - Spatial position within the periodic box (of the particle with the minimum gravitational potential energy)
    self.group_m_mean200 = np.fromfile(f, dtype=np.float32, count=ngroups)
    self.group_m_crit200 = np.fromfile(f, dtype=np.float32, count=ngroups)
    self.group_m_tophat200 = np.fromfile(f, dtype=np.float32, count=ngroups)
    self.group_veldisp = np.fromfile(f, dtype=np.float32, count=ngroups)

    if have_veldisp:
      self.group_veldisp_mean200 = np.fromfile(f, dtype=np.float32, count=ngroups)
      self.group_veldisp_crit200 = np.fromfile(f, dtype=np.float32, count=ngroups)
      self.group_veldisp_tophat200 = np.fromfile(f, dtype=np.float32, count=ngroups)

    self.group_nsubs = np.fromfile(f, dtype=np.uint32, count=ngroups)
    self.group_firstsub = np.fromfile(f, dtype=np.uint32, count=ngroups)

    if nsubs > 0:
      self.sub_len = np.fromfile(f, dtype=np.uint32, count=nsubs)
      self.sub_offset = np.fromfile(f, dtype=np.uint32, count=nsubs)
      self.sub_grnr = np.fromfile(f, dtype=np.uint64, count=nsubs)
      self.sub_nr = np.fromfile(f, dtype=np.uint64, count=nsubs)
      self.sub_pos = np.fromfile(f, dtype=np.dtype((np.float32,3)), count=nsubs)
      self.sub_vel = np.fromfile(f, dtype=np.dtype((np.float32,3)), count=nsubs)
      self.sub_cm = np.fromfile(f, dtype=np.dtype((np.float32,3)), count=nsubs)
      # group_cm  (N,3) - Center of mass of the subgroup
      # group_vel (N,3) - Velocity of the subgroup 
      # group_pos (N,3) - Spatial position within the periodic box 
      self.sub_spin = np.fromfile(f, dtype=np.dtype((np.float32,3)), count=nsubs)
      self.sub_veldisp = np.fromfile(f, dtype=np.float32, count=nsubs)
      # sub_veldisp N   - One-dimensional velocity dispersion of all the member particles/cells (the 3D dispersion divided by sqrt(3) ).
      self.sub_vmax = np.fromfile(f, dtype=np.float32, count=nsubs)
      self.sub_vmaxrad = np.fromfile(f, dtype=np.float32, count=nsubs)
      self.sub_halfmassrad = np.fromfile(f, dtype=np.float32, count=nsubs)
      #self.sub_shape = np.fromfile(f, dtype=np.dtype((np.float32,6)), count=nsubs)
      self.sub_ebind = np.fromfile(f, dtype=np.float32, count=nsubs)
      self.sub_pot = np.fromfile(f, dtype=np.float32, count=nsubs)
      #self.sub_profile = np.fromfile(f, dtype=np.dtype((np.float32,9)), count=nsubs)
      self.sub_parent = np.fromfile(f, dtype=np.uint32, count=nsubs)
      self.sub_idbm = np.fromfile(f, dtype=np.uint64, count=nsubs)

      curpos = f.tell()
      f.seek(0,os.SEEK_END)
      if curpos != f.tell():
        print("Warning: finished reading before EOF for file",filenum)
      f.close()  

def read_groups_allfiles_mpi(basedir, snapnum, blocks): 
    '''
    read halo information from all subfind group files with MPI
    Wirtten by Qingyang Li, Yizhou Gu

    Parameters
    ----------
    basefir: path of base files
    snapnum: snapshot names
    blocks: name of blocks

    Returns
    -------
    DATAALL: datasets of blocks. for example, group_pos = DATAALL['group_pos'][:]
    DATATYPE: data type of blocks
    
    Note also
    ---------
    Host halo and subhalo information can be obtained together.
    '''

    import time 
    import glob 
    import numpy as np
    import os
    from mpi4py import MPI

    start_time  = time.time()                    
    mp = 0.03722953 #dark matter particle mass

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    DATAALL = {}
    DATATYPE = {}
    for block in blocks:
        DATAALL[block] = []
        DATATYPE[block] = []

    curfile_=  basedir +  \
                  "/groups_" + str(snapnum).zfill(3) +  \
                  "/subhalo_tab_" + str(snapnum).zfill(3) + "."

    filenum = np.arange( len( glob.glob(curfile_ + '*') )  ) 
    #filenum = filenum[:10] # for a test
    filenum = np.array_split(filenum, size ) 
    filenum = filenum[rank]
    #-------------------------------------------------------------------
    NSUB = []
    NGROUP = []
    for ii in filenum: 
        if ii%1000 == 0: print('rank %s: reading num. %s'%(rank, ii) )
        curfile = curfile_ + str(ii)
        grp_= subfind_catalog(curfile)
        nh_ = np.shape(grp_.group_nr)[0] #number of groups
        ns_ = np.shape(grp_.sub_nr)[0] #number of subhalos
        NSUB.extend([ns_])
        NGROUP.extend([nh_])

        #group information
        for block in blocks:
            if block == 'group_mass':
                DATAALL[block].extend(np.log10(grp_.group_len*mp)+10)
                DATATYPE[block] = 'float32'
            elif block == 'sub_mass':
                DATAALL[block].extend(np.log10(grp_.sub_len*mp)+10)
                DATATYPE[block] = 'float32'
            else:
                DATAALL[block].extend(getattr(grp_, block))
                DATATYPE[block] = getattr(grp_, block).dtype
        
    for block in blocks:
        DATAALL[block] = np.array(DATAALL[block], dtype = DATATYPE[block])

    NSUB = np.sum(NSUB, dtype = np.int64)
    NGROUP = np.sum(NGROUP, dtype = np.int64)


    #shift in Gatherv and total number of subhalos in each process
    num_subhalo = None
    num_subhalo = comm.gather(NSUB, root=0)
    num_subhalo = comm.bcast(num_subhalo, root=0)
    num_subhalo = np.int64(num_subhalo)
    total_num_subhalo = np.sum(num_subhalo)
    disp_subhalo = np.insert(np.cumsum(num_subhalo), 0, 0)[:-1]

    num_group = None
    num_group = comm.gather(NGROUP, root=0)
    num_group = comm.bcast(num_group, root=0)
    num_group = np.int64(num_group)
    total_num_group = np.sum(num_group)
    disp_group = np.insert(np.cumsum(num_group), 0, 0)[:-1]
    
    blockin3D = ['group_cm', 'group_vel', 'group_pos', 'sub_pos', 'sub_vel', 'sub_cm', 'sub_spin']
    
    #collect dataset
    DATACOLLECT = {}
    comm.Barrier()
    for block in blocks:
        DATATEMP = None

        #groups or subhalos
        if block[:5] == 'group':
            total_num_object = total_num_group
            num_object = num_group
            disp_object = disp_group
        else:
            total_num_object = total_num_subhalo
            num_object = num_subhalo
            disp_object = disp_subhalo

        #3D or 1D
        if rank == 0:
            if block in blockin3D:
                DATATEMP = np.empty((total_num_object, 3), dtype = DATATYPE[block])
            else:
                DATATEMP = np.empty((total_num_object), dtype = DATATYPE[block])
        
        #set MPI data type
        if DATATYPE[block] == 'uint32':
            mpidtype = MPI.UNSIGNED_INT
        elif DATATYPE[block] == 'uint64':
            mpidtype = MPI.UNSIGNED_LONG
        else:
            mpidtype = MPI.FLOAT
        
        #3D or 1D data
        if block in blockin3D:
            comm.Gatherv(DATAALL[block], [DATATEMP, num_object*3, disp_object*3, mpidtype], root = 0)
        else:
            comm.Gatherv(DATAALL[block], [DATATEMP, num_object, disp_object, mpidtype], root = 0)

        if rank == 0:
            DATACOLLECT[block] = DATATEMP 

    # MPI.Finalize()

    end_time = time.time()
#-------------------------------
#    comm.Barrier()
#    MPI.Finalize()
    # print('rank = %s' % rank)
#    if rank != 0 :  os._exit(0)
# added by yizhou
#--------------------------------
    if rank == 0:
        print('running time is: %.3f s' % (end_time - start_time))
        # print(DATACOLLECT['group_pos'][:])
        # print(DATATYPE['group_pos'])
        return DATACOLLECT, DATATYPE


# basedir = "/home/cossim/Jiutian/M1000/groups"
# snapnum = 127 
# blocks = ['group_pos', 'group_mass', 'group_nr', 'sub_vel']
# DATACOLLECT, DATATYPE = read_groups_allfiles_mpi(basedir, snapnum, blocks)


def read_hbt_head(filename):
    '''
    read parameters in hbt files. 
    Cosmology in all files are same at same snapshot.
    '''
    HEAD = {}
    data = h5py.File(filename,'r')
    HEAD['HubbleParam'] = data['Cosmology/HubbleParam'][0]
    HEAD['OmegaLambda0'] = data['Cosmology/OmegaLambda0'][0]
    HEAD['OmegaM0'] = data['Cosmology/OmegaM0'][0]
    HEAD['ScaleFactor'] = data['Cosmology/ScaleFactor'][0]
    HEAD['NumberOfSubhalosInAllFiles'] = data['NumberOfSubhalosInAllFiles'][0]
    HEAD['NumberOfFiles'] = data['NumberOfFiles'][0]

    #we supplyment information here
    HEAD['NumberOfSubhalosThisFile'] = data['Subhalos'].fields('TrackId')[:].shape[0]
    data.close()

    return HEAD


def read_hbt_subhalos(filename, blocks):
    '''
    read subhalo information from one SubSnap file
    Wirtten by Qingyang Li

    Parameters
    ----------
    filename: path of file
    blocks: name of blocks

    Returns
    -------
    Nsub: subhalo number
    DATA: datasets of blocks
    DATATYPE: data type of blocks 
    '''
    
    import numpy as np
    import h5py

    #all names of blocks
    names = ['TrackId', 'Nbound', 'Mbound', 'HostHaloId', 'Rank', 'Depth', 'LastMaxMass', 'SnapshotIndexOfLastMaxMass', 'SnapshotIndexOfLastIsolation', 'SnapshotIndexOfBirth', 'SnapshotIndexOfDeath', 'SnapshotIndexOfSink', 'RmaxComoving', 'VmaxPhysical', 'LastMaxVmaxPhysical', 'SnapshotIndexOfLastMaxVmax', 'R2SigmaComoving', 'RHalfComoving', 'BoundR200CritComoving', 'BoundM200Crit', 'SpecificSelfPotentialEnergy', 'SpecificSelfKineticEnergy', 'SpecificAngularMomentum', 'InertialEigenVector', 'InertialEigenVectorWeighted', 'InertialTensor', 'InertialTensorWeighted', 'ComovingAveragePosition', 'PhysicalAverageVelocity', 'ComovingMostBoundPosition', 'PhysicalMostBoundVelocity', 'MostBoundParticleId', 'SinkTrackId']
    types = ['<i8', '<i8', '<f4', '<i8', '<i8', '<i4', '<f4', '<i4', '<i4', '<i4', '<i4', '<i4', '<f4', '<f4', '<f4', '<i4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<i8', '<i8']

    data = h5py.File(filename,'r')
    nbound= data['Subhalos'].fields('Nbound')[:]
    idx = np.where(nbound > 1)[0]
    Nsub = idx.shape[0]
    
    DATA = {}
    DATATYPE = {}
    for block in blocks:
        readdata = data['Subhalos'].fields(block)[:]
        readdata = readdata[idx] #select Nbound != 1
        nblock = names.index(block)
        datatype = types[nblock]
        # datatype = data['Subhalos'].dtype[nblock]
        # datatype = data['Subhalos'].dtype.descr[nblock][1]
        # print('Read subhalo blocks: %s,' %block + ' dtype:', datatype)

        DATA[block] = np.array(readdata, dtype = datatype)
        DATATYPE[block] = datatype
    data.close()
    return Nsub, DATA, DATATYPE

def read_hbt_allfiles_mpi(basedir, snapnum, blocks):
    '''
    read subhalo information from all SubSnap file with MPI
    Wirtten by Qingyang Li

    Parameters
    ----------
    basefir: path of base files
    snapnum: snapshot names
    blocks: name of blocks

    Returns
    -------
    DATAALL: datasets of blocks. for example, Nbound = DATAALL['Nbound'][:]
    DATATYPE: data type of blocks
    
    Note also
    ---------
    Subhalos with Nbound<=1 are not included. 
    '''
    import glob
    import time
    from mpi4py import MPI
    import numpy as np
    start_time = time.time()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    #initialize data
    DATAALL = {}
    for block in blocks:
        DATAALL[block] = []

    #reading data
    # fileallname = np.sort(glob.glob(basedir + '%s/SubSnap_%s.*.hdf5' %(snapnum,snapnum)))
    nfile = len(glob.glob(basedir + '%s/SubSnap_%s.*.hdf5' %(snapnum,snapnum)))
    fileallname = [basedir + '%s/SubSnap_%s.%s.hdf5' %(snapnum,snapnum,ii) for ii in range(nfile)]
    #fileallname = fileallname[:10] #for a test
    fileallname = np.array_split(fileallname, size)
    NSUB = []
    for filename in fileallname[rank]:
        # print('reading %s ...' %filename)
        Nsub, DATA, DATATYPE = read_hbt_subhalos(filename, blocks)
        for block in blocks:
            DATAALL[block].extend(DATA[block])
        NSUB.extend([Nsub])

    #from list to array
    for block in blocks:
        DATAALL[block] = np.array(DATAALL[block], dtype = DATATYPE[block])
    NSUB = np.sum(NSUB, dtype = np.int64)

    #shift in Gatherv and total number of subhalos in each process
    num_subhalo = None
    num_subhalo = comm.gather(NSUB, root=0)
    num_subhalo = comm.bcast(num_subhalo, root=0)
    num_subhalo = np.int64(num_subhalo)
    total_num = np.sum(num_subhalo)
    disp_subhalo = np.insert(np.cumsum(num_subhalo), 0, 0)[:-1]
    
    blockin3D = ['ComovingAveragePosition', 'PhysicalAverageVelocity', 'ComovingMostBoundPosition', 'PhysicalMostBoundVelocity']
    
    #collect dataset
    DATACOLLECT = {}
    comm.Barrier()
    for block in blocks:
        DATATEMP = None
        if rank == 0:
            if block in blockin3D:
                DATATEMP = np.empty((total_num, 3), dtype = DATATYPE[block])
            else:
                DATATEMP = np.empty((total_num), dtype = DATATYPE[block])
        
        #set MPI data type
        if DATATYPE[block] == '<i4':
            mpidtype = MPI.INT
        elif DATATYPE[block] == '<i8':
            mpidtype = MPI.LONG
        else:
            mpidtype = MPI.FLOAT
        
        #3D or 1D data
        if block in blockin3D:
            comm.Gatherv(DATAALL[block], [DATATEMP, num_subhalo*3, disp_subhalo*3, mpidtype], root = 0)
        else:
            comm.Gatherv(DATAALL[block], [DATATEMP, num_subhalo, disp_subhalo, mpidtype], root = 0)

        if rank == 0:
            DATACOLLECT[block] = DATATEMP

    # MPI.Finalize()

    end_time = time.time()
#-----------------------------
#    comm.Barrier()
#    MPI.Finalize()
    # print('rank = %s' % rank)
#   if rank != 0 :  os._exit(0)
# added by yizhou
#--------------------------------
    if rank == 0:
        print('running time is: %.3f s' % (end_time - start_time))
        # print('running as single core: %.3f s' %(end_time - start_time))
#         print('expected time for all sample: %.3f s' %((end_time - start_time) / 20 *2400))
        
#         print(DATACOLLECT['TrackId'][0], DATACOLLECT['TrackId'][1],  DATACOLLECT['TrackId'][-2],DATACOLLECT['TrackId'][-1])
#         print(DATACOLLECT['ComovingAveragePosition'][0], DATACOLLECT['ComovingAveragePosition'][1], DATACOLLECT['ComovingAveragePosition'][-2], DATACOLLECT['ComovingAveragePosition'][-1])
        return DATACOLLECT, DATATYPE

# a test
# basedir = '/home/cossim/Jiutian/M1000/hbt/'
# snapnum = 127
# blocks = ['TrackId','ComovingAveragePosition']
# DATACOLLECT, DATATYPE = read_hbt_allfiles_mpi(basedir, snapnum, blocks)
    


def collect_input_mockdata(snapnum, basedir_hbt, basedir_groups):
    '''
    collect input data used to construct CSSTmock from Jiutian 

    Parameters:
    ----------
    snapnum: snapshot
    basedir_hbt: file path for hbt data
    basedir_groups: file path for subfind data (FoF)

    Return:
    -------
    INPUT_FORT: include all input data information, see details below

    subhalo information:
    (mass in unit of Msun/h using log values, cordinates in unit of Mpc/h, velocity in unit of km/s),

    1. 'sub_id': subhalo ID 
    2. 'sub_order': subhalo order (rank==1: central halo (most massive); rank ==0: subhalos) 
    3. 'sub_pos': position (x, y, z) 
    4. 'sub_velocity': velocity (vx, vy, vz): 
    5. 'sub_mass': subhalo mass:  
    6. 'sub_accmass': subhalo mass at accretion
    7. 'sub_velmax': maximum circular velocity of subhalo
    8. 'host_id': host halo ID in present snapshot
    9. 'host_mass': host halo mass
    10. 'host_nextid': host halo ID in next snapshot (-99 means not found)
    '''

    #import read_jiutian_hbt_mpi as read_hbt
    #import read_jiutian_groups_mpi as read_groups
    import numpy as np
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    import time

    start_time = time.time()

    #subhalo information from HBT
    
    blocks = ['ComovingAveragePosition', 'PhysicalAverageVelocity', 'LastMaxMass','VmaxPhysical','Mbound', 'Nbound', 'Rank', 'TrackId', 'HostHaloId']
    DATACOLLECT_hbt, DATATYPE_hbt = read_hbt_allfiles_mpi(basedir_hbt, snapnum, blocks)

    #host halo information from subfind
    blocks = ['group_nr', 'group_mass']
    DATACOLLECT_groups, DATATYPE_groups = read_groups_allfiles_mpi(basedir_groups, snapnum, blocks)
    
    blocks = ['group_nr']
    DATACOLLECT_groups_next, DATATYPE_groups_next = read_groups_allfiles_mpi(basedir_groups, snapnum-1, blocks)

    blocks = ['TrackId', 'HostHaloId']
    DATACOLLECT_hbt_next, DATATYPE_hbt_next = read_hbt_allfiles_mpi(basedir_hbt, snapnum-1, blocks)

    INPUT_FORT = {}
    idx_host = None; idx_host_next = None
    next_groups_id = None 
    idxinnext = None
    matchsubid = None; matchhosthaloid = None
    groupsid_next = None; idxin = None
    end_time = None
#-------------------------------
#    comm.Barrier()
#    MPI.Finalize()
#    if rank != 0 :  os._exit(0)
# added by yizhou
#-------------------------------
    if rank == 0:
        print('Finish reading data set')
        # print(np.where(DATACOLLECT_hbt['Nbound'][:] == 0)[0].shape)
        INPUT_FORT['sub_id'] = DATACOLLECT_hbt['TrackId'][:]
        INPUT_FORT['sub_order'] = DATACOLLECT_hbt['Rank'][:]
        INPUT_FORT['sub_pos'] = DATACOLLECT_hbt['ComovingAveragePosition'][:]
        INPUT_FORT['sub_vel'] = DATACOLLECT_hbt['PhysicalAverageVelocity'][:]
        INPUT_FORT['sub_mass'] = np.log10(DATACOLLECT_hbt['Mbound'][:] * 1e10) #need to convert log Mh [Msun/h]
        INPUT_FORT['sub_accmass'] = np.log10(DATACOLLECT_hbt['LastMaxMass'][:] * 1e10)
        INPUT_FORT['sub_velmax'] = DATACOLLECT_hbt['VmaxPhysical'][:]

        idx_host = DATACOLLECT_hbt['HostHaloId'][:]
        INPUT_FORT['host_id'] = DATACOLLECT_groups['group_nr'][idx_host] #hbt gives the index of host halos
        INPUT_FORT['host_mass'] = DATACOLLECT_groups['group_mass'][idx_host] #hbt gives the index of host halos log Mh [Msub/h]

        # #===========================================
        # #determine the host halo id in next snapshot

        idx_host_next = DATACOLLECT_hbt_next['HostHaloId'][:]
        next_groups_id = DATACOLLECT_groups_next['group_nr'][idx_host_next]

        print('finding subhalo id in next snapshot')
        # #take host halo id from next snapshot
        idxinnext = np.in1d(DATACOLLECT_hbt_next['TrackId'][:], DATACOLLECT_hbt['TrackId'][:], assume_unique = True)
        matchsubid = DATACOLLECT_hbt_next['TrackId'][idxinnext]
        matchhosthaloid = next_groups_id[idxinnext]

        print('inserting host halo id in next snapshot to the present subhalos')
        # #find the corresonding id in present halo
        groupsid_next = np.zeros(DATACOLLECT_hbt['TrackId'][:].shape[0], dtype = np.int64)-99
        idxin = np.nonzero(np.in1d(DATACOLLECT_hbt['TrackId'][:], matchsubid, assume_unique = True))[0]
        groupsid_next[idxin] = matchhosthaloid    
        INPUT_FORT['host_nextid'] = groupsid_next

        end_time = time.time()
        print('reading a snapshot time is %.3f s' %(end_time - start_time)) 
        return INPUT_FORT


if __name__ == '__main__':
    snapnum = 127
    basedir_hbt = '/home/cossim/Jiutian/M1000/hbt/'
    basedir_groups = '/home/cossim/Jiutian/M1000/groups/'
    INPUT_FORT = collect_input_mockdata(snapnum, basedir_hbt, basedir_groups)
