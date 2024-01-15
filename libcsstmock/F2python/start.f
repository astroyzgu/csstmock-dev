      PROGRAM makemock
c************************************************************************
c  Construct a mock galaxy redshift survey.
c  written by Xiaohu Yang: xyang@sjtu.edu.cn
c************************************************************************
      use omp_lib      
      include 'common.inc'
c================================================================      

      REAL, allocatable:: gaxsim(:,:),gaxobs(:,:)
      REAL, allocatable:: rGLF(:,:),gaxmock(:,:)
      INTEGER*8, allocatable:: indexgax(:,:),idgax(:)

      INTEGER LMtype,Nd1,Nd2

      REAL, external:: get_r
      REAL x1,Veff
      INTEGER*8 Nsel,Nmock

      INTEGER, EXTERNAL:: lblnk

c==================================================================

      call read_global_param

c--- the relation between comoving distance and redshift on a grid

      call init_comoving_distance

c--- check the parameters:

      xlum_cut=6.0  !! the minimum luminosity. 
      z_cut0=z_cut1
      IF(z_cut2.ge.zzmax) print*, 'initial zzmax too small'
      x1=get_r(z_cut1)
      rLmax=get_r(z_cut2)
      print*,'Redshift range',z_cut1,z_cut2,x1,rLmax


      IF(match_only.eq.'Y') GOTO 999   !!! only make the match

c==================================================================
c--- This part is used to test the I/O and debug the code...
c      call test_io
c==================================================================
c---  get the observed cumulative luminosity function ...

      outLF='MF/GLF.dat'

      IF(CLFfile(1:lblnk(CLFfile)).ne.'N') THEN
        call getGLFfromLFdata
      ELSE

c--- Input the DESI imaging legacy surveys data for calibation...

        IF(obs_name(1:lblnk(obs_name)).eq.'DESIDR9sv3') THEN
          Ngaldim=5
          Ngalmax=5000000
          allocate(rGLF(Ngaldim,Ngalmax))
          call inputDESIsv3(rGLF,Ngaldim,Ngalmax,Nsel)
          LMtype=1  !!! get the Lumonisity Function=1
          call getGLF_MF(rGLF,Ngaldim,Ngalmax,Nsel,LMtype)
          deallocate(rGLF)
        ENDIF

        IF(obs_name(1:lblnk(obs_name)).eq.'DESIDR9') THEN
          Ngaldim=8
          Ngalmax=70000000
          allocate(gaxobs(Ngaldim,Ngalmax))
          allocate(idgax(Ngalmax))
          call inputDESI(gaxobs,idgax,Ngaldim,Ngalmax)
          Ngaldim1=5
          allocate(rGLF(Ngaldim1,Ngalmax))
          call selectDESI(gaxobs,rGLF,Ngaldim,Ngalmax,Ngaldim1,Ngalmax)
          Nsel=Ngalsel
          LMtype=1  !!! get the Lumonisity Function=1
          call getGLF_MF(rGLF,Ngaldim1,Ngalmax,Nsel,LMtype)
          deallocate(rGLF)
        ENDIF

      ENDIF

c==================================================================
c---  input galaxies in the simulation box
     
      sig_vel=35.  !!! velocity error of spectro
 
      IF(sim_name(1:lblnk(sim_name)).eq.'Uchuu') THEN
        Nsimmax=880000000
        allocate(gaxsim(Nsimdim,Nsimmax))
        allocate(indexgax(3,Nsimmax))
        call inputSIM_GuoHong(gaxsim,indexgax,Nsimdim,Nsimmax)
        call generate_cubic(gaxsim,indexgax,Nsimdim,Nsimmax)
        call assignL(gaxsim,indexgax,Nsimdim,Nsimmax)
        call geometry(gaxsim,indexgax,Nsimdim,Nsimmax,Nmock)
        deallocate(gaxsim)
        deallocate(indexgax)
      ENDIF

      IF(sim_name(1:lblnk(sim_name)).eq.'Mill') THEN
        Nsimmax=500000000
        allocate(gaxsim(Nsimdim,Nsimmax))
        allocate(indexgax(3,Nsimmax))
        call inputSIM_Millennium(gaxsim,indexgax,Nsimdim,Nsimmax)
        call generate_cubic(gaxsim,indexgax,Nsimdim,Nsimmax)
        call assignL(gaxsim,indexgax,Nsimdim,Nsimmax)
        call geometry(gaxsim,indexgax,Nsimdim,Nsimmax,Nmock)
        deallocate(gaxsim)
        deallocate(indexgax)
      ENDIF

      IF(sim_name(1:lblnk(sim_name)).eq.'9tian') THEN
        Nsimmax=700000000
        call generate_MGRS_9tian(Nsimdim,Nsimmax)
      ENDIF

c---  match the observed galaxies according to group properties.

999   continue 

      Nmockdim=8
      Nmock=140000000
      allocate(gaxmock(Nmockdim,Nmock))
      call inputmock(gaxmock,Nmockdim,Nmock)

      IF(obs_name(1:lblnk(obs_name)).NE.'DESIDR9') THEN
        Ngaldim=8
        Ngalmax=70000000
        allocate(gaxobs(Ngaldim,Ngalmax))
        allocate(idgax(Ngalmax))
        call inputDESI(gaxobs,idgax,Ngaldim,Ngalmax)
      ENDIF

      call match_obs(gaxmock,Nmockdim,Nmock,gaxobs,idgax,
     & Ngaldim,Ngalmax)

      deallocate(gaxobs)
      deallocate(gaxmock)
      deallocate(idgax)     
 
      STOP
      END

c**********************************************************************

      INTEGER FUNCTION lblnk(char)
c--------------------------------------------------------------------
c Function gives NUMBER of characters in a character variable `char'
c--------------------------------------------------------------------
      character char*(*)
c---
      lblnk=index(char,' ')-1

      RETURN
      END
