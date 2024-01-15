c*****************************************************************
      SUBROUTINE generate_MGRS_9tian(Nd1,Nd2max)
c================================================================
c This file tries to generate the mock galaxy redshift surveys
c from Jiutian simulations. 
c The first step is to read the snapshots in simulation box.
c Next, we generate the lightcone subhalo catalogs.
c Third, we generate mock galaxies in the lightcone.
c 
c output rg format  i: 1-3=x,y,z; 4-6=vx,vy,vz; 7=log Mst, 
c 8=log Mh, indexrg (1: halo id, 2: subhalo id/cs, 3, track ID)
c================================================================
      USE omp_lib
      use, intrinsic :: iso_c_binding !
      include 'common.inc'
c================================================================
      INTEGER Nd1
      INTEGER*8 Nd2max

      REAL z9tian(0:128),xx,compl
      REAL x(3),v(3),xMst,xMh
      REAL zsmin,zsmax,boxsize,z_obs,z_snap
      INTEGER lblnk,isnap_min,isnap_max,isnap,iztype
      INTEGER*8 i1,Nin,Nout,Novlap,i2,i3
      INTEGER N_ass
      EXTERNAL lblnk

      REAL*8, allocatable::sub_ovlap(:,:)
      INTEGER*8,allocatable:: id_ovlap(:,:)
      REAL*8, allocatable:: rg(:,:)
      INTEGER*8,allocatable:: indexrg(:,:)
      real*8, allocatable :: sub_in(:,:)           
      integer*8, allocatable :: idsub_in(:,:)

      skycover=(180./pi)**2*4.*pi  !! all skycov
      call initzbins

      indir='/home/cossim/Jiutian/M1000/simu/'
      infile='zlist_ns_128.txt'
      isnap_min=127
      isnap_max=1
      z9tian(128)=0.
      OPEN(10,file=indir(1:lblnk(indir))//infile(1:lblnk(infile)),
     &     status='OLD')
      do isnap=0,127
        READ(10,*,end=200) xx
        z9tian(isnap)=1./xx -1.
      enddo
      do isnap=0,127
        zsmax=z9tian(isnap)
        if(z_cut1.lt.zsmax) isnap_min=isnap+1
        if(z_cut2.le.zsmax) isnap_max=isnap
      enddo
200   CLOSE(10)
      Print*,'9tian snapshots:',isnap_min,isnap_max,
     & z9tian(isnap_min),z9tian(isnap_max)

      allocate(sub_in(Nd1,Nd2max))
      allocate(sub_ovlap(Nd1,Nd2max))
      allocate(rg(Nd1,Nd2max))
      allocate(indexrg(3,Nd2max))
      allocate(idsub_in(3,Nd2max))
      allocate(id_ovlap(3,Nd2max))
      id_ovlap=0

      boxsize=1000.
      Nin=0
      Nout=0
      Novlap=0

      IF(test_only.eq.'Y') isnap_min=isnap_max !! for test

      Nsimsel=0

        isnap=isnap_min
        call readsnapshot_jiutian(isnap,sub_in,idsub_in,Nd1,Nd2max,Nin)
        WRITE(*,*)' '
        print*,'Snapshot',isnap,Nin
        do i1=Nin,Nin,-1
          write(*,*) idsub_in(1:3,i1), sub_in(1:Nd1,i1)
        enddo
        z_snap=z9tian(isnap)
        call assignL_cube(sub_in,Nd1,Nd2max,idsub_in,Nin,z_snap)

!!! output the galaxy for test----------------
        IF(test_only.eq.'Y')  THEN
          call outputgax(sub_in,Nd1,Nd2max,idsub_in,Nin,z_snap)
        ENDIF
!!! ------------------------------------------

      do isnap=isnap_min-1,isnap_max, -1
         
!!!   store information from previous snapshot...
        sub_ovlap=sub_in
        id_ovlap=idsub_in        
        Novlap=Nin

        zsmin=z9tian(isnap+1)
        zsmax=z9tian(isnap)
        zsmin=max(z_cut1,zsmin)
        zsmax=min(z_cut2,zsmax)
 
        call readsnapshot_jiutian(isnap,sub_in,idsub_in,Nd1,Nd2max,Nin)
        WRITE(*,*)' '
        print*,'Snapshot',isnap,zsmin,zsmax,Nin
        do i1=Nin,Nin,-1
          write(*,*) idsub_in(1:3,i1), sub_in(1:Nd1,i1)
        enddo
        z_snap=z9tian(isnap)
        call assignL_cube(sub_in,Nd1,Nd2max,idsub_in,Nin,z_snap)

        call subhalo_lightcone2(sub_in,idsub_in,Nd1,Nd2max,Nin,
     & zsmin,zsmax,boxsize,sub_ovlap,id_ovlap,Novlap,
     & rg,indexrg,Nout)


c ----- apply mask -----

       call skycov(rg(1,1:Nout),rg(2,1:Nout),
     &  rg(6,1:Nout),'desidr9', Nout, Nout)


       do i1=1,Nout
         z_obs=rg(5,i1)
         compl=rg(6,i1)
         IF(compl.le.0.5) GOTO 222
         IF(z_obs.le.z_cut2 .and. z_obs.ge.z_cut1) THEN
          Nsimsel=Nsimsel+1
          rg(1:Nd1,Nsimsel)=rg(1:Nd1,i1)
          rg(6,Nsimsel)=isnap
          indexrg(1:3,Nsimsel)=indexrg(1:3,i1)
         ENDIF
222    enddo
       print*,'Total subhalos selected in the redshift range:', Nsimsel
       print*, ' '
      enddo

      WRITE(*,*) 'From the snapshots selected ',Nsimsel,' subhalos.'
      WRITE(*,*)' '

      deallocate(sub_in)
      deallocate(sub_ovlap)
      deallocate(idsub_in)
      deallocate(id_ovlap)


      OPEN(13,file=outdir(1:lblnk(outdir))//outfile(1:lblnk(outfile))
     &,     status='replace')
      DO i1=1,Nsimsel
         WRITE(13,133)indexrg(1:3,i1),rg(1:10,i1)
      END DO
 133  FORMAT(I10,1X,I6,1X,I10,1x,2(F12.7,1X),F11.3,7(1X,F9.5))

      deallocate(rg)
      deallocate(indexrg)


      RETURN
      END

     

c****************************************************************

      subroutine outputgax(rg,Nd1,Nd2max,idrg,Nin,redshift)
c================================================================
c Output galaxy in a snapshot for testing purpose... 
c================================================================
      include 'common.inc'
c================================================================
      INTEGER Nd1,LMtype
      INTEGER*8 Nd2max,Nin,i1

      REAL*8 rg(Nd1,Nd2max)
      INTEGER*8 idrg(3,Nd2max)
      REAL redshift

      integer, external:: lblnk  

      OPEN(13,file=outdir(1:lblnk(outdir))//outfile(1:lblnk(outfile))
     & //'_snap',     status='replace')
      write(13,*) redshift
      DO i1=1,Nin
         IF(rg(10,i1).ge.9.0) THEN
           WRITE(13,163)idrg(1:3,i1),rg(7:10,i1),rg(1:6,i1)
         ENDIF
      END DO
 163  FORMAT(I10,1X,I6,1X,I10,1x,4(1X,F9.5),6(1x,f11.4))

      STOP

      RETURN
      END

