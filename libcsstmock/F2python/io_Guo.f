c*****************************************************************
      SUBROUTINE inputSIM_GuoHong(rg,indexrg,Nd1,Nd2max)
c================================================================
c This file input the mock galaxy data in simulation box.
c with rg format  i: 1-3=x,y,z; 4-6=vx,vy,vz; 7=log Mst, 
c 8=log Mh, indexrg (1: galaxy id, 2: halo id, 3, cen/sat type)
c================================================================
      include 'common.inc'
c================================================================
      INTEGER Nd1
      INTEGER*8 Nd2max

      REAL rg(Nd1,Nd2max)
      INTEGER*8 indexrg(3,Nd2max)
      INTEGER,EXTERNAL:: lblnk
      INTEGER*8 i1,ngal,i0
      REAL xMst,xMh,xMmax,xMmin,zmedian

      integer*8,allocatable:: idhalo(:),ihost(:)
      real,allocatable:: xh(:),yh(:),zh(:),vx(:),vy(:),vz(:)
     & ,halomass(:),halomass_acc(:),stellarmass(:),hostmass(:)


      indir='/home/yczhang/MOCK/'
      infile='gallist_z0p00.binary'
      zmedian=(z_cut1+z_cut2)/2.
      IF(zmedian.gt.0.49) infile='gallist_z0p49binary'
      IF(zmedian.gt.1.03) infile='gallist_z1p03binary'
      IF(zmedian.gt.2.03) infile='gallist_z2p03binary'
      IF(zmedian.gt.3.93) infile='gallist_z3p93binary'

      rLbox=2000.

      OPEN(10,file=indir(1:lblnk(indir))//infile(1:lblnk(infile)),
     & form='unformatted',status='OLD')
      read(10)ngal                           !  number of galaxies

      print*, infile, ngal, 'galaxies...'

      allocate(idhalo(ngal),ihost(ngal),stellarmass(ngal))
      allocate(halomass(ngal),halomass_acc(ngal))
      allocate(xh(ngal),yh(ngal),zh(ngal))
      allocate(vx(ngal),vy(ngal),vz(ngal))

      read(10)idhalo                         !   Subhalo ID
      read(10)ihost                          !   host halo ID, -1: for distinct halos, and others for subhalos
      read(10)halomass                       !   Current Subhalo Mass, Msun/h
      read(10)stellarmass                    !   stellar mass, Msun/h
      read(10)halomass_acc                   !   Subhalo Mass at accretion time, Msun/h
      read(10)xh                             !   x Position of this subhalo, Mpc/h
      read(10)yh                             !   y Position of this subhalo
      read(10)zh                             !   z Position of this subhalo
      read(10)vx                             !   x velocity of this subhalo, km/s
      read(10)vy                             !   y velocity of this subhalo, km/s
      read(10)vz                             !   z velocity of this subhalo, km/s
      close(10)

      print*, 'File reading OK.'
      i1=0
      do i0=1,ngal
        i1=max(idhalo(i0),i1)
      enddo
      print*, ihost(1:10),i1


      infile=infile(1:lblnk(infile))//'.hostmass'
      OPEN(10,file=indir(1:lblnk(indir))//infile(1:lblnk(infile)),
     & form='unformatted',status='OLD')
      read(10)ngal                        !  number of galaxies
      allocate(hostmass(ngal))            ! host halo mass
      read(10) hostmass
      close(10)

      i1=0
      DO i0=1,Ngal
        xMst=stellarmass(i0)
        IF(xMst.le.7.0) GOTO 101
        i1=i1+1
        rg(1,i1)=xh(i0)
        rg(2,i1)=yh(i0)
        rg(3,i1)=zh(i0)
        rg(4,i1)=vx(i0)
        rg(5,i1)=vy(i0)
        rg(6,i1)=vz(i0)
        xMh=hostmass(i0)
        rg(7,i1)= xMst
        rg(8,i1)= xMh
        indexrg(1,i1)=idhalo(i0)
        indexrg(2,i1)=ihost(i0)
101   END DO
  
      Nsimsel = i1

      WRITE(*,*) infile, ' consists of ',Nsimsel,' galaxies.', Ngal
      WRITE(*,*)' '

      deallocate(idhalo,ihost)
      deallocate(stellarmass,halomass,halomass_acc)
      deallocate(xh,yh,zh)
      deallocate(vx,vy,vz)
      deallocate(hostmass)

      RETURN
      END



