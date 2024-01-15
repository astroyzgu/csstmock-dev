c*****************************************************************
      SUBROUTINE inputDESI(rg,idrg,Nd1,Nd2max)
c================================================================
c This file input the DESI image legacy surveys data
c with format  i: 1=ra; 2=dec; 3=z; 4=z_error; 5=absmag_z;
c 6=z_max
c================================================================
      include 'common.inc'
c================================================================
      INTEGER Nd1
      INTEGER*8 Nd2max,Nd2sp

      REAL rg(Nd1,Nd2max)
      INTEGER*8 idrg(Nd2max)
      REAL x1,x2,x3,x4,xmag1,xmag2,abm1,abm2,xkcor1,xkcor2
      INTEGER lblnk,itype,error,iselect
      INTEGER*8 id,idb1,idb2,Ndesimax,Nin,iskip,i1,i2
      REAL get_d,getz_max,dist,compl,stellar
      REAL*8 x10,x20,xxx(20)
      EXTERNAL lblnk,get_d,getz_max
      INTEGER*8, allocatable:: id_sp(:),idrg_sp(:)
      REAL*8, allocatable:: rg_sp(:,:)
      INTEGER, allocatable:: ics(:)
      REAL, allocatable:: weight(:),halomass(:)
      Character*3 region

      print*,'reading DESI galaxy data...'
      NDESImax=150000000

      indir='/home/xhyang/work/Gfinder/DESIDR9/data/'
      infile='DESIDR9Y1_NGC_group'
      allocate(halomass(0:NDESImax))
      halomass=0.
      allocate(weight(0:NDESImax))
      weight=0.

      OPEN(66,FILE=indir(1:lblnk(indir))//
     &  infile(1:lblnk(infile)),STATUS='OLD')
        do idb1=1, NDESImax
         READ(66,*,end=41) i1,i2,x1,x2,x3,x4
         weight(idb1)=x4 !! log halo mass
       enddo
41    CLOSE(66)
      Print*,'Input DESIDR9Y1_NGC groups:', idb1-1

      infile='iDESIDR9Y1_NGC_1'
      NDESImax=150000000
      allocate(ics(0:NDESImax))
      ics=-1
      OPEN(66,FILE=indir(1:lblnk(indir))//
     &  infile(1:lblnk(infile)),STATUS='OLD')
        do id=1, NDESImax
         READ(66,*,end=81) idb1,idb2,itype
         ics(idb1)=itype
         halomass(idb1)=weight(idb2)
       enddo
81    CLOSE(66)
      print*,'Read c/s type OK.'

      weight=1.
!      indir='/home/yzgu/SJTU/desi-dr9-latest/seedcat/'
!      infile2='DESIDR9_galaxy_wht2d1.txt'
!      OPEN(68,FILE=indir(1:lblnk(indir))//
!     &  infile2(1:lblnk(infile2)),STATUS='OLD')
!        do id=1, NDESImax
!         READ(68,*,iostat=error) idble,compl
!         if(error>0)GOTO 91
!         if(error<0)GOTO 92
!         weight(idble)=compl
!91      enddo
!92    CLOSE(68)
!      print*,'Read weight OK.'

      indir='/home/yrwang/work_wyr/DESI_DR9/new_catalogue/'
      infile='DESIDR9_galaxy_comb'
      region='NGC'

      skycover=9622.
      V_fact=skycover/(180./pi)**2/3
      call initzbins

ccc input photometric data...

      iskip=-1
      IF(region.eq.'SGC') iskip=72687991

      id=0
      OPEN(66,FILE=indir(1:lblnk(indir))//
     &  infile(1:lblnk(infile)),STATUS='OLD')

        do idb1=0, 72687992+iskip

         READ(66,*,end=121) idble,x10,x20,x3,xxx(1:4),
     & x4,stellar,xkcor2,xmag2,xkcor1,xmag1,iselect

         IF(iselect.le.0) GOTO 101
         IF(idb1.le.iskip) GOTO 101

         IF(x3.ge.z_cut2) GOTO 101
         IF(x3.le.z_cut1) GOTO 101

         dist=get_d(x3)
         abm2=xmag2-5.*alog10(dist)-25.-xkcor2

         IF(abm2.ge.-10.) THEN
!           print*,idb1,abm2,xmag2,xkcor2,x3,dist
           GOTO 101
         ENDIF

         IF(abm2.le.-26.) THEN
!           print*,idb1,abm2,xmag2,xkcor2,x3,dist
           abm2=-26.0
         ENDIF

         id=id+1
         rg(1,id)=x10 !! ra
         rg(2,id)=x20 !! dec
         rg(3,id)=x3 !! redshift
         rg(5,id)=halomass(idble) !!  halo mass 

         compl=weight(idble)
         IF(compl.le.0.05) print*,id,i,compl,xmag2
         compl=max(compl,0.01)

         rg(4,id)=(Msunx-abm2)/2.5        !! z-band luminosity
         rg(6,id)=getz_max(x3,xmag2)     !! max redshift
         rg(7,id)=compl
         rg(8,id)=ics(idble)

         idrg(id)=idble

101   enddo
121   CLOSE(66)
      CLOSE(68)
      Ngalsel=id
      
      Print*,'Input DESI galaxies:', Ngalsel
      Print*,rg(1:8,Ngalsel)

      call rangeLM(rg(5,1:Ngalsel),Ngalsel)


      deallocate(ics)
      deallocate(weight)
      deallocate(halomass)

      RETURN
      END

c*****************************************************************
      SUBROUTINE selectDESI(rg,rgax,Nd1,Nd2max,Nd1s,Nd2s)
c================================================================
c This file selects a few values from the DESI input data
c with format  i:1=redshift bin, 2=log L; 3=v_eff; 4=completeness
c================================================================
      include 'common.inc'
c================================================================
      INTEGER Nd1,Nd1s
      INTEGER*8 Nd2max,Nd2s

      REAL rg(Nd1,Nd2max),rgax(Nd1s,Nd2s)
      REAL x1,x2,x3,veff,zmaxi,zmin,xLmin,xLmax,Vmaxi
      REAL get_r
      EXTERNAL get_r
      INTEGER*8 id1,i1,Nsel
      INTEGER iztype

      id1=0
      do i1=1, Ngalsel
         id1 = id1+1
         x3=rg(3,i1)             !! redshift
         iztype=int((x3-z_cut1)/dzbin)+1
         IF(iztype.gt.Nzbin) iztype=Nzbin
         IF(iztype.lt.1) iztype=1
         rgax(1,id1) = float(iztype)+0.1 !! redshift bin
         rgax(2,id1) = rg(4,i1) !! log L
         zmaxi = rg(6,i1)
         IF(zmaxi.lt.z_max(iztype)) THEN
           Vmaxi=V_fact*(get_r(zmaxi))**3
           veff=Vmaxi-v_min(iztype)
         ELSE
           veff=v_eff(iztype)
         ENDIF

         IF(veff.le.1.0) print*,veff,x3,zmaxi,iztype

         rgax(3,id1) = veff
         rgax(4,id1) = rg(7,i1)
         rgax(5,id1) = rg(8,i1)
101   enddo

      Nsel=id1
     
      Print*,'Select DESI galaxies:', Nsel

      RETURN
      END
      
     
c*****************************************************************
      SUBROUTINE inputDESIgroup(rg,Nd1,Nd2max)
c================================================================
c This file input the DESI galaxy and group data
c with format  i: 1=ra; 2=dec; 3=z; 4=log L 
c 5=log group mass
c================================================================
      include 'common.inc'
c================================================================
      INTEGER Nd1
      INTEGER*8 Nd2max

      REAL rg(Nd1,Nd2max)
      REAL x1,x2,x3,x4,xmag1,xmag2,abm1,abm2,xkcor1,xkcor2
      INTEGER*8 id,i1,i2,Ngrp
      INTEGER lblnk
      REAL get_d,getz_max,dist,xmin
      EXTERNAL lblnk,get_d,getz_max

      print*,'reading DESI galaxy and group data...'
      indir='/home/xhyang/work/Gfinder/DESIDR9/data/'

      infile='DESIDR9_NGC_galaxy'
      OPEN(66,FILE=indir(1:lblnk(indir))//
     &  infile(1:lblnk(infile)),STATUS='OLD')
        id=0
        do i=1, Nd2max
         READ(66,*,end=91) x1,x2,x3,x4,xmag1
         id=id+1
         rg(1,id)=x1 !! ra
         rg(2,id)=x2 !! dec
         rg(3,id)=x3 !! redshift
         rg(4,id)=log10(x4)+10. !! log Luminosity
       enddo
91    CLOSE(66)
      Ngalsel=id
      Print*,'Input DESI galaxies:', Ngalsel

      xmin=100.
      infile='DESIDR9_NGC_group'
      OPEN(66,FILE=indir(1:lblnk(indir))//
     &  infile(1:lblnk(infile)),STATUS='OLD')
        id=0
        do i=1, Nd2max
         READ(66,*,end=141) i1,i2,x1,x2,x3,x4
         id=id+1
         rg(6,id)=x4 !! log halo mass
         xmin=min(xmin,x4)
       enddo
141   CLOSE(66)
      Ngrp=id
      Print*,'Input DESI groups:', Ngrp
      print*, 'Minimum halo mass', xmin

      infile='iDESIDR9_NGC_1'
      OPEN(66,FILE=indir(1:lblnk(indir))//
     &  infile(1:lblnk(infile)),STATUS='OLD')
        id=0
        do i=1, Nd2max
         READ(66,*,end=161) i1,i2
         id=id+1
         rg(5,id)=rg(6,i2) !!! log halo mass
       enddo
161   CLOSE(66)
      Ngalsel=id
      Print*,'Get DESI groups mass:', Ngalsel

      RETURN
      END 


c*****************************************************************
      SUBROUTINE inputmock(rg,Nd1,Nmock)
c================================================================
c This file input mock DESI galaxy data
c with format  i: 1=z; 2=log Mh; 3=absmag 
c================================================================
      include 'common.inc'
c================================================================
      INTEGER Nd1
      INTEGER*8 Nmock

      REAL rg(Nd1,Nmock),xx(10)
      INTEGER*8 id(10),i1
      INTEGER lblnk
      EXTERNAL lblnk

      print*,'reading mock galaxy data ...'
      OPEN(66,FILE=outdir(1:lblnk(outdir))//
     &  outfile(1:lblnk(outfile)),STATUS='OLD')
        do i1=1, Nmock
         READ(66,*,end=121) id(1:3),xx(1:10)
         rg(1,i1)=xx(5) !! redshift
         rg(2,i1)=xx(10) !! log L
         rg(3,i1)=max(xx(7),xx(9)) !! log Mh
         IF(mod(i1,10000000).eq.0) print*, i1
       enddo
121   CLOSE(66)
      Nsimsel=i1-1
      Print*,'Input mock galaxies:', Nsimsel, Nmock

      RETURN
      END


