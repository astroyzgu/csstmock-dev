c*****************************************************************
      SUBROUTINE inputDESIsv3(rgax,Nd1,Nd2max,Ngal)
c================================================================
c This file input the DESI image legacy surveys data
c with format  i: 1=ra; 2=dec; 3=z; 4=z_error; 5=absmag_z;
c 6=z_max
c================================================================
      include 'common.inc'
c================================================================
      INTEGER Nd1
      INTEGER*8 Nd2max,Ngal

      REAL rgax(Nd1,Nd2max)
      REAL x1,x2,x3,x4,xmag1,xmag2,abm1,abm2,xkcor1,xkcor2
      REAL amag_cut0,xLmin,xLmax,zmaxi,vmaxi,veff
      INTEGER lblnk,itype,iztype,error
      INTEGER*8 id,NDESImax,id1,id2,i1
      REAL get_d,getz_max,get_r,dist,compl
      EXTERNAL lblnk,get_d,getz_max,get_r
      INTEGER, allocatable:: ics(:)
      REAL, allocatable:: weight(:)
      REAL xxx(20)
      REAL*8 total

      amag_cut0=amag_cut  !! copy amag_cut
      amag_cut=19.0   !! use this to calculate LFs.

      print*,'reading DESI sv3 galaxy data...'
      indir='/home/xhyang/work/Gfinder/DESIDR9/data/'
      infile='iDESIDR9v2_NGC_1'
      NDESImax=150000000
      allocate(ics(0:NDESImax))
      ics=-1
      OPEN(66,FILE=indir(1:lblnk(indir))//
     &  infile(1:lblnk(infile)),STATUS='OLD')
        do id=1, NDESImax
         READ(66,*,end=121) id1,id2,itype
         ics(id1)=itype
       enddo
121   CLOSE(66)
      print*,'Read c/s type OK.'

      allocate(weight(0:NDESImax))
      total=0.
      weight=1.
c      indir='/home/yzgu/SJTU/desi-dr9-latest/seedcat/'
c      infile2='DESIDR9_galaxy_wht2d1.txt'
c      OPEN(68,FILE=indir(1:lblnk(indir))//
c     &  infile2(1:lblnk(infile2)),STATUS='OLD')
c        do id=1, NDESImax
c         READ(68,*,iostat=error) idble,compl
c         if(error>0)GOTO 91
c         if(error<0)GOTO 92
c         weight(idble)=compl
c         total=total+1./compl
c91      enddo
c92    CLOSE(68)
c      print*,'Read weight OK.',idble,total

      indir='/home/yrwang/work_wyr/DESI_DR9/new_catalogue/'
      infile='DESIDR9_sv3_galaxyv2'

      total=0.
      skycover=124. !!! 132.
      V_fact=skycover/(180./pi)**2/3.
      call initzbins

      id1=0
      OPEN(68,FILE=indir(1:lblnk(indir))// 
     &  infile(1:lblnk(infile)),STATUS='OLD')
        do i1=1, Nd2max
         READ(68,*,end=141) id2,xxx(1:11)

         x3=xxx(3)
         xmag2=xxx(11)
         xkcor2=xxx(10)
         IF(x3.gt.z_cut2) GOTO 131
         IF(x3.lt.z_cut1) GOTO 131
         IF(xmag2.gt.amag_cut) GOTO 131
         itype=ics(id2)
         IF(itype.eq.-1) GOTO 131

         id1=id1+1

         iztype=int((x3-z_cut1)/dzbin)+1
         IF(iztype.gt.Nzbin) iztype=Nzbin
         IF(iztype.lt.1) iztype=1
         rgax(1,id1) = float(iztype)+0.1 !! redshift bin

         dist=get_d(x3)
         abm2=xmag2-5.*alog10(dist)-25.-xkcor2
         IF(abm2.le.-26.) THEN
!           print*,idb1,abm2,xmag2,xkcor2,x3,dist
           abm2=-26.0
         ENDIF
         rgax(2,id1) = (Msunx-abm2)/2.5 !! log L

         zmaxi = getz_max(x3,xmag2)
         IF(zmaxi.lt.z_max(iztype)) THEN
           Vmaxi=V_fact*(get_r(zmaxi))**3
           veff=Vmaxi-v_min(iztype)
         ELSE
           veff=v_eff(iztype)
         ENDIF
         rgax(3,id1) = veff

         compl=weight(idble)
         IF(compl.le.0.05) print*,id,i,compl,xmag2
         compl=max(compl,0.01)
         rgax(4,id1)=compl
         rgax(5,id1) =itype 
         total=total+1./compl

c         print*,rgax(1:5,id1)

131   enddo
141   close(68)

      Ngal=id1
     
      Print*,'Select DESI galaxies:', Ngal,i1,total

      call rangeLM(rgax(2,1:Ngal),Ngal)

      deallocate(ics)
      deallocate(weight)

      amag_cut=amag_cut0  !! get back amag_cut

      RETURN
      END
      
     


