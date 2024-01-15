c*****************************************************************
      SUBROUTINE inputSIM_Millennium(rg,indexrg,Nd1,Nd2max)
c================================================================
c This file input the mock galaxy data in simulation box.
c with rg format  i: 1-3=x,y,z; 4-6=vx,vy,vz; 7=log Mst, 
c 8=log Mh, indexrg (1: galaxy id, 2, cen=0, 3: halo id)
c================================================================
      include 'common.inc'
c================================================================
      INTEGER Nd1
      INTEGER*8 Nd2max

      REAL rg(Nd1,Nd2max)
      INTEGER*8 indexrg(3,Nd2max)
      REAL x(3),v(3),xMst,xMh,xMmax,xMmin
      INTEGER lblnk
      INTEGER*8 i1
      EXTERNAL lblnk

      indir='/home/xhyang/work/DESI/MGRS/data/'
      infile='Millennium'

      i1=0

      OPEN(10,file=indir(1:lblnk(indir))//infile(1:lblnk(infile)),
     &     status='OLD')
      DO WHILE (.TRUE.)
         
        i1=i1+1
         
        READ(10,*,end=200) indexrg(1:3,i1),x(1:3),v(1:3)
     &, xMst,xMh

        rg(1:3,i1)=x(1:3)
        rg(4:6,i1)=v(1:3)
        rg(7,i1)= xMst
        rg(8,i1)= xMh
        
        IF (i1.GE.Nsimmax) Print*, 'Too many galaxies'

      END DO
 200  CONTINUE
      CLOSE(10)
  
      Nsimsel = i1-1

      WRITE(*,*) infile, ' consists of ',Nsimsel,' galaxies.'
      WRITE(*,*)' '

      RETURN
      END

     

