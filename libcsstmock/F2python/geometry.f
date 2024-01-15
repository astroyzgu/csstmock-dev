c**********************************************************************
      subroutine geometry(galpar,indexgal,Nd1,Nd2max,Nmock)
c================================================================
c This subroutine is used to select galaxies that are within the 
c magnititude limit and survey geometry.
c================================================================      
      include 'common.inc'
c================================================================
      INTEGER Nd1
      INTEGER*8 Nd2max,Nmock

      REAL galpar(Nd1,Nd2max)
      INTEGER*8 indexgal(3,Nd2max)

      INTEGER ix,iy,iz,N_ass
      INTEGER*8 igal,Nblock,i0,i1,i2,i3
      REAL    x,y,z,vx,vy,vz,r,v,phi,theta
      REAL    ct,st,cp,sp,vr,appmag,compl
      REAL    ra,dec,z_cos,z_obs,kecorr
      REAL      D_L

      REAL*8,allocatable:: raxx(:),decxx(:),amagxx(:),wwxx(:)

      INTEGER,EXTERNAL::  lblnk

      REAL, EXTERNAL::ekcor_z05

      Nblock=Nsimsel
      print*, Nblock
      allocate(raxx(Nblock))
      allocate(decxx(Nblock))
      allocate(amagxx(Nblock))
      allocate(wwxx(Nblock))

      i0=0
      do igal=1,Nsimsel
         i0=i0+1
         raxx(i0)=galpar(1,igal)
         decxx(i0)=galpar(2,igal)
         amagxx(i0)=galpar(10,igal)
      enddo

       call skycov(raxx,decxx,wwxx,
     &  'desidr9', Nblock, Nblock)


      OPEN(13,file=outdir(1:lblnk(outdir))//outfile(1:lblnk(outfile))
     &,     status='replace')

      Nmock=0
      DO i1=1,i0
        galpar(6,i1)=wwxx(i1)
c---check whether galaxy falls within survey limits
        IF(wwxx(i1).ge.0.5) THEN
          Nmock = Nmock + 1
          WRITE(13,133)Nmock,indexgal(1:3,i1),galpar(1:10,i1)
        ENDIF
      END DO
      print*,'The number of galaxies written...',Nmock
c===============================================================

      CLOSE(13)


      WRITE(*,*)' '
      WRITE(*,*)' Total of ',Nmock,' galaxies in mock survey.'

      WRITE(*,*)' '

c---FORMATS

 133  FORMAT(I9,1X,I15,1X,I7,1X,I12,1X,2(F12.7,1X),8(F9.5,1X))

      deallocate(raxx)
      deallocate(decxx)
      deallocate(amagxx)
      deallocate(wwxx)

      RETURN
      END



