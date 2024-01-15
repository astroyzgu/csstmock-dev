c**********************************************************************
      subroutine match_obs(gaxmock,Nd1m,Nd2m,gaxobs,idgax,Nd1o,Nd2o)
c================================================================
c This file use the observed galaxy data to match with the 
c mock galaxy data, and assign each galaxy the cordinates of
c the matched observed galaxy...
c Here we use redshift, halo mass and luminosity of the galaxy
c to make the match... 
c================================================================      
      include 'common.inc'
c================================================================
      INTEGER Nd1m,Nd1o
      INTEGER*8 Nd2m,Nd2o

      REAL gaxmock(Nd1m,Nd2m), gaxobs(Nd1o,Nd2o)
      INTEGER*8 idgax(Nd2o)

      INTEGER*8 Nmock,igal
      REAL    x,y,z,diff0,diff,dx1,dy1,dz1
      INTEGER ix,iy,iz,iy1,iz1,L,izlow,izup,iylow,iyup,irange
      REAL    rxx3_min(3),rxx3_max(3),rxx3_2L(3)

      REAL, allocatable:: rxx3(:,:)
      INTEGER*8,allocatable:: ihoc(:,:,:),ll(:),idmatch(:)

      INTEGER, EXTERNAL::  lblnk

      L=128
      Nmock=Nsimsel

      allocate(ihoc(L,L,L))
      allocate(ll(Ngalsel))
      allocate(rxx3(3,Ngalsel))
      allocate(idmatch(Nmock))

      rxx3_min=100.
      rxx3_max=-100.
      do igal=1,Ngalsel
        do j=1,3
         IF(j.eq.1) x=gaxobs(3,igal) !! redshift
         IF(j.eq.2) x=gaxobs(4,igal) !! log L
         IF(j.eq.3) x=gaxobs(5,igal) !! log Mh
         rxx3(j,igal)=x
         rxx3_min(j)=min(rxx3_min(j),x)
         rxx3_max(j)=max(rxx3_max(j),x)
        enddo
      enddo

      print*,'Min properties', rxx3_min(1:3)
      print*,'Max properties', rxx3_max(1:3)

      do j=1,3
        rxx3_2L(j)=float(L)/(rxx3_max(j)-rxx3_min(j))
      enddo
      print*,'Convert', rxx3_2L(1:3)

      do igal=1,Ngalsel
        do j=1,3
         rxx3(j,igal)=(rxx3(j,igal)-rxx3_min(j))*rxx3_2L(j)
        enddo
      enddo

!      do igal=1,20
!        print*,rxx3(1:3,igal),gaxobs(3:5,igal)
!      enddo

      call linklist(rxx3,ll,Ngalsel,ihoc,L)


!      open(18,file='checklinklist.txt',status='replace')
!      do i=1,L
!        do j=1,L
!          do k=1,L
!            WRITE(18,*) i,j,k, ihoc(i,j,k)
!          enddo
!        enddo
!      enddo
!      close(18)


!$OMP PARALLEL DO default(none)
!$OMP& private(igal,x,y,z,ix,iy,iz,iy1,iz1,j)
!$OMP& private(diff0,diff,dx1,dy1,dz1)             
!$OMP& private(izlow,izup,iylow,iyup,irange)
!$OMP& shared(gaxmock,rxx3_min,rxx3_2L,ihoc,rxx3,L,ll)
!$OMP& shared(gaxobs,Nmock,idmatch,idgax)      
       
      DO igal=1,Nmock

          diff0=float(L)

          x  = (gaxmock(1,igal)-rxx3_min(1))*rxx3_2L(1)
          y  = (gaxmock(2,igal)-rxx3_min(2))*rxx3_2L(2)
          z  = (gaxmock(3,igal)-rxx3_min(3))*rxx3_2L(3)

          ix = int(x)+1
          iy = int(y)+1
          iz = int(z)+1

          ix=min(max(1,ix),L)
          iy=min(max(1,iy),L)
          iz=min(max(1,iz),L)
      
          irange=0
          j=ihoc(ix,iy,iz)
          IF(j.gt.0) GOTO 67

50        irange=irange+1
          iylow=min(max(1,iy-irange),L)
          iyup=min(max(1,iy+irange),L)
          izlow=min(max(1,iz-irange),L)
          izup=min(max(1,iz+irange),L)

          do iz1=izlow,izup
            j=ihoc(ix,iy,iz1) 
            IF(j.gt.0) GOTO 67
          enddo

          do iy1=iylow,iyup
            j=ihoc(ix,iy1,iz)
            IF(j.gt.0) GOTO 67
          enddo

          do iy1=iylow,iyup
            do iz1=izlow,izup
              j=ihoc(ix,iy1,iz1)
              IF(j.gt.0) GOTO 67
            enddo
          enddo

          goto 50
      !!    print*, ix,iy,iz1, gaxmock(1:3,igal), 'check...'

67        if(j.gt.0)then
            dx1=abs(x-rxx3(1,j))
            dy1=abs(y-rxx3(2,j))
            dz1=abs(z-rxx3(3,j))
            diff=dx1+dy1+dz1
            IF(diff.lt.diff0) THEN
              gaxmock(4:8,igal)=gaxobs(1:5,j) !! ra, dec
              idmatch(igal)=idgax(j)
              diff0=diff
            ENDIF 
68         j=ll(j)
           goto 67
          end if
              
      END DO

!$OMP END PARALLEL DO
      
c===============================================================

      print*, 'Writing the matched galaxy information...'
      outfile=outfile(1:lblnk(outfile))//'_obs'
      OPEN(13,file=outdir(1:lblnk(outdir))//outfile(1:lblnk(outfile))
     &,     status='replace')
      DO i=1,Nmock
         WRITE(13,133)idmatch(i),gaxmock(1:3,i),gaxmock(6:8,i)
      END DO
      print*,'The number of galaxies written',Nmock
      CLOSE(13)

c===============================================================

      WRITE(*,*)' '
      WRITE(*,*)' Total of ',Nmock,'galaxies matched with observation.'

      WRITE(*,*)' '

c---FORMATS

 132  FORMAT(i10)
 133  FORMAT(i10,1x,3(F9.5,1X),3(F9.5,1X))

      RETURN
      END

c****************************************************************
      subroutine Linklist(r,ll,nobj,ihoc,L)
      INTEGER*8 nobj
      real r(3,nobj)
      integer L,Ls,Ls2,Ls3
      integer*8 ihoc(L,L,L),ll(nobj)
c
      Ls=L
      Ls2=L*L
      Ls3=Ls2*L
      N=nobj

      do i=1,L
         do j=1,L
            do k=1,L
               ihoc(i,j,k)=0
            enddo
         enddo
      enddo

      do i=1,nobj
      ll(i)=0
      enddo
cc      print*,nobj

      rL=float(L)
      rL1=0.

      do 2 i=1,Nobj
         rrx=r(1,i)+rL1+1.
         rry=r(2,i)+rL1+1.
         rrz=r(3,i)+rL1+1.
         ix= int(rrx)
         iy= int(rry)
         iz= int(rrz)
         if (ix==L+1) ix=L
         if (iy==L+1) iy=L
         if (iz==L+1) iz=L
         inst=ihoc(ix,iy,iz)
         ll(i)=inst
         ihoc(ix,iy,iz)=i
2     continue

      return
      end


