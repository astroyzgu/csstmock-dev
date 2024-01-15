c****************************************************************

      subroutine assignL_cube(rg,Nd1,Nd2max,idrg,Nin,redshift)
c================================================================
c First, make the abundance matching. Then
c get the lumonisities of galaxies in the whole simulation box.
c Then remove galaxies that are below the minimum luminosity.
c================================================================
      include 'common.inc'
c================================================================
      INTEGER Nd1,LMtype
      INTEGER*8 Nd2max,Nin,i1

      REAL*8 rg(Nd1,Nd2max)
      INTEGER*8 idrg(3,Nd2max)
      REAL redshift,v_factc
      REAL  xMst,xL,z_obs,xMmax,xMmin
      REAL, EXTERNAL:: get_Ld,get_L,ekcor_z05
      REAL, allocatable::rGLF(:,:)
      integer, external:: lblnk  
      real, external:: gasdev
      INTEGER itype
 
      z_obs=redshift 
      Ngaldim1=5
      v_factc=rLbox**3

c-- check if more than 1 sample is generated: to shift and rovate --
      do i=1,isample
        call rotate8(rg,Nd1,Nd2max,Nin)
      enddo

      allocate(rGLF(Ngaldim1,Nin))
      rGLF(1,1:Nin)=1.1  !! redshift bin
      rGLF(2,1:Nin)=rg(9,1:Nin)
      rGLF(3,1:Nin)=v_factc !! volume
      rGLF(4,1:Nin)=1.0
      rGLF(5,1:Nin)=float(idrg(2,1:Nin)+1)
      call rangeLM(rGLF(2,1:Nin),Nin)
      outLF='MF/GMF.dat'
      LMtype=2  !!! get the stellar mass Function=2
      call getGLF_MF(rGLF,Ngaldim1,Nin,Nin,LMtype)

      outLF='MF/MLgax.dat'
      call setMLtable

      do idble=1,Nin
         xMst=rg(9,idble)
         itype=3  !!! satellite
         IF(idrg(2,idble).eq.0) itype=2  !! central
         itype=1 !!! use all relation, will test later
         xL=get_L(xMst,z_obs,itype)+gasdev(iseed)*0.15 
         rg(10,idble)=xL
c         IF(abs(xL).gt.100.) print*,idble,xMst,xL,z_obs,itype
!!! scatter added
      enddo
      print*,'Luminosity assignment in cubic box with scatter.'

!!! ... and rank again... 

      rGLF(2,1:Nin)=rg(10,1:Nin)
      call rangeLM(rGLF(2,1:Nin),Nin)
      outLF='MF/GMF_scatter.dat'
      LMtype=2  !!! get the stellar mass Function=2
      call getGLF_MF(rGLF,Ngaldim1,Nin,Nin,LMtype)

      outLF='MF/MLgax_scatter.dat'
      call setMLtable

      do idble=1,Nin
         xMst=rg(10,idble)
         itype=3  !!! satellite
         IF(idrg(2,idble).eq.0) itype=2  !! central
         itype=1
         rg(10,idble)=get_L(xMst,z_obs,itype)
      enddo
      print*,'Luminosity assignment in cubic box done.'

      deallocate(rGLF)

      RETURN
      END


      SUBROUTINE setMLtable
c================================================================
      include 'common.inc'
c================================================================
      REAL, EXTERNAL:: get_Ld
      integer, external:: lblnk

      do i=1,N_GLF
         ML_gax(1,1,i)=GLFa(1,i)
         ML_gax(2,1,i)=GLFc(1,i)
         ML_gax(3,1,i)=GLFs(1,i)
         do j=1,Nzbin
          ML_gax(1,1+j,i)=get_Ld(1,GLFa(1+j,i),1) !! only use first bin
          ML_gax(2,1+j,i)=get_Ld(2,GLFc(1+j,i),1) !! only use first bin
          ML_gax(3,1+j,i)=get_Ld(3,GLFs(1+j,i),1) !! only use first bin
         enddo
      enddo

      open(3,FILE=outLF(1:lblnk(outLF))
     x   ,status='replace')
      do j=1,3
        do i=1,N_GLF
          write(3,152) ML_gax(j,1:Nzbin+1,i)
        enddo
      enddo
      close(3)
152   format(6(1x,f9.5))

      RETURN 
      END      
