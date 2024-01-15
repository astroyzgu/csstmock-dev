c*************************************************************
      subroutine getGLF_MF(rGLF,Nd1,Nd2max,Nsel,LMtype)
c================================================================
c  Get the ACCUMULATIVE galaxy luminosity (stellar mass) functions
c  which are used to get the stellar-luminosity relations.
c  rGLF: 1: z_G (redshift) 2: L_G/M_F  (luminosity,
c    or stellar mass)  3: v_eff  
c Writen by: Xiaohu Yang
c================================================================ 
      include 'common.inc'
c================================================================
      integer Nboot,Nd1,LMtype
      INTEGER*8 Nd2max,Nsel,i1,j1,ibt

      real rGLF(Nd1,Nd2max)
      real*8 xlmbin(N_GLF),tinyx
      real*8 acc(N_GLF,Nzbin),ccc(N_GLF,Nzbin),scc(N_GLF,Nzbin)
      real*8 av(N_GLF,Nzbin),cv(N_GLF,Nzbin),sv(N_GLF,Nzbin)
      INTEGER izt,index,kindex,itype
      real v_xx,fact1,compl
      real xLgrp,xtmpa,xtmpc,xtmps
      integer, external:: lblnk
 
      tinyx=1.e-4
      acc=tinyx
      ccc=tinyx
      scc=tinyx
      av=0.
      do i=1,N_GLF
         xlmbin(i)=xlum_1+dxlum*(i-1)
      enddo
      
!!! cccc now count ...

      DO i1=1,Nsel
         xLgrp=rGLF(2,i1)
         izt=int(rGLF(1,i1))
         v_xx=rGLF(3,i1)/1000000.
         compl=rGLF(4,i1)
         itype=int(rGLF(5,i1)+0.1)
         fact1=1./v_xx/compl

         index=int((xLgrp-xlum_1)/dxlum+1.0)  !!! for accumulative +1.0,
         IF(index.le.0 .or. index.gt.N_GLF) GOTO 131
            acc(index,izt)=acc(index,izt)+fact1
            IF(itype.eq.1) THEN
               ccc(index,izt)=ccc(index,izt)+fact1
            ELSE
               scc(index,izt)=scc(index,izt)+fact1
            ENDIF

131      index=int((xLgrp-xlum_1)/dxlum+1.5)  !!! for LF
         IF(index.le.0 .or. index.gt.N_GLF) GOTO 141
            av(index,izt)=av(index,izt)+fact1/dxlum
            IF(itype.eq.1) THEN
               cv(index,izt)=cv(index,izt)+fact1/dxlum
            ELSE
               sv(index,izt)=sv(index,izt)+fact1/dxlum
            ENDIF

141   enddo 


cccc here we try to correct the inpact of low luminosity cut on the GLFa
cccc measurements, by replacing the values using the ones in the lower 
cccc redshift bin.

      IF(LMtype.eq.1) THEN
        do i=2,Nzbin
         do j=1,N_GLF-1
          IF(av(j,i).le.tinyx .or. av(j,i).lt.av(j+1,i))THEN
            acc(j,i)=max(acc(j,i),acc(j,i-1))
            ccc(j,i)=max(ccc(j,i),ccc(j,i-1))
            scc(j,i)=max(scc(j,i),scc(j,i-1))
          ENDIF
         enddo
        enddo
      ENDIF


      do izt=1,Nzbin
        do k=N_GLF-1,1,-1   !!! accumulative 
          acc(k,izt)=acc(k,izt)+acc(k+1,izt)
          ccc(k,izt)=ccc(k,izt)+ccc(k+1,izt)
          scc(k,izt)=scc(k,izt)+scc(k+1,izt)
        enddo
      enddo

      open(3,FILE=outLF(1:lblnk(outLF))
     x   ,status='replace')
      do i=1,N_GLF
         write(3,142) xlmbin(i),(acc(i,j),ccc(i,j),scc(i,j),j=1,Nzbin)
      enddo
      close(3)
142   format(1x,f9.5,24(1x,e11.4))


      open(3,FILE=outLF(1:lblnk(outLF))//'LF'
     x   ,status='replace')
      do i=1,N_GLF
         write(3,152) xlmbin(i),(av(i,j),j=1,Nzbin)
      enddo
      close(3)
152   format(1x,f9.5,8(1x,e11.4))

      IF(LMtype.eq.2) THEN
        do i=1,N_GLF
            GMFa(1,i)=xlmbin(i)
            GMFc(1,i)=xlmbin(i)
            GMFs(1,i)=xlmbin(i)
            GMFa(2:Nzbin+1,i)=acc(i,1:Nzbin)**(-1./3.)*100.  !!! average distance
            GMFc(2:Nzbin+1,i)=ccc(i,1:Nzbin)**(-1./3.)*100.  !!! average
            GMFs(2:Nzbin+1,i)=scc(i,1:Nzbin)**(-1./3.)*100.  !!! average
        enddo
      ENDIF         

      IF(LMtype.eq.1) THEN
        do i=1,N_GLF         
            GLFa(1,i)=xlmbin(i)
            GLFc(1,i)=xlmbin(i)
            GLFs(1,i)=xlmbin(i)
            GLFa(2:Nzbin+1,i)=acc(i,1:Nzbin)**(-1./3.)*100.  !!! average distance
            GLFc(2:Nzbin+1,i)=ccc(i,1:Nzbin)**(-1./3.)*100.  !!! average
            GLFs(2:Nzbin+1,i)=scc(i,1:Nzbin)**(-1./3.)*100.  !!! average
        enddo
      ENDIF

      return
      end
      

c****************************************************************

      subroutine assignL(rg,indexrg,Nd1,Nd2max)
c================================================================
c First, make the abundance matching. Then
c get the lumonisities of galaxies in the whole simulation box.
c Then remove galaxies that are below the minimum luminosity.
c================================================================
      include 'common.inc'
c================================================================
      INTEGER Nd1,LMtype
      INTEGER*8 Nd2max,i1,Nin

      REAL rg(Nd1,Nd2max)
      INTEGER*8 indexrg(3,Nd2max)
      REAL  xMst,xL,z_obs,r,absmag,D_L,kecorr,appmag
      REAl  xMmax,xMmin
      integer, external:: lblnk
      REAL, EXTERNAL:: get_Ld,get_L,ekcor_z05
      REAL, allocatable:: rGLF(:,:)
      INTEGER itype

      Nin=Nsimsel
      call rangeLM(rg(7,1:Nin),Nin)

      Ngaldim1=5
      allocate(rGLF(Ngaldim1,Nsimsel))
      rGLF(1,1:Nsimsel)=rg(9,1:Nsimsel) !! redshift bin
      rGLF(2,1:Nsimsel)=rg(7,1:Nsimsel)
      rGLF(3,1:Nsimsel)=rg(10,1:Nsimsel) !! volume
      rGLF(4,1:Nsimsel)=1.0
      outLF='MF/GMF.dat'
      LMtype=2  !!! get the stellar mass Function=2
      call getGLF_MF(rGLF,Ngaldim1,Nsimsel,Nsimsel,LMtype)
      deallocate(rGLF)

c---  get the luminosity for galaxies in the simulation box

      outLF='MF/MLgax.dat'
      call setMLtable

      jdble=0
      do idble=1,Nsimsel
         r=rg(3,idble)
         z_obs=rg(5,idble)
         xMst=rg(7,idble)

         itype=3
         IF(indexrg(2,idble).eq.0) itype=2 !!! central
         xL=get_L(xMst,z_obs,itype)

         IF(xL.lt.GLFa(1,1)) GOTO 544 

         absmag=Msunx-xL*2.5
         D_L = (1.0 + z_obs) * r
         appmag = absmag + 5.0*ALOG10(D_L) + 25.0
         kecorr = ekcor_z05(z_obs)
         appmag = appmag + kecorr !!-1.62*z_obs

         IF (appmag.GT.amag_cut) GOTO 544
         jdble=jdble+1
         rg(1:8,jdble)=rg(1:8,idble)
         rg(9,jdble)=absmag  !!! absolute magnitude
         indexrg(1:3,jdble)=indexrg(1:3,idble)
544   enddo
      Nsimsel=jdble
      print*,'After luminosity assignment, galaxies remained:', Nsimsel

      RETURN
      END

c****************************************************************
      REAL FUNCTION get_L(xMst,zz,itype)
c================================================================
c Calculate the galaxy luminosity according to its stellar mass
c both in log space.   
c================================================================
      include 'common.inc'
c================================================================
      REAL xMst,zz
      REAL r1,r2,z1,z2,xM,tt1
      INTEGER iztype,jj,itype
      REAL MLtmp(N_GLF)
      
      CALL locate(z_eff,Nzbin,zz,jj)
      IF (jj.LE.0) jj=1
      IF (jj.GE.Nzbin) jj=Nzbin-1
      z1 = z_eff(jj)
      z2 = z_eff(jj+1)
      do i=1,N_GLF
        r1 = ML_gax(itype,jj+1,i)
        r2 = ML_gax(itype,jj+2,i)
        xM = r1 + ((zz-z1)/(z2-z1)) * (r2-r1)
        MLtmp(i)=xM
      enddo

      tt1=xMst
      CALL locate(MLtmp,N_GLF,tt1,jj)

      IF (jj.LE.0) THEN
         get_L=ML_gax(itype,1,1)
         return
      END IF
      
      IF (jj.GE.N_GLF) THEN
         get_L=ML_gax(itype,1,N_GLF)
         return
      END IF

      z1 = MLtmp(jj)
      z2 = MLtmp(jj+1)

      r1 = ML_gax(itype,1,jj)
      r2 = ML_gax(itype,1,jj+1)
      get_L = r1 + ((tt1-z1)/(z2-z1)) * (r2-r1)
   
    
      RETURN
      END      
      
c****************************************************************
      REAL FUNCTION get_Ld(ics,d,iz)
c================================================================
c Calculate the galaxy/group luminosity according to its mean 
c separation (abundance).
c================================================================
      include 'common.inc'
c================================================================
      REAL d
      REAL r1,r2,z1,z2,xM,tt1
      INTEGER iz,ics


      tt1=d

      IF(ics.eq.1) THEN
        CALL locate(GMFa(1+iz,1:N_GLF),N_GLF,tt1,j)
        IF (j.EQ.0) THEN
          get_Ld=GMFa(1,1)
          return
        ENDIF
        IF (j.EQ.N_GLF) THEN
          get_Ld=GMFa(1,N_GLF)
          return
        ENDIF
        r1 = GMFa(1+iz,j)
        r2 = GMFa(1+iz,j+1)
        z1 = GMFa(1,j)
        z2 = GMFa(1,j+1)
        get_Ld = z1 + ((tt1-r1)/(r2-r1)) * (z2-z1)
        return
      ENDIF

      IF(ics.eq.2) THEN
        CALL locate(GMFc(1+iz,1:N_GLF),N_GLF,tt1,j)
        IF (j.EQ.0) THEN
          get_Ld=GMFc(1,1)
          return
        ENDIF
        IF (j.EQ.N_GLF) THEN
          get_Ld=GMFc(1,N_GLF)
          return
        ENDIF
        r1 = GMFc(1+iz,j)
        r2 = GMFc(1+iz,j+1)
        z1 = GMFc(1,j)
        z2 = GMFc(1,j+1)
        get_Ld = z1 + ((tt1-r1)/(r2-r1)) * (z2-z1)
        return
      ENDIF

      IF(ics.eq.3) THEN
        CALL locate(GMFs(1+iz,1:N_GLF),N_GLF,tt1,j)
        IF (j.EQ.0) THEN
          get_Ld=GMFs(1,1)
          return
        ENDIF
        IF (j.EQ.N_GLF) THEN
          get_Ld=GMFs(1,N_GLF)
          return
        ENDIF
        r1 = GMFs(1+iz,j)
        r2 = GMFs(1+iz,j+1)
        z1 = GMFs(1,j)
        z2 = GMFs(1,j+1)
        get_Ld = z1 + ((tt1-r1)/(r2-r1)) * (z2-z1)
        return
      ENDIF

      RETURN
      END      
