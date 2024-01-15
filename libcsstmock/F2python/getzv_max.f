            
c*************************************************************
      real function getv_eff(zmin,zmax)
c================================================================
c     Get the effective volume that the galaxy can be observed:
c     1: zmin, 2: zmax   (redshift range)  
c     Writen by: Xiaohu Yang
c================================================================ 
      include 'common.inc'
c================================================================
    
      real zmin,zmax
      real get_r
      external get_r
      real Vmin, Vmax

      IF(zmax.le.zmin) print*,'Check zmin,zmax...'
      Vmin=V_fact*(get_r(zmin))**3      
      Vmax=V_fact*(get_r(zmax))**3
      getv_eff=Vmax-Vmin

      return
      end

      
c*************************************************************
      real function getz_max(zg,apmag)
c================================================================
c     Get the maximum redshift that the galaxy can be observed:
c     1: zg, 2: apmag   (redshift and apparent magnitude)  
c     Writen by: Xiaohu Yang
c================================================================ 
      include 'common.inc'
c================================================================
    
      real zg,apmag
      real get_d,get_zd,ekcor_z05
      external get_d,get_zd,ekcor_z05
      real dmag,dL1,dL2,z1

      dmag=amag_cut-apmag
      dL1=5.0*alog10(get_d(zg))  !! distance module
      dL2=10.**((dmag+dL1)/5.)   !! first try of the D_L
      z1=get_zd(dL2)

      dL1=5.0*alog10(get_d(zg))+ekcor_z05(zg)-ekcor_z05(z1)
      dL2=10.**((dmag+dL1)/5.)   !! second try of the D_L
      z1=get_zd(dL2)      

      dL1=5.0*alog10(get_d(zg))+ekcor_z05(zg)-ekcor_z05(z1)
      dL2=10.**((dmag+dL1)/5.)   !! third try of the D_L
      z1=get_zd(dL2)
      
      getz_max=z1
      

      return
      end


c****************************************************************

      REAL FUNCTION get_zd(r)
c================================================================
c  Compute the redshift z from a luminosity distance d (in Mpc/h)
c  We use linear interpolation on the vectors ddd and zzz which we
c  have initialized with the subroutine init_comoving_distance.
c================================================================
      include 'common.inc'
c================================================================

      REAL     r1,r2,z1,z2,r
c---
      CALL locate(ddd,Nz,r,j)

      IF (j.EQ.0) THEN
         get_zd=zzz(1)
         return
      END IF

      IF (j.EQ.Nz) THEN
         get_zd=zzz(Nz)
         return
      END IF
      
      r1 = ddd(j)
      r2 = ddd(j+1)

      z1 = zzz(j)
      z2 = zzz(j+1)

      get_zd = z1 + ((r-r1)/(r2-r1)) * (z2-z1)

      RETURN
      END

c****************************************************************

      REAL FUNCTION get_d(z)
c================================================================
c  Compute the luminosity distance d (in Mpc/h) for a given z
c  We use linear interpolation on the vectors ddd and zzz which we
c  have initialized with the subroutine init_comoving_distance.
c================================================================
      include 'common.inc'
c================================================================

      REAL  r1,r2,z1,z2,z
      
c---
      CALL locate(zzz,Nz,z,j)

      IF (j.EQ.0) THEN
        get_d=ddd(1)
        return
      ENDIF

      IF (j.EQ.Nz) THEN
        get_d=ddd(Nz)
        return
      END IF

      r1 = ddd(j)
      r2 = ddd(j+1)

      z1 = zzz(j)
      z2 = zzz(j+1)

      get_d = r1 + ((z-z1)/(z2-z1)) * (r2-r1)

      END

c****************************************************************

      REAL FUNCTION get_r(z)
c================================================================
c  Compute the comoving distance d (in Mpc/h) for a given z
c  We use linear interpolation on the vectors rrr and zzz which we
c  have initialized with the subroutine init_comoving_distance.
c================================================================
      include 'common.inc'
c================================================================

      REAL     r1,r2,z1,z2,z
c---
      CALL locate(zzz,Nz,z,j)

      IF (j.EQ.0) THEN
        get_r=rrr(1)
        return
      ENDIF

      IF (j.EQ.Nz) THEN
        get_r=rrr(Nz)
        return        
      END IF

      r1 = rrr(j)
      r2 = rrr(j+1)

      z1 = zzz(j)
      z2 = zzz(j+1)

      get_r = r1 + ((z-z1)/(z2-z1)) * (r2-r1)

      END

      

c****************************************************************
      
      real function ekcor_z05(z)
        xkcor_r=-0.33-0.65*z+0.88*z*z
        ekcor_z05=xkcor_r  !!+ecorr_r
      return
      end


c**********************************************************************

      REAL FUNCTION get_z(r)
c---------------------------------------------------------------------- 
c  Compute the redshift for a given comoving distance r (in Mpc/h)
c  We use linear interpolation on the vectors rrr and zzz which we
c  have initialized with the subroutine init_comoving_distance.              
c----------------------------------------------------------------------
      include 'common.inc'
c----------------------------------------------------------------------

      REAL     r,r1,r2,z1,z2

c---

      CALL locate(rrr,Nz,r,j)

      IF (j.EQ.0) THEN
        get_z=zzz(1)
        return
      ENDIF

      IF (j.EQ.Nz) THEN
        get_z=zzz(Nz)
        return
      END IF
      
      r1 = rrr(j)
      r2 = rrr(j+1)

      z1 = zzz(j)
      z2 = zzz(j+1)

      get_z = z1 + ((r-r1)/(r2-r1)) * (z2-z1)

      END

c**********************************************************************

      SUBROUTINE init_comoving_distance
c---------------------------------------------------------------------- 
c  Compute the comoving distances on a grid of redshift.
c----------------------------------------------------------------------
      include 'common.inc'
c----------------------------------------------------------------------

      INTEGER  iz
      REAL     SS1

c---

      DO iz = 1,Nz
        zzz(iz) = FLOAT(iz-1)/FLOAT(Nz-1) * zzmax
        IF (iz.EQ.1) THEN
           rrr(iz) = 0.0
           ddd(iz) = 0.0
        ELSE
          CALL QROMO_distance(0.0,zzz(iz),SS1)
          rrr(iz) = speed_of_light/100.0 * SS1
          ddd(iz) = rrr(iz)*(1.+zzz(iz))
        END IF
      END DO

      END

c**********************************************************************

      REAL FUNCTION toint8(zz)
c---------------------------------------------------------------------- 
c  The integral of this function from zero to z gives the comoving
c  distance for an object at redshift z.
c----------------------------------------------------------------------
      include 'common.inc'
c----------------------------------------------------------------------
      
      REAL    zz,f1
c---

      f1 = omega_m * (1.0+zz)**3 +  omega_lambda
      f1 = f1+ omega_k * (1.0+zz)**2
      toint8 = 1.0/SQRT(f1)
      
      END

      SUBROUTINE qromo_distance(a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL a,b,ss,EPS,toint8
      EXTERNAL toint8
      PARAMETER (EPS=1.e-6, JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint
cu    use midpnt
      INTEGER j
      REAL dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call midpnt(toint8,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.
11    continue
      print*, 'too many steps in qromo'
      END
