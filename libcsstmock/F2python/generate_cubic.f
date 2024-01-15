c**********************************************************************
      subroutine generate_cubic(galpar,indexgal,Nd1,Nd2max)
c================================================================
c This file use the mock galaxy data in simulation box to
c generate mock galaxy redshift surveys.
c with rg format  i: 1-3=x,y,z; 4-6=vx,vy,vz; 7=log Mst, 
c 8=log Mh, indexrg (1: galaxy id, 2: halo id, 3, cen/sat type)
c================================================================      
      include 'common.inc'
c================================================================
      INTEGER Nd1
      INTEGER*8 Nd2max,Nmock

      REAL galpar(Nd1,Nd2max)
      INTEGER*8 indexgal(3,Nd2max)
      INTEGER*8,allocatable:: idgaxblock(:,:)
      REAL,allocatable::  gaxblock(:,:)

      INTEGER ix,iy,iz,N_ass,iztype
      INTEGER*8 igal,Nblock,i0,i1,i2,i3
      REAL    x,y,z,vx,vy,vz,r,v,phi,theta
      REAL    ct,st,cp,sp,vr
      REAL    ra,dec,z_cos,z_obs

      REAL     gasdev,ran1,get_z
      EXTERNAL gasdev,ran1,get_z

      skycover=(180./pi)**2*4.*pi  !! all skycov
      call initzbins

      allocate(idgaxblock(3,Nd2max))
      allocate(gaxblock(Nd1,Nd2max))

c---  set up of the mock redshift suvey...

      sig_vel=35.  !!! velocity error of spectro
      Lduplicate=int(rLmax*2./rLbox)+1    !!! duplication of boxes along 1-direction
      rLcent=rLbox*float(Lduplicate)*0.5  !!! move to center
      print*,Lduplicate,rLcent,rLbox

c---  'Generate random sample (Y/N)?'

      if(samp_ran.EQ.'Y')  call randomize(galpar,Nd1,Nd2max,Nsimsel)

c---  rotate boxes if you want to generate more---

      do i=1,isample
         call rotate(galpar,Nd1,Nd2max,Nsimsel)
      enddo

      Nmock = 0
      i0=0

      DO ix = 1,Lduplicate
        DO iy = 1,Lduplicate
          DO iz = 1,Lduplicate


!$OMP PARALLEL DO default(none)
!$OMP& private(igal,x,y,z,vx,vy,vz,r,z_cos,z_obs)
!$OMP& private(vr,ct,cp,theta,phi)             
!$OMP& private(ra,dec,iztype)
!$OMP& shared(galpar,ix,iy,iz,rLbox,rLcent,Nsimsel,Nd1)
!$OMP& shared(iseed,z_cut1,z_cut2,sig_vel,Nmock,v_eff,dzbin)
!$OMP& shared(i0,indexgal,idgaxblock,gaxblock)
             
            DO igal=1,Nsimsel

              x = galpar(1,igal) + FLOAT(ix-1) *rLbox - rLcent
              y = galpar(2,igal) + FLOAT(iy-1) *rLbox - rLcent
              z = galpar(3,igal) + FLOAT(iz-1) *rLbox - rLcent
 
              vx = galpar(4,igal)
              vy = galpar(5,igal)
              vz = galpar(6,igal)

c---  compute comoving distance

              r = SQRT(x**2 + y**2 + z**2)
              
c---  compute velocity in radial direction
              
              vr = (vx*x + vy*y + vz*z)/r

c---  add random error to redshift, and compute final los velocity

              vr= vr+sig_vel*gasdev(iseed)

c---  compute cosmological redshift and add contribution from peculiar
C     velocity and random error
              
              z_cos = get_z(r)
              z_obs = z_cos + (vr/speed_of_light) * (1.0+z_cos)
              
              IF (z_obs.LT.z_cut1.OR.z_obs.GT.z_cut2) GOTO 544

c---  compute spherical coordinates theta and phi
              
              ct = z/r
              ct = min(1.0,ct)
              ct = max(-1.0,ct)
              cp = x / SQRT(x**2 + y**2)
              cp = min(1.0,cp)
              cp = max(-1.0,cp)

              theta = ACOS(ct)
              IF (y.GT.0.0) THEN
                phi = ACOS(cp)
              ELSE
                phi = 2.0*pi - ACOS(cp)
              END IF
         
c---convert theta and phi to ra and dec (both in degrees)

              ra = phi * rads_degs
              dec = ((pi/2.0) - theta) * rads_degs  

c---compute coordinates for plotting and write to file

!$OMP CRITICAL(UPDATE_SHARED)  
              i0=i0+1
              idgaxblock(2:4,i0)=indexgal(1:3,igal)
              gaxblock(1,i0)=ra
              gaxblock(2,i0)=dec
              gaxblock(3,i0)=r
              gaxblock(4,i0)=z_cos
              gaxblock(5,i0)=z_obs
              gaxblock(6,i0)=vr
              gaxblock(7:8,i0)=galpar(7:8,igal)
              iztype=int((z_obs-z_cut1)/dzbin)+1
              gaxblock(9,i0)= float(iztype)+0.1
              gaxblock(10,i0)=v_eff(iztype)
!$OMP END CRITICAL(UPDATE_SHARED)
              
544          CONTINUE
            END DO

!$OMP END PARALLEL DO

          END DO
        END DO
      END DO

      Nmock=i0
      Nsimsel=Nmock

      do i1=1,Nmock
        galpar(1:Nd1,i1)=gaxblock(1:Nd1,i1)
        indexgal(1:3,i1)=idgaxblock(1:3,i1)
      enddo

      WRITE(*,*)' '
      WRITE(*,*)' Total of ',Nmock,' galaxies in redshift range.'
      WRITE(*,*)' '

      deallocate(idgaxblock)
      deallocate(gaxblock)

      RETURN
      END




