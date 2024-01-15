c**********************************************************************
      subroutine generate(galpar,indexgal,Nd1,Nd2max,Nmock)
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

      INTEGER ix,iy,iz,N_ass
      INTEGER*8 igal,Nblock,i0,i1,i2,i3
      REAL    x,y,z,vx,vy,vz,r,v,phi,theta
      REAL    ct,st,cp,sp,vr,appmag,compl
      REAL    ra,dec,z_cos,z_obs,kecorr
      REAL      D_L
      REAL,allocatable::  gaxblock(:,:)

      REAL*8,allocatable:: raxx(:),decxx(:),wwxx(:)

      INTEGER  lblnk
      REAL     gasdev,ran1,get_z
      EXTERNAL lblnk,gasdev,ran1,get_z

      REAL, EXTERNAL::ekcor_z05

      Nblock=50000000
      print*, Nblock
      allocate(idgaxblock(4,Nblock))
      allocate(gaxblock(11,Nblock))
      allocate(raxx(Nblock))
      allocate(decxx(Nblock))
      allocate(wwxx(Nblock))


      OPEN(13,file=outdir(1:lblnk(outdir))//outfile(1:lblnk(outfile))
     &,     status='replace')


c---decide whether or not to use the L100 simulation

      Nmock = 0
      i0=0

      DO ix = 1,Lduplicate
        DO iy = 1,Lduplicate
          DO iz = 1,Lduplicate


!$OMP PARALLEL DO default(none)
!$OMP& private(igal,x,y,z,vx,vy,vz,r,z_cos,z_obs)
!$OMP& private(vr,D_L,appmag,kecorr,ct,cp,theta,phi)             
!$OMP& private(ra,dec,compl)
!$OMP& shared(galpar,ix,iy,iz,rLbox,rLcent,Nsimsel)
!$OMP& shared(iseed,z_cut1,z_cut2,amag_cut,Nmock,Nblock)
!$OMP& shared(i0,indexgal,idgaxblock,gaxblock)
!$OMP& shared(raxx,decxx,wwxx)
             
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

c---  compute cosmological redshift and add contribution from peculiar
C     velocity
              
              z_cos = get_z(r)
              z_obs = z_cos + (vr/speed_of_light) * (1.0+z_cos)
              
c---  add random error to redshift, and compute final los velocity

!!!              z_obs = z_obs + (sig_vel*gasdev(iseed))/speed_of_light   

              IF (z_obs.LT.z_cut1.OR.z_obs.GT.z_cut2) GOTO 544

c---compute luminosity distance (Only valid for flat cosmology)

              D_L = (1.0 + z_cos) * r
              
c---compute apparent magnitude, using luminosity distance

              appmag = galpar(9,igal) + 5.0*ALOG10(D_L) + 25.0

c---apply negative k+e-correction 

              kecorr = ekcor_z05(z_obs)
              appmag = appmag + kecorr !!-1.62*z_obs
            
c---apply apparent magnitude limit

              IF (appmag.GT.amag_cut) GOTO 544

c---  compute spherical coordinates theta and phi
              
              ct = z/r
              ct = min(1.0,ct)
              ct = max(-1.0,ct)
             !! st = SQRT(1.0 - ct**2)
              cp = x / SQRT(x**2 + y**2)
             !! sp = y / SQRT(x**2 + y**2)
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
              gaxblock(3,i0)=z_cos
              gaxblock(4,i0)=z_obs
              gaxblock(5,i0)=appmag
              gaxblock(6,i0)=amag_cut
              gaxblock(8:10,i0)=galpar(7:9,igal)
              raxx(i0)=ra
              decxx(i0)=dec
!$OMP END CRITICAL(UPDATE_SHARED)
              
544          CONTINUE
            END DO

!$OMP END PARALLEL DO

 
c===============================================================
      print*,ix,iy,iz,i0
      N_ass=(i0/5000000)+1
      do i=1,N_ass
       i1=i0*(i-1)/N_ass+1
       i2=i0*(i)/N_ass
       i3=i2-i1+1
!!       call python_assign(raxx(i1:i2),decxx(i1:i2),wwxx(i1:i2),i3)
      enddo

      DO i1=1,i0
        gaxblock(7,i1)=wwxx(i1)
c---check whether galaxy falls within survey limits
!!      IF(wwxx(i1).ge.0.5) THEN
         Nmock = Nmock + 1
         idgaxblock(1,i1)=Nmock
         WRITE(13,133)idgaxblock(1:4,i1),gaxblock(1:10,i1)
!!       ENDIF
      END DO
      i0=0
      print*,'The number of galaxies written...',Nmock
c===============================================================

          END DO
        END DO
      END DO

      CLOSE(13)


      WRITE(*,*)' '
      WRITE(*,*)' Total of ',Nmock,' galaxies in mock survey.'

      WRITE(*,*)' '

c---FORMATS

 133  FORMAT(I9,1X,2(I15,1X),I8,1x,2(F12.7,1X),8(F9.5,1X))


      deallocate(idgaxblock)
      deallocate(gaxblock)
      deallocate(raxx)
      deallocate(decxx)
      deallocate(wwxx)

      RETURN
      END


c**********************************************************************

      SUBROUTINE set_Survey_mask
c----------------------------------------------------------------------
c This subroutine sets the ra/dec boundaries of the survey 
c================================================================      
      include 'common.inc'
c================================================================
      
      ra_min(1)= 300.0 * degs_rads
      ra_max(1)= 80.0 * degs_rads     
      dec_min(1)= -12.0 * degs_rads    
      dec_max(1)= 17.0 * degs_rads    

      ra_min(2)= 108.0 * degs_rads
      ra_max(2)= 268.0 * degs_rads
      dec_min(2)= -6.5 * degs_rads
      dec_max(2)= 71.5 * degs_rads

      RETURN 
      END


