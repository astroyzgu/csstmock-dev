**********************************************************************
      subroutine subhalo_lightcone2(sub_in,idsub_in,Nd1,Nd2max,Nin,
     & zsmin,zsmax,boxsize,sub_ovlap,idsub_ovlap,Novlap,
     & sub_out,idsub_out,Nout)
c================================================================
c This file input the subhalos in different snapshots to
c generate light-cone subhalo catalogs. 
c Note: subhalos in the same host halo should be put together,
c starting with subhalo id=0 (main subhalo).
c
c sub_in: input subhalo catalog; format: 1-3=x,y,z; 4-6=vx,vy,vz;
c    7-10=halo mass, subhalo mass, max mass, max circular velocity.
c
c idsub_in: input IDs; format: 1: halo id, 2: subhalo id, 3: subhalo trackID
c
c Nd1=10: the number of columns.
c Nd2max: the maximum row of data.
c
c Nin: the total number of subhalos inputted
c
c The redshift range and boxsize shall be provided in generating the light-cone.
c
c zsmin, zsmax: the minmum and maximum redshift.
c boxsize: the boxsize of the simulation in unit of Mpc/h
c
c ================================================================
c We use two snapshots to interpolate the related information
c xxx_ovlap: information in previous snapshot
c
c ================================================================
c In order to generate multi-mock galaxy redshift surveys, 
c we rotate and shift the box at begining.
c
c It is specified by isample.
c
c ================================================================
c The output subhalo catalog format: 
c
c sub_out: 1-3=ra,dec,r, 4-6=z_cos, z_obs,vr, 
c   others are the same as input sub_in.
c
c idsub_out: the same format as idsub_in
c
c Nout: the number of subhalos outputted.
c
c Written by: Xiaohu Yang, 2022, 10, 20 
c
c================================================================      
      include 'common.inc'
c================================================================
      INTEGER Nd1
      INTEGER*8 Nd2max
      INTEGER*8 Nin,Novlap,Nout

      REAL*8 sub_in(Nd1,Nd2max)
      REAL*8 sub_out(Nd1,Nd2max)
      REAL*8 sub_ovlap(Nd1,Nd2max)
      REAL zsmin,zsmax,boxsize
      INTEGER*8 idsub_in(3,Nd2max),idsub_out(3,Nd2max)
      INTEGER*8 idsub_ovlap(3,Nd2max)

      INTEGER*8, allocatable::idtrack(:)

      INTEGER ix,iy,iz
      INTEGER*8 nobj,inst,i_sub
      INTEGER*8 igal,i0,i1
      INTEGER*8 isubmax,isubmin
      REAL    x1,y1,z1,r1,x2,y2,z2,r2
      REAL    x,y,z,vx,vy,vz,r,phi,theta,r_host,r_ovlap
      REAL    ct,cp,vr,rLmin,w_1,w_2
      REAL    ra,dec,z_cos,z_obs,dd
      REAL    xL,absmag,D_L,kecorr,appmag

      INTEGER  lblnk
      REAL     gasdev,ran1,get_z,get_r
      EXTERNAL lblnk,gasdev,ran1,get_z,get_r
      REAL, EXTERNAL:: ekcor_z05

      rLbox=boxsize

      isubmax=0
      isubmin=100
      do i1=1,Novlap
        isubmax=max(isubmax,idsub_ovlap(3,i1))
        isubmin=min(isubmin,idsub_ovlap(3,i1))
      enddo
      write(*,*) 'TrankID range: ', isubmax,isubmin
      allocate(idtrack(isubmin:isubmax))
      idtrack=0
      do i1=1,Nin
        i_sub=idsub_in(3,i1)
        idtrack(i_sub)=i1
      enddo

c--  'Generate random sample (Y/N)?'
      IF(samp_ran.EQ.'Y') THEN
        call randomize8(sub_in,Nd1,Nd2max,Nin)
        call randomize8(sub_ovlap,Nd1,Nd2max,Novlap)
      ENDIF

c-- check the redshift range of this snapshot -----------------

      IF(zsmax.ge.zzmax) print*, 'initial zzmax too small'
      rLmin=get_r(zsmin)
      rLmax=get_r(zsmax)
      rLbox=boxsize
      print*,'Redshift range',zsmin,zsmax,rLmin,rLmax

      Lduplicate=int(rLmax/rLbox+0.51)    !!! duplication of boxes 
      rLcent=rLbox*0.5  !!! move to center
c      print*,Lduplicate,rLcent     

c-- start to select subhalos ------------------

      i0=0
      DO ix = -Lduplicate,Lduplicate
        DO iy = -Lduplicate,Lduplicate
          DO iz = -Lduplicate,Lduplicate

!$OMP PARALLEL DO default(none)
!$OMP& private(igal,x,y,z,vx,vy,vz,r,z_cos,z_obs,vr)
!$OMP& private(x1,y1,z1,r1,x2,y2,z2,r2)
!$OMP& private(ct,cp,theta,phi,dd,appmag,absmag,xL,d_L)            
!$OMP& private(ra,dec,i_sub,kecorr,i1,w_1,w_2)
!$OMP& shared(sub_in,sub_out,ix,iy,iz,rLbox,rLcent,Nsimsel)
!$OMP& shared(zsmin,zsmax,rLmax,rLmin)
!$OMP& shared(i0,idsub_ovlap,idsub_out,sub_ovlap,idtrack)
!$OMP& shared(iseed,sig_vel,amag_cut,Nin,Novlap)
             
            DO igal=1,Novlap

c--- compute comoving distance of the subhalo in previous snapshot 
              x1 = sub_ovlap(1,igal) + FLOAT(ix) *rLbox - rLcent
              y1 = sub_ovlap(2,igal) + FLOAT(iy) *rLbox - rLcent
              z1 = sub_ovlap(3,igal) + FLOAT(iz) *rLbox - rLcent
              r1 =SQRT(x1**2 + y1**2 + z1**2)

c--- compute comoving distance and select subhalos
              i_sub=idsub_ovlap(3,igal)
              i1=idtrack(i_sub)
              IF(i1.eq.0) THEN
                w_1=1.
                w_2=0.
                i1=1
                GOTO 111
              ENDIF
              x2 = sub_in(1,i1) + FLOAT(ix) *rLbox - rLcent
              y2 = sub_in(2,i1) + FLOAT(iy) *rLbox - rLcent
              z2 = sub_in(3,i1) + FLOAT(iz) *rLbox - rLcent
              r2 = SQRT(x2**2 + y2**2 + z2**2)


              IF(abs(r1-r2).ge.20.) print*,'wrong',r1,r2,igal

              IF(min(r1,r2).gt.rLmax) GOTO 544
              IF(max(r1,r2).le.rLmin) GOTO 544

              r=(r1+r2)/2.
              w_1=(rLmax-r)/(rLmax-rLmin)
              w_2=(r-rLmin)/(rLmax-rLmin)
111           r=r1*w_1+r2*w_2 !! more accurate distance

              IF(r.gt.rLmax) GOTO 544
              IF(r.le.rLmin) GOTO 544

              x=x1*w_1+x2*w_2
              y=y1*w_1+y2*w_2
              z=z1*w_1+z2*w_2

c---  compute velocity in radial direction
             
              vx = sub_ovlap(4,igal)*w_1+sub_in(4,i1)*w_2
              vy = sub_ovlap(5,igal)*w_1+sub_in(5,i1)*w_2
              vz = sub_ovlap(6,igal)*w_1+sub_in(6,i1)*w_2
              vr = (vx*x + vy*y + vz*z)/r

c---  add random error to redshift, and compute final los velocity

              vr= vr+sig_vel*gasdev(iseed)

c---  compute cosmological redshift and add contribution from peculiar
c     velocity
              
              z_cos = get_z(r)
              z_obs = z_cos + (vr/speed_of_light) * (1.0+z_cos)

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

c--- calculate the apparent magnitude of the galaxies
c--- select those that fulfill the mag cut.

             xL=sub_ovlap(10,igal)*w_1+sub_in(10,i1)*w_2
             absmag=Msunx-xL*2.5
             D_L = (1.0 + z_obs) * r
             appmag = absmag + 5.0*ALOG10(D_L) + 25.0
             kecorr = ekcor_z05(z_obs)
             appmag = appmag + kecorr !!-1.62*z_obs

             IF (appmag.GT.amag_cut) GOTO 544


!$OMP CRITICAL(UPDATE_SHARED)  
              i0=i0+1
              idsub_out(1:3,i0)=idsub_ovlap(1:3,igal)
              sub_out(1,i0)=ra
              sub_out(2,i0)=dec
              sub_out(3,i0)=r
              sub_out(4,i0)=z_cos
              sub_out(5,i0)=z_obs
              sub_out(7:10,i0)=w_1*sub_ovlap(7:10,igal)
     &       +w_2*sub_in(7:10,i1)
              sub_out(8,i0)=appmag
!$OMP END CRITICAL(UPDATE_SHARED)
              
544          CONTINUE
            END DO

!$OMP END PARALLEL DO

c          print*,'The number of subhalo selected...',ix,iy,iz,i0
c===============================================================

          END DO
        END DO
      END DO

      Nout=i0     

      WRITE(*,*)'Total',Nout,' subhalos in the snapshot redshift range.'

      deallocate(idtrack)

      RETURN
      END






