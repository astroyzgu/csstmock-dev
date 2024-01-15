**********************************************************************
      subroutine subhalo_lightcone(sub_in,idsub_in,Nd1,Nd2max,Nin,
     & zsmin,zsmax,boxsize,sub_ovlap,id_ovlap,Novlap,
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
c In order to add possible missing subhalos at previous snapshot
c due to the structure evolution, we also need the host halo ID 
c and ra, dec at the outskirt of previous snapshot with \Delta
c r = 3Mpc/h.
c
c sub_ovlap: 1=ra, 2=dec
c
c id_ovlap: host halo ID in previous snapshot
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
      REAL*8 sub_ovlap(2,Nd2max)
      REAL zsmin,zsmax,boxsize
      INTEGER*8 idsub_in(3,Nd2max),idsub_out(3,Nd2max)
      INTEGER*8 id_ovlap(Nd2max)

      INTEGER*8, allocatable::ihoc(:), ll(:),ihalomain(:),imainsub(:)

      INTEGER ix,iy,iz
      INTEGER*8 nobj,inst,i_sub
      INTEGER*8 igal,i0,i1,ihalo,ihalomax,ihalomin
      INTEGER*8 isubmax,isubmin
      REAL    x,y,z,vx,vy,vz,r,phi,theta,r_host
      REAL    ct,cp,vr,rLovlap,rLmin
      REAL    ra,dec,z_cos,z_obs,dd
      REAL    xL,absmag,D_L,kecorr,appmag

      INTEGER  lblnk
      REAL     gasdev,ran1,get_z,get_r
      EXTERNAL lblnk,gasdev,ran1,get_z,get_r
      REAL, EXTERNAL:: ekcor_z05

      rLbox=boxsize

cccc main subhalo ID 

      ihalomax=0
      ihalomin=100
      isubmax=0
      isubmin=100
      do i1=1,Nin
        ihalomax=max(ihalomax,idsub_in(1,i1))
        ihalomin=min(ihalomin,idsub_in(1,i1))
        isubmax=max(isubmax,idsub_in(3,i1))
        isubmin=min(isubmin,idsub_in(3,i1))
      enddo
      write(*,*) 'Halo ID range: ', ihalomax,ihalomin
      write(*,*) 'TrankID range: ', isubmax,isubmin

      allocate(ihalomain(Nin))
      allocate(imainsub(ihalomin:ihalomax))

      do i1=1,Nin
        ihalo=idsub_in(1,i1)
        IF(idsub_in(2,i1).eq.0) imainsub(ihalo)=i1 !!! main subhalo
      enddo

      do i1=1,Nin
        ihalo=idsub_in(1,i1)
        IF(ihalo.gt.0) THEN
          ihalomain(i1)=imainsub(ihalo) 
        ELSE
          ihalomain(i1)=i1
        ENDIF
      enddo

      deallocate(imainsub)


c-- check if more than 1 sample is generated: to shift and rovate --

      do i=1,isample
        call rotate8(sub_in,Nd1,Nd2max,Nin)
      enddo

c--  'Generate random sample (Y/N)?'

      IF(samp_ran.EQ.'Y') THEN
        call randomize8(sub_in,Nd1,Nd2max,Nin)
        do i1=1,Nin
         ihalomain(i1)=i1   !!! main
        enddo
      ENDIF

c-- prepare for the overlap region ---------------------------
c-- here we use TrackID to check and remove the duplications 

      rLovlap=0.0

!      nobj=Novlap
!      ihalomax=0
!      ihalomin=100
!      do igal=1,nobj
!        ihalomax=max(ihalomax,id_ovlap(igal))
!        ihalomin=min(ihalomin,id_ovlap(igal))
!      enddo
!      write(*,*) 'TrackID in overlap region:'
!     & ,ihalomax,ihalomin,nobj
!
!
!      ihalomax=max(ihalomax,isubmax)
!      ihalomin=min(ihalomin,isubmin)
!      allocate(ihoc(ihalomin:ihalomax))
!      allocate(ll(nobj))
!
!      ll=0
!      ihoc=0
!
!      do igal=1,nobj
!         ihalo=id_ovlap(igal)
!         IF(ihalo.lt.0) GOTO 88
!         inst=ihoc(ihalo)
!         ll(igal)=inst
!         ihoc(ihalo)=igal
!88    enddo

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
!$OMP& private(ct,cp,theta,phi,dd,appmag,absmag,xL,d_L)            
!$OMP& private(ra,dec,r_host,ihalo,i_sub,kecorr)
!$OMP& shared(sub_in,sub_out,ix,iy,iz,rLbox,rLcent,Nsimsel)
!$OMP& shared(zsmin,zsmax,rLmax,rLmin,rLovlap,ll,ihoc)
!$OMP& shared(i0,idsub_in,idsub_out,sub_ovlap,ihalomain)
!$OMP& shared(iseed,sig_vel,amag_cut,Nin)
             
            DO igal=1,Nin

c--- compute comoving distance of the host halo 
c--- here we keep subhalos according to the distance of the host halo!!!
              ihalo=ihalomain(igal)
              x = sub_in(1,ihalo) + FLOAT(ix) *rLbox - rLcent
              y = sub_in(2,ihalo) + FLOAT(iy) *rLbox - rLcent
              z = sub_in(3,ihalo) + FLOAT(iz) *rLbox - rLcent
              r_host=SQRT(x**2 + y**2 + z**2)

              IF(r_host.gt.rLmax) GOTO 544
              IF(r_host.le.rLmin-rLovlap*0.5) GOTO 544

c--- compute comoving distance and select subhalos

              x = sub_in(1,igal) + FLOAT(ix) *rLbox - rLcent
              y = sub_in(2,igal) + FLOAT(iy) *rLbox - rLcent
              z = sub_in(3,igal) + FLOAT(iz) *rLbox - rLcent
              r = SQRT(x**2 + y**2 + z**2)

c---  compute velocity in radial direction
             
              vx = sub_in(4,igal)
              vy = sub_in(5,igal)
              vz = sub_in(6,igal)
              vr = (vx*x + vy*y + vz*z)/r

c---  add random error to redshift, and compute final los velocity

c              vr= vr+sig_vel*gasdev(iseed)

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

c---check if the host halo is already there from the previous snapshot

!              IF(r_host.le.rLmin) THEN
!                ihalo=idsub_in(3,igal)
!                i_sub=ihoc(ihalo)
!                do while (i_sub.ne.0) 
!                  dd=(ra-sub_ovlap(1,i_sub))**2
!                  dd=dd+(dec-sub_ovlap(2,i_sub))**2 
!                  IF(dd.le.10.) GOTO 544
!                  i_sub=ll(i_sub)
!                enddo
!              ENDIF

c--- calculate the apparent magnitude of the galaxies
c--- select those that fulfill the mag cut.

             xL=sub_in(10,igal)
             absmag=Msunx-xL*2.5
             D_L = (1.0 + z_obs) * r
             appmag = absmag + 5.0*ALOG10(D_L) + 25.0
             kecorr = ekcor_z05(z_obs)
             appmag = appmag + kecorr !!-1.62*z_obs

             IF (appmag.GT.amag_cut) GOTO 544


!$OMP CRITICAL(UPDATE_SHARED)  
              i0=i0+1
              idsub_out(1:3,i0)=idsub_in(1:3,igal)
              sub_out(1,i0)=ra
              sub_out(2,i0)=dec
              sub_out(3,i0)=r
              sub_out(4,i0)=z_cos
              sub_out(5,i0)=z_obs
              sub_out(6,i0)=r_host
              sub_out(7:10,i0)=sub_in(7:10,igal)
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

!      deallocate(ll)
!      deallocate(ihoc)
      deallocate(ihalomain)

ccc here we update the subhalo in the next snapshot overlap region...

!      i0=0
!      i1=0
!      do igal=1,Nout
!        IF(sub_out(6,igal).le.rLmin) i1=i1+1
!
!        IF(sub_out(6,igal).ge.(rLmax-rLovlap))THEN
!          i0=i0+1
!          sub_ovlap(1:2,i0)=sub_out(1:2,igal)
!          id_ovlap(i0)=idsub_out(3,igal)  !! keep trackID
!        ENDIF
!      enddo
!      Novlap=i0
!      print*,i1,'halos in the overlap region added.'
!      print*,i0,'halos in the new overlap region updated.'


      RETURN
      END






