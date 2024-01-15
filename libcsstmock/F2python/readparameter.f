c**********************************************************************
c  Subroutine used to read the global parameters
c**********************************************************************

      SUBROUTINE read_global_param
c---------------------------------------------------------------------------
c  Subroutine to read in global parameters
c---------------------------------------------------------------------------
      include 'common.inc'
      
      CHARACTER name_1,z_range1*3,z_range2*3
      INTEGER, EXTERNAL:: lblnk

      OPEN(9,FILE='parameters', STATUS='OLD')

c--- read observational data to be reproduced ----------

      READ(9,'(A)')obs_name
      WRITE(*,*)' Give observation name : ',obs_name

      READ(9,*)amag_cut
      WRITE(*,*)' Apparent magnitude limit : ',amag_cut

      READ(9,*)z_cut1
      WRITE(*,*)' Give the minimum redshift : ', z_cut1

      READ(9,*)z_cut2
      WRITE(*,*)' Give the maximum redshift : ', z_cut2

      READ(9,'(A)')CLFfile
      WRITE(*,*)' Give cumulative LF name (N) : ', CLFfile

c--- read halo/subhalo and output catalog information ------------

      READ(9,'(A)')sim_name
      WRITE(*,*)' Give simulation catalog name : ',sim_name

      READ(9,*)isample
      WRITE(*,*)' Rotate simulation box times : ',isample

      READ(9,*)Nsimmax
      WRITE(*,*)' Maximum number of simulation data : ',Nsimmax

      READ(9,*)Nsimdim
      WRITE(*,*)' Number of columns of the simulation data : ',Nsimdim

      READ(9,'(A)')outdir
      WRITE(*,*)' Give output catalog directory : ',outdir

      READ(9,*)mask_name
      WRITE(*,*)' Give sky-coverage (mask file) type : ',mask_name

      READ(9,*)samp_ran
      WRITE(*,*)' Generate random sample (Y/N) : ',samp_ran

      READ(9,*)iseed
      WRITE(*,*)' Give random seed: ',iseed

      READ(9,*)match_only
      WRITE(*,*)' Only need make the match: ',match_only

      READ(9,*)test_only
      WRITE(*,*)' Only for test snapshot: ',test_only

      CLOSE(9)
c--- end read parameters --- 

c--- simulation box and cosmology information ---      

      IF(sim_name(1:lblnk(sim_name)).eq.'9tian') THEN
        rLbox=1000.        
        omega_m=0.3111
        omega_lambda=0.6889
      ENDIF

      IF(sim_name(1:lblnk(sim_name)).eq.'Uchuu') THEN
        rLbox=2000.
        omega_m=0.3089
        omega_lambda=1.0-omega_m
      ENDIF


      IF(sim_name(1:lblnk(sim_name)).eq.'Mill') THEN
        rLbox=500.
        omega_m=0.25
        omega_lambda=1.0-omega_m
      ENDIF

      omega_k = 1.0 - omega_m - omega_lambda

      WRITE(*,*)' Simulation box size : ',rLbox

      WRITE(*,*)' Omega_m : ',omega_m

      WRITE(*,*)' Omega_Lambda : ',omega_lambda


c--- set output file name ----

      write(name_1,'(i1)') isample
      write(z_range1,'(f3.1)') z_cut1
      write(z_range2,'(f3.1)') z_cut2

      outfile='mock'//obs_name(1:lblnk(obs_name))//'_'//  
     &  sim_name(1:lblnk(sim_name))//name_1//
     &  '_z'//z_range1//'_'//z_range2
      IF(samp_ran.EQ.'Y') THEN
      outfile='rand'//obs_name(1:lblnk(obs_name))//'_'//
     &  sim_name(1:lblnk(sim_name))//name_1//
     &  '_z'//z_range1//'_'//z_range2
      ENDIF
      
      WRITE(*,*)' Output file name : ',outfile


      RETURN
      END




