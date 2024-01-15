c=======================================================================
c
c                      INCLUDE FILE common.inc
c
c=======================================================================
c  Written by: Xiaohu Yang  
c=======================================================================

        use, intrinsic :: iso_c_binding ! Fortran内部定义的和C语言类型互通的模块iso_c_binding

        implicit none

        REAL  pi,pival,degs_rads,tpi,rads_degs
        PARAMETER (pi=3.141592654)
        PARAMETER (pival=pi)
        PARAMETER (degs_rads = pival/180.0)
        PARAMETER (rads_degs = 180./pival)
        PARAMETER (tpi = 2.*pival)

c---  a small number

        REAL eps
        PARAMETER (eps=1.e-5)
 
c---  The absolute magnitude of the sun in x band
 
        REAL   Msunx
        PARAMETER (Msunx = 4.80)

c---  The speed of light in km/s
  
        REAL   speed_of_light
        PARAMETER (speed_of_light = 3.0E+5)

c---  number of grid points for comoving_radius(redshift)
c---  maxmum redshift in this table: zzmax
        
        INTEGER Nz
        PARAMETER (Nz = 3000)
        REAL zzmax
        PARAMETER (zzmax=3.0)

c---  cosmology

        REAL  omega_m,omega_lambda,omega_k,zzz(Nz),rrr(Nz),ddd(Nz)
        COMMON/cosmo/omega_m,omega_lambda,omega_k,zzz,rrr,ddd

c---  boundaries (ra,dec) of SDSS regions 1, 2 and 3 ...

        REAL  ra_min(4),ra_max(4),dec_min(4),dec_max(4)
        COMMON/boundaries/ra_min,ra_max,dec_min,dec_max

c---  I/O related parameters ...
        
        INTEGER*8 Ngaldim,Ngalmax,Ngalsel,Nsimdim,Nsimmax,Nsimsel
     &   ,Nmockdim,Nmockmax,Nmocksel,Ngaldim1   
        COMMON/inputnumbers/Ngaldim,Ngalmax,Ngalsel,Ngaldim1
     &  ,Nsimdim,Nsimmax,Nsimsel,Nmockdim,Nmockmax,Nmocksel


c---  IO directory and file names

        CHARACTER*60 indir,infile,infile2,outLF
        CHARACTER*40 outdir,CLFfile,outfile,file_sp
        CHARACTER*20 obs_name,sim_name,mask_name
        CHARACTER samp_ran,cat_type,match_only,test_only
        COMMON/dirinfo/indir,outdir,infile,outfile,outLF,infile2
     &   ,obs_name,sim_name,mask_name,CLFfile,samp_ran,cat_type
     &   ,file_sp,match_only,test_only

c---  M/L ratio
      
        INTEGER N_GLF,Nzbin
        PARAMETER (N_GLF=100,Nzbin=5)
      
        REAL GMFa(Nzbin+1,N_GLF),GLFa(Nzbin+1,N_GLF)
        REAL GMFc(Nzbin+1,N_GLF),GLFc(Nzbin+1,N_GLF)
        REAL GMFs(Nzbin+1,N_GLF),GLFs(Nzbin+1,N_GLF)
        REAL ML_gax(3,Nzbin+1,N_GLF)
        COMMON/CLF/GMFa,GLFa,GMFc,GLFc,GMFs,GLFs,ML_gax

c---  luminosity function parameters

        REAL amag_cut,xlum_cut,xlum_1,xlum_2,dxlum,V_fact
        COMMON/LFDESI/amag_cut,xlum_cut,xlum_1,xlum_2,dxlum,V_fact

c--- redshift bins

        real z_cut1,z_cut2,z_cut0,dzbin,skycover,fact
        real z_min(Nzbin),z_max(Nzbin),z_eff(Nzbin)
        real v_min(Nzbin),v_max(Nzbin),v_eff(Nzbin)
        COMMON /dzinfo1/z_cut1,z_cut2,z_cut0,dzbin,skycover,fact
        COMMON /dzinfo2/z_min,z_max,z_eff,v_min,v_max,v_eff

c---  boxsize, duplication, redshift range, etc.

        REAL rLbox, rLcent,rLmax,sig_vel
        INTEGER Lduplicate
        COMMON/boxsize/Lduplicate,rLbox,rLcent,rLmax,sig_vel
        
        INTEGER i,j,k,iseed,isample
        INTEGER*8 idble,jdble
        COMMON/random/iseed,isample

c=======================================================================
c                             END
c=======================================================================


