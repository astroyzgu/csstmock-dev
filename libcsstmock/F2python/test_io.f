      subroutine test_io
c************************************************************************
c  Construct a mock galaxy redshift survey.
c  written by Xiaohu Yang: xyang@sjtu.edu.cn
c************************************************************************
      use omp_lib      
      include 'common.inc'
c================================================================      
      REAL, allocatable:: gaxsim(:,:),gaxobs(:,:)
      REAL, allocatable:: gaxmock(:,:)
      INTEGER*8, allocatable:: indexgax(:,:)
      INTEGER*8 Nsel,Nmock
c==================================================================
c--- This part is used to test the I/O and debug the code...
      
      outdir='../../data/'
      outfile='mock1'
      Nmock=13475020
      Nmockdim=8
      allocate(gaxmock(Nmockdim,Nmock))
      call inputmock(gaxmock,Nmockdim,Nmock)

      Ngaldim=6
      Ngalmax=70000000
      allocate(gaxobs(Ngaldim,Ngalmax))
      call inputDESIgroup(gaxobs)
      call match_obs(gaxmock,Nmock,gaxobs)

      deallocate(gaxmock)
      deallocate(gaxobs)

      Nsimdim=9
      Nsimmax=812693223
      allocate(gaxsim(Nsimdim,Nsimmax))
      allocate(indexgax(3,Nsimmax))
      call inputSIM_GuoHong(gaxsim,indexgax)
      deallocate(gaxsim)
      deallocate(indexgax)
      
      RETURN
      END



