# 测试:
结论： 目前，开启mpi后，python程序能成功运行多个子进程，并给出结果。 
  但是！！！，若在python段关闭mpi， 则无法再次开启其他含有mpi4py的程序。 
            若在python段不关闭mpi，则产生多个平行的fortran程序。 
            当然，可以通过返回rank在fortran端关闭mpi, 但是不利于后续程序的编写。 
            
cd tset  
gfortran -o testomp testomp.f  -lcsstmock  -fopenmp 
gfortran -o test    test.f     -lcsstmock  

mpirun -np 2 ./test    # test the data 
mpirun -np 2 ./testomp # test the speed using both mpi in python and openmp in fortran    

c=============================================================================
c
c--- the details of test directory  
c    1. pureomp.f    # 测试openmp和fortran是否运行正常 
c       ifort -o pureomp pureomp.f -fopenmp; ./pureomp # 编译和执行 
c
c    2. sharematrix.f # 测试fortran和python的各类数组(long, double, char)是否正确共享内存
c       ifort -o sharematrix sharematrix.f -lcsstmock; sharematrix
c 
c       call dummydata(larr, darr, carr, Nl1, Nd1, Nc1, N2)
c       subroutine dummydata(larr, darr, carr,Nl1, Nd1, Nc1, N2) bind (c)
c       extern  void dummydata_c(long *, double *, char *, long *, long *, long *, long *);
c       def dummydata_c(larr, darr, carr, Nl1, Nd1, Nc1, N2):
c 
c    3. sharematrix_mpi.f # 测试在python部分使用mpi加速后下是否能够在fortran部分仍正常正确共享内存
c       ifort -o sharematrix_mpi sharematrix_mpi.f -lcsstmock; sharematrix_mpi
c
c    4. omp_plus_mpi.f    # 测试混合mpi和openmp后程序是否仍然可以高效运行(使用mpi是否会限制后续omp资源的请求)
c       ifort -o test_omp_plus_mpi test_omp_plus_mpi.f -lcsstmock -fopenmp 
c       mpirun -np 2 omp_plus_mpi
c
c--- the subroutine name in each part of this test 
c 
c main program    io_dummy.f/io_jiutian.f  mymodule.py/plugin.h  myfunc.py/haloinfo.py
c exec in test    <-   subroutine      <-     cffi wrapper      <- origin  python 
c  test_share    call  dummy_data          dummy_data_c           dummy_data_py 
c  test_mpi      call  dummy_data_mpi      dummy_data_mpi_c       dummy_data_mpi_py 
c  test_mpiomp         ... 
c  test          call  read_jiutian        read_jiutian_mpi_c     read_jiutian_mpi   
c 
c=============================================================================
c 
c run_test.sh 
