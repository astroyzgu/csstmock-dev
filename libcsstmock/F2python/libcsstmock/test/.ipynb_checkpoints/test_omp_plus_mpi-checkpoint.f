        PROGRAM mainomp
        USE omp_lib
        use, intrinsic :: iso_c_binding ! 和C语言类型互通的模块iso_c_binding
        implicit none

        INTEGER(c_long) :: i,j 
        real :: g(1000, 1000)
        real :: time_begin, time_end, time_rate
        integer(c_long) :: Nin, Noutmax
        real(c_double), allocatable :: x_in(:), x_out(:)

        INTEGER(8) :: snap 
        integer(8) :: nd1, nl1, ngal 
        character*256 :: basedir 
        integer(8), allocatable :: larr(:,:)
        real(8), allocatable :: darr(:,:)

        basedir = '/home/cossim/Jiutian/M1000/'
        snap    = 127 
        nd1     =  3
        nl1     =  2
        ngal    = 200
        allocate( larr(ngal, nl1), darr(ngal,nd1) )
        call dummpydata_mpi(snap, basedir, 
     &                   darr, larr, nd1, nl1, ngal)

c        write(*,"(A10, 3F10.6)")'double:', darr(:, 1)
c        write(*,"(A10, 3F10.6)")'double:', darr(:, 1)
c        write(*,"(A10,  3I10)") 'long:', larr(:, 1)

c        CALL omp_set_num_threads(2)
c        !$OMP PARALLEL
c        WRITE(*,*)"Hello World"
c        !$OMP END PARALLEL

c        time_begin = omp_get_wtime()
c        !$omp parrallel private(i,j)  
c        !$omp do 
c        DO i = 1, 1000
c            DO j = 1, 1000
c               g(i,j) = exp(1.0*i*j) 
c            END DO
c        END DO
c        !$omp end do 
c        !$omp end parrallel  
c        time_end = omp_get_wtime()
c        WRITE(*,*)'The time wasted on serial computing is: ', 
c     &    1.0*(time_end - time_begin), 's'
c 
        write(*,*) 'end test (omp + mpi) '
        stop
        END

