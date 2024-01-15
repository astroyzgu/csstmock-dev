        PROGRAM mainomp
        USE omp_lib
        use, intrinsic :: iso_c_binding ! 和C语言类型互通的模块iso_c_binding
        implicit none

        INTEGER(8) :: i,j,N 
        INTEGER(8) :: x, y
        INTEGER :: time_begin, time_end, time_rate

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
        allocate( larr(ngal, nl1), darr(ngal, nd1) )

        write(*,*) 'start test (mpi) '
        write(*,"(A20, 3i10)")'Input  dim:', nd1, nl1, ngal
        call dummydata_mpi(snap, basedir, 
     &                   darr, larr, nd1, nl1, ngal)
        
        write(*,"(A20, 3F10.6)")'1st line in fortran:', darr(1, :)
        write(*,"(A20, 3F10.6)")'2st line in fortran:', darr(2, :)
        write(*,"(A20, 3i10)")'Output dim:', nd1, nl1, ngal

        write(*,*) 'start test (omp) '
        call omp_set_num_threads(3) 
        !$OMP PARALLEL
                PRINT *, "Hello World", OMP_GET_THREAD_NUM()
        !$OMP END PARALLEL 

        N = 20000000
        call SYSTEM_CLOCK(time_begin)
        x = 0 
        do i=1,N
           x = x + i
           j = OMP_get_num_threads()
        end do 
        call SYSTEM_CLOCK(time_end)
        print *, 'x = ', x
        WRITE(*,*)'The time wasted on serial computing is: ', 
     &    (time_end - time_begin), 's'

c        time_begin = omp_get_wtime()
        y = 0 
        i = 0 
        call SYSTEM_CLOCK(time_begin)
        !$OMP PARALLEL DO REDUCTION(+:y) 
             do i=1,N
                y = y + i
                j = OMP_get_num_threads()
             end do 
        !$OMP end parallel do 

c        time_end   = omp_get_wtime()
        call SYSTEM_CLOCK(time_end)
        WRITE(*,*) 'y = ', y
        WRITE(*,*)'The time wasted on serial computing is: ', 
     &    (time_end - time_begin), 's'

        write(*,*) 'end test (omp + mpi) '
        stop
        END

