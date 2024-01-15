        PROGRAM mainomp

        USE omp_lib

        integer :: rank, numprocs 
        INTEGER(8) :: i,j,N 
        INTEGER(8) :: x, y
        INTEGER :: time_begin, time_end, time_rate

        call omp_set_num_threads(4)

        numprocs = omp_get_num_procs() ! 核心数

        PRINT *, 'number of procs', numprocs 
        !PRINT *, 'n threads', omp_get_num_threads() 

        !$OMP PARALLEL
                PRINT *, "Hello World", OMP_GET_THREAD_NUM() 
        !$OMP END PARALLEL

        write(*,*) 'start test (omp) '

        N = 20000000
        call omp_set_num_threads(1)
        call SYSTEM_CLOCK(time_begin)
        x = 0 
        do i=1,N
           x = x + i
           j = OMP_get_num_threads()
        end do 
        call SYSTEM_CLOCK(time_end)
        WRITE(*,*) 'x = ', x
        WRITE(*,*)'The time wasted on serial computing is: ', 
     &    (time_end - time_begin), 's'

        call omp_set_num_threads(2)
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
        end 
