        PROGRAM test 
        USE omp_lib
        use, intrinsic :: iso_c_binding ! 和C语言类型互通的模块iso_c_binding
        implicit none

        INTEGER(8) :: i,j,N 
        INTEGER(8) :: x, y
        INTEGER :: time_begin, time_end, time_rate

        INTEGER(4) :: snap, nl1, nd1
        integer(8) :: ngal, ngalmax 
        character*256 :: basedir 
        integer(c_long), allocatable :: larr(:,:)
        real(c_double), allocatable  :: darr(:,:)

        basedir = '/home/cossim/Jiutian/M1000/'
        snap    = 127 
        nd1     =  10
        nl1     =   3
        ngalmax = 500000000
        allocate( larr(nl1, ngalmax), darr(nd1, ngalmax) )

        write(*,*) 'start test (jiutian) '
        write(*,"(A20, 3i10)")'Input  dim:', nd1, nl1, ngalmax
        call readsnapshot_jiutian( snap, darr
     &                   , larr, nd1, ngalmax, ngal)
        write(*,"(A20, 3i10)")'Output  dim:', nd1, nl1, ngal
        write(*,"(A20, 3i10)")'1st line in fortran:', larr(:, 1)
        write(*,"(A20, 3i10)")'2st line in fortran:', larr(:, 2)
        write(*,"(A20, 3i10)")'3st line in fortran:', larr(:, 3)
        write(*,"(A20, 3i10)")'last line in fortran:', larr(:, ngal)

        write(*,*) darr(:, 1)
        write(*,*) darr(:, 2)
        write(*,*) darr(:, 3)
        write(*,*) darr(:, ngal)

        write(*,"(A20, 3i10)")'Output dim:', nd1, nl1, ngal
        
        END 
