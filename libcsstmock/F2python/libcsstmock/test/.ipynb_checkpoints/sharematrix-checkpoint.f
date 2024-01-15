        PROGRAM main
        USE omp_lib
        use, intrinsic :: iso_c_binding ! 和C语言类型互通的模块iso_c_binding
        implicit none

        integer(c_long) :: Nin, Noutmax, snap 
        real(c_double), allocatable :: x_in(:), x_out(:)
        integer(c_long)  :: Nl1, Nd1, Nc1, N2
        integer(c_long) :: Ngal 
        integer(c_long),  allocatable :: larr(:,:)
        real(c_double),   allocatable :: darr(:,:)
        character*256,    allocatable :: carr(:)
        Ngal  = 4 
        Nd1   = 3 
        Nl1   = 3 
        Nc1   = 4 

        allocate( larr(Nl1, Ngal), darr(Nd1, Ngal) ) 
        allocate( carr(Nc1) ) 

!c---   evaluation array of long, double and char types

        larr(1,:) = (/1, 2, 3, 4/) 
        larr(2,:) = (/1, 2, 3, 4/) 
        larr(3,:) = (/1, 2, 3, 4/) 
        darr(1,:) = (/1.1, 2.2, 3.3, 4.4/) 
        darr(2,:) = (/1.1, 2.2, 3.3, 4.4/) 
        darr(3,:) = (/1.1, 2.2, 3.3, 4.4/) 
        carr(1)  = 'this is a test1'
        carr(2)  = 'this is a test2'
        write(*,"(A10, 3F10.6)")'double:', darr(:, 1)
        write(*,"(A10,  3I10)") 'long:', larr(:, 1)
        write(*,"(A10  A20)")"char(1):",  carr(1) 
        write(*,"(A10  A20)")"char(2):",  carr(2) 

        call dummydata(larr, darr, carr, Nl1, Nd1, Nc1, Ngal ) 

        !write(*,*) 'readin ', N2, ' sub halos'
        !write(*,*) Nl1, Nd1 
        !write(*,*) larr(:,1), larr(:,2)
        !write(*,*) darr(:,1), darr(:,2)
	stop   
        END

