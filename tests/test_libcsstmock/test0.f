        PROGRAM test 
        use, intrinsic :: iso_c_binding ! 和C语言类型互通的模块iso_c_binding
        implicit none
        integer(c_long)  :: Nl1, Nd1, Nc1, N2
        character*256,   allocatable :: carr(:)
        integer(c_long), allocatable :: larr(:,:)
        real(c_double),  allocatable :: darr(:,:)
        Nl1= 3
        Nd1= 5
        Nc1= 3
        N2 = 4
        allocate( larr(nl1, N2) )
        allocate( darr(nd1, N2) )
        allocate( carr(nc1) )
        carr(1) = 'Hello world'
        carr(2) = 'test for dataio'
        carr(3) = 'char array item 3'
        call dummydata(larr, darr, carr, Nl1, Nd1, Nc1, N2)
        write(*,*)(larr(:,1)) 
        write(*,*)(darr(:,1)) 
        write(*,*)('#------------------------ end test0')
 
        END 
