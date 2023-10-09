	    program test3 
        use, intrinsic :: iso_c_binding ! Fortran内部定义的和C语言类型互通的模块iso_c_binding
        implicit none 
        
        integer(8) :: ngal, ngalmax 
        character*256 :: survey  
        real(c_double), allocatable :: x(:), y(:)
        real(c_double), allocatable  :: w(:)
        ngalmax = 6; ngal = 3
        survey  = 'desidr9'
        allocate( x(ngalmax), y(ngalmax), w(ngalmax) )

        x = [0.0, 30.0, 240.0, 120.0, 90.0, 0.0] ! ra的fortran数组
        y = [0.0, 50.0,  44.0,  60.0,  5.0, 0.0] ! dec的fortran数组
        call skycov(x,y,w,survey,ngalmax, ngal) 
        write(*, '(6f12.6)') w

        end program test3 
