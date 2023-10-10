	program test3 
        use, intrinsic :: iso_c_binding 
        implicit none 
        
        integer(8) :: ngal 
        character*256 :: survey  
        real(4), allocatable :: x(:), y(:)
        real(4), allocatable :: w(:)

        ngal    = 6; 
        survey  = 'lsdr9'
        allocate( x(ngal), y(ngal), w(ngal) )

        x = [0.0, 30.0, 240.0, 120.0, 90.0, 0.0] ! ra的fortran数组
        y = [0.0, 50.0,  44.0,  60.0,  5.0, 0.0] ! dec的fortran数组
        call isin_survey(x,y,ngal,survey,w) 
        write(*, '(6f12.6)') w

        end program test3 
