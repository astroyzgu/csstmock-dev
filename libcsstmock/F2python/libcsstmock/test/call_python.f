        subroutine loaddata(Nin, Nout)
        !use, intrinsic :: iso_c_binding ! 和C语言类型互通的模块iso_c_binding
        implicit none 
        integer(4) :: Nin, Nout  
        Nin  = 5
        Nout = 5 
        !real(8) :: x_in(Nin), x_out(Nout)
        !x_out(2) = 2
        !x_out(1) = 10 
        return 
	end subroutine loaddata
