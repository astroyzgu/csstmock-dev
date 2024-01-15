        subroutine loaddata(rg, Nd1, Nd2max) 
        !use, intrinsic :: iso_c_binding ! 和C语言类型互通的模块iso_c_binding
        include 'common.inc'
        INTEGER Nd1
        INTEGER Nd2max
        REAL rg(Nd1,Nd2max)
        rg = 0.0
        rg(1,1)=1 !! ra
        rg(2,1)=2 !! dec
        rg(3,1)=3 !! z 
        return 
	end subroutine loaddata

        PROGRAM main         
        use omp_lib
        use, intrinsic :: iso_c_binding ! 和C语言类型互通的模块iso_c_binding
        include 'common.inc'
        REAL, allocatable:: gaxobs(:,:)
        Ngaldim = 4
        Ngalmax = 3
        allocate(gaxobs(Ngaldim,Ngalmax))
        call loaddata(gaxobs) !, Ngaldim,Ngalmax) 
        write(*,*) gaxobs(:,1) 
        !write(*,*) rg(:,0) 
        END

