c------------------------------------------------------------------------------
	subroutine dummydata(larr, darr, carr, Nl1, Nd1, Nc1, N2) 
c------------------------------------------------------------------------------

        use iso_c_binding 
        implicit none 
        interface 
		subroutine dummydata_c(
     &                       larr, darr, carr,  
     &                       Nl1, Nd1, Nc1, N2) bind (c)  
        	 use iso_c_binding 
                 implicit none 
       		 integer(c_long) :: Nl1(1), Nd1(1), Nc1(1)
       		 integer(c_long) :: N2(1)
        	 integer(c_long) :: larr(Nl1(1),N2(1))
        	 real(c_double)  :: darr(Nd1(1),N2(1))
                 character(c_char) :: carr(Nc1(1))
		end subroutine dummydata_c
	end interface 

        integer(c_long)  :: Nl1, Nd1, Nc1, N2 
        integer(c_long)  :: Nl1arr(1), Nd1arr(1), Nc1arr(1)
        integer(c_long)  :: N2arr(1) 
        integer(c_long)  :: larr(Nl1,N2)
        real(c_double)   :: darr(Nd1,N2)
        character*256 :: carr(Nc1)
	Nl1arr(1) = Nl1
	Nd1arr(1) = Nd1
	Nc1arr(1) = Nc1
	N2arr(1)  = N2 
	call dummydata_c(larr, darr, carr, 
     &                   Nl1arr, Nd1arr, Nc1arr, N2arr)
	Nl1 = Nl1arr(1)
	Nd1 = Nd1arr(1)
	Nc1 = Nc1arr(1)
	N2  =  N2arr(1) 
c------------------------------------------------------------------------------
	end subroutine dummydata 
c------------------------------------------------------------------------------

c------------------------------------------------------------------------------
	subroutine dummydata_mpi(snap, basedir, 
     &                           darr, larr, nd1, nl1, ngal) 
c------------------------------------------------------------------------------
        use iso_c_binding 
        implicit none 
    
        interface 
		subroutine dummydata_mpi_c(snap, basedir, 
     &              darr, larr, nd1, nl1, ngal) bind (c)
        	use iso_c_binding 
       		integer(c_long)   :: snap(1) 
                character(c_char) :: basedir 
       		integer(c_long) :: nd1(1), nl1(1), ngal(1) 
        	integer(c_long) :: larr( ngal(1), nl1(1) )
        	real(c_double)  :: darr( ngal(1), nd1(1) )
		end subroutine dummydata_mpi_c
	end interface 

        integer(c_long)   :: snaparr(1), snap 
        character(c_char) :: basedir 
        integer(c_long) :: nd1arr(1), nl1arr(1), ngalarr(1) 
        integer(c_long) :: nd1, nl1, ngal 
        integer(c_long) :: larr( ngal, nl1 )
        real(c_double)  :: darr( ngal, nd1 )

	snaparr(1) = snap
	ngalarr(1) = ngal
	nl1arr(1)  = nl1
	nd1arr(1)  = nd1

	call dummydata_mpi_c(snaparr(1), basedir, 
     &        darr, larr, nd1arr, nl1arr, ngalarr)
	snap = snaparr(1)
	ngal = ngalarr(1)
	nl1  = nl1arr(1)
	nd1  = nd1arr(1)
c------------------------------------------------------------------------------
	end subroutine dummydata_mpi 
c------------------------------------------------------------------------------
