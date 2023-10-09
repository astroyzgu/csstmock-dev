c------------------------------------------------------------------------------
	subroutine isin_survey(ra, dec, n, survey, veto)
c------------------------------------------------------------------------------
        !> Description 
        !> -----------
        !> is in survey, or not?  
        !> 
        !> Syntax
        !> -----------
        !> call isin_survey(ra, dec, n, 'lsdr9-ngc', veto)
        !>  
        !> Argument(s)
        !> -----------
        !> ra(in):  real(4), RA
	!> 
        !> dec(in): real(4), DEC
	!> 
        !> n(in):  integer(8), the number of input coordinates 
        !>   
        !> survey(in): name of the embedded survey 
        !>       []
	!> 
	!> veto(out): real(4)
	!>        if 1, in survey. if 0, not in survey. 

        use iso_c_binding 
        implicit none 
    
        interface
                subroutine isin_survey_c(          
     &              ra, dec, n, survey, veto) bind (c)
                use iso_c_binding 
                integer(c_long) :: n(1)
                real(c_float)  :: ra(n(1)), dec(n(1)) 
                real(c_float)  :: veto(n(1))
                character(c_char) :: survey 
                end subroutine isin_survey_c
        end interface 


        integer(c_long) :: n
        integer(c_long) :: narr(1)
        real(c_float)  :: ra(n), dec(n) 
        real(c_float)  :: veto(n)
        character(c_char) :: survey 
        
        narr(1)    = n 

        call isin_survey_c(ra, dec, narr, survey, veto)
        n       = narr(1)
c------------------------------------------------------------------------------
        end subroutine isin_survey 
c------------------------------------------------------------------------------
