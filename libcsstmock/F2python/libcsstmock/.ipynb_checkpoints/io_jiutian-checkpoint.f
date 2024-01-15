c------------------------------------------------------------------------------
        subroutine readsnapshot_jiutian(
     &             isnap,sub_in,idsub_in,Nd1,Nd2max,Nin) 
c------------------------------------------------------------------------------
        use iso_c_binding 
        implicit none 
    
        interface 
                subroutine readsnapshot_jiutian_c(isnap,          
     &              sub_in,idsub_in, Nd1,Nd2max,Nin) bind (c)
                use iso_c_binding 
                integer(c_long) :: Nin(1), Nd2max(1)
                integer(c_int)  :: Nd1(1), isnap(1)
                integer(c_long) :: idsub_in( 3,   Nd2max(1) )
                real(c_double)  :: sub_in(Nd1(1), Nd2max(1) )
                end subroutine readsnapshot_jiutian_c
        end interface 

        integer(c_int)  :: Nd1arr(1), isnaparr(1)
        integer(c_long) :: Ninarr(1), Nd2maxarr(1) 
        
        integer(c_long) :: Nin,Nd2max 
        integer(c_int)  :: Nd1, isnap 
        integer(c_long) :: idsub_in( 3, Nd2max )
        real(c_double)  :: sub_in( Nd1, Nd2max )
        
        isnaparr(1) = isnap
        Nd1arr(1)   = Nd1
        Nd2maxarr(1)= Nd2max
        Ninarr(1)   = Nin 
        call readsnapshot_jiutian_c(isnaparr,          
     &              sub_in, idsub_in, Nd1arr,Nd2maxarr,Ninarr)
     
        isnap = isnaparr(1)
        Nd1   = Nd1arr(1)
        Nd2max  = Nd2maxarr(1)
        Nin  = Ninarr(1)
c------------------------------------------------------------------------------
        end subroutine readsnapshot_jiutian 
c------------------------------------------------------------------------------

c------------------------------------------------------------------------------
        subroutine readlightcone_jiutian(
     &             sub_in,idsub_in, Nd1,Nd2max,Nin) 
c------------------------------------------------------------------------------
        use iso_c_binding 
        implicit none 
    
        interface 
                subroutine readlightcone_jiutian_c(          
     &              sub_in,idsub_in, Nd1,Nd2max,Nin) bind (c)
                use iso_c_binding 
                integer(c_long) :: Nin(1), Nd2max(1)
                integer(c_int)  :: Nd1(1) 
                integer(c_long) :: idsub_in( 2,   Nd2max(1) )
                real(c_double)  :: sub_in(Nd1(1), Nd2max(1) )
                end subroutine readlightcone_jiutian_c
        end interface 

        integer(c_int)  :: Nd1arr(1), isnaparr(1)
        integer(c_long) :: Ninarr(1), Nd2maxarr(1) 
        
        integer(c_long) :: Nin,Nd2max 
        integer(c_int)  :: Nd1, isnap 
        integer(c_long) :: idsub_in( 2, Nd2max )
        real(c_double)  :: sub_in( Nd1, Nd2max )
        
        Nd1arr(1)   = Nd1
        Nd2maxarr(1)= Nd2max
        Ninarr(1)   = Nin 
        call readlightcone_jiutian_c(          
     &              sub_in, idsub_in, Nd1arr,Nd2maxarr,Ninarr)
     
        Nd1   = Nd1arr(1)
        Nd2max  = Nd2maxarr(1)
        Nin  = Ninarr(1)
c------------------------------------------------------------------------------
        end subroutine readlightcone_jiutian 
c------------------------------------------------------------------------------