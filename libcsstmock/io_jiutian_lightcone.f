c------------------------------------------------------------------------------
        subroutine fullsky_z2_filenum_jiutian(
     &             snapnum, filenum)
c------------------------------------------------------------------------------
        !> Description 
        !> -----------
        !> the file number of fullsky_z2 lightcone at the given snapshot 
        !> 
        !> Syntax
        !> -----------
        !> call fullsky_z2_filenum_jiutian(127, filenum)
        !> 
        !> Argument(s)
        !> -----------
        !> snapnum(in):  integer(4), the given snapshot number [127-65]
        !> 
        !> filenum(out): integer(4), the number of divided files
        !> 

        use iso_c_binding
        implicit none
        interface
                subroutine fullsky_z2_filenum_jiutian_c(
     &              snapnum, filenum) bind (c)
                use iso_c_binding
                integer(c_int) :: snapnum(1)
                integer(c_int) :: filenum(1)
                end subroutine fullsky_z2_filenum_jiutian_c
        end interface 
        integer(c_int) :: snapnum, snapnumarr(1)
        integer(c_int) :: filenum, filenumarr(1)
        snapnumarr(1) = snapnum 
        call fullsky_z2_filenum_jiutian_c(snapnumarr, filenumarr)
        filenum = filenumarr(1)
c------------------------------------------------------------------------------
        end subroutine fullsky_z2_filenum_jiutian
c------------------------------------------------------------------------------

c------------------------------------------------------------------------------
        subroutine fullsky_z2_jiutian(
     &             snapnum, ifile, larr, darr, nmax)
c------------------------------------------------------------------------------
        !> Description 
        !> -----------
        !> read the fullsky_z2 lightcone of jiutian simulation
        !> 
        !> Syntax
        !> -----------
        !> call fullsky_z2_jiutian(127, larr, darr, nsub)
        !> 
        !> Argument(s)
        !> -----------
        !> snapnum(in): integer(4), the given snapshot number [127-65]
        !> 
        !> larr(out): integer(8), long type data
        !>      (1-5) 'trackID', 'hosthaloID', 'rank', 'snapNum','group_nr'
        !>
        !> darr(out): real(4), float type data
        !>      ( 1-7 ) 'posx', 'posy', 'posz', 'velx', 'vely', 'velz', 'v_lfs', 
        !>      ( 8-14) 'shMbound', 'd_comoving', 'ra', 'dec', 'vmax', 'PeakMass', 'PeakVmax',
        !>      (15-18) 'shBoundM200Crit', 'redshift_true', 'redshift_obs', 'shMbound_at_ac',
        !>      (19-20) 'group_mass', 'group_m_mean200'
        !> 
        !> nmax(inout): integer(8), max number of subhalos
        !>      
        
        use iso_c_binding 
        implicit none 
        interface 
                subroutine fullsky_z2_jiutian_c(         
     &              snapnum, ifile, larr, darr, nmax) bind (c)
                use iso_c_binding 
                integer(c_int)  :: snapnum(1), ifile(1)
                integer(c_long) :: nmax(1)
                integer(c_long) :: larr( 5, nmax(1) )
                real(c_float)  :: darr(20, nmax(1) )
                end subroutine fullsky_z2_jiutian_c 
        end interface 

        integer(c_int)  :: snapnum, snapnumarr(1)
        integer(c_int)  :: ifile,   ifilearr(1)
        integer(c_long) :: nmax   , nmaxarr(1)
        integer(c_long) :: Nl1, Nd1, Nc1, N2
        integer(c_long) :: larr( 5,nmax)
        real(c_float)  :: darr(20,nmax)
        snapnumarr(1) = snapnum
        ifilearr(1)   = ifile
        nmaxarr(1)    = nmax
        call fullsky_z2_jiutian_c(snapnumarr, ifilearr,
     &                 larr, darr, nmaxarr) 
        snapnum = snapnumarr(1)
        nmax    = nmaxarr(1)
c------------------------------------------------------------------------------
        end subroutine fullsky_z2_jiutian 
c------------------------------------------------------------------------------

