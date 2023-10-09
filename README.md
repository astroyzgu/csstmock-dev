##### INSTALL ON GRAVITY 

``` 
# install package

module load anaconda 
module load compiler/intel-2018 # load ifort 
conda activate python3          # activate defualt python3
make all
```

```
# if use libcsstmock, need to set environment variable, 
#
# >> e.g., 

echo PYTHONPATH=\$PYTHONPATH:$(pwd)/libcsstmock >> ~/.bash_profile  
echo LIBRARY_PATH=\$LIBRARY_PATH:$(pwd)/libcsstmock >> ~/.bash_profile
echo LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$(pwd)/libcsstmock >> ~/.bash_profile
echo DYLD_LIBRARY_PATH=\$DYLD_LIBRARY_PATH:$(pwd)/libcsstmock >> ~/.bash_profile 
source ~/.bash_profile 

# then, -lcsstmock works. 
# >> e.g, to read the fullsky_z2 lightcone  
cd tests/test_libcsstmock/
ifort -o test4 test4.f -lcsstmock 
test4
``` 

###### subroutine Description 

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


