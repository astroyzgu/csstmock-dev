Gravity location: /home/yzgu/csstmock/libcsstmock/readme.md

The usage of csstmock library, which can call python program in fortan. 
In current stage, the following subroutine is embedded:   
   
   1. readsnap_jiutian(snap, darr, larr, nd1, ngalmax, ngal):
   nd1 = 10, ngalmax = 500000000, (5亿， 460267974 subhaloσ in factor) 
   
   larr (int*8)
    1. 'host_id': host halo ID in present snapshot
    2. 'sub_id':  subhalo ID 
    3. 'host_nextid': host halo ID in next snapshot (-99 means not found)
    
   darr (real*8) 
    1-3. 'sub_pos': position (x, y, z) 
    4-6. 'sub_velocity': velocity (vx, vy, vz): 
    8. 'host_mass': host halo mass
    8. 'sub_mass': subhalo mass:  
    9. 'PeakMass': peak mass  
    10. 'sub_velmax': maximum circular velocity of subhalo


   2. readlightcone_jiutian(darr, larr, nd1, ngalmax, ngal)
   nd1 = 10, ngalmax = 5000000, (5百万, 2018662 in factor)
   larr (int*8)
    1. 'group_nr' halo ID from fof, 'group_nr' 
    2. 'TrackID': subhalo id

   darr (real*8)
    1.   'ra': Right Ascension
    2.   'dec': Declination 
    3.   'd_comoving': comoving distance
    4.   'redshift_true': cosmology redshift, z_cos
    5.   'redshift_obs': observed redshift, z_obs 
    6.   'v_lfs': line-of-sight velocity 
    7.   'group_mass':  host halo mass from fof, 'group_mass'
    8.   'shMbound': subhalo mass 
    9.   'PeakMass': peak mass  
    10.   'vmax':     max circular velocity 
    
Writen by: Li, Qingyang & Yizhou Gu 
c=============================================================================
c see compile_on_gravity.sh 

# 编译动态链接库
module load compiler/intel-2018 # 加载ifort(gravity上)
python build.py                 # 生成静态链接库 libcsstplugin.a 
gfortran -fPIC -shared -o libcsstmock.so io_dummy.f io_jiutian.f libcsstplugin.a # 合并为动态链接库 

# 添加路径:
export PYTHONPATH=$PYTHONPATH:$(pwd)
export LIBRARY_PATH=$LIBRARY_PATH:$(pwd)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$(pwd)

# 编译程序 readsnapshot_jiutian的例子
ifort -o ./test1 test1.f -lcsstmock
./test1 

# 编译程序 readlightcone_jiutian的例子
ifort -o ./test2 test2.f -lcsstmock
./test2 


c the subroutine name in each part of the test1.f 
c 
c main program          io_jiutian.f    mymodule.py/plugin.h       io_jiutian.py (编辑这里就行了)
c exec in test    <-   subroutine      <-     cffi wrapper      <- origin  python 
c  test.f             readsnapshot      readsnapshot_jiutian_c    collect4csstmock 


