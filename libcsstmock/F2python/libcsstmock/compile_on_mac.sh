python build.py # 生成库文件libplugin.so(linux)和libplugin.dylib(Macos)

export PYTHONPATH=$PYTHONPATH:$(pwd)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$(pwd)
export LIBRARY_PATH=$LIBRARY_PATH:$(pwd)

gfortran -fPIC -shared -o libcsstmock.so io_*.f libcsstplugin.a # 合并为动态链接库

#install_name_tool -add_rpath /opt/anaconda2/envs/py37/lib ./test 
#gfortran -o test0 test.f  -lcsstmock 
#./test0 
#gfortran -o testomp0 testomp.f  -lcsstmock -fopenmp 
#./testomp0 

