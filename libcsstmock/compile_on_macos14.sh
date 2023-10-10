python build.py # 生成库文件libplugin.so(linux)和libplugin.dylib(Macos)

export PYTHONPATH=$PYTHONPATH:$(pwd)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$(pwd)
export LIBRARY_PATH=$LIBRARY_PATH:$(pwd)

gfortran -fPIC -shared -o libcsstmock.so isin_*.f io_*.f libcsstplugin.a # 合并为动态链接库

#gfortran -o test/test0 test/test0.f  -lcsstmock 
#./test/test0 
#gfortran -o testomp0 testomp.f  -lcsstmock -fopenmp 
#./testomp0 

