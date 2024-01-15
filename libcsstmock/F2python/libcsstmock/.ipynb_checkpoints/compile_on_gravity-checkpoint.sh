# running in Gravity
conda activate csstmock
module load compiler/intel-2018
python build.py # 生成静态链接库 
ifort -fPIC -shared -o libcsstmock.so io_jiutian.f libcsstplugin.a # 合并为动态链接库
# 添加路径 
export PYTHONPATH=$PYTHONPATH:$(pwd)
export LIBRARY_PATH=$LIBRARY_PATH:$(pwd)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$(pwd)

#ifort  -o test2 test2.f  -lcsstmock 
#./test2

ifort  -o test1 test1.f  -lcsstmock 
./test1
