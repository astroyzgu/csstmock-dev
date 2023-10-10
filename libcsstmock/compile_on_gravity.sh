# if running on Gravity: 
module load anaconda 
module load compiler/intel-2018
source activate python3  

python build.py # 生成静态链接库 
ifort -fPIC -shared -o libcsstmock.so i*.f libcsstplugin.a # 合并为动态链接库

export PYTHONPATH=$PYTHONPATH:$(pwd)
export LIBRARY_PATH=$LIBRARY_PATH:$(pwd)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$(pwd)

# 添加路径 
echo 
echo 
echo 
echo "# >>> paste the following into ~/.bash_profile to add PATH"
echo 
echo "export PYTHONPATH=\$PYTHONPATH:$(pwd)" 
echo "export LIBRARY_PATH=\$LIBRARY_PATH:$(pwd)"
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$(pwd)"
echo "export DYLD_LIBRARY_PATH=\$DYLD_LIBRARY_PATH:$(pwd)"
echo 
echo "# >>> source ~/.bash_profile" 
echo 
echo "END"

#ifort  -o test0 test0.f  -lcsstmock 
#./test0

#ifort  -o test2 test2.f  -lcsstmock 
#./test2
#ifort  -o test1 test1.f  -lcsstmock 
#./test1

#ifort  -o test3 test3.f  -lcsstmock 
#./test3

#ifort  -o test4 test4.f  -lcsstmock 
#./test4
