python build.py # 生成库文件libplugin.so(linux)和libplugin.dylib(Macos)
gfortran -fPIC -shared -o libcsstmock.so isin_*.f io_*.f libcsstplugin.a # 合并为动态链接库

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


