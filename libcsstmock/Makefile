FFFLAG= -fPIC -shared 
LIB_DIR=$(PWD)
build-gravity: io_*.f 
	FF=ifort
	python build.py # 生成静态链接库 
	$(FF) $(FFFLAG) -o libcsstmock.so $^  libcsstplugin.a # 合并为动态链接库

build-macos: io_*.f 
	FF=gfortran
	python build.py # 生成静态链接库 
	$(FF) $(FFFLAG) -o libcsstmock.so $^  libcsstplugin.a # 合并为动态链接库

