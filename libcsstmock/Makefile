FF=ifort
FFFLAG= -fPIC -shared 
LIB_DIR=$(PWD)

all: io_*.f 
	python build.py # 生成静态链接库 
	$(FF) $(FFFLAG) -o libcsstmock.so $^  libcsstplugin.a # 合并为动态链接库

