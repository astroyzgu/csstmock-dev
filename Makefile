.PHONY: test build

build-gravity: build build-ifort test 

build-macos: build build-gfortran test 

build:
	python -m pip install -e ./

build-ifort: 
	make -C ./libcsstmock

build-gfortran:
	make -C ./libcsstmock

test: 
	pytest 

