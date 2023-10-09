.PHONY: test build

all: build build-lib test 

build:
	python -m pip install -e ./

build-lib: 
	make -C ./libcsstmock
test: 
	pytest 

