.PHONY: test build

# echo "if running on Gravity"
# echo "module load anaconda"
# echo "module load compiler/intel-2018 # load ifort"
# echo "conda activate python3          # activate defualt python3 env" 

all: build test 

build:
	python -m pip install -e ./

test: 
	pytest 

