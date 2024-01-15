# running in Gravity
module purge 
source activate
conda deactivate 
module load anaconda 
module load compiler/intel-2018
source activate python3 
python build.py
ifort -fPIC -shared -o libcsstmock.so io_jiutian.f libcsstplugin.a
export PYTHONPATH=$PYTHONPATH:$(pwd)
export LIBRARY_PATH=$LIBRARY_PATH:$(pwd)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$(pwd)
ifort  -o test1 test1.f  -lcsstmock  
./test1

