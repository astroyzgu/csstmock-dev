module load compiler/intel-2018
module load compiler/intel-2018u3
module load mpi 
export OMP_NUM_THREADS=10 
ifort -o mainomp mainomp.f -fopenmp 
./mainomp 


ifort -o test test.f -lcsstmock  
./test
mpirun -np 2 ./test 

ifort -o testomp testomp.f -lcsstmock -fopenmp
./testomp
mpirun -np 2 ./testomp 


