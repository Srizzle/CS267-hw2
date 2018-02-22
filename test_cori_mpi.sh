srun -n 1 -N 1 ./serial -n 500 -no -s mpi.txt
srun -n 1 -N 1 ./mpi -p 1 -n 500 -no -s mpi.txt
srun -n 2 -N 2 ./mpi -p 2 -n 500 -no -s mpi.txt
srun -n 4 -N 4 ./mpi -p 4 -n 500 -no -s mpi.txt
srun -n 6 -N 6 ./mpi -p 6 -n 500 -no -s mpi.txt
srun -n 2 -N 2 ./mpi -p 2 -n 1000 -no -s mpi.txt
srun -n 4 -N 4 ./mpi -p 4 -n 2000 -no -s mpi.txt
srun -n 6 -N 6 ./mpi -p 6 -n 3000 -no -s mpi.txt