rm -rf *.txt
rm -rf log
# rm -rf output*

./serial -n 500
export OMP_NUM_THREADS=1
./openmp -n 500  -no  -s openmp.txt
export OMP_NUM_THREADS=2
./openmp -n 500  -no  -s openmp.txt
export OMP_NUM_THREADS=4
./openmp -n 500  -no  -s openmp.txt
export OMP_NUM_THREADS=8
./openmp -n 500  -no  -s openmp.txt
export OMP_NUM_THREADS=16
./openmp -n 500  -no  -s openmp.txt
export OMP_NUM_THREADS=24
./openmp -n 500  -no  -s openmp.txt
export OMP_NUM_THREADS=32
./openmp -n 500  -no  -s openmp.txt


#./autograder -v openmp -s openmp.txt


# rm openmp.txt
# ./openmp -n 500 -no -s openmp.txt
# ./openmp -n 1000 -no -s openmp.txt
# ./openmp -n 2000 -no -s openmp.txt
# ./openmp -n 4000 -no -s openmp.txt
# ./openmp -n 8000 -no -s openmp.txt
# ./autograder -v openmp -s openmp.txt
~                                             