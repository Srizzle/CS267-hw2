rm -rf *.txt
rm -rf log
# rm -rf output*

./serial -n 50   -s serial.txt
./serial -n 100  -s serial.txt
./serial -n 150  -s serial.txt
./serial -n 250  -s serial.txt
./serial -n 300  -s serial.txt
./serial -n 350  -s serial.txt
./serial -n 450  -s serial.txt
./serial -n 500  -s serial.txt
./serial -n 650  -s serial.txt
./serial -n 750  -s serial.txt
./serial -n 850  -s serial.txt
./serial -n 900  -s serial.txt
./serial -n 1050  -s serial.txt
./serial -n 1100  -s serial.txt
./serial -n 1150  -s serial.txt
./serial -n 1250  -s serial.txt
./serial -n 1300  -s serial.txt
./serial -n 1350  -s serial.txt
./serial -n 2100  -s serial.txt
./serial -n 2150  -s serial.txt
./serial -n 2250  -s serial.txt
./serial -n 2300  -s serial.txt
./serial -n 2350  -s serial.txt
# ./autograder -v serial -s serial.txt


# rm serial.txt
# ./serial -n 500 -no -s serial.txt
# ./serial -n 1000 -no -s serial.txt
# ./serial -n 2000 -no -s serial.txt
# ./serial -n 4000 -no -s serial.txt
# ./serial -n 8000 -no -s serial.txt
# ./autograder -v serial -s serial.txt
