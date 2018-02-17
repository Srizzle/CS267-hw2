#
# Edison - NERSC
#
# Intel Compilers are loaded by default; for other compilers please check the module list
#
CC = g++
MPCC = mpic++
OPENMP = -fopenmp #Note: this is the flag for Intel compilers. Change this to -fopenmp for GNU compilers. See http://www.nersc.gov/users/computational-systems/edison/programming/using-openmp/
CFLAGS = -O3 -w
LIBS = -lm


TARGETS = serial openmp mpi autograder

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) serial.o common.o -lm
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) autograder.o common.o -lm
openmp: openmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o -lm
mpi: mpi.o common.o
	$(MPCC) -o $@ $(LIBS) mpi.o common.o -lm

autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp -lm
openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp -lm
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp -lm
mpi.o: mpi.cpp mpi.h common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp -lmsss
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) -lm common.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
