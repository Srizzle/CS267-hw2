#
# Edison - NERSC 
#
# Intel Compilers are loaded by default; for other compilers please check the module list
#
CC = CC
MPCC = CC
OPENMP = -openmp #Note: this is the flag for Intel compilers. Change this to -fopenmp for GNU compilers. See http://www.nersc.gov/users/computational-systems/edison/programming/using-openmp/
CFLAGS = -O3
LIBS =


TARGETS = serial openmp mpi autograder

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) serial.o common.o
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) autograder.o common.o
openmp: openmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o
mpi: mpi.o common.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) $(OPENMP) mpi.o common.o

autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) $(OPENMP) autograder.cpp
openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) $(OPENMP) openmp.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) $(OPENMP) serial.cpp
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) $(OPENMP) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS)  $(OPENMP) common.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.error *.txt
