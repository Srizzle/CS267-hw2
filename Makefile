# Load CUDA using the following command
# module load cuda
#
CC = nvcc  -Wno-deprecated-gpu-targets
CFLAGS = -O3 -arch=compute_37 -code=sm_37
NVCCFLAGS = -O3 -arch=compute_37 -code=sm_37
WALL = -Wno-deprecated-gpu-targets
LIBS = 

TARGETS = serial gpu autograder

all:	$(TARGETS)

serial: 
	$(CC) -o serial $(CFLAGS) serial.cu common.cu
gpu: gpu.o common.o
	$(CC) -o $@ $(NVCCLIBS) $(WALL) gpu.o common.o
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) $(WALL) autograder.o common.o

autograder.o: autograder.cu common.h
	$(CC) -c $(CFLAGS) $(WALL) autograder.cu
gpu.o: gpu.cu gpu.h common.h
	$(CC) -c $(NVCCFLAGS) $(WALL) gpu.cu
common.o: common.cu common.h
	$(CC) -c $(CFLAGS) $(WALL) common.cu

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt