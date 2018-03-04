# Load CUDA using the following command
# module load cuda
#
CC = nvcc  -Wno-deprecated-gpu-targets
CFLAGS = -O3 -arch=compute_37 -code=sm_37
NVCCFLAGS = -O3 -arch=compute_37 -code=sm_37
WALL = -Wno-deprecated-gpu-targets
LIBS = 

TARGETS = serial serial_dummy gpu gpu_dummy autograder

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) $(WALL) serial.o common.o
serial_dummy: serial_dummy.o common.o
	$(CC) -o $@ $(LIBS) $(WALL) serial_dummy.o common.o
gpu: gpu.o common.o
	$(CC) -o $@ $(NVCCLIBS) $(WALL) gpu.o common.o
gpu_dummy: gpu_dummy.o common.o
	$(CC) -o $@ $(NVCCLIBS) $(WALL) gpu_dummy.o common.o
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) $(WALL) autograder.o common.o

serial.o: serial.cu common.h
	$(CC) -c $(CFLAGS) $(WALL) serial.cu
serial_dummy.o: serial_dummy.cu common.h
	$(CC) -c $(CFLAGS) $(WALL) serial_dummy.cu
autograder.o: autograder.cu common.h
	$(CC) -c $(CFLAGS) $(WALL) autograder.cu
gpu.o: gpu.cu gpu.h common.h
	$(CC) -c $(NVCCFLAGS) $(WALL) gpu.cu
gpu_dummy.o: gpu_dummy.cu common.h
	$(CC) -c $(NVCCFLAGS) $(WALL) gpu_dummy.cu
common.o: common.cu common.h
	$(CC) -c $(CFLAGS) $(WALL) common.cu

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
