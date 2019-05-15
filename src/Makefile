ROOT_DIR = $(PWD)
SERIAL_DIR = $(ROOT_DIR)/serial
OMP_DIR = $(ROOT_DIR)/openmp
MPI_DIR = $(ROOT_DIR)/mpi

TARGET_SERIAL = simpar
TARGET_OMP = simpar-omp
TARGET_MPI = simpar-mpi

CC =gcc 
MPICC =mpicc
COMPP =kinst-ompp 
CFLAGS =-c -g3 -Wall -pedantic -Werror
CFLAGSLINK =-lm
COMPFLAGS =-fopenmp

TEST_DIR = $(ROOT_DIR)/tests
Tests:=$(shell cd $(TEST_DIR); ls  | egrep '^test[0-9]+$$' | sort -n  )
program = serial
test = test01
ompp = false

all: serial omp mpi

serial.o: 
	$(CC) $(CFLAGS) $(SERIAL_DIR)/$(TARGET_SERIAL).c -o $(SERIAL_DIR)/$(TARGET_SERIAL).o 

serial: serial.o
	$(CC) $(SERIAL_DIR)/$(TARGET_SERIAL).o $(CFLAGSLINK) -o $(SERIAL_DIR)/$(TARGET_SERIAL)

omp.o: 
ifeq ($(ompp), true)
	$(COMPP) $(CC) $(CFLAGS) $(COMPFLAGS) $(OMP_DIR)/$(TARGET_OMP).c -o $(OMP_DIR)/$(TARGET_OMP).o 
else
	$(CC) $(CFLAGS) $(COMPFLAGS) $(OMP_DIR)/$(TARGET_OMP).c -o $(OMP_DIR)/$(TARGET_OMP).o
endif

omp: omp.o
ifeq ($(ompp), true)
	$(COMPP) $(CC) $(OMP_DIR)/$(TARGET_OMP).o $(CFLAGSLINK) $(COMPFLAGS) -o $(OMP_DIR)/$(TARGET_OMP)
else
	$(CC) $(OMP_DIR)/$(TARGET_OMP).o $(CFLAGSLINK) $(COMPFLAGS) -o $(OMP_DIR)/$(TARGET_OMP)
endif

mpi.o: 
	$(MPICC) $(CFLAGS) $(MPI_DIR)/$(TARGET_MPI).c -o $(MPI_DIR)/$(TARGET_MPI).o 

mpi: mpi.o
	$(MPICC) $(MPI_DIR)/$(TARGET_MPI).o $(CFLAGSLINK) -o $(MPI_DIR)/$(TARGET_MPI)

# To run a test have 2 parameters: test is the name of the test and program is the program to use to test
# The program can be: serial, omp, mpi
test: 
ifeq ($(program), serial)
	$(SERIAL_DIR)/$(TARGET_SERIAL) $(shell cat $(TEST_DIR)/$(test)) > $(TEST_DIR)/output/$(test).myout
endif
ifeq ($(program), omp)
	$(OMP_DIR)/$(TARGET_OMP) $(shell cat $(TEST_DIR)/$(test)) > $(TEST_DIR)/output/$(test).myout
endif
ifeq ($(program), mpi)
	$(MPI_DIR)/$(TARGET_MPI) $(shell cat $(TEST_DIR)/$(test)) > $(TEST_DIR)/output/$(test).myout
endif

	@if diff -s $(TEST_DIR)/output/$(test).out $(TEST_DIR)/output/$(test).myout; \
		then echo $(test): PASSED; \
		else echo $(test): FAILED; \
	fi
	rm $(TEST_DIR)/output/$(test).myout

testall: 
	$(foreach x, $(Tests), $(MAKE) test program=$(program) test=$x;)

rompp:
	find . -type f -name 'opari.rc' -delete
	$(MAKE) test

clean:
	find . -type f -name '*.o' -delete
	find . -type f -name '*.myout' -delete
	find . -type f -name '$(TARGET_SERIAL)' -delete
	find . -type f -name '$(TARGET_OMP)' -delete
	find . -type f -name '$(TARGET_OMP).c.*' -delete
	find . -type f -name '$(TARGET_OMP).mod.c' -delete
	find . -type f -name 'opari*' -delete
	find . -type f -name '*.txt' -delete
	find . -type f -name '$(TARGET_MPI)' -delete