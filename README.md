# CPD - 2018/2019
Parallel and Distributed Computing Project - Particles Simulator

## Introduction
The project is a sequential and two parallel implementations of a simulator of particles moving in free 2D space. It uses parallel programming on both UMA and multicomputer systems, using OpenMP and MPI, respectively.

## Requirements
You must have installed the following tools:
- **gcc** (A version that supports **_#pragma_** directive with OpenMP pack)

## How to Compile
The compilation of the implementations can be made using the Makefile in **/src**. 

To compile, you just need to issue the following command:
```
    make [implementation]
```

Where **\<implementation> := serial | omp | mpi** representing Sequential, OpenMP and MPI respectively. If not present, all implementations will be compiled.

To clean the object files and the binaries, just run:
```
    make clean
```

## Run the Tests
There are a set of 4 tests that were provided in the statement of the project.

To run all tests:
```
    make testall [program=implementation]
```

To run a specific test:
```
    make test [test=testfile] [program=implementation]
```

## How to use the profiler
In order to profile the parallel versions of the project you need to install [ompP](http://ompp-tool.com/), read the [manual](http://www.ompp-tool.com/downloads/ompp-manual.pdf).
Don't forget to change the file **Makefile.defs**. The makefile expects the existence of the command **_kinst-ompp_**.

The Makefile integrates the use of the profiler. So, in order to get this feature, you need first to compile the code with the followiing command:
```
    make [program=implementation] ompp=true
```
Then, you just need to run the test you want to profile with the command:
```
    make rompp [test=testfile] [program=implementation]
```

## Notes
- By default **program=_serial_** and **test=_test01_**