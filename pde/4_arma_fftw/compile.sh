#!/bin/sh

# Pre compilation
rm main

# Personal lapack libraries
#g++ --std=c++11 -Wall -g -pg main.cpp -o main -O2 -DARMA_DONT_USE_WRAPPER -DARMA_USE_BLAS -DARMA_USE_LAPACK  -lopenblas -llapack -lfftw3 

##g++ example.cpp -o example -O3 -larmadillo -DARMA_DONT_USE_WRAPPER -DARMA_USE_BLAS -DARMA_USE_LAPACK -llapack -lblas
##g++ example.cpp -o example -O3 -larmadillo -DARMA_DONT_USE_WRAPPER -DARMA_USE_BLAS -DARMA_USE_LAPACK -lblas

# Default wrapper
#g++ example.cpp -o example -O3 -larmadillo 


# Personal libraries and mpi 
#mpic++ example.cpp -o example -O3 -DARMA_DONT_USE_WRAPPER -DARMA_USE_BLAS -DARMA_USE_LAPACK -larmadillo  -lopenblas -llapack

#   Cmake
cmake .
make
rm -Rf CMakeFiles cmake_install.cmake CMakeCache.txt


./main > data
gnuplot plot.gnu
#evince sol.eps


# Analyze
#gprof main gmon.out > analysis.txt
#valgrind --tool=callgrind ./main
#valgrind --tool=massif ./main
#ms_print massif.out.???? > heap.txt

