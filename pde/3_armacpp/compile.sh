#!/bin/sh

# Personal lapack libraries
#g++ example.cpp -o example -O3 -DARMA_DONT_USE_WRAPPER -DARMA_USE_BLAS -DARMA_USE_LAPACK -larmadillo  -lopenblas -llapack
##g++ example.cpp -o example -O3 -larmadillo -DARMA_DONT_USE_WRAPPER -DARMA_USE_BLAS -DARMA_USE_LAPACK -llapack -lblas
##g++ example.cpp -o example -O3 -larmadillo -DARMA_DONT_USE_WRAPPER -DARMA_USE_BLAS -DARMA_USE_LAPACK -lblas

# Default wrapper
#g++ example.cpp -o example -O3 -larmadillo 


# Personal libraries and mpi 
#mpic++ example.cpp -o example -O3 -DARMA_DONT_USE_WRAPPER -DARMA_USE_BLAS -DARMA_USE_LAPACK -larmadillo  -lopenblas -llapack

rm main
cmake .
make
rm -Rf CMakeFiles cmake_install.cmake CMakeCache.txt
./main > data
gnuplot plot.gnu
evince sol.eps

