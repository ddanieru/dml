cmake_minimum_required(VERSION 3.0)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

SET (CMAKE_CXX_COMPILER "/usr/bin/c++")
SET (CMAKE_CXX_FLAGS                     "-Wall -g -pg -std=c++11")
SET (CMAKE_CXX_FLAGS_INIT                "-Wall -g -pg -std=c++11")
SET (CMAKE_CXX_FLAGS_DEBUG_INIT          "-g")
SET (CMAKE_CXX_FLAGS_MINSIZEREL_INIT     "-Os -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELEASE_INIT        "-O4 -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g -pg")

# Compiler options: I need to be sure that the the wraper is not used
#SET(GCC_COVERAGE_COMPILE_FLAGS "-O2 -DARMA_DONT_USE_WRAPPER -DARMA_USE_BLAS -DARMA_USE_LAPACK")
#SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
##SET(GCC_COVERAGE_LINK_FLAGS    "-larmdillo")
##SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )


####### (1) with wrapper
#find_package(Armadillo REQUIRED)
#include_directories(${ARMADILLO_INCLUDE_DIRS})
#add_executable(example example.cpp)
#target_link_libraries(example ${ARMADILLO_LIBRARIES})

####### (2) with libraries
#set(ENV{BLA_VENDOR} "OpenBLAS" )
#find_package(BLAS REQUIRED)
set(OpenBLAS_LIBRARIES "/usr/lib/libopenblas.so")
#set(LAPACK_DIR /usr/lib/lapack/)
#find_package(LAPACK REQUIRED)
#set(LAPACK_LIBRARIES " " )
#MESSAGE( STATUS "LAPACK LIBRARIES: " ${LAPACK_LIBRARIES} )

#FIND_LIBRARY(FFTW3_LIBRARIES NAMES fftw3 libfftw3 PATHS /usr/lib )
find_package(FFTW REQUIRED)

set(ARMADILLO_INCLUDE_DIR /usr/local/include )
set(ARMADILLO_LIBRARY /usr/local/lib/libarmadillo.so)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
add_executable(main main.cpp)
#target_link_libraries(main ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${OpenBLAS_LIBRARIES} ${ARMADILLO_LIBRARIES})
target_link_libraries(main ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${OpenBLAS_LIBRARIES} ${FFTW_LIBRARIES})
#target_link_libraries(example ${OpenBLAS_LIBRARIES})



