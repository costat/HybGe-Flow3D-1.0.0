CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7 FATAL_ERROR)

PROJECT(hgf)

SET(HGF_BUILD_EXECUTABLE ON CACHE BOOL "Build HGF Driver Executable.")

### FIND PACKAGES ###
## OpenMP ##
FIND_PACKAGE(OpenMP REQUIRED)
SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

## CUDA ##
FIND_PACKAGE(CUDA REQUIRED)
INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})
INCLUDE(FindCUDA)

## Boost ##
FIND_PACKAGE(Boost REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

## MAGMA ##

## BLAS/LAPACK ##

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

LIST(APPEND CMAKE_CXX_FLAGS "-std=c++11 -O3")

SET(LIBRARY_SRCS ../src/hgfM.cpp
                 ../src/hgfStokesM.cpp
                 ../src/hgfPoreNetworkM.cpp
                 ../src/hgfArrays.cpp
                 ../src/hgfBC.cpp
                 ../src/hgfIB.cpp
                 ../src/hgfPPM.cpp)

SET(LIBRARYCU_SRCS ../src/hgfMeshCu.cu)

SET(EXECUTABLE_SRCS ../src/driver.cpp)

SET(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)

CUDA_ADD_LIBRARY(hgflibcu STATIC ${LIBRARYCU_SRCS})

ADD_LIBRARY(hgflib STATIC ${LIBRARY_SRCS})

CUDA_ADD_EXECUTABLE(hgf ${EXECUTABLE_SRCS})
TARGET_LINK_LIBRARIES(hgf ${CUDA_LIBRARIES}
                          ${CUDA_CUBLAS_LIBRARIES}
                          ${CUDA_cusparse_LIBRARY}
                          ${Boost_LIBRARIES}
                          hgflibcu
                          hgflib)
