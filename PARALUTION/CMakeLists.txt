CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7 FATAL_ERROR)

PROJECT(hgf)

SET(HGF_BUILD_EXECUTABLE ON CACHE BOOL "Build HGF Driver Executable.")

SET(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)

### FIND PACKAGES ###
## OpenMP ##
FIND_PACKAGE(OpenMP REQUIRED)

## CUDA ##
FIND_PACKAGE(CUDA REQUIRED)
INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})
INCLUDE(FindCUDA)

## Boost ##
FIND_PACKAGE(Boost COMPONENTS serialization filesystem system REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

## Paralution ##
FIND_PACKAGE(PARALUTION REQUIRED)
INCLUDE_DIRECTORIES(${PARALUTION_INCLUDE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${CUDA_CXX_FLAGS}")

LIST(APPEND CUDA_NVCC_FLAGS "-std=c++11 -O3 -Xcompiler -DSOLVER=1")

SET(LIBRARY_SRCS ../src/hgf.cpp
                 ../src/hgfUtils.cpp
                 ../src/hgfMeshCu.cu
                 ../src/hgfStokesP.cpp
                 ../src/hgfPoreNetworkP.cpp
                 ../src/hgfArrays.cpp
                 ../src/hgfBC.cpp
                 ../src/hgfIB.cpp
                 ../src/hgfPP.cpp)

SET(EXECUTABLE_SRCS ../src/driver.cpp)

SET(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)

CUDA_ADD_LIBRARY(hgflib STATIC ${LIBRARY_SRCS})

CUDA_ADD_EXECUTABLE(hgf ${EXECUTABLE_SRCS})

TARGET_LINK_LIBRARIES( hgf ${CUDA_LIBRARIES}
                           ${CUDA_CUBLAS_LIBRARIES}
                           ${CUDA_cusparse_LIBRARY}
                           ${MAGMA_LIBRARY_DIRS}
                           hgflib
                           ${PARALUTION_LIBRARY}
                           ${Boost_LIBRARIES})
