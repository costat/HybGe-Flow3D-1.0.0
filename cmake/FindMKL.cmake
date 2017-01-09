# - Find Intel MKL
# Find the MKL libraries
#
# Options:
#
#   MKL_SHARED        :   use shared library
#   MKL_SEQUENTIAL    :   use sequential library
#
# This module defines the following variables:
#
#   MKL_FOUND            : True if MKL_INCLUDE_DIR are found
#   MKL_INCLUDE_DIR      : where to find mkl.h, etc.
#   MKL_INCLUDE_DIRS     : set when MKL_INCLUDE_DIR found
#   MKL_LIBRARIES        : the library to link against.


include(FindPackageHandleStandardArgs)

set(MKL_ROOT $ENV{MKLROOT} CACHE PATH "Folder contains MKL")
if (WIN64)
  set(IOMP_ROOT $ENV{IOMPROOT} CACHE PATH "Folder contains IOMP")
endif()

# Find include dir
find_path(MKL_INCLUDE_DIR mkl.h
    PATHS ${MKL_ROOT}/include)

# Find libraries
# suffix
set(_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

if(WIN32)
  if(MKL_SHARED)
    set(CMAKE_FIND_LIBRARY_SUFFIXES _dll.lib)
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
  endif()
else()
  if(MKL_SHARED)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  endif()
endif()

## INTERFACE ##
set(MKL_INTERFACE_LIBNAME mkl_intel_lp64)
## THREADING ##
if (MKL_SEQUENTIAL)
  set(MKL_THREADING_LIBNAME mkl_sequential)
else()
  set(MKL_THREADING_LIBNAME mkl_intel_thread)
endif()
## COMPUTATION ##
set(MKL_CORE_LIBNAME mkl_core)

find_library(MKL_INTERFACE_LIBRARY ${MKL_INTERFACE_LIBNAME}
    PATHS ${MKL_ROOT}/lib/intel64/)

find_library(MKL_THREADING_LIBRARY ${MKL_THREADING_LIBNAME}
    PATHS ${MKL_ROOT}/lib/intel64/)

find_library(MKL_CORE_LIBRARY ${MKL_CORE_LIBNAME}
    PATHS ${MKL_ROOT}/lib/intel64/)

if(WIN32)
  set(IOMP_LIBNAME iomp5md)
  find_library(INTEL_OMP_LIBRARY ${IOMP_LIBNAME}
  PATHS ${IOMP_ROOT}/lib/intel64/)
  set(MKL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${INTEL_OMP_LIBRARY})

  find_package_handle_standard_args(MKL DEFAULT_MSG
  MKL_INCLUDE_DIR MKL_LIBRARY MKL_INTERFACE_LIBRARY MKL_THREADING_LIBRARY MKL_CORE_LIBRARY INTEL_OMP_LIBRARY)
else()
  set(MKL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY})

  find_package_handle_standard_args(MKL DEFAULT_MSG
    MKL_INCLUDE_DIR MKL_LIBRARY MKL_INTERFACE_LIBRARY MKL_THREADING_LIBRARY MKL_CORE_LIBRARY)
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES ${_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

if(MKL_FOUND)
    set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
    set(MKL_LIBRARIES ${MKL_LIBRARY})
endif()
