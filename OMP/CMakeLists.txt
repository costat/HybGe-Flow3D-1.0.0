CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7 FATAL_ERROR)

PROJECT(hgf)

SET(HGF_BUILD_PYTHON_HOOKS ON CACHE BOOL "Build HGF Python Hooks Using Swig.")
SET(HGF_BUILD_EXECUTABLE ON CACHE BOOL "Build HGF Driver Executable.")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(LIBRARY_SRCS ../src/hgf.cpp
                 ../src/hgfMesh.cpp
                 ../src/hgfStokes.cpp
                 ../src/hgfPoreNetwork.cpp
                 ../src/hgfArrays.cpp
                 ../src/hgfBC.cpp
                 ../src/hgfIB.cpp
                 ../src/hgfPP.cpp)

SET(EXECUTABLE_SRCS ../src/driver.cpp)

SET(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)

FIND_PACKAGE(PARALUTION REQUIRED)
INCLUDE_DIRECTORIES(${PARALUTION_INCLUDE_DIR})

FIND_PACKAGE(OpenMP REQUIRED)
SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

FIND_PACKAGE(OpenCL)
IF(OpenCL_FOUND)
  INCLUDE_DIRECTORIES(${OpenCL_INCLUDE_DIRS})
ENDIF()

FIND_PACKAGE(CUDA)
IF(CUDA_FOUND)
  INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})
ENDIF()

ADD_LIBRARY(hgflib STATIC ${LIBRARY_SRCS})

IF(HGF_BUILD_PYTHON_HOOKS)
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})

  FIND_PACKAGE(PythonInterp REQUIRED)
  FIND_PACKAGE(PythonLibs REQUIRED)
  FIND_PACKAGE(NumPy REQUIRED)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIRS})

  SET(CMAKE_SWIG_FLAGS "")
  SET_SOURCE_FILES_PROPERTIES(hgf.i PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES(hgf.i PROPERTIES SWIG_FLAGS "-includeall")
  SET_SOURCE_FILES_PROPERTIES(hgf.i PROPERTIES SWIG_FLAGS "-cpperraswarn")
  SWIG_ADD_MODULE(hgf python hgf.i ${LIBRARY_SRCS})
  IF(CUDA_FOUND)
    IF(OPENCL_FOUND)
      SWIG_LINK_LIBRARIES(hgf ${PARALUTION_LIBRARY}
                              ${OpenCL_LIBRARIES}
                              ${CUDA_LIBRARIES}
                              ${CUDA_CUBLAS_LIBRARIES}
                              ${CUDA_cusparse_LIBRARY}
                              ${PYTHON_LIBRARIES}
                              hgflib)
    ELSE()
      SWIG_LINK_LIBRARIES(hgf ${PARALUTION_LIBRARY}
                              ${CUDA_LIBRARIES}
                              ${CUDA_CUBLAS_LIBRARIES}
                              ${CUDA_cusparse_LIBRARY}
                              ${PYTHON_LIBRARIES}
                              hgflib)
    ENDIF()
  ELSEIF(OPENCL_FOUND)
    SWIG_LINK_LIBRARIES(hgf ${PARALUTION_LIBRARY}
                            ${OpenCL_LIBRARIES}
                            ${PYTHON_LIBRARIES}
                            hgflib)
  ELSE()
    SWIG_LINK_LIBRARIES(hgf ${PARALUTION_LIBRARY}
                            ${PYTHON_LIBRARIES}
                            hgflib)
  ENDIF()
ENDIF()

IF(HGF_BUILD_EXECUTABLE)
  ADD_EXECUTABLE(hgf ${EXECUTABLE_SRCS})
  IF(CUDA_FOUND)
    IF(OPENCL_FOUND)
      TARGET_LINK_LIBRARIES(hgf ${OpenCL_LIBRARIES}
                                ${CUDA_LIBRARIES}
                                ${CUDA_CUBLAS_LIBRARIES}
                                ${CUDA_cusparse_LIBRARY}
                                hgflib
                                ${PARALUTION_LIBRARY})
    ELSE()
      TARGET_LINK_LIBRARIES(hgf ${CUDA_LIBRARIES}
                                ${CUDA_CUBLAS_LIBRARIES}
                                ${CUDA_cusparse_LIBRARY}
                                hgflib
                                ${PARALUTION_LIBRARY})
    ENDIF()
  ELSEIF(OPENCL_FOUND)
    TARGET_LINK_LIBRARIES(hgf ${OpenCL_LIBRARIES}
                              hgflib
                              ${PARALUTION_LIBRARY})
  ELSE()
    TARGET_LINK_LIBRARIES(hgf hgflib
                          ${PARALUTION_LIBRARY})
  ENDIF()
ENDIF()
