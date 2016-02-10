CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7 FATAL_ERROR)

PROJECT(hgf)

SET(HGF_BUILD_PYTHON_HOOKS OFF CACHE BOOL "Build HGF Python Hooks using Swig")

macro (add_sources)
    file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND SRCS "${_relPath}/${_src}")
        else()
            list (APPEND SRCS "${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SRCS to parent directory
        set (SRCS ${SRCS} PARENT_SCOPE)
    endif()
endmacro()

add_sources(hgf.cpp)
add_sources(hgfMesh.cpp)
add_sources(hgfStokes.cpp)
add_sources(hgfArrays.cpp)
add_sources(hgfBC.cpp)
add_sources(hgfIB.cpp)
add_sources(hgfPP.cpp)

SET(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)

FIND_PACKAGE(PARALUTION REQUIRED)
INCLUDE_DIRECTORIES(${PARALUTION_INCLUDE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${PARALUTION_INCLUDE_DIR})

IF(HGF_BUILD_PYTHON_HOOKS)
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})

  FIND_PACKAGE(PythonLibs)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

  SET(CMAKE_SWIG_FLAGS "")

  SET_SOURCE_FILES_PROPERTIES(hgf.i PROPERTIES CPLUSPLUS ON)
  SWIG_ADD_MODULE(hgf python hgf.i ${SRCS})
  SWIG_LINK_LIBRARIES(hgf ${PYTHON_LIBRARIES})
ENDIF()

ADD_EXECUTABLE(hgf ${SRCS})
