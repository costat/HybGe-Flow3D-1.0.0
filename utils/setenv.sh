#!/bin/bash

# bash script to set paths for boost and linear algebra libraries
# set export appropriate root directories before sourcing

## BOOST ##
if [ -z "${BOOST_ROOT+xxx}" ]; then
  echo BOOST_ROOT is not set
elif [ -z "${BOOST_ROOT}" ] && [ "${BOOST_ROOT+xxx}" = "xxx" ]; then
  echo BOOST_ROOT is set but empty
else
  export LIBRARY_PATH=$LIBRARY_PATH:$BOOST_ROOT/lib
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BOOST_ROOT/lib
  export CMAKE_LIBRARY_PATH=$CMAKE_LIBRARY_PATH:$BOOST_ROOT/lib

  export INCLUDE_PATH=$INCLUDE_PATH:$BOOST_ROOT/include
  export C_INCLUDE_PATH=$C_INCLUDE_PATH:$BOOST_ROOT/include
  export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$BOOST_ROOT/include
  export CMAKE_INCLUDE_PATH=$CMAKE_INCLUDE_PATH:$BOOST_ROOT/include

  echo BOOST paths set
fi

## PARALUTION, OPTIONAL ##
if [ -z "${PARALUTION_ROOT+xxx}" ]; then
  echo PARALUTION_ROOT is not set
elif [ -z "${PARALUTION_ROOT}" ] && [ "{PARALUTION_ROOT+xxx}" = "xxx" ]; then
  echo PARALUTION_ROOT is set but empty
else
  export LIBRARY_PATH=$LIBRARY_PATH:$PARALUTION_ROOT/lib
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PARALUTION_ROOT/lib
  export CMAKE_LIBRARY_PATH=$CMAKE_LIBRARY_PATH:$PARALUTION_ROOT/lib

  export INCLUDE_PATH=$INCLUDE_PATH:$PARALUTION_ROOT/inc
  export C_INCLUDE_PATH=$C_INCLUDE_PATH:$PARALUTION_ROOT/inc
  export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$PARALUTION_ROOT/inc
  export CMAKE_INCLUDE_PATH=$CMAKE_INCLUDE_PATH:$PARALUTION_ROOT/inc

  echo PARALUTION paths set
fi

## MKL, OPTIONAL ##
if [ -z "${MKL_ROOT+xxx}" ]; then
  echo MKL_ROOT is not set
elif [ -z "${MKL_ROOT}" ] && [ "{PARALUTION_ROOT+xxx}" = "xxx" ]; then
  echo MKL_ROOT is set but empty
else
  export LIBRARY_PATH=$LIBRARY_PATH:$MKL_ROOT/lib/intel64
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKL_ROOT/lib/intel64
  export CMAKE_LIBRARY_PATH=$CMAKE_LIBRARY_PATH:$MKL_ROOT/lib/intel64

  export INCLUDE_PATH=$INCLUDE_PATH:$MKL_ROOT/include
  export C_INCLUDE_PATH=$C_INCLUDE_PATH:$MKL_ROOT/include
  export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$MKL_ROOT/include
  export CMAKE_INCLUDE_PATH=$CMAKE_INCLUDE_PATH:$MKL_ROOT/include

  export MKLROOT=$MKL_ROOT
  echo MKL paths set
fi

