#!/bin/bash

# bash script to set paths for installed hgf library
# set HGF_ROOT before sourcing

## HGF ##
if [ -z "${HGF_ROOT+xxx}" ]; then
  echo HGF_ROOT is not set
elif [ -z "${HGF_ROOT}" ] && [ "${HGF_ROOT+xxx}" = "xxx" ]; then
  echo HGF_ROOT is set but empty
else
  export LIBRARY_PATH=$LIBRARY_PATH:$HGF_ROOT/lib
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HGF_ROOT/lib
  export CMAKE_LIBRARY_PATH=$CMAKE_LIBRARY_PATH:$HGF_ROOT/lib

  export INCLUDE_PATH=$INCLUDE_PATH:$HGF_ROOT/include
  export C_INCLUDE_PATH=$C_INCLUDE_PATH:$HGF_ROOT/include
  export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$HGF_ROOT/include
  export CMAKE_INCLUDE_PATH=$CMAKE_INCLUDE_PATH:$HGF_ROOT/include

  echo HGF paths set
fi


