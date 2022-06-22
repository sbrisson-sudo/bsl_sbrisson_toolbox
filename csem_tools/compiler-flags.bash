#! /bin/bash

if [ $# -eq 0 ]
then
    BUILD_DIR=build 
else
    BUILD_DIR=$1
fi

if [ ! -e "${BUILD_DIR}/CMakeCache.txt" ]
then 
    echo "Error : ${BUILD_DIR}, not a cmake build directory."
fi 

cat "${BUILD_DIR}/CMakeCache.txt" | grep -E "CMAKE_Fortran_FLAGS:STRING=|CMAKE_C_FLAGS:STRING="