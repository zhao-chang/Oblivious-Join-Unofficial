#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: ./run_odbj_load mode scale"
    exit -1
fi

make

echo "===================================================================================================="

./odbj_load $1 $2

echo "===================================================================================================="