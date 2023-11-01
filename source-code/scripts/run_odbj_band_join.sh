#!/bin/bash

if [ $# -ne 5 ]; then
    echo "Usage: ./run_odbj_band_join server_name sort_type mode scale query_file"
    exit -1
fi

make

echo "===================================================================================================="

./odbj_band_join $2 $3 $4 $5
