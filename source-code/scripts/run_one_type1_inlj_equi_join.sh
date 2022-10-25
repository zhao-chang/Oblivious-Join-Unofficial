#!/bin/bash

if [ $# -ne 5 ]; then
    echo "Usage: ./run_one_type1_inlj_equi_join server_name mode scale outsourced_height query_file"
    exit -1
fi

make

echo "===================================================================================================="

./one_type1_inlj_equi_join $2 $3 $4 $5
