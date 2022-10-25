#!/bin/bash

if [ $# -ne 4 ]; then
    echo "Usage: ./run_obli_database_equi_join server_name mode scale query_file"
    exit -1
fi

make

echo "===================================================================================================="

./obli_database_join $2 $3 $4
