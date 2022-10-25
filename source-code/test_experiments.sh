#!/bin/bash

mkdir -p results

make clean && make

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Parameters:

# server_name: the server name

# mode = 1: tpc-h dataset
# mode = 2: twitter dataset

# scale:
# for tpc-h dataset, the raw data size = scale factor * 1 GB, e.g., scale factor = 0.05x means the raw data size is 50 MB;
# for twitter dataset, the scale factor is the number of total users in all three tables, e.g., scale factor = 20k means the number of total users is 20,000.

# outsourced_height = 100: no index caching optimization
# outsourced_height = 1: caching the index above the leaf level

# query_file: the sql query file

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Examples:

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load Data

scripts/run_oij_type1_load.sh 1 0.05x 1 |  tee results/tpch_0.05x_binary_oij_type1_load_1.txt
scripts/run_oij_type1_load.sh 1 0.05x 100 |  tee results/tpch_0.05x_binary_oij_type1_load_100.txt
scripts/run_oij_type1_multi_load.sh 1 0.05x 1 |  tee results/tpch_0.05x_multi_oij_type1_load_1.txt
scripts/run_oij_type1_multi_load.sh 1 0.05x 100 |  tee results/tpch_0.05x_multi_oij_type1_load_100.txt

scripts/run_one_type1_load.sh 1 0.05x 1 |  tee results/tpch_0.05x_binary_one_type1_load_1.txt
scripts/run_one_type1_load.sh 1 0.05x 100 |  tee results/tpch_0.05x_binary_one_type1_load_100.txt
scripts/run_one_type1_multi_load.sh 1 0.05x 1 |  tee results/tpch_0.05x_multi_one_type1_load_1.txt
scripts/run_one_type1_multi_load.sh 1 0.05x 100 |  tee results/tpch_0.05x_multi_one_type1_load_100.txt

scripts/run_oij_type1_load.sh 2 20k 1 |  tee results/twitter_20k_binary_oij_type1_load_1.txt
scripts/run_oij_type1_load.sh 2 20k 100 |  tee results/twitter_20k_binary_oij_type1_load_100.txt
scripts/run_oij_type1_multi_load.sh 2 20k 1 |  tee results/twitter_20k_multi_oij_type1_load_1.txt
scripts/run_oij_type1_multi_load.sh 2 20k 100 |  tee results/twitter_20k_multi_oij_type1_load_100.txt

scripts/run_one_type1_load.sh 2 20k 1 |  tee results/twitter_20k_binary_one_type1_load_1.txt
scripts/run_one_type1_load.sh 2 20k 100 |  tee results/twitter_20k_binary_one_type1_load_100.txt
scripts/run_one_type1_multi_load.sh 2 20k 1 |  tee results/twitter_20k_multi_one_type1_load_1.txt
scripts/run_one_type1_multi_load.sh 2 20k 100 |  tee results/twitter_20k_multi_one_type1_load_100.txt

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Binary Equi-Join

scripts/run_oij_type1_inlj_equi_join.sh localhost 1 0.05x 100 queries/tpch_binary_equi_q2.txt | tee results/tpch_0.05x_equi_oij_type1_bnlj_join_100_q2.txt
scripts/run_oij_type1_inlj_equi_join.sh localhost 1 0.05x 1 queries/tpch_binary_equi_q2.txt | tee results/tpch_0.05x_equi_oij_type1_bnlj_join_1_q2.txt
scripts/run_oij_type1_smj_equi_join.sh localhost 1 0.05x 1 queries/tpch_binary_equi_q2.txt | tee results/tpch_0.05x_equi_oij_type1_smj_join_1_q2.txt

scripts/run_one_type1_inlj_equi_join.sh localhost 1 0.05x 100 queries/tpch_binary_equi_q2.txt | tee results/tpch_0.05x_equi_one_type1_bnlj_join_100_q2.txt
scripts/run_one_type1_inlj_equi_join.sh localhost 1 0.05x 1 queries/tpch_binary_equi_q2.txt | tee results/tpch_0.05x_equi_one_type1_bnlj_join_1_q2.txt
scripts/run_one_type1_smj_equi_join.sh localhost 1 0.05x 1 queries/tpch_binary_equi_q2.txt | tee results/tpch_0.05x_equi_one_type1_smj_join_1_q2.txt

scripts/run_oij_type1_inlj_equi_join.sh localhost 2 20k 100 queries/twitter_binary_equi_q2.txt | tee results/twitter_20k_equi_oij_type1_bnlj_join_100_q2.txt
scripts/run_oij_type1_inlj_equi_join.sh localhost 2 20k 1 queries/twitter_binary_equi_q2.txt | tee results/twitter_20k_equi_oij_type1_bnlj_join_1_q2.txt
scripts/run_oij_type1_smj_equi_join.sh localhost 2 20k 1 queries/twitter_binary_equi_q2.txt | tee results/twitter_20k_equi_oij_type1_smj_join_1_q2.txt

scripts/run_one_type1_inlj_equi_join.sh localhost 2 20k 100 queries/twitter_binary_equi_q2.txt | tee results/twitter_20k_equi_one_type1_bnlj_join_100_q2.txt
scripts/run_one_type1_inlj_equi_join.sh localhost 2 20k 1 queries/twitter_binary_equi_q2.txt | tee results/twitter_20k_equi_one_type1_bnlj_join_1_q2.txt
scripts/run_one_type1_smj_equi_join.sh localhost 2 20k 1 queries/twitter_binary_equi_q2.txt | tee results/twitter_20k_equi_one_type1_smj_join_1_q2.txt

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Multiway Equi-Join

scripts/run_oij_type1_multi_join.sh localhost 1 0.05x 100 queries/tpch_obl_multi_q2.txt | tee results/tpch_0.05x_multi_oij_type1_join_100_q2.txt
scripts/run_oij_type1_multi_join.sh localhost 1 0.05x 1 queries/tpch_obl_multi_q2.txt | tee results/tpch_0.05x_multi_oij_type1_join_1_q2.txt

scripts/run_one_type1_multi_join.sh localhost 1 0.05x 100 queries/tpch_obl_multi_q2.txt | tee results/tpch_0.05x_multi_one_type1_join_100_q2.txt
scripts/run_one_type1_multi_join.sh localhost 1 0.05x 1 queries/tpch_obl_multi_q2.txt | tee results/tpch_0.05x_multi_one_type1_join_1_q2.txt

scripts/run_oij_type1_multi_join.sh localhost 2 20k 100 queries/twitter_obl_multi_q2.txt | tee results/twitter_20k_multi_oij_type1_join_100_q2.txt
scripts/run_oij_type1_multi_join.sh localhost 2 20k 1 queries/twitter_obl_multi_q2.txt | tee results/twitter_20k_multi_oij_type1_join_1_q2.txt

scripts/run_one_type1_multi_join.sh localhost 2 20k 100 queries/twitter_obl_multi_q2.txt | tee results/twitter_20k_multi_one_type1_join_100_q2.txt
scripts/run_one_type1_multi_join.sh localhost 2 20k 1 queries/twitter_obl_multi_q2.txt | tee results/twitter_20k_multi_one_type1_join_1_q2.txt

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Binary Band Join

scripts/run_oij_type1_inlj_band_join.sh localhost 1 0.05x 100 queries/tpch_band_q1.txt | tee results/tpch_0.05x_band_oij_type1_inlj_join_100_q1.txt
scripts/run_oij_type1_inlj_band_join.sh localhost 1 0.05x 1 queries/tpch_band_q1.txt | tee results/tpch_0.05x_band_oij_type1_inlj_join_1_q1.txt
scripts/run_oij_type1_smj_band_join.sh localhost 1 0.05x 1 queries/tpch_band_q1.txt | tee results/tpch_0.05x_band_oij_type1_smj_join_1_q1.txt

scripts/run_one_type1_inlj_band_join.sh localhost 1 0.05x 100 queries/tpch_band_q1.txt | tee results/tpch_0.05x_band_one_type1_inlj_join_100_q1.txt
scripts/run_one_type1_inlj_band_join.sh localhost 1 0.05x 1 queries/tpch_band_q1.txt | tee results/tpch_0.05x_band_one_type1_inlj_join_1_q1.txt
scripts/run_one_type1_smj_band_join.sh localhost 1 0.05x 1 queries/tpch_band_q1.txt | tee results/tpch_0.05x_band_one_type1_smj_join_1_q1.txt
