#!/bin/bash

mkdir -p ../output/PISCES_training_$1_nodiscard

time ../code/roundtrip -base "../output/PISCES_training_$1_nodiscard/$1" -top_directory "../../PDB_files/entire_PDB/files/$1" -discard 0 -training_top_directory "../../PDB_files/entire_PDB/files" -list_file "../../PDB_files/PISCES/140401/cullpdb_pc20_res1.6_R0.25_d140401_chains2733"

