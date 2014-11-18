#!/bin/bash

mkdir -p ../output/PISCES_training_entire_nodiscard

time ../code/roundtrip -base "../output/PISCES_training_entire_nodiscard/entire" -top_directory "../../PDB_files/entire_PDB/files" -discard 0 -training_top_directory "../../PDB_files/entire_PDB/files" -list_file "../../PDB_files/PISCES/140401/cullpdb_pc20_res1.6_R0.25_d140401_chains2733"

