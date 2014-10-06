#!/bin/bash

mkdir -p ../output/$1_nodiscard

time ../code/roundtrip -base "../output/$1_nodiscard/$1" -top_directory "../../PDB_files/entire_PDB/files/$1" -discard 0 

