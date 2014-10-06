#!/bin/bash

mkdir -p ../output/entire_nodiscard

time ../code/roundtrip -base "../output/entire_nodiscard/entire" -top_directory "../../PDB_files/entire_PDB/files" -discard 0 

