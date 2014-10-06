#!/bin/bash

mkdir -p ../output/$1_discard

time ../code/roundtrip -base "../output/$1_discard/$1" -top_directory "../../PDB_files/entire_PDB/files/$1" 

