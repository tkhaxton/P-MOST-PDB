#!/bin/bash

discard=1

directory="../molecular_files/$1"
filebase="$1"
pdbdir=$(echo $1 | head -c 3 | tail -c 2)
inputpdb="../../PDB_files/entire_PDB/files/$pdbdir/pdb$1.ent"
inputbase="../output/entire_nodiscard/entire"

mkdir -p $directory

time ../code/create_molecular_files -directory $directory -filebase $filebase -pdbfile $inputpdb -input_base $inputbase

