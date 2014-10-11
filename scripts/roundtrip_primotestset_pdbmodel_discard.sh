#!/bin/bash

mkdir -p ../output/primotestset_pdbmodel_discard

time ../code/roundtrip -base "../output/primotestset_pdbmodel_discard/primotestset" -top_directory "../../PDB_files/primotestset" -input_base "../output/entire_discard/entire" -discard 1 -dbref 0

