#!/bin/bash

mkdir -p ../output/primotestset_pdbmodel_nodiscard

time ../code/roundtrip -base "../output/primotestset_pdbmodel_nodiscard/primotestset" -top_directory "../../PDB_files/primotestset" -input_base "../output/entire_nodiscard/entire" -discard 0 -dbref 0

