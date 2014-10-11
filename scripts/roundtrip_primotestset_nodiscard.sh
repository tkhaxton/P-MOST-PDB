#!/bin/bash

mkdir -p ../output/primotestset_nodiscard

time ../code/roundtrip -base "../output/primotestset_nodiscard/primotestset" -top_directory "../../PDB_files/primotestset" -discard 0 -dbref 0

