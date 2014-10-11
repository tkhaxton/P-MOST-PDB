#!/bin/bash

mkdir -p ../output/primotestset_discard

time ../code/roundtrip -base "../output/primotestset_discard/primotestset" -top_directory "../../PDB_files/primotestset" -discard 1 -dbref 0

