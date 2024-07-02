#!/bin/bash

all_ops=("FM0" "FM1" "FM2" "FM3" "FM4" "FM5" "FM7" "FS02" "FS1" "FT0" "FT1" "FT2" "FT5" "FT6" "FT7")
processes=("WmZ" "WpZ" "ZZ")

for process in "${processes[@]}"; do
    for op in "${all_ops[@]}"; do
        python run_rivet.py --evtMax 10000 --conf "MGPy8EG_aQGC${op}_QUAD_1_${process}_llqq" --DOCUT "YES" --redoRivet "yes" --redoPlots "yes" --keep "True" --name "All_op_plot"
    done
done