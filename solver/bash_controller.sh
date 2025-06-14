#!/bin/bash
dw=0.04

for w in $(seq 0.0 0.05 0.75)
do
#echo $w &
matlab -batch "mainSolver($w)" &

done


