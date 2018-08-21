#!/bin/bash

for l in {0..10}
do
	echo "$l    $k"
	echo $(date) > log_$1_$l.out
        ./comb_fit input/z1_k$1_l$l.dat >> log_$1_$l.out
	echo $(date) >> log_$1_$l.out
done
