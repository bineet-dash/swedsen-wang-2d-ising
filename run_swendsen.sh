#!/bin/sh
mkdir -p results
g++ -o sw_mc.out sw_mc.cpp
read -p "Enter lattice size: " size
read -p "Enter number of sweeps: " no_sweeps
./sw_mc.out $size $no_sweeps
gnuplot -e "latticeSize=$size" plot_results.gnuplot
rm sw_mc.out
