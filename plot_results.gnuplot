#!/usr/bin/gnuplot

set terminal pdfcairo

set title "Swendsen-Wang Simulation for  latticeSize x latticeSize  lattice"
set xl "Temperature"
set yl "<M>"
set output "results/magnetisation.pdf"
plot "results.dat" u 1:2 w l t "Magnetisation" lt rgb 'black'
reset

set title "Swendsen-Wang Simulation for  latticeSize x latticeSize  lattice"
set xl "Temperature"
set yl "sigma_M"
set output "results/error_magnetisation.pdf"
plot "results.dat" u 1:3 w l t "Error in Magnetisation" lt rgb 'black'
reset

set title "Swendsen-Wang Simulation for  latticeSize x latticeSize  lattice"
set xl "Temperature"
set yl "<E>"
set output "results/energy.pdf"
plot "results.dat" u 1:4 w l t "Internal Energy" lt rgb 'black'
reset

set title "Swendsen-Wang Simulation for  latticeSize x latticeSize  lattice"
set xl "Temperature"
set yl "sigma_E"
set output "results/error_energy.pdf"
plot "results.dat" u 1:5 w l t "Error in Energy" lt rgb 'black'
reset

set title "Swendsen-Wang Simulation for  latticeSize x latticeSize  lattice"
set xl "Temperature"
set yl "Xi"
set output "results/magnetic_susceptibility.pdf"
plot "results.dat" u 1:7 w l t "Magnetic Susceptibility" lt rgb 'black'
reset

set title "Swendsen-Wang Simulation for  latticeSize x latticeSize  lattice"
set xl "Temperature"
set yl "C_v"
set output "results/specific_heat.pdf"
plot "results.dat" u 1:6 w l t "Specific Heat" lt rgb 'black'
reset

set title "Swendsen-Wang Simulation for  latticeSize x latticeSize  lattice"
set xl "Temperature"
set yl "C_v"
set output "results/empty.pdf"
plot "results.dat" u 1:6 w l t "Specific Heat" lt rgb 'black'
reset
