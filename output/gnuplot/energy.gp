set terminal png enhanced
set terminal postscript enhanced
set terminal postscript eps enhanced
set term postscript enhanced color
set output "colorindex.ps"
plot "../data/densityandenergy.dat" using 1:3 title "" with line  lw 4 lc 1
set ylabel "Energy (eV)" 
set xlabel "Number of Iterations"
set grid
set title "Potential Energy per SiO_{2} atom"
set terminal png font arial 23 size 1024,768
set output "../graph/Potential Energy versus Iterations.png"
replot

