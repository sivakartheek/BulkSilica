set terminal png enhanced
set terminal postscript enhanced
set terminal postscript eps enhanced
set term postscript enhanced color
set output "colorindex.ps"
plot "../data/densityandenergy.dat" using 1:2 title "" with line  lw 4 lc 165
set ylabel "Density (g/cm^3)" 
set xlabel "Number of Iterations"
set grid
set title "Density"
set terminal png font arial 23 size 1024,768
set output "../graph/Density versus Iterations.png"
replot
