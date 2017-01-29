set terminal png enhanced
set terminal postscript enhanced
set terminal postscript eps enhanced
set term postscript enhanced color
set output "colorindex.ps"
set style line 1 lt 1 lw 2
set style line 2 lt 2 lw 2
set style line 3 lt 3 lw 2
plot 	"../data/RDF.dat" using 1:4 title "" with line ls 4
set autoscale
set xrange [0:6]
set xlabel "r (Å)" 
set ylabel "RDF G_{O-O}(r)"
set title "Oxygen - Oxygen"
set terminal png font arial 22 size 1024,768
set output "../graph/Partial RDF O-O.png"
replot

