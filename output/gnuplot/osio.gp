set terminal png enhanced
set terminal postscript enhanced
set terminal postscript eps enhanced
set term postscript enhanced color
set output "colorindex.ps"
set key box
set key left
set style line 1 lt 1 lw 2
set style line 2 lt 2 lw 2
set style line 3 lt 1 lw 2
set style line 4 lt 2 lw 2
plot 	"../data/BAD.dat" using 1:3 title "" with line ls 1
set autoscale
set xrange [60:200]
set title "O-Si-O angle"
set xlabel "Bond Angle (Degrees)" 
set ylabel "Angular Distribution"
set terminal png font arial 22 size 1024,768
set output "../graph/osio.png"
replot



