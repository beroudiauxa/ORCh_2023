#File to print ...

set encoding iso_8859_1
set terminal postscript portrait color noenhanced "Times-Roman" 12

set style line 2000 linetype 2 linecolor rgb "#000000" linewidth 1.5 pointtype 5 pointsize 0.5
set style line 1000 linetype 2 linecolor rgb "#00c000" linewidth 1.5 pointtype 5 pointsize 0.5


set macros
TMARGIN = "set tmargin at screen 0.90; set bmargin at screen 0.60"
LMARGIN = "set lmargin at screen 0.25; set rmargin at screen 0.80"

set key at graph 0.5,-0.2
set key font "0,8"
set key spacing 0.8

set output "fitnessEvolution.eps"

set xlabel "generation" 

@TMARGIN
@LMARGIN

set ylabel "fitness"
plot "fitness.dat" using 1:2 with linespoints ls 1000 title 'mean fitness' ,\
     "fitness.dat" using 1:3 with linespoints ls 2000 title 'best fitness'
 

