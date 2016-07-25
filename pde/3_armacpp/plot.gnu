set terminal postscript eps color enhanced "Helvetica" 20
#set terminal pdfcairo color enhanced #font "Helvetica, 20"
set output 'sol.eps'
#set encoding iso_8859_1
#save_encoding = GPVAL_ENCODING
set encoding utf8

## Styles for lines
set style line 1  lt 1 lw 1.5 lc rgb "#000000" pt 5  ps 0.4 # black
set style line 2  lt 1 lw 1.5 lc rgb "#0000FF" pt 9  ps 0.5 # blue
set style line 3  lt 2 lw 1.5 lc rgb "#FF0000" pt 7  ps 0.4 # red
set style line 4  lt 2 lw 1.5 lc rgb "#33A02C" pt 13 ps 0.6 # dark green
set style line 5  lt 3 lw 1.5 lc rgb "#5d478b" pt 15 ps 0.5
set style line 6  lt 3 lw 1.5 lc rgb "#FF7F00" pt 12 ps 0.5 # dark orange
set style line 7  lt 4 lw 1.5 lc rgb "#a0522d" pt 5  ps 0.4
set style line 8  lt 4 lw 1.5 lc rgb "#E7298A" pt 2  ps 0.4 # dark magenta
set style line 9  lt 1 lw 1.5 lc rgb "#AAAA00" pt 3  ps 0.5
set style line 10 lt 4 lw 1.5 lc rgb "#666666" pt 7  ps 0.4 # dark gray
set style line 11 lt 4 lw 1.5 lc rgb "#66A61E" pt 7  ps 0.5 # dark lime green 
set style line 15 lt 4 lw 1.5 lc rgb "#999999" pt 7  ps 0.5 # light gray
#set style line 8  lt 4 lw 1.5 lc rgb "#FB9A99" pt 7  ps 1.0 # light red
set style line 99 lt 1 lw 1.5 lc rgb "#000000"
set style data linespoints

## Styles for points
# set style line 1 pt 1  lc rgb "#000000" ps 0.7
# set style line 2 pt 2  lc rgb "#0000FF" ps 0.7
# set style line 3 pt 3  lc rgb "#FF0000" ps 0.7
# set style line 4 pt 5  lc rgb "#6e9b3d" ps 0.7
# set style line 5 pt 7  lc rgb "#5d478b" ps 0.7
# set style line 6 pt 9  lc rgb "#daa520" ps 0.7
# set style line 7 pt 11 lc rgb "#a0522d" ps 0.7
# set style line 8 pt 13 lc rgb "#999999" ps 0.7
# set style line 8 pt 15 lc rgb "#AAAA00" ps 0.7
# set style line 9 pt 66 lc rgb "#000000" ps 0.7
# set style line 10 pt 68 lc rgb "#888888" ps 0.7
#

set xtics pi/2.
#set format x '%.0P π'
set format x ''
#set xtics (0, 'π' pi, '2π' 2*pi)

p 'data' u 1:2 w p t "Numerical",\
      '' u 1:3 w l t "Analytic"
