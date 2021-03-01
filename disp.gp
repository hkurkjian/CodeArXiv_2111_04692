set terminal postscript color size 12cm,8cm lw 2
set encoding utf8
set key off

system  "rm disp.eps"
set output 'disp.eps'

disp ="DONNEES/solautocortest.dat"

x0=4
eBCS(x)=sqrt((x**2-x0)**2+1)

set xrange [2.86:2.9]
set yrange [2.99:3.0]

plot eBCS(x) with line,\
     disp index 0 u 1:2 with points,\
     disp index 1 u 1:2 with points,\
     disp index 2 u 1:2 with points,\
     disp index 3 u 1:2 with points,\
     disp index 4 u 1:2 with points,\
     disp u 1:3 with line,\
     disp u 1:4 with line
