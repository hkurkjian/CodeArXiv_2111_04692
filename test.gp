set terminal postscript color size 12cm,8cm lw 2
set encoding utf8
set key off

system  "rm test.eps"
set output 'test.eps'

dispLU="DONNEES/dispLU.dat"

#set xrange [2.053:2.060]
#set yrange [0:0.0001]

plot dispLU u 1:2 with linespoints
