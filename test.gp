set terminal epslatex color size 12cm,8cm lw 2
set encoding utf8
set key off

system  "rm test.eps"
set output 'test.tex'

intq ="intq.dat"
intq2="intq1.dat"

system "rm test.eps"

plot intq  u 1:2 with points,\
     intq2 u 1:2 with points
