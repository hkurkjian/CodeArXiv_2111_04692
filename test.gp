set terminal postscript color size 12cm,8cm lw 2
set encoding utf8
set key off

system  "rm test.eps"
set output 'test.eps'

intq ="intqgrdk1.dat"
intq2="intqgrdk2.dat"

intq ="intqavecbestMres.dat"
intq2="intqexact2.dat"
#intq3="intqtest4.dat"

set xrange [0:5]
#set yrange [0:0.0001]

plot intq2   u 1:($2) with points,\
     intq  u 1:($2) with points,\
#     intq3  u 1:($3) with points,\
#     intq   u 1:($3) with points,\
