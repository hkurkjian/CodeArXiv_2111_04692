set terminal postscript color size 12cm,8cm lw 2
set encoding utf8
set key off

system  "rm carte.eps"
set output 'carte.eps'

carte  ="DONNEES/cartebove3.dat"
carte2 ="DONNEES/detbelow3.dat"

x0=4

set cbrange [0:5]
#set xrange [2.86:2.9]
#set yrange [2.99:3.0]

set pm3d map

splot carte  u 1:2:(1/($3**2+$4**2)) with pm3d,\
      carte2 u 1:2:(1/($3**2+$4**2)) with pm3d

