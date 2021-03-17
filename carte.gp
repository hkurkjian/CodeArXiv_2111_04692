set terminal postscript color size 12cm,8cm lw 2
set encoding utf8
set key off

system  "rm carte.eps"
set output 'carte.eps'

carte  ="DONNEES/cartebove3.dat"
carte2 ="DONNEES/detbelow3.dat"
carte3 ="DONNEES/dettgk2.dat"
cartegk  ="DONNEES/detgk4.dat"
cartetgk ="DONNEES/dettgk4.dat"
cartettgk ="DONNEES/detttgk4.dat"
tab=      "Tab.dat"

x0=4

epsBCS(x)=sqrt((x**2-x0)**2+1)
set pm3d interpolate 2,2

#set cbrange [0:0.01]
set xrange [2.8:6.0]
set yrange [3.0:25.0]

set pm3d map

splot cartegk    u 1:2:(log(1/($3**2+$4**2))) with pm3d,\
      cartetgk   u 1:2:(log(1/($3**2+$4**2))) with pm3d,\
      cartettgk  u 1:2:(log(1/($3**2+$4**2))) with pm3d,\
      carte2     u 1:2:(log(1/($3**2+$4**2))) with pm3d,\
      tab       u 1:(epsBCS($1)):(0) with line lc rgb "green",\

#carte3 u 1:2:(1/($3**2+$4**2)) with pm3d,\
#      carte  u 1:2:(1/($3**2+$4**2)) with pm3d,\
#      carte2 u 1:2:(1/($3**2+$4**2)) with pm3d,\
