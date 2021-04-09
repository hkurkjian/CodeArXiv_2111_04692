set terminal postscript color size 12cm,8cm lw 2
set encoding utf8
set key bottom right

system  "rm carte.eps"
set output 'carte.eps'

carte  ="DONNEES/cartebove3.dat"
carte2 ="DONNEES/detbelow3.dat"
carte3 ="DONNEES/dettgk2.dat"
cartegk  ="DONNEES/detgk4.dat"
cartepk  ="DONNEES/detpk4.dat"
cartetgk ="DONNEES/dettgk4.dat"
cartettgk ="DONNEES/detttgk4.dat"
cartet3gk ="DONNEES/dett3gk4.dat"
coupegk ="DONNEES/coupegk.dat"
coupepk ="DONNEES/coupepk.dat"
coupefr ="DONNEES/coupefr.dat"
tab=      "Tab.dat"

x0=4

epsBCS(x)=sqrt((x**2-x0)**2+1)

#set cbrange [0:0.01]

set xlabel "omega/Delta"
set ylabel "Im(G_{22}/det G(k,omega+i0+))"

plot coupepk index 0  u 2:((1/($3**2+$4**2))) with linespoints lc rgb "red"    title "k=",\
     coupepk index 3  u 2:((1/($3**2+$4**2))) with linespoints lc rgb "orange" title "k=",\
#     coupepk index 1  u 2:((1/($3**2+$4**2))) with linespoints lc rgb "blue"   title "k=",\
#     coupepk index 2  u 2:((1/($3**2+$4**2))) with linespoints lc rgb "green"  title "k=",\
#     coupepk index 4  u 2:((1/($3**2+$4**2))) with lines lc rgb "brown"  title "k=",\
#     coupepk index 5  u 2:((1/($3**2+$4**2))) with lines lc rgb "black"  title "k=",\
#     coupepk index 5  u 2:((1/($3**2+$4**2))) with lines lc rgb "black"  title "k=",\
#     coupepk index 6  u 2:((1/($3**2+$4**2))) with lines lc rgb "black"  title "k=",\
#     coupepk index 7  u 2:((1/($3**2+$4**2))) with lines lc rgb "black"  title "k=",\
#     coupepk index 8  u 2:((1/($3**2+$4**2))) with lines lc rgb "black"  title "k=",\
#     coupepk index 9  u 2:((1/($3**2+$4**2))) with lines lc rgb "black"  title "k=",\
#     coupepk index 10 u 2:((1/($3**2+$4**2))) with lines lc rgb "black"  title "k=",\

plot coupefr index 0 u 2:6 with lines lc rgb "red"    title "k=6.04",\
     coupefr index 1 u 2:6 with lines lc rgb "blue"   title "k=6.54",\
     coupefr index 2 u 2:6 with lines lc rgb "green"  title "k=7.48",\
     coupefr index 3 u 2:6 with lines lc rgb "orange" title "k=7.97",\
     coupefr index 4 u 2:6 with lines lc rgb "brown"  title "k=8.47",\
     coupefr index 5 u 2:6 with lines lc rgb "black"  title "k=8.66",\

set ylabel "abs(1/det G(k,omega+i0+))"
set key top right

plot coupegk index 0 u 2:((1/($3**2+$4**2))) with lines lc rgb "red"    title "k=6.04",\
     coupegk index 1 u 2:((1/($3**2+$4**2))) with lines lc rgb "blue"   title "k=6.54",\
     coupegk index 2 u 2:((1/($3**2+$4**2))) with lines lc rgb "green"  title "k=7.48",\
     coupegk index 3 u 2:((1/($3**2+$4**2))) with lines lc rgb "orange" title "k=7.97",\
     coupegk index 4 u 2:((1/($3**2+$4**2))) with lines lc rgb "brown"  title "k=8.47",\
     coupegk index 5 u 2:((1/($3**2+$4**2))) with lines lc rgb "black"  title "k=8.66",\

set pm3d map
set xrange [2.0:4.0]
#set yrange [1.0:9.0]

set pm3d interpolate 2,2
#splot cartet3gk    u 1:2:((1/($3**2+$4**2))) with pm3d,\
#
set key off
set ylabel "omega/Delta"
set xlabel "k/k_0"

#splot cartegk    u 1:2:((1/($3**2+$4**2))) with pm3d,\
#      cartepk    u 1:2:((1/($3**2+$4**2))) with pm3d,\
#      cartetgk   u 1:2:((1/($3**2+$4**2))) with pm3d,\
#      cartettgk  u 1:2:((1/($3**2+$4**2))) with pm3d,\
#      cartet3gk  u 1:2:((1/($3**2+$4**2))) with pm3d,\
#      carte2     u 1:2:((1/($3**2+$4**2))) with pm3d,\
#      tab       u 1:(epsBCS($1)):(0) with line lc rgb "green",\

splot cartegk    u 1:2:((1/($3**2+$4**2))) with pm3d,\
      cartetgk   u 1:2:((1/($3**2+$4**2))) with pm3d,\
      cartepk    u 1:2:((1/($3**2+$4**2))) with pm3d,\
      cartettgk  u 1:2:((1/($3**2+$4**2))) with pm3d,\
      cartet3gk  u 1:2:((1/($3**2+$4**2))) with pm3d,\
      carte2     u 1:2:((1/($3**2+$4**2))) with pm3d,\
      tab        u 1:(epsBCS($1)):(0) with line lc rgb "green",\

#carte3 u 1:2:(1/($3**2+$4**2)) with pm3d,\
#      carte  u 1:2:(1/($3**2+$4**2)) with pm3d,\
#      carte2 u 1:2:(1/($3**2+$4**2)) with pm3d,\
