set terminal postscript color size 12cm,8cm lw 2
set encoding utf8
set key off

system  "rm disp.eps"
set output 'disp.eps'

disp  ="DONNEES/solautocortest.dat"
disp2 ="DONNEES/selfEtotabove4.dat"
disp2 ="DONNEES/detabove3.dat"
disp2 ="DONNEES/coupe.dat"

x0=4
eBCS(x)=sqrt((x**2-x0)**2+1)

#set xrange [2.86:2.9]
#set yrange [2.99:3.0]

#detGres=(-zk+xik-sig(1))*(-zk-xik-sig(2))-(1.0_qp-sig(3))**2

#plot eBCS(x) with line,\
#     disp index 0 u 1:2 with points,\
#     disp index 1 u 1:2 with points,\
#     disp index 2 u 1:2 with points,\
#     disp index 3 u 1:2 with points,\
#     disp index 4 u 1:2 with points,\
#     disp u 1:3 with line,\
#     disp u 1:4 with line

plot disp2 u 2:(1/($3**2+$4**2)) with linespoints,\

#plot disp2 index 00 u 2:(1/($3**2+$4**2)) with linespoints,\
#     disp2 index 10 u 2:(1/($3**2+$4**2)) with linespoints,\
#     disp2 index 20 u 2:(1/($3**2+$4**2)) with linespoints,\
#     disp2 index 30 u 2:(1/($3**2+$4**2)) with linespoints,\
#     disp2 index 40 u 2:(1/($3**2+$4**2)) with linespoints,\
#     disp2 index 43 u 2:(1/($3**2+$4**2)) with linespoints,\
#     disp2 index 46 u 2:(1/($3**2+$4**2)) with linespoints,\
#     disp2 index 49 u 2:(1/($3**2+$4**2)) with linespoints,\
#     disp2 index 80 u 2:(1/($3**2+$4**2)) with linespoints,\
