set terminal postscript color size 12cm,8cm lw 2
set encoding utf8
set key on

system  "rm test.eps"
set output 'test.eps'

dispLU="DONNEES/dispLU.dat"
dispLU ="DONNEES/selfEtot100.dat"
dispLU2="DONNEES/selfEtot100bis.dat"

dispLU1="DONNEES/selfELU1.dat"
dispLU2="DONNEES/selfELU2.dat"

intq1  ="intq100_k10/intqpasres_0_7.dat"
intq2  ="intq100_k10/intqpasres_0_6.dat"
intq3  ="intq100_k10/intqres_0_8.dat"
intq4  ="intq100_k10/intqres_0_9.dat"
intq5  ="intq100_k10/intqres_0_12.dat"

intq1  ="intqpasres_0_0.dat"
intq2  ="intqres_0_0.dat"

intom1  ="intpasresom.dat"
intom2  ="intomEPS35.dat"
intom3  ="intom.dat"

#set xrange [0:0.1]
#set xrange [9:10]
#set yrange [0:0.0001]

#set logscale x
set xrange [250:400]

plot intom1 u 2:3,\
     intom2 u 2:3,\
     intom3 u 2:3,\

set xrange [0:25]

plot intq1 u 1:2,\
     intq2 u 1:2,\

#     intq3 u 1:2,\
#     intq4 u 1:2,\

plot intq4 u 1:2,\

unset xrange

plot dispLU1 u 2:3 with linespoints,\
     dispLU2 u 2:3 with linespoints,\
