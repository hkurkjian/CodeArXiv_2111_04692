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

intq1  ="intqpasres_0_20.dat"
intq2  ="intqpasres_0_24.dat"
intq3  ="intqres_0_25.dat"
intq4  ="intqres_0_27.dat"
intq5  ="intqres_0_30.dat"
intq6  ="intqres_0_28.dat"

set xrange [0:3]
#set yrange [0:0.0001]

plot intq4 u 1:2,\

plot intq2 u 1:2,\
     intq3 u 1:2,\


unset xrange

plot dispLU1 u 2:3 with linespoints,\
     dispLU2 u 2:3 with linespoints,\
