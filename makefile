COMP=gfortran
MOD=./libperso
LIB_PERSO=-L./libperso -lmalib

%.eps : %.gp
	gnuplot $<

%.o : %.f90
	$(COMP) -I$(MOD) -fopenmp -fcheck=all -fbacktrace -fstack-arrays -fmax-array-constructor=300000000 -O3 -c $< -o $@

encorr: dspec.o intldc.o encorr.o
	$(COMP) -fcheck=bounds -O3 $^ $(LIB_PERSO) -o encorr

HKtest: eqdetat.o dspec.o Zerom.o estM.o bestM.o intldc.o intpole.o HKtest.o
	$(COMP) -O3 -fcheck=all -fbacktrace -fstack-arrays -fmax-array-constructor=300000000 $^ $(LIB_PERSO) -o HKtest

SVLtest: eqdetat.o dspec.o Zerom.o bestM.o intldc.o intpole.o SVLtest.o
	$(COMP) -O3 $^ $(LIB_PERSO) -o SVLtest

pointsM: dspec.o pointsM.o
	$(COMP) -O3 $^ $(LIB_PERSO) -o pointsM

pointsMparra: dspec.o pointsMparra.o
	$(COMP) -O3 -fopenmp $^ $(LIB_PERSO) -o pointsMparra

suppointsM: dspec.o suppointsM.o
	$(COMP) -O3 -fopenmp $^ $(LIB_PERSO) -o suppointsM

spectre: eqdetat.o dspec.o Zerom.o spectre.o
	$(COMP) -O3 $^ $(LIB_PERSO) -o spectre

mkInfoFile: eqdetat.o mkInfoFile.o
	$(COMP) -O3 $^ $(LIB_PERSO) -o mkInfoFile

propre : 
	rm -f *.o 
	rm -f *.o95 
	rm -f *.mod
	rm -f intldc
	rm -f HKtest
	rm -f SVLtest
	rm -f encorr
	rm -f *.aux
	rm -f *.syntec.gz
	rm -f spectre
	rm -f pointsM
	rm -f pointsMparra
	rm -f suppointsM
	rm -f *.exe
