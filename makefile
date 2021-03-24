COMP=gfortran
MOD=./libperso
LIB_PERSO=-L./libperso -lmalib

%.eps : %.gp
	gnuplot $<

%.o : %.f90
	$(COMP) -I$(MOD) -fopenmp -fcheck=all -fbacktrace -O3 -c $< -o $@

selfcons: eqdetat.o dspec.o bestM.o angularint.o intldc.o intpole.o selftot.o selfcons.o
	$(COMP) -fopenmp -fcheck=all -fbacktrace -O3 $^ $(LIB_PERSO) -o selfcons

encorr: eqdetat.o dspec.o bestM.o angularint.o intldc.o intpole.o selftot.o encorr.o
	$(COMP) -fopenmp -fcheck=all -fbacktrace -O3 $^ $(LIB_PERSO) -o encorr

hktest: eqdetat.o dspec.o Zerom.o bestM.o angularint.o intldc.o intpole.o selftot.o hktest.o
	$(COMP) -O3 -fopenmp -fcheck=all -fbacktrace $^ $(LIB_PERSO) -o hktest

hktest2: eqdetat.o dspec.o Zerom.o bestM.o angularint.o intldc.o intpole.o hktest2.o
	$(COMP) -O3 -fcheck=all -fbacktrace $^ $(LIB_PERSO) -o hktest2

SVLtest: eqdetat.o dspec.o Zerom.o angularint.o intpole.o SVLtest.o
	$(COMP) -fopenmp -fcheck=all -fbacktrace -O3 $^ $(LIB_PERSO) -o SVLtest

pointsM: dspec.o pointsM.o
	$(COMP) -O3 -fopenmp $^ $(LIB_PERSO) -o pointsM

suppointsM: eqdetat.o dspec.o suppointsM.o
	$(COMP) -O3 -fopenmp $^ $(LIB_PERSO) -o suppointsM

spectre: eqdetat.o dspec.o Zerom.o spectre.o
	$(COMP) -O3 $^ $(LIB_PERSO) -o spectre

mkInfoFile: eqdetat.o mkInfoFile.o
	$(COMP) -O3 $^ $(LIB_PERSO) -o mkInfoFile

traitement: traitement.o
	$(COMP) -O3 $^ $(LIB_PERSO) -o traitement

propre : 
	rm -f *.o 
	rm -f *.o95 
	rm -f *.mod
	rm -f intldc
	rm -f HKtest
	rm -f SVLtest
	rm -f encorr
	rm -f mkInfoFile
	rm -f *.aux
	rm -f *.syntec.gz
	rm -f spectre
	rm -f pointsM
	rm -f pointsMparra
	rm -f suppointsM
	rm -f *.exe
