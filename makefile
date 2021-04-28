COMP=gfortran
COMP=ifort
MOD=./libperso
LIB_PERSO=-L./libperso -lmalib
BALISES=-fopenmp -fcheck=all -fbacktrace -O3

%.eps : %.gp
	gnuplot $<

%.o : %.f90
	$(COMP) -I$(MOD) $(BALISES) -c $< -o $@

selfcons: eqdetat.o dspec.o estM.o angularint.o intldc.o intpole.o selftot.o selfcons.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o selfcons

encorr: eqdetat.o dspec.o estM.o angularint.o intldc.o intpole.o selftot.o encorr.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o encorr

hktest: eqdetat.o dspec.o Zerom.o bestM.o angularint.o intldc.o intpole.o selftot.o hktest.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o hktest

testdspec: eqdetat.o dspec.o testdspec.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o testdspec

test: eqdetat.o dspec.o estM.o test.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o test

SVLtest: eqdetat.o dspec.o Zerom.o angularint.o intpole.o SVLtest.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o SVLtest

pointsM: dspec.o pointsM.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o pointsM

suppointsM: eqdetat.o dspec.o suppointsM.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o suppointsM

spectre: eqdetat.o dspec.o Zerom.o spectre.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o spectre

mkInfoFile: eqdetat.o mkInfoFile.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o mkInfoFile

traitement: eqdetat.o dspec.o angularint.o intpole.o traitement.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o traitement

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
