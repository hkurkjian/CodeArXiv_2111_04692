MOD=-I./libperso
LIB_PERSO=-L./libperso -lmalib

COMP=ifort
BALISES=-qopenmp -check all -backtrace -O3

COMP=gfortran
BALISES=-fopenmp -fcheck=all -fbacktrace -O3

%.eps : %.gp
	gnuplot $<

%.o : %.f90
	$(COMP) $(MOD) $(BALISES) -c $< -o $@

selfcons: eqdetat.o dspec.o estM.o angularint.o intldc.o intpole.o selftot.o selfcons.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o selfcons

encorr: eqdetat.o dspec.o estM.o angularint.o intldc.o intpole.o selftot.o encorr.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o encorr

suppointsM: eqdetat.o dspec.o suppointsM.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o suppointsM

spectre: eqdetat.o dspec.o Zerom.o spectre.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o spectre

mkInfoFile: eqdetat.o mkInfoFile.o
	$(COMP) $(BALISES) $^ $(LIB_PERSO) -o mkInfoFile


propre : 
	rm -f *.o 
	rm -f *.mod
	rm -f encorr
	rm -f mkInfoFile
	rm -f spectre
	rm -f suppointsM
