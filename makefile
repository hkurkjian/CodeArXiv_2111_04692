COMP=gfortran
MOD=./libperso
LIB_PERSO=-L./libperso -lmalib

%.eps : %.gp
	gnuplot $<

%.o : %.f90
	$(COMP) -I$(MOD) -O3 -c $< -o $@

encorr: dspec.o intldc.o encorr.o
	$(COMP) -fcheck=bounds -O3 $^ $(LIB_PERSO) -o encorr

test: dspec.o Zerom.o intldc.o intpole.o test.o
	$(COMP) -O3 $^ $(LIB_PERSO) -o test

pointsM: dspec.o pointsM.o
	$(COMP) -O3 $^ $(LIB_PERSO) -o pointsM

spectre: dspec.o Zerom.o spectre.o
	$(COMP) -O3 $^ $(LIB_PERSO) -o spectre

propre : 
	rm -f *.o 
	rm -f *.mod
	rm -f intldc
	rm -f test
	rm -f encorr
	rm -f *.aux
	rm -f *.syntec.gz
	rm -f spectre
	rm -f pointsM
