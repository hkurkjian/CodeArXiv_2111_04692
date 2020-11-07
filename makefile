COMP=gfortran
MOD=./libperso
LIB_PERSO=-L./libperso -lmalib

%.o : %.f90
	$(COMP) -I$(MOD) -O3 -c $< -o $@

encorr: dspec.o intldc.o encorr.o
	$(COMP) -fcheck=bounds -O3 $^ $(LIB_PERSO) -o encorr

test: dspec.o test.o
	$(COMP) -O3 $^ $(LIB_PERSO) -o test

test2: dspec2.o test2.o
	$(COMP) -O3 $^ $(LIB_PERSO) -o test2

propre : 
	rm -f *.o 
	rm -f *.mod
	rm -f intldc
	rm -f test
	rm -f encorr
	rm -f *.aux
	rm -f *.syntec.gz
