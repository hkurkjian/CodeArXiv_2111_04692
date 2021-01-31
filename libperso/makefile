COMP=gfortran

%.o : %.f90
	$(COMP) -fopenmp -O3 -c $< -o $@

malib: nrtype.o nrutil.o modpol.o modsim.o recettes.o
	ar rcv libmalib.a $^ 

propre : 
	rm -f *.o 
	rm -f *.exe
	rm -f *.mod
	rm -f libmalib.a
