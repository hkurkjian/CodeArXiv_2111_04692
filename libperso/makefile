COMP=ifort
BALISES= -qopenmp -O3

COMP=gfortran
BALISES= -fopenmp -O3

%.o : %.f90
	$(COMP) $(BALISES) -c $< -o $@

malib: nrtype.o nrutil.o modpol.o modsim.o recettes.o
	ar rcv libmalib.a $^ 

propre : 
	rm -f *.o 
	rm -f *.exe
	rm -f *.mod
	rm -f libmalib.a
