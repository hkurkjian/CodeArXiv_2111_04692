This version of the code uses a library for modules containing 
basic numerical operations (nrtype, nrutil, recettes and modsim)
Library modules are in the "libperso" directory. 
The library file "libmalib.a" is produced by the command

> make malib

(executed in the libperso subdirectory)

This compiles all the module files (gfortran -O3 -c module.f90 -o module.o)
and compresses them in a library file 
(ar rcv libmalib.a nrtype.o nrutil.o modsim.o recettes.o).

Then the programs test and intldc are compiled using

> make test
> make intldc

(in the main directory). The -I flag tells the compiler where the 
modules (except vars2.mod and dspec2.mod) are when producing the .o files 
(for example gfortran -I./modules -O3 -c dspec2.f90 -o dspec2.o).
Finally the main program is compiled using dspec2.o and the L flag
to include the library (gfortran -O3 dspec2.o test.o -L./modules -lmalib -o test)
