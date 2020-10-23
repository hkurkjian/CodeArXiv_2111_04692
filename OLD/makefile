COMP=gfortran
LIB_TH=-pthread -L/opt/intel/mkl/10.0.1.014/lib/em64t -lguide -lmkl_blacs_lp64 -lmkl_cdft -lmkl_cdft_core -lmkl_core -lmkl_em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_lapack -lmkl_solver_lp64

%.o : %.f90
	$(COMP) -O3 -c $< -o $@

example: nrtype.o nrutil.o modsim.o recettes.o intsuromega.o example.f90
	$(COMP) -O3 $^ -o example

propre : 
	rm -f *.o 
	rm -f *.exe
	rm -f *.mod
	rm -f example
