This repository contains the program used in the numerical calculation of

[1] S. Van Loon, J. Tempere and H. Kurkjian, Quasiparticle disintegration in superfluid Fermi gases, SciPost ?? (arXiv:2111.04692)

The main program is written in fortran 90. Clean executables and .o files using

> make propre

***********************************LIBRARY

The subdirectory "libperso" contains a library of basic operations
(such as integration or interpolation routines). It should be be built by the command

> make malib

(executed in the libperso subdirectory)
 
***********************************EXECUTABLES

Build an executable <programme> using 

> make <programme>

(in the main directory)

** spectre.f90 USES eqdetat, dspec, Zerom

Calls Zerom repeatedly to store the dispersion of the Anderson-Bogoliubov branch q |-> omega_q in a datafile. See the comments in the file for content of spectre.inp

** mkInfofile.f90 USES eqdetat

Completes the info file produced by spectre (which should be specified in mkInfofile.inp)

** suppointsM.f90 USES dspec

Produce a table of value of M for varying q and om. This is used by encorr, selfcons through the routine intldc/intres. 
Two data files should be created: one on an interval [0,q1] and one on [q1,qmax] with qmax>>1 

** encorr.f90 USES selftot

Repeatedly calls detG or detGres to produce a grid of values of Sigma for varying k,zk. See the comments in the file for content of encorr.inp

** selfcons.f90 USES selftot

Find the real root of detG in the off-resonant case. See the comments in the file for content of selfcons.inp. If detGres is not used
the production of the datafile by suppointsM can be avoided by providing "bidon" as input fichiers(1)

***********************************MODULES

** eqdetat.f90

Contains subroutines and functions that evaluate the
BCS equation-of-state (for instance the correspondance between 1/kF*a and mu/Delta).

** dspec.f90

Evaluates the inverse pair-propagator \Gamma^{-1}=M (see equation (12) in [1]).
Contains implicit input variables described in the MODULE vars.
Front-end routines: *mat_pairfield         (computes M in the general case)
                    *mat_pairfield_pttq    (computes M at low q using analytic formulas)
                    *mat_pairfield_gom0    (computes M at large omega)
                    *mat_pairfield_gom0_gq (computes M at large omega and large q)
                    *oangpp                (computes the angular points of the pair-breaking continuum and writes them in opp)

** Zerom.f90 USES dspec

Computes the energy of the Anderson-Bogoliubov mode for a given value of q.
                    *Zero (Finds the real root (\omega_q) of det M below the pair-breaking continuum)
    
** angularint.f90

Analytic formulas for the angular integration in the self-energy (see appendix B in [1])

** estM.f90 USES dspec

Estimates \Gamma by interpolating over tabulated data. Calls mat_pairfield from dspec if unsuccessful.
The data file (produced by suppointsM) should loaded by load_data before estmat_pairfield can be used.
                   *estmat_pairfield      (interpolates M from the data contained in donnees and donneesgq, loaded by load_data(fich))

** intldc.f90   USES dspec estM angularint. Data files produced by suppointsM and recursively by intldc/intpasres.

Evaluates the branch contribution to the self-energy (eq (20) of [1]). Either by
                  *intpasres              (Evaluates Sigma^bc in the off-resonant case (no imaginary part). Reads (lecture=T) or write (ecriture=T) its own datafile. Initialisation by ini_intpasres is needed)
                  *intres                 (Evaluates Sigma^bc in the resonant case. Uses estM to evaluate M if interpolation=T)
For intres a careful (too careful?) splitting of the integration over q,om is done as explained in app A. of [1]
                  *bornesk,lignesenergie  (Compute the inner boundaries of the 1->3 disintegration continuum see Fig. 13 in [1])
                  *bornesq,bornesom       (Compute the value of q and om where the resonance configuration changes, see Tables I-II-III and IV in [1])

** intpole.f90  USES dspec angularint. Data file produced by spectre. Info file produced by spectre, completed by mkInfofile

Evaluates the pole contribution to the self-energy (eq (18) of [1]) using
                  *selfEpole              (Evaluates Sigma^p. rdInfo(fich) should be used to load the data file produced by spectre)

** selftot USES dspec, estM, intpol, intldc

Front-end module aimed at rationalize the calls of the modules.
                 *initialisation          (Initialize the variables of dspec, loads the data for estM, intldc and intpol)
                 *detG                    (Computes the determinant of the Green's function, and the self-energy matrix (sigcomb) with both the pole and branch cut contributions. Uses intpasres for Sigma^bc)
                 *detGres                 (Computes the determinant of the Green's function, and the self-energy matrix (sigcomb) with both the pole and branch cut contributions. Uses intres    for Sigma^bc)
                 *thresholds              (Computes the boundaries of the 1->2 and 1->3 disintegration continua)
