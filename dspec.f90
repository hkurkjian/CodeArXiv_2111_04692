!Bibliography
![1] Scientific Reports 10, 11591
![2] https://doi.org/10.5802/crphys.1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE vars
 USE nrtype
 IMPLICIT NONE
 !$OMP THREADPRIVATE(xq,opp)
 !!$OMP THREADPRIVATE(x0,xq,beta,xqjoin,opp,opt,ptbranchmtpp,ptbranchmtpt)
! Global variables

! Input variables
 REAL(QP) xq,x0,beta !momentum q/sqrt(2m*Delta), interaction parameter mu/Delta, temperature parameter Delta/kT
 REAL(QP) EPSpp,EPSpt,EPSrpp,EPSrpt !precisions used respectively for the real (EPSpp for quasiparticle-quasiparticle part,EPSpt for the quasiparticle-quasihole part) and imaginary parts (EPSrpp,EPSrpt) of the matrix elements
!EPSrpp|EPSrpt should be much lower than EPSpp|EPSpt
!Suggested configuration EPSpp=1.0e-6_qp, EPSrpp=1.0e-10_qp, EPSpt=1.0e-6_qp, EPSrpt=1.e-7_qp
 LOGICAL bla1,bla2   !verbose level
 LOGICAL temperaturenulle ! .TRUE. to compute the matrix at T=0 (beta=+oo), using analytical formulas to speed-up the calculation. .FALSE. otherwise 

! Internal variables (calculated by oangpp and oangpt subroutines)
 INTEGER ptbranchmtpp,ptbranchmtpt !Number of angular points on the qp-qp and qp-qh branch cuts, written by oangpp and oangpt
 REAL(QP) opp(1:4)   !angular points om1,om2,om3 of the qp-qp branch cut and k such that om2=E(k+q/2)+E(k-q/2), written by oangpp
 REAL(QP) opt(1:2)   !angular points om_ph of the qp-qh branch cut and k such that om_ph=E(k+q/2)-E(k-q/2), written by oangpt
 REAL(QP) xqjoin ! value of xq after which opp(2) and opp(3) coincide (for x0>0)
 REAL(QP) x0crit     !In the far BCS regime (x0>x0crit) use an optimized routine for the omega integral

! Parameters
 REAL(QP), PARAMETER, DIMENSION(1:6) :: eta=(/1.0,1.0, 1.0,1.0,-1.0,-1.0/) !Matrice eta   (B3) dans [1]: eta=signe devant Sigma dans (36)
 REAL(QP), PARAMETER, DIMENSION(1:6) :: sig=(/1.0,1.0,-1.0,1.0,-1.0, 1.0/) !Matrice sigma (B5) dans [1]: sig=+1 pour S^epsilon et -1 pour S^omega
END MODULE vars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE dspec
 USE nrtype
 USE nrutil, ONLY : ecrit
 USE vars
 USE modsim
 IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE mat(om0,e,det,M,fr)
!mat computes the 3*3 fluctuation matrix in the phase-modulus-density basis

!input as global variables (using the "vars" module):  
!xq,x0,beta,temperaturenulle,EPSpp,EPSpt,EPSrpp,EPSrpt,bla1

!input as arguments:
!complex energy z=om0+I*e WARNING: current version of the program limited to the real axis: use e=0.0 

!output: 
!M(1:3,1:3) a 3*3 symmetric matrix in the basis 1:phase, 2:modulus, 3:density
!det=M(1,1)*M(2,2)-M(1,2)**2, the determinant of the 2*2 order-parameter matrix
!fr a 3*3 symmetric matrix of response functions in the basis 1:phase, 2:modulus, 3:density (see equation (43) in [1])

IMPLICIT NONE
COMPLEX(QPC), DIMENSION(1:3,1:3), INTENT(OUT) :: M,fr
COMPLEX(QPC), INTENT(OUT) :: det
REAL(QP), INTENT(IN)  :: om0,e

! Compute the angular points for the present value of xq and x0
call oangpp
call oangpt

M(:,:)=cmplx(0.0_qp,0.0_qp,kind=qpc)
M(1,1)=elementmat(1.5_qp,om0,e)
M(2,2)=elementmat(2.5_qp,om0,e)
M(1,2)=elementmat(3.5_qp,om0,e)
M(3,3)=elementmat(4.5_qp,om0,e)
M(1,3)=elementmat(5.5_qp,om0,e)
M(2,3)=elementmat(6.5_qp,om0,e)
M(2,1)=M(1,2)
M(3,1)=M(1,3)
M(3,2)=M(2,3)

det=M(1,1)*M(2,2)-M(1,2)**2.0_qp

fr(1,1)=M(2,2)/det
fr(2,2)=M(1,1)/det
fr(1,2)=-M(1,2)/det
fr(1,3)=(M(1,2)*M(2,3)-M(1,3)*M(2,2))/det
fr(2,3)=(M(1,2)*M(1,3)-M(2,3)*M(1,1))/det
fr(3,3)=M(3,3)-(M(1,3)**2.0_qp*M(2,2)+M(2,3)**2.0_qp*M(1,1)-2.0_qp*M(1,3)*M(2,3)*M(1,2))/det
fr(2,1)=fr(1,2)
fr(3,1)=fr(1,3)
fr(3,2)=fr(2,3)

END SUBROUTINE mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE mat_pairfield(om0,e,det,M,fr)
!mat_pairfield computes the 2*2 fluctuation matrix in the phase-modulus basis
!input: same for mat

!output: 
!M(1:2,1:2) a 2*2 symmetric matrix in the basis 1:phase, 2:modulus
!det=M(1,1)*M(2,2)-M(1,2)**2, the determinant of the matrix
!fr=1/M the pair propagator

IMPLICIT NONE
COMPLEX(QPC), DIMENSION(1:2,1:2), INTENT(OUT) :: M,fr
COMPLEX(QPC), INTENT(OUT) :: det
REAL(QP), INTENT(IN)  :: om0,e

call oangpp
call oangpt

M(:,:)=cmplx(0.0_qp,0.0_qp,kind=qpc)
M(1,1)=elementmat(1.5_qp,om0,e)
M(2,2)=elementmat(2.5_qp,om0,e)
M(1,2)=elementmat(3.5_qp,om0,e)
M(2,1)=M(1,2)

det=M(1,1)*M(2,2)-M(1,2)**2.0_qp

fr(1,1)=M(2,2)/det
fr(2,2)=M(1,1)/det
fr(1,2)=-M(1,2)/det
fr(2,1)=fr(1,2)

END SUBROUTINE mat_pairfield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE der_mat(om0,e,M,dM)

IMPLICIT NONE
REAL(QP), DIMENSION(1:2,1:2), INTENT(OUT) :: M,dM
REAL(QP), INTENT(IN)  :: om0,e

call oangpp
call oangpt

M(:,:)=0.0_qp
M(1,1)=real(elementmat(1.5_qp,om0,e))
M(2,2)=real(elementmat(2.5_qp,om0,e))
M(1,2)=real(elementmat(3.5_qp,om0,e))
M(2,1)=M(1,2)

dM(:,:)=0.0_qp
dM(1,1)=deriveepp(1.5_qp,om0,e,-1.0_qp)
dM(2,2)=deriveepp(2.5_qp,om0,e,-1.0_qp)
dM(1,2)=deriveepp(3.5_qp,om0,e,-1.0_qp)

dM(2,1)=dM(1,2)

END SUBROUTINE der_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION elementmat(num,om0,e)
!Computes one of the 3*3 matrix element
!input:
!complex energy z=om0+I*e  WARNING: current version of the program limited to the real axis: use e=0.0  
!num: a REAL selector 
! 1.5 -> phase-phase
! 2.5 -> modulus-modulus
! 3.5 -> modulus-phase
! 4.5 -> density-density
! 5.5 -> density-phase
! 6.5 -> density-modulus
! input as global variables same as for mat
!input:
IMPLICIT NONE
REAL(QP), INTENT(IN)  :: num,om0,e
COMPLEX(QPC) elementmat
REAL(QP) re,im 

if(temperaturenulle)then
 re=eta(floor(num))*intpp(num,om0,e,-1.0_qp)
 im=-PI*eta(floor(num))*rhopp(num,om0,-1.0_qp)
else
! attention: intpp(T>0) est l'intégrale de (e_+ ou om)*ab(n+ + n-)/(om^2-e+^2), d'où le signe moins
! intpt est la partie réelle de -S_ab -> signe + 
 re=eta(floor(num))*(intpp(num,om0,e,-1.0_qp)-intpp(num,om0,e,+1.0_qp))+intpt(num,om0,e) 
 im=-PI*(eta(floor(num))*(rhopp(num,om0,-1.0_qp)-rhopp(num,om0,+1.0_qp))+rhopt(num,om0)) !même remarque que pour la partie réelle
endif
elementmat=cmplx(re,im,kind=qpc)

END FUNCTION elementmat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION intpp(num,om0,e,T)
!Computes the real part of the quasip-quasip term of a matrix element as an integral over energies I(om0)=int_0^+oo d(om) rhopp(om)/(om0-om)
!When rhopp(om0)\=0 the integrand has a singularity in om0. The principal part is calculated by subtracting rhopp(om0)/(om0-om) to the integrand.

!input:
!complex energy z=om0+I*e  WARNING: current version of the program limited to the real axis: use e=0.0  
!num: a REAL selector (see elementmat)
!T: a REAL selector
!T=-1. intpp returns the "+1" part of the matrix element: sum_k (e_+ or om)*ab/(om^2-e+^2) (see eqs. (36-37) in [1])
!T=+1. intpp returns the "n_+ + n_-" part of the matrix element: sum_k (e_+ or om)*ab*(n+ + n-)/(om^2-e+^2)

!input as global variables 
!xq,x0,beta,temperaturenulle,bla1
!ptbranchmtpp must have been defined by oangpp

USE recettes
IMPLICIT NONE
REAL(QP), INTENT(IN)  :: num,om0,e,T
REAL(QP) intpp
REAL(QP) r0 !spectral density in om0
REAL(QP) bmax,db,omcrit!bmax: energy cutoff, omcrit: critical value of om after which racinfq should be used to integrate the tail
REAL(QP) Icomp,Iinf!partial integrals
REAL(QP) arg(1:1,1:10),bornes(1:10)
REAL(QP), PARAMETER :: singu=1.0_qp,nonsingu=-1.0_qp !singu/nonsingu: call inter with or without the regularizing term r0/(om0-om)
REAL(QP), PARAMETER :: su=1.0_qp,nsu=-1.0_qp !short for singu/nonsingu
REAL(QP) ag(1:4),mil
INTEGER choix(1:10)
LOGICAL axer !calculation on or outside the real axis (currently not active)

axer=.TRUE.

Iinf=0.0_qp
Icomp=0.0_qp

if(bla1)then
 write(6,*)
 write(6,*)"intpp avec num, om0, T=",floor(num),om0,T
 write(6,*)
endif

omcrit=100.0_qp !at om0>omcrit switch from midpntq to racinfq to integration the ~1/om^(3/2) tail between om(3) and 2.0_qp*om0-opp(3)

if(axer)then
 if(T<0.0_qp)then ! Intégrale du "+1" de 1+nP+nM

  if(floor(num)<4)then !energy cutoff in the integrals
   db=100.0_qp
   bmax=max(2.0e6_qp,2*om0-opp(3)+db,2*opp(3)+db,4*x0+db,2*om0-2*x0+db,2*om0+db) !In case 2*om0-opp or 2*opp overflows bmax. This formula works whatever the value of ptbranchmtpp
  else
   bmax=1.0e55_qp
  endif

  Iinf=intinfini(bmax) !analytic integration of the 1/om^(3/2) tail from bmax to +oo for M11, M22 and M12
  call ecrit(bla1,'Iinf=',Iinf)

  r0=rhopp(num,om0,T)
  if(x0<x0crit)then
    if(ptbranchmtpp==1)then
     if(om0<opp(1))then !cut in opp(1)
      arg(1,1:2)=         (/nsu,    nsu/)
      choix(1:2)=         (/mpnt,   rinf/)
      intpp=decoupe(inter,(/0.0_qp, opp(1),bmax/),arg(1:1,1:2),choix(1:2),EPSpp,bla1)
     else!cut in opp(1), om0 and 2.0_qp*om0-opp(1)
      arg(1,1:4)=         (/nsu,    su,     su,    nsu/)
      choix(1:4)=         (/mpnt,   mpnt,   mpnt,  rinf/)
      if(om0<omcrit) choix(2:3)=rinf
      intpp=decoupe(inter,(/0.0_qp, opp(1), om0,   2*om0-opp(1),bmax/),arg(1:1,1:4),choix(1:4),EPSpp,bla1)
     endif


    elseif(ptbranchmtpp==2)then
     if(om0<opp(1))then!cut in opp(1), opp(2) and 2*opp(2)
      arg(1,1:4)=         (/nsu,    nsu,    nsu,    nsu/)
      choix(1:4)=         (/mpnt,   mpnt,   msql,   rinf/)
      intpp=decoupe(inter,(/0.0_qp, opp(1), opp(2), 2*opp(2),bmax/),arg(1:1,1:4),choix(1:4),EPSpp,bla1)
     elseif(om0<opp(2))then!cut in opp(1), om0, opp(2) and 2*opp(2)
      arg(1,1:5)=         (/nsu,    su,     su,     nsu,    nsu/)
      choix(1:5)=         (/mpnt,   mpnt,   mpnt,   msql,   rinf/)
      intpp=decoupe(inter,(/0.0_qp, opp(1), om0,    opp(2), 2*opp(2),bmax/),arg(1:1,1:5),choix(1:5),EPSpp,bla1)
  
      Icomp=-r0*log((opp(2)-om0)/(om0-opp(1)))!remove int_opp(1)^opp(2) r0/(om0-om)
      call ecrit(bla1,'Icomp=',Icomp)
     else!cut in opp(1), opp(2), om0 and 2.0_qp*om0-opp(2)
      arg(1,1:5)=         (/nsu,    nsu,    su,    su,    nsu/)
      choix(1:5)=         (/mpnt,   mpnt,   msql,  mpnt,  rinf/)
      if(om0<omcrit)choix(3:4)=rinf
      intpp=decoupe(inter,(/0.0_qp, opp(1), opp(2),om0,   2*om0-opp(2),bmax/),arg(1:1,1:5),choix(1:5),EPSpp,bla1)
     endif


    elseif(ptbranchmtpp==3)then
     if(om0<opp(1))then !cut in opp(1), opp(2), opp(3)
      arg(1,1:4)=         (/nsu,    nsu,    nsu,    nsu/)
      choix(1:4)=         (/mpnt,   mpnt,   msql,   rinf/)
      intpp=decoupe(inter,(/0.0_qp, opp(1), opp(2), opp(3),bmax/),arg(1:1,1:4),choix(1:4),EPSpp,bla1)
     elseif(om0<opp(2))then!cut in opp(1), om0, opp(2), opp(3)
      arg(1,1:5)=         (/nsu,    su,     su,     nsu,    nsu/)
      choix(1:5)=         (/mpnt,   mpnt,   mpnt,   msql,   rinf/)
      intpp=decoupe(inter,(/0.0_qp, opp(1), om0,    opp(2), opp(3),bmax/),arg(1:1,1:5),choix(1:5),EPSpp,bla1)
  
      Icomp=-r0*log((opp(2)-om0)/(om0-opp(1)))!remove int_opp(1)^opp(2) r0/(om0-om)
      call ecrit(bla1,'Icomp=',Icomp)
     elseif(om0<opp(3))then!cut in opp(1), opp(2), om0 and opp(3)
      arg(1,1:5)=         (/nsu,    nsu,    su,     su,     nsu/)
      choix(1:5)=         (/mpnt,   mpnt,   msql,   mpnt,   rinf/)
      intpp=decoupe(inter,(/0.0_qp, opp(1), opp(2), om0,    opp(3),bmax/),arg(1:1,1:5),choix(1:5),EPSpp,bla1)
  
      Icomp=-r0*log((opp(3)-om0)/(om0-opp(2)))!remove int_opp(2)^opp(3) r0/(om0-om)
      call ecrit(bla1,'Icomp=',Icomp)
     else!cut in opp(1), opp(2), opp(3), om0 and 2.0_qp*om0-opp(3)
      arg(1,1:6)=         (/nsu,    nsu,    nsu,   su,    su,    nsu/)
      choix(1:6)=         (/mpnt,   mpnt,   msql,  mpnt,  mpnt,  rinf/)
      if(om0<omcrit) choix(4:5)=rinf
      intpp=decoupe(inter,(/0.0_qp, opp(1), opp(2),opp(3),om0,   2*om0-opp(3),bmax/),arg(1:1,1:6),choix(1:6),EPSpp,bla1)
     endif
    endif


  else!x0>x0crit


    if(ptbranchmtpp==1)then
     ag=(/opp(1),2*x0,-bidon,-bidon/)
     call tri(ag)
     if(om0<ag(1))then
      arg(1,1:5)=         (/nsu,    nsu,   nsu,                  nsu,  nsu/)
      choix(1:5)=         (/msqu,   msql,  msqu,                 msql, rinf/)
      intpp=decoupe(inter,(/0.0_qp, ag(1), (ag(1)+ag(2))/2.0_qp, ag(2),2*ag(2),bmax/),arg(1:1,1:5),choix(1:5),EPSpp,bla1)
     elseif(om0<ag(2))then
      arg(1,1:5)=         (/nsu,    su,    su,    nsu,  nsu/)
      choix(1:5)=         (/msqu,   msql,  msqu,  msql, rinf/)
      intpp=decoupe(inter,(/0.0_qp, ag(1), (ag(1)+ag(2))/2.0_qp,   ag(2),2*ag(2),bmax/),arg(1:1,1:5),choix(1:5),EPSpp,bla1)

      Icomp=-r0*log((ag(2)-om0)/(om0-ag(1)))
     else
      arg(1,1:6)=         (/nsu,    nsu,   nsu,                 su,   su,  nsu/)
      choix(1:6)=         (/msqu,   msql,  msqu,                msql, rinf,rinf/)
      intpp=decoupe(inter,(/0.0_qp, ag(1), (ag(1)+ag(2))/2.0_qp,ag(2),om0, 2*om0-ag(2),bmax/),arg(1:1,1:6),choix(1:6),EPSpp,bla1)
     endif
     call ecrit(bla1,'Icomp=',Icomp)


    elseif(ptbranchmtpp==2)then
     ag=(/opp(1),opp(2),2*x0,-bidon/)
     call tri(ag)
     if(om0<ag(1))then
      arg(1,1:6)=         (/nsu,    nsu,  nsu,  nsu,                 nsu,  nsu/)
      choix(1:6)=         (/mpnt,   msqu, msql, msqu,                msql, rinf/)
      intpp=decoupe(inter,(/0.0_qp, ag(1),ag(2),(ag(2)+ag(3))/2.0_qp,ag(3),2*ag(3),bmax/),arg(1:1,1:6),choix(1:6),EPSpp,bla1)
     elseif(om0<ag(2))then
       arg(1,1:7)=         (/nsu,   su,   su,  nsu,  nsu,                 nsu,  nsu/)
       choix(1:7)=         (/mpnt,  mpnt, msqu,msql, msqu,                msql, rinf/)
       bornes(1:8)=        (/0.0_qp,ag(1),(ag(1)+ag(2))/2.1_qp,ag(2),(ag(2)+ag(3))/2.0_qp,ag(3),2*ag(3),bmax/)
       bornes(1:8)=        (/0.0_qp,ag(1),om0                 ,ag(2),(ag(2)+ag(3))/2.0_qp,ag(3),2*ag(3),bmax/)
       intpp=decoupe(inter,bornes(1:8),arg(1:1,1:7),choix(1:7),EPSpp,bla1)

       Icomp=-r0*log((ag(2)-om0)/(om0-ag(1)))
     elseif(om0<ag(3))then
      arg(1,1:6)=         (/nsu,   nsu,   su,     su,   nsu,   nsu/)
      choix(1:6)=         (/mpnt,  msqu,  msql,   msqu, msql,  rinf/)
      intpp=decoupe(inter,(/0.0_qp,ag(1), ag(2),  (ag(2)+ag(3))/2.0_qp,  ag(3), 2*ag(3),bmax/),arg(1:1,1:6),choix(1:6),EPSpp,bla1)

      Icomp=-r0*log((ag(3)-om0)/(om0-ag(2)))
     else
      arg(1,1:7)=         (/nsu,   nsu,  nsu,  nsu,                 su,   su,  nsu/)
      choix(1:7)=         (/mpnt,  msqu, msql, msqu,                msql, rinf,rinf/)
      intpp=decoupe(inter,(/0.0_qp,ag(1),ag(2),(ag(2)+ag(3))/2.0_qp,ag(3),om0, 2*om0,bmax/),arg(1:1,1:7),choix(1:7),EPSpp,bla1)

      Icomp=-r0*log((2*om0-om0)/(om0-ag(3)))
     endif
     call ecrit(bla1,'Icomp=',Icomp)


    elseif(ptbranchmtpp==3)then
     ag=(/opp(1),opp(2),opp(3),2*x0/)
     call tri(ag)
     if(om0<ag(1))then
      arg(1,1:9)=   (/nsu,    nsu,  nsu,                 nsu,  nsu,                 nsu,  nsu,                 nsu,  nsu/)
      choix(1:9)=   (/mpnt,   mpnt, msqu,                msql, msqu,                msql, msqu,                msql, rinf/)
      if(floor(num)==5) choix(2)=msql
      bornes(1:10)= (/0.0_qp, ag(1),(ag(1)+ag(2))/2.0_qp,ag(2),(ag(2)+ag(3))/2.0_qp,ag(3),(ag(3)+ag(4))/2.0_qp,ag(4),2*ag(4),bmax/)
      intpp=decoupe(inter,bornes(1:10),arg(1:1,1:9),choix(1:9),EPSpp,bla1)
     elseif(om0<ag(2))then
      arg(1,1:9)  =       (/nsu,   su,   su,  nsu,  nsu,                 nsu,  nsu,                 nsu,  nsu/)
      choix(1:9)  =       (/mpnt,  mpnt, msqu,msql, msqu,                msql, msqu,                msql, rinf/)
      if(floor(num)==5) choix(3)=rinf
      bornes(1:10)=       (/0.0_qp,ag(1),(ag(1)+ag(2))/2,ag(2),(ag(2)+ag(3))/2.0_qp,ag(3),(ag(3)+ag(4))/2.0_qp,ag(4),2*ag(4),bmax/)
      intpp=decoupe(inter,bornes(1:10),arg(1:1,1:9),choix(1:9),EPSpp,bla1)
  
      Icomp=-r0*log((ag(2)-om0)/(om0-ag(1)))
     elseif(om0<ag(3))then
      mil=max(om0,(ag(2)+ag(3))/2.0_qp)
      arg(1,1:8)=         (/nsu,    nsu,     su,     su,   nsu,   nsu,                 nsu,   nsu/)
      choix(1:8)=         (/mpnt,   msqu,    msql,   msqu, msql,  msqu,                msql,  rinf/)
      bornes(1:9)=        (/0.0_qp, ag(1),   ag(2),  mil,  ag(3), (ag(3)+ag(4))/2.0_qp,ag(4), 2*ag(4),bmax/)
      intpp=decoupe(inter,bornes(1:9),arg(1:1,1:8),choix(1:8),EPSpp,bla1)

      Icomp=-r0*log((ag(3)-om0)/(om0-ag(2)))
     elseif(om0<ag(4))then
      arg(1,1:8)=         (/nsu,    nsu,    nsu,   nsu,                 su,   su,   nsu,   nsu/)
      choix(1:8)=         (/mpnt,   msqu,   msql,  msqu,                msql, msqu, msql,  rinf/)
      bornes(1:9)=        (/0.0_qp, ag(1),  ag(2), (ag(2)+ag(3))/2.0_qp,ag(3),om0,  ag(4), 2*ag(4),bmax/)
      intpp=decoupe(inter,bornes(1:9),arg(1:1,1:8),choix(1:8),EPSpp,bla1)

      Icomp=-r0*log((ag(4)-om0)/(om0-ag(3)))
     else
      arg(1,1:9)=         (/nsu,    nsu,    nsu,   nsu,                 nsu,  nsu,                 su,   su,  nsu/)
      choix(1:9)=         (/mpnt,   msqu,   msql,  msqu,                msql, msqu,                msql, rinf,rinf/)
      bornes(1:10)=       (/0.0_qp, ag(1),  ag(2), (ag(2)+ag(3))/2.0_qp,ag(3),(ag(3)+ag(4))/2.0_qp,ag(4),om0, 2*om0,bmax/)
      intpp=decoupe(inter,bornes(1:10),arg(1:1,1:9),choix(1:9),EPSpp,bla1)
      
      Icomp=-r0*log((2*om0-om0)/(om0-ag(4)))
     endif
     call ecrit(bla1,'Icomp=',Icomp)
    endif
  endif

  intpp=intpp+Icomp+Iinf
  call ecrit(bla1,'intpp=',intpp)
 elseif(T>0.0_qp)then ! Intégrale du terme en nP+nM de 1+nP+nM
  bmax=1111840.0_qp
  r0=rhopp(num,om0,T)
  if(ptbranchmtpp==1)then
   if(om0<opp(1))then
    arg(1,1:2)=         (/nsu,    nsu/)
    choix(1:2)=         (/mpnt,   minf/)
    intpp=decoupe(inter,(/0.0_qp, opp(1),bmax/),arg(1:1,1:2),choix(1:2),EPSpp,bla1)
   else
    arg(1,1:4)=         (/nsu,    su,     su,    nsu/)
    choix(1:4)=         (/mpnt,   mpnt,   mpnt,  minf/)
    intpp=decoupe(inter,(/0.0_qp, opp(1), om0,   2.0_qp*om0-opp(1),bmax/),arg(1:1,1:4),choix(1:4),EPSpp,bla1)
   endif
  elseif(ptbranchmtpp==2)then
   if(om0<opp(1))then
    arg(1,1:3)=         (/nsu,    nsu,    nsu/)
    choix(1:3)=         (/mpnt,   mpnt,   minf/)
    intpp=decoupe(inter,(/0.0_qp, opp(1), opp(2),bmax/),arg(1:1,1:3),choix(1:3),EPSpp,bla1)
   elseif(om0<opp(2))then
    arg(1,1:5)=         (/nsu,    su,     su,     nsu,    nsu/)
    choix(1:5)=         (/mpnt,   mpnt,   mpnt,   msql,   minf/)
    intpp=decoupe(inter,(/0.0_qp, opp(1), om0,    opp(2), 2*opp(2),bmax/),arg(1:1,1:5),choix(1:5),EPSpp,bla1)

    Icomp=-r0*log((opp(2)-om0)/(om0-opp(1)))
    call ecrit(bla1,'Icomp=',Icomp)
   else
    arg(1,1:5)=         (/nsu,    nsu,    su,     su,    nsu/)
    choix(1:5)=         (/mpnt,   mpnt,   msql,   mpnt,  minf/)
    intpp=decoupe(inter,(/0.0_qp, opp(1), opp(2), om0,   2*om0-opp(2),bmax/),arg(1:1,1:5),choix(1:5),EPSpp,bla1)
   endif
  elseif(ptbranchmtpp==3)then
   if(om0<opp(1))then
    arg(1,1:4)=         (/nsu,    nsu,    nsu,    nsu/)
    choix(1:4)=         (/mpnt,   mpnt,   msql,   minf/)
    intpp=decoupe(inter,(/0.0_qp, opp(1), opp(2), opp(3), bmax/),     arg(1:1,1:4),choix(1:4),EPSpp,bla1)
   elseif(om0<opp(2))then
    arg(1,1:5)=         (/nsu,    su,     su,   nsu,    nsu/)
    choix(1:5)=         (/mpnt,   mpnt,   mpnt, msql,   minf/)
    intpp=decoupe(inter,(/0.0_qp, opp(1), om0,  opp(2), opp(3),bmax/),arg(1:1,1:5),choix(1:5),EPSpp,bla1)

    Icomp=-r0*log((opp(2)-om0)/(om0-opp(1)))
    call ecrit(bla1,'Icomp=',Icomp)
   elseif(om0<opp(3))then
    arg(1,1:5)=         (/nsu,    nsu,    su,     su,   nsu/)
    choix(1:5)=         (/mpnt,   mpnt,   msql,   mpnt, minf/)
    intpp=decoupe(inter,(/0.0_qp, opp(1), opp(2), om0,  opp(3),bmax/),arg(1:1,1:5),choix(1:5),EPSpp,bla1)

    Icomp=-r0*log((opp(3)-om0)/(om0-opp(2)))
    call ecrit(bla1,'Icomp=',Icomp)
   else
    arg(1,1:6)=         (/nsu,    nsu,    nsu,    su,    su,   nsu/)
    choix(1:6)=         (/mpnt,   mpnt,   msql,   mpnt,  mpnt, minf/)
    intpp=decoupe(inter,(/0.0_qp, opp(1), opp(2), opp(3),om0,  2*om0-opp(3),bmax/),arg(1:1,1:6),choix(1:6),EPSpp,bla1)
   endif
  else
   STOP "Erreur de ptbranchmntpp"
  endif
 intpp=intpp+Icomp
 call ecrit(bla1,'intpp=',intpp)
 endif
endif

CONTAINS
  FUNCTION intinfini(bm) !analytically integrates the 1/om^(3/2) tail from bm to +oo for M11, M22 and M12
  REAL(QP), INTENT(IN)  :: bm
  REAL(QP) :: intinfini

  if((floor(num)==1).OR.(floor(num)==2))then
   intinfini= PI*(x0+xq**2.0_qp/4.0_qp)*sqrt(2.0_qp/bm)
  elseif(floor(num)==3)then
   intinfini= -PI*om0*sqrt(2.0_qp/bm)*(1.0_qp+(x0-xq**2.0_qp/4.0_qp)/(3.0_qp*bm))
  else 
   intinfini=0.0_qp
  endif
  END FUNCTION intinfini

  FUNCTION inter(om,arg)
! energy integrand
  IMPLICIT NONE
  
  REAL(QP), INTENT(IN), DIMENSION(:)  ::  om,arg
  REAL(QP), DIMENSION(size(om))       ::  inter
  REAL(QP) r,s,ome,cterm
  INTEGER is
  
  s=arg(1)
  
  inter(:)=0.0_qp

  do is=1,size(om)
   ome=om(is)

   r =rhopp(num,ome,T)

   if(num<3.0_qp)then !compute the counter-term for M++ and M--
    if(T<0.0_qp)then 
     cterm=PI*sqrt(ome)/sqrt(8.0_qp*(1.0_qp+(ome/2.0_qp-x0)**2.0_qp)) ! Contre-terme 1/2eps_k pour M++ et M-- à T=0
    else
     cterm=-PI*(tanh(beta*sqrt(1.0_qp+(ome/2.0_qp-x0)**2.0_qp)/2.0_qp)-1.0_qp)*sqrt(ome)&
     /sqrt(8.0_qp*(1.0_qp+(ome/2.0_qp-x0)**2.0_qp)) ! Contre-terme 1/2eps_k pour M++ et M-- à T\=0
    endif
   else
    cterm=0.0_qp
   endif

   if(s>0.0_qp)then !compensate the singularity in om0
    inter(is)=-sig(floor(num))*r/(om0+ome)+(r-r0)/(om0-ome)+cterm
   else !no singularity to compensate
    inter(is)=-sig(floor(num))*r/(om0+ome)+r     /(om0-ome)+cterm
   endif

   if(isnan(inter(is)))then
    write(6,*) 'isnan(inter(is))'
    write(6,*)"inter(is)=",inter(is)
    write(6,*)"ome,om0,num,T,r,r0,inter(is),opp=",ome,om0,num,T,r,r0,inter(is),opp
    stop
   endif
   if(bla2)then
    write(6,*)"om0,ome,num,inter(is)=",om0,ome,floor(num),inter(is)*ome**(1.5_qp)
   endif

  enddo
  END FUNCTION inter
END FUNCTION intpp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION deriveepp(num,om0,e,T)
USE recettes
IMPLICIT NONE
REAL(QP), INTENT(IN)  :: num,om0,e,T
REAL(QP) deriveepp
REAL(QP) r0 !spectral density in om0
REAL(QP) bmax,db,omcrit!bmax: energy cutoff, omcrit: critical value of om after which racinfq should be used to integrate the tail
REAL(QP) Iinf!partial integrals
REAL(QP) arg(1:1,1:10),bornes(1:10)
REAL(QP) ag(1:4)
INTEGER choix(1:10)


if(bla1)then
 write(6,*)
 write(6,*)"deriveepp avec num, om0, T=",floor(num),om0,T
 write(6,*)
endif

arg(:,:)=bidon
omcrit=100.0_qp !at om0>omcrit switch from midpntq to racinfq to integration the ~1/om^(3/2) tail between om(3) and 2.0_qp*om0-opp(3)

if(T<0.0_qp)then ! Intégrale du "+1" de 1+nP+nM

 if(floor(num)<4)then !energy cutoff in the integrals
  db=100.0_qp
  bmax=max(2.0e7_qp,2*om0-opp(3)+db,2*opp(3)+db,4*x0+db,2*om0-2*x0+db,2*om0+db) !In case 2*om0-opp or 2*opp overflows bmax. This formula works whatever the value of ptbranchmtpp
 else
  bmax=1.0e55_qp
 endif

 Iinf=0.0_qp
 Iinf=intinfini(bmax) !analytic integration of the 1/om^(5/2) tail from bmax to +oo for dM11, dM22 and 1/om^(3/2) for M12

 call ecrit(bla1,'Iinf=',Iinf)

 r0=rhopp(num,om0,T)
 if(x0<x0crit)then
   if(ptbranchmtpp==1)then
    choix(1:1)=             (/rinf/)
    deriveepp=decoupe(inter,(/opp(1),bmax/),arg(1:1,1:1),choix(1:1),EPSpp,bla1)
   elseif(ptbranchmtpp==2)then
    choix(1:3)=             (/mpnt,   msql,   rinf/)
    deriveepp=decoupe(inter,(/opp(1), opp(2), 2*opp(2),bmax/),arg(1:1,1:3),choix(1:3),EPSpp,bla1)
   elseif(ptbranchmtpp==3)then
    choix(1:3)=             (/mpnt,   msql,   rinf/)
    deriveepp=decoupe(inter,(/opp(1), opp(2), opp(3),bmax/),arg(1:1,1:3),choix(1:3),EPSpp,bla1)
   endif
 else!x0>x0crit
   if(ptbranchmtpp==1)then
    ag=(/opp(1),2*x0,-bidon,-bidon/)
    call tri(ag)
    choix(1:4)=             (/msql,  msqu,                 msql, rinf/)
    deriveepp=decoupe(inter,(/ag(1), (ag(1)+ag(2))/2.0_qp, ag(2),2*ag(2),bmax/),arg(1:1,1:4),choix(1:4),EPSpp,bla1)
   elseif(ptbranchmtpp==2)then
    ag=(/opp(1),opp(2),2*x0,-bidon/)
    call tri(ag)
    choix(1:5)=             (/msqu, msql, msqu,                msql, rinf/)
    deriveepp=decoupe(inter,(/ag(1),ag(2),(ag(2)+ag(3))/2.0_qp,ag(3),2*ag(3),bmax/),arg(1:1,1:5),choix(1:5),EPSpp,bla1)
   elseif(ptbranchmtpp==3)then
    ag=(/opp(1),opp(2),opp(3),2*x0/)
    call tri(ag)
    choix(1:8)=  (/mpnt, msqu,                msql, msqu,                msql, msqu,                msql, rinf/)
    if(floor(num)==5) choix(2)=msql
    bornes(1:9)= (/ag(1),(ag(1)+ag(2))/2.0_qp,ag(2),(ag(2)+ag(3))/2.0_qp,ag(3),(ag(3)+ag(4))/2.0_qp,ag(4),2*ag(4),bmax/)
    deriveepp=decoupe(inter,bornes(1:9),arg(1:1,1:8),choix(1:8),EPSpp,bla1)
   endif
 endif
 deriveepp=deriveepp+Iinf
 call ecrit(bla1,'deriveepp=',deriveepp)
elseif(T>0.0_qp)then ! Intégrale du terme en nP+nM de 1+nP+nM
 bmax=1111840.0_qp
 r0=rhopp(num,om0,T)
 if(ptbranchmtpp==1)then
  choix(1:2)=         (/mpnt,   minf/)
  deriveepp=decoupe(inter,(/0.0_qp, opp(1),bmax/),arg(1:1,1:2),choix(1:2),EPSpp,bla1)
 elseif(ptbranchmtpp==2)then
  choix(1:3)=         (/mpnt,   mpnt,   minf/)
  deriveepp=decoupe(inter,(/0.0_qp, opp(1), opp(2),bmax/),arg(1:1,1:3),choix(1:3),EPSpp,bla1)
 elseif(ptbranchmtpp==3)then
  choix(1:4)=         (/mpnt,   mpnt,   msql,   minf/)
  deriveepp=decoupe(inter,(/0.0_qp, opp(1), opp(2), opp(3), bmax/),     arg(1:1,1:4),choix(1:4),EPSpp,bla1)
 endif
 call ecrit(bla1,'deriveepp=',deriveepp)
endif

CONTAINS
  FUNCTION intinfini(bm) !analytically integrates the 1/om^(3/2) tail from bm to +oo for M11, M22 and M12
  REAL(QP), INTENT(IN)  :: bm
  REAL(QP) :: intinfini

  if((floor(num)==1).OR.(floor(num)==2))then
   intinfini= -4.0_qp*PI*om0/(3.0_qp*sqrt(2.0_qp)*bm**(3.0_qp/2.0_qp))
  elseif(floor(num)==3)then
!   intinfini= -2.0_qp*PI/(sqrt(2.0_qp*bm))-2.0_qp*PI*(x0-xq**2/4.0_qp)/(3.0_qp*sqrt(2.0_qp*bm**3)) !Subleading order in der_+-
   intinfini= -2.0_qp*PI/(sqrt(2.0_qp*bm))
  else
   stop "Mauvais num dans deriveepp"
  endif
  END FUNCTION intinfini

  FUNCTION inter(om,arg)
! energy integrand
  IMPLICIT NONE
  
  REAL(QP), INTENT(IN), DIMENSION(:)  ::  om,arg
  REAL(QP), DIMENSION(size(om))       ::  inter
  REAL(QP) r,ome
  INTEGER is
  
  inter(:)=0.0_qp

  do is=1,size(om)
   ome=om(is)
   r =rhopp(num,ome,T)
   inter(is)=sig(floor(num))*r/(om0+ome)**2-r     /(om0-ome)**2
  enddo
  END FUNCTION inter
END FUNCTION deriveepp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION rhopp(num,om,T)
!Computes quasip-quasip spectral densities on the real axis
!See appendix B in [1] for numerical integration and Eqs.(66-68) and table 2 in [2] for the T=0 analytic formulas 
!input:
!eigenfrequency om, otherwise same as intpp
USE recettes
IMPLICIT NONE
REAL(QP), INTENT(IN)  :: num,om,T
REAL(QP) rhopp

REAL(QP) rho
REAL(QP) cht,sht,tht,th2t
REAL(QP) xipp(1:3),arg(1:1)
REAL(QP) s1,s2

INTEGER secteur

rhopp=0.0_qp
rho=0.0_qp
arg(1)=bidon

if(om<opp(1))then
 rhopp=0.0_qp
else
 !t variable
 cht=om/2.0_qp
 tht=sqrt(om**2.0_qp-4.0_qp)/om
 sht=tht*cht
 !avoid square
 th2t=tht**2.0_qp

 if(om<opp(2))then
  secteur=1
 elseif(om<opp(3))then!Numerically solve equation (71) in [2] for s1 and s2, WARNING s1 here is -s1 from table 2 in [2]
  secteur=2
  call bornesbrutalpp(om,secteur,xipp) 
  s1=-xipp(1)/sht
  s2= xipp(2)/sht
 else
  secteur=3
  call bornesbrutalpp(om,secteur,xipp)
  s2=xipp(1)/sht
 endif
 s2=min(s2,1.0_qp) !In case s2 overflows 1 when om becomes very large

 if(T<0.0_qp)then
  rho=rhoanalytic(secteur,floor(num))
 else
  if(secteur==1)then
   if((floor(num)==1).OR.(floor(num)==2).OR.(floor(num)==4).OR.(floor(num)==5))then
    rho=2.0_qp*qromoq(integrande,0.0_qp,tht,arg,midsquq,EPSrpp)
   else
    rho=0.0_qp
   endif
  elseif(secteur==2)then
   rho=    qromoq(integrande,-tht,2.0_qp*xipp(1)/om,arg,midsqlq,EPSrpp)
   rho=rho+qromoq(integrande,2.0_qp*xipp(2)/om,tht,arg,midsquq,EPSrpp)
  else
   rho=qromoq(integrande,2.0_qp*xipp(1)/om,tht,arg,midsquq,EPSrpp)
  endif
 endif
endif
rhopp=rho*PI/(2.0_qp*xq*om)
!if(rhopp>1.0e20)then
! write(6,*) 'rhopp=Infini'
! write(6,*)"om,num,T,opp,secteur=",om,num,T,opp,secteur
! stop
!endif
!!!!!!!!!!!!!!!!
CONTAINS
 FUNCTION integrande(y,arg)
 !Numerical integration over y for the spectral densities at nonzero temperature
 !see Eqs.(B6--B14) in [1]
 IMPLICIT NONE
 REAL(QP), INTENT(IN), DIMENSION(:)  :: y
 REAL(QP), INTENT(IN), DIMENSION(:)  :: arg
 REAL(QP), DIMENSION(size(y)) :: integrande
 REAL(QP), DIMENSION(size(y)) :: h,xii,np,nm,denocommun,nume

 integrande(:)=0.0_qp

 xii=y*cht
 h=sqrt((th2t-y**2.0_qp)/(1.0_qp-y**2.0_qp))
 if(T<0.0_qp)then
  denocommun=1.0_qp/sqrt((th2t-y**2.0_qp)*(1.0_qp-y**2.0_qp))
 else
  np=1.0_qp/(1.0_qp+exp(beta*sqrt((xii+cht*h)**2.0_qp+1.0_qp)))
  nm=1.0_qp/(1.0_qp+exp(beta*sqrt((xii-cht*h)**2.0_qp+1.0_qp)))
  denocommun=(np+nm)/sqrt((th2t-y**2.0_qp)*(1.0_qp-y**2.0_qp))
 endif

!See Eqs. (B8--B13) in [1]
 if(num<2.0_qp)then
  nume=2.0_qp/(1.0_qp-y**2.0_qp)
 elseif(num<3.0_qp)then
  nume=2.0_qp*y**2.0_qp/(1.0_qp-y**2.0_qp)
 elseif(num<4.0_qp)then
  nume=2.0_qp*y/(1.0_qp-y**2.0_qp)
 elseif(num<5.0_qp)then
  nume=2.0_qp*cht**2.0_qp*(1.0_qp-y**2.0_qp)
 elseif(num<6.0_qp)then
  nume=2.0_qp*cht
 elseif(num<7.0_qp)then
  nume=2.0_qp*cht*y
 endif

 integrande=nume*denocommun
 END FUNCTION integrande
!!!!!!!!!!!!!!!!
 FUNCTION rhoanalytic(sec,num)
 !Analytic formulas for the spectral densities is sector sec and for the matrix element number num
 !For num=1..3 see 66--68 and table 2 in [2]
 !cht,sht,tht and th2t must have been defined
 !s1 and s2 must have been defined (WARNING s1 here is -s1 from table 2 in [2])
 INTEGER, INTENT(IN) :: sec,num
 REAL(QP) rhoanalytic
 REAL(QP) rrC(1:2)
 
 rhoanalytic=0.0_qp
 if    (num==1)then
  rhoanalytic=2.0_qp*rBp(sec)
 elseif(num==2)then
  rhoanalytic=2.0_qp*(rBp(sec)-rA(sec))
 elseif(num==3)then
  rrC=rCC(sec)
  rhoanalytic=2.0_qp*tht*(rrC(1)+rrC(2)) 
 elseif(num==4)then
  rhoanalytic=2.0_qp*rB(sec)
 elseif(num==5)then
  rhoanalytic=2.0_qp*cht*rA(sec)
 elseif(num==6)then
  rhoanalytic=2.0_qp*sht*rDD(sec)
 endif
 END FUNCTION rhoanalytic
!!!!!!!!!!!!!!!!!
! !Auxiliary functions
 FUNCTION rA(sec) !Compare with (66--68) in [2] using ellf(PI/2.0_qp,tht)-ellf(asin(s1),tht)=ellf(PI/2.0_qp-asin(s1),iiq*sht)*cht
 INTEGER, INTENT(IN) :: sec
 REAL(QP) rA

 if    (sec==1)then
  rA=   2.0_qp*ellf(PI/2.0_qp,tht)
 elseif(sec==2)then
  rA=   ellf(PI/2.0_qp,tht)-ellf(asin(s1),tht)
  rA=rA+ellf(PI/2.0_qp,tht)-ellf(asin(s2),tht)
 elseif(sec==3)then
  rA=   ellf(PI/2.0_qp,tht)-ellf(asin(s2),tht)
 else 
  stop "Mauvais secteur dans rA"
 endif

 END FUNCTION rA
!!!!!!!!!!!!!!!!
 FUNCTION rB(sec)
 INTEGER, INTENT(IN) :: sec
 REAL(QP) rB

 if    (sec==1)then
  rB=    2.0_qp*elle(PI/2.0_qp,tht)/(1.0_qp-th2t)
 elseif(sec==2)then
  rB=   (elle(PI/2.0_qp,tht)-elle(asin(s1),tht))/(1.0_qp-th2t)
  rB=rB+(elle(PI/2.0_qp,tht)-elle(asin(s2),tht))/(1.0_qp-th2t)
 elseif(sec==3)then
  rB=   (elle(PI/2.0_qp,tht)-elle(asin(s2),tht))/(1.0_qp-th2t)
 else 
  stop "Mauvais secteur dans rB"
 endif

 END FUNCTION rB
!!!!!!!!!!!!!!!!
! FUNCTION rAbis(sec) !Alternate version of the elliptic integrals. Disactivated to avoid calling elle,ellf with complex arguments
! INTEGER, INTENT(IN) :: sec
! REAL(QP) rAbis
! 
! if    (sec==1)then 
!  rAbis=   2.0_qp*ellf(PI/2.0_qp,tht)
! elseif(sec==2)then
!  rAbis=   ellf(PI/2.0_qp-asin(s1),iiq*sht)*cht
!  rAbis=rAbis+ellf(PI/2.0_qp-asin(s2),iiq*sht)*cht
! elseif(sec==3)then
!  rAbis=   ellf(PI/2.0_qp-asin(s2),iiq*sht)*cht
! endif
! END FUNCTION rAbis
!!!!!!!!!!!!!!!!!
! FUNCTION rBbis(sec)
! INTEGER, INTENT(IN) :: sec
! REAL(QP) rBbis
! 
! if    (sec==1)then 
!  rBbis=   2.0_qp*elle(PI/2.0_qp,tht)/(1.0_qp-th2t)
! elseif(sec==2)then
!  rBbis=   elle(PI/2.0_qp-asin(s1),iiq*sht)*cht
!  rBbis=rBbis+elle(PI/2.0_qp-asin(s2),iiq*sht)*cht
! elseif(sec==3)then
!  rBbis=   elle(PI/2.0_qp-asin(s2),iiq*sht)*cht
! endif
! END FUNCTION rBbis
!!!!!!!!!!!!!!!!
 FUNCTION rCC(sec)
 INTEGER, INTENT(IN) :: sec
 REAL(QP) rCC(1:2)
 
 if    (sec==1)then 
  rCC(:)=   0.0_qp
 elseif(sec==2)then
  rCC(1)=(-sqrt((1.0_qp-s1**2.0_qp)/(1.0_qp-th2t*s1**2.0_qp)))/(1.0_qp-th2t)
  rCC(2)=( sqrt((1.0_qp-s2**2.0_qp)/(1.0_qp-th2t*s2**2.0_qp)))/(1.0_qp-th2t)
 elseif(sec==3)then
  rCC(1)=0.0_qp
  rCC(2)=( sqrt((1.0_qp-s2**2.0_qp)/(1.0_qp-th2t*s2**2.0_qp)))/(1.0_qp-th2t)
 else 
  stop "Mauvais secteur dans rCC"
 endif
 END FUNCTION rCC
!!!!!!!!!!!!!!!!
 FUNCTION rBp(sec)
 INTEGER, INTENT(IN) :: sec
 REAL(QP) rBp
 REAL(QP) rrC(1:2)
 
 if    (sec==1)then
  rBp=rB(sec)
 elseif(sec==2)then
  rrC=rCC(sec)
  rBp=rB(sec)-th2t*(s1*rrC(1)-s2*rrC(2))   
 elseif(sec==3)then
  rrC=rCC(sec)
  rBp=rB(sec)+th2t*s2*rrC(2)
 else 
  stop "Mauvais secteur dans rBp"
 endif
 END FUNCTION rBp
!!!!!!!!!!!!!!!!
 FUNCTION rDD(sec)
 INTEGER, INTENT(IN) :: sec
 REAL(QP) rDD
 
 if    (sec==1)then 
  rDD=   0.0_qp
 elseif(sec==2)then
  rDD=  -log((sqrt(1.0_qp-th2t*s1**2.0_qp)+sqrt(th2t*(1.0_qp-s1**2.0_qp)))/sqrt(1.0_qp-th2t))/tht
  rDD=rDD+log((sqrt(1.0_qp-th2t*s2**2.0_qp)+sqrt(th2t*(1.0_qp-s2**2.0_qp)))/sqrt(1.0_qp-th2t))/tht
 elseif(sec==3)then
  rDD=   log((sqrt(1.0_qp-th2t*s2**2.0_qp)+sqrt(th2t*(1.0_qp-s2**2.0_qp)))/sqrt(1.0_qp-th2t))/tht
 else 
  stop "Mauvais secteur dans rDD"
 endif
 END FUNCTION rDD
!!!!!!!!!!!!!!!!
END FUNCTION rhopp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION intpt(num,om0,e)
!Computes the real part of the quasip-quasihole term of a matrix element as an integral over energies I(om0)=int_0^+oo d(om) rhopt(om)/(om0-om)
!When rhopt(om0)\=0 the integrand has a singularity in om0. The principal part is calculated by subtracting rhopt(om0)/(om0-om) to the integrand.

!input:
!complex energy z=om0+I*e  WARNING: current version of the program limited to the real axis: use e=0.0  
!num: a REAL selector (see elementmat)

!input as global variables 
!xq,x0,beta,temperaturenulle,bla1
!ptbranchmtpt must have been defined by oangpt

USE recettes
IMPLICIT NONE
REAL(QP), INTENT(IN)  :: num,om0,e
REAL(QP) intpt
REAL(QP) bmax,db,Icomp
REAL(QP) r0,arg(1,1:10)
REAL(QP), PARAMETER :: singu=1.0_qp,nonsingu=-1.0_qp
REAL(QP), PARAMETER :: su=1.0_qp,nsu=-1.0_qp
LOGICAL axer
INTEGER choix(1:10)

axer=.TRUE.

r0=rhopt(num,om0)
db=100.0_qp
bmax=max(2320.0_qp,2.0_qp*om0-opt(1)+db,2.0_qp*om0+db) !energy cutoff of the (exponentially small) high-energy tail

Icomp=0.0_qp

if(bla1)then
 write(6,*)
 write(6,*)"intpt avec num, om0=",floor(num),om0
 write(6,*)
endif

if(axer)then
 if(ptbranchmtpt==1)then
  if(om0<opt(1))then!cut in om0,opt(1)
   arg(1,1:4)=         (/su,    su,    nsu,   nsu/)
   choix(1:4)=         (/mpnt,  mpnt,  msql,  mexp/)
   intpt=decoupe(inter,(/0.0_qp,om0,   opt(1),2*opt(1),bmax/),arg(1:1,1:4),choix(1:4),EPSpt,bla1)

   Icomp=-r0*log((opt(1)-om0)/om0)
   call ecrit(bla1,'Icomp=',Icomp)
  else!cut in opt(1),om0,2.0_qp*om0-opt(1)
   arg(1,1:4)=         (/nsu,   su,    su,     nsu/)
   choix(1:4)=         (/mpnt,  mpnt,  mpnt,   mexp/)
   intpt=decoupe(inter,(/0.0_qp,opt(1),om0,    2.0_qp*om0-opt(1),bmax/),arg(1:1,1:4),choix(1:4),EPSpt,bla1)
  endif
 else!cut in om0,2.0_qp*om0
  arg(1,1:3)=         (/su,    su,    nsu/)
  choix(1:3)=         (/mpnt,  mpnt,  mexp/)
  intpt=decoupe(inter,(/0.0_qp,om0,   2.0_qp*om0,bmax/),arg(1:1,1:3),choix(1:3),EPSpt,bla1)

 endif
 intpt=intpt+Icomp
 call ecrit(bla1,'intpt=',intpt)
endif

CONTAINS
  FUNCTION inter(om,arg)
  IMPLICIT NONE
  
  REAL(QP), INTENT(IN), DIMENSION(:)  ::  om,arg
  REAL(QP), DIMENSION(size(om))       ::  inter
  REAL(QP) r,s,ome
  INTEGER is
  
  s=arg(1)
  
  do is=1,size(om)
   ome=om(is)
   r =rhopt(num,ome)
   if(s>0.0_qp)then
     inter(is)=(r-r0)/(om0-ome)-sig(floor(num))*r/(om0+ome)
   else
     inter(is)=r/(om0-ome)-sig(floor(num))*r/(om0+ome)
   endif
   if(isnan(inter(is)))STOP 'isnan(inter(is))'
  enddo
  END FUNCTION inter
END FUNCTION intpt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION rhopt(num,om) 
!Computes quasip-quasihole spectral densities on the real axis -S_ab^{epsilon or omega}
!see equation (B15) in [1]
USE recettes
IMPLICIT NONE
REAL(QP), INTENT(IN)  :: num,om
REAL(QP) rhopt,rho1,rho2
REAL(QP) xipt(1:3),arg(1:2)
REAL(QP) cost,sin2t,tan2t


cost=om/2.0_qp
sin2t=1.0_qp-cost**2.0_qp
tan2t=sin2t/cost**2.0_qp

call bornesbrutalpt(om,xipt)

rhopt=0.0_qp
rho1=0.0_qp
rho2=0.0_qp

rho1=qromoq(integrande,0.0_qp,om/xipt(1)/2.0_qp,arg,midsqlq,EPSrpt)
if((ptbranchmtpt==1).AND.(om<opt(1)))then
 rho2=qromoq(integrande,om/xipt(3)/2.0_qp,om/xipt(2)/2.0_qp,arg,midpntq,EPSrpt)
endif
 
rhopt=-(rho1-rho2)*PI/(2.0_qp*xq*om)

CONTAINS
 FUNCTION integrande(y,arg)
 IMPLICIT NONE
 REAL(QP), INTENT(IN), DIMENSION(:)  :: y
 REAL(QP), INTENT(IN), DIMENSION(:)  :: arg
 REAL(QP), DIMENSION(size(y)) :: integrande
 REAL(QP), DIMENSION(size(y)) :: h,xii,np,nm,denocommun,nume

 xii=cost/y
 h=sqrt((1.0_qp+y**2.0_qp*tan2t)/(1.0_qp-y**2.0_qp))
 np=1.0_qp/(1.0_qp+exp(beta*sqrt((xii+cost*h)**2.0_qp+1.0_qp)))
 nm=1.0_qp/(1.0_qp+exp(beta*sqrt((xii-cost*h)**2.0_qp+1.0_qp)))
 denocommun=(np-nm)/sqrt((1.0_qp+y**2.0_qp*tan2t)*(1.0_qp-y**2.0_qp))

 if(num<2.0_qp)then
  nume=2.0_qp*y**2.0_qp/(1.0_qp-y**2.0_qp)
 elseif(num<3.0_qp)then
  nume=2.0_qp/(1.0_qp-y**2.0_qp)
 elseif(num<4.0_qp)then
  nume=2.0_qp*y/(1.0_qp-y**2.0_qp)
 elseif(num<5.0_qp)then
  nume=2.0_qp*cost**2.0_qp*(1.0_qp-y**2.0_qp)/y**2.0_qp
 elseif(num<6.0_qp)then
  nume=2.0_qp*cost
 elseif(num<7.0_qp)then
  nume=2.0_qp*cost/y
 endif
 integrande=nume*denocommun
 END FUNCTION integrande
 
END FUNCTION rhopt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE bornesbrutalpp(om,sect,xi)
!Trouve les bornes d'intégration sur k des densités spec qp-qp
!ap integrale angulaire dans les secteurs 2 et 3
USE recettes
IMPLICIT NONE
REAL(QP), INTENT(IN)  :: om
INTEGER, INTENT(IN) :: sect
REAL(QP), INTENT(OUT) :: xi(1:2)
REAL(QP) rac

xi=(/-1.e100_qp,-1.e100_qp/)

if(sect==2)then
 rac=rtsafeq(emaxpp,(/om/),1.e-20_qp,opp(4),1.e-18_qp*om)
 xi(1)=rac**2.0_qp-x0+xq**2.0_qp/4.0_qp
 rac=rtsafeq(emaxpp,(/om/),opp(4),1.e9_qp,1.e-18_qp*om)
 xi(2)=rac**2.0_qp-x0+xq**2.0_qp/4.0_qp

elseif(sect==3)then
 rac=rtsafe(emaxpp,(/om/),1.e-20_qp,1.e13_qp,1.e-18_qp*om)
 xi(1)=rac**2.0_qp-x0+xq**2.0_qp/4.0_qp
 xi(2)=-1.e200_qp

else
 STOP 'Mauvais secteur'
endif
END SUBROUTINE bornesbrutalpp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE bornesbrutalpt(om,xi)
!Trouve les bornes d'intégration sur k des densités spec qp-qt
!ap integrale angulaire
USE recettes
IMPLICIT NONE
REAL(QP), INTENT(IN)  :: om
REAL(QP), INTENT(OUT) :: xi(1:3)
REAL(QP) rac

xi=(/-1.e100_qp,-1.e100_qp,-1.e100_qp/)
if(x0-xq**2.0_qp/4.0_qp>0.0_qp)then
 rac=rtsafeq(emaxpt,(/om/),sqrt(x0-xq**2.0_qp/4.0_qp),1.e+20_qp,1.e-18_qp*om)
else
 rac=rtsafeq(emaxpt,(/om/),0.0_qp,1.e+20_qp,1.e-18_qp*om)
endif
xi(1)=rac**2.0_qp-x0+xq**2.0_qp/4.0_qp

if(om<opt(1))then
 rac=rtsafeq(emaxpt,(/-om/),1.e-20_qp,opt(2),1.e-18_qp*om)
 xi(2)=rac**2.0_qp-x0+xq**2.0_qp/4.0_qp
 rac=rtsafeq(emaxpt,(/-om/),opt(2),sqrt(x0-xq**2.0_qp/4.0_qp),1.e-18_qp*om)
 xi(3)=rac**2.0_qp-x0+xq**2.0_qp/4.0_qp
endif
END SUBROUTINE bornesbrutalpt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE emaxpp(x,arg,ene,den)
!Renvoie ekq+(x,u=1)-om et sa dérivée par rapport à x

IMPLICIT NONE
REAL(QP), INTENT(IN) :: x
REAL(QP), DIMENSION(:), INTENT(IN) :: arg
REAL(QP), INTENT(OUT) :: ene,den
REAL(QP)  xiP,xiM,eP,eM,om

om=arg(1)

xiP=x**2.0_qp+x*xq+xq**2.0_qp/4.0_qp-x0
xiM=x**2.0_qp-x*xq+xq**2.0_qp/4.0_qp-x0

eP=sqrt(xiP**2.0_qp+1.0_qp)
eM=sqrt(xiM**2.0_qp+1.0_qp)
ene=eP+eM-om
den=xiP/eP*(2.0_qp*x+xq)+xiM/eM*(2.0_qp*x-xq)
END SUBROUTINE emaxpp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE emaxpt(x,arg,ene,den)
!Renvoie ekq-(x,u=1)-om et sa dérivée par rapport à x

IMPLICIT NONE
REAL(QP), INTENT(IN) :: x
REAL(QP), DIMENSION(:), INTENT(IN) :: arg
REAL(QP), INTENT(OUT) :: ene,den
REAL(QP)  xiP,xiM,eP,eM,om

om=arg(1)

xiP=x**2.0_qp+x*xq+xq**2.0_qp/4.0_qp-x0
xiM=x**2.0_qp-x*xq+xq**2.0_qp/4.0_qp-x0

eP=sqrt(xiP**2.0_qp+1.0_qp)
eM=sqrt(xiM**2.0_qp+1.0_qp)
ene=eP-eM-om
den=xiP/eP*(2.0_qp*x+xq)-xiM/eM*(2.0_qp*x-xq)
END SUBROUTINE emaxpt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION Enpp(x,u)
!Renvoie ekq+

IMPLICIT NONE
REAL(QP), INTENT(IN) :: x,u
REAL(QP)  Enpp
REAL(QP)  xiP,xiM,eP,eM

xiP=x**2.0_qp+x*xq*u+xq**2.0_qp/4.0_qp-x0
xiM=x**2.0_qp-x*xq*u+xq**2.0_qp/4.0_qp-x0

eP=sqrt(xiP**2.0_qp+1.0_qp)
eM=sqrt(xiM**2.0_qp+1.0_qp)
Enpp=eP+eM
END FUNCTION Enpp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION Enpt(x,u)
!Renvoie ekq-

IMPLICIT NONE
REAL(QP), INTENT(IN) :: x,u
REAL(QP)  Enpt
REAL(QP)  xiP,xiM,eP,eM

xiP=x**2.0_qp+x*xq*u+xq**2.0_qp/4.0_qp-x0
xiM=x**2.0_qp-x*xq*u+xq**2.0_qp/4.0_qp-x0

eP=sqrt(xiP**2.0_qp+1.0_qp)
eM=sqrt(xiM**2.0_qp+1.0_qp)
Enpt=eP-eM
END FUNCTION Enpt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE derivpp(x,arg,dEdx,ddEddx)
!Renvoie dE(x,u=1)/dx et d2E(x,u=1)/dx2

IMPLICIT NONE
REAL(QP), INTENT(IN)  :: x
REAL(QP), DIMENSION(:), INTENT(IN)  :: arg
REAL(QP), INTENT(OUT) :: dEdx,ddEddx
REAL(QP) xiP,xiM,eP,eM

xiP=x**2.0_qp+x*xq+xq**2.0_qp/4.0_qp-x0
xiM=x**2.0_qp-x*xq+xq**2.0_qp/4.0_qp-x0

eP=sqrt(xiP**2.0_qp+1.0_qp)
eM=sqrt(xiM**2.0_qp+1.0_qp)

dEdx=xiP/eP*(2.0_qp*x+xq)+xiM/eM*(2.0_qp*x-xq)

ddEddx=(2.0_qp*x+xq)**2.0_qp/eP+(2.0_qp*x-xq)**2.0_qp/eM &
&       +2.0_qp*xiP/eP+2.0_qp*xiM/eM &
&       -xiP**2.0_qp/eP**3.0_qp*(2.0_qp*x+xq)**2.0_qp &
&       -xiM**2.0_qp/eM**3.0_qp*(2.0_qp*x-xq)**2.0_qp
END SUBROUTINE derivpp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE derivpt(x,arg,dEdx,ddEddx)
!Renvoie dE(x,u=1)/dx et d2E(x,u=1)/dx2

IMPLICIT NONE
REAL(QP), INTENT(IN)  :: x
REAL(QP), DIMENSION(:), INTENT(IN)  :: arg
REAL(QP), INTENT(OUT) :: dEdx,ddEddx
REAL(QP) xiP,xiM,eP,eM

xiP=x**2.0_qp+x*xq+xq**2.0_qp/4.0_qp-x0
xiM=x**2.0_qp-x*xq+xq**2.0_qp/4.0_qp-x0

eP=sqrt(xiP**2.0_qp+1.0_qp)
eM=sqrt(xiM**2.0_qp+1.0_qp)

dEdx=xiP/eP*(2.0_qp*x+xq)-xiM/eM*(2.0_qp*x-xq)

ddEddx=(2.0_qp*x+xq)**2.0_qp/eP-(2.0_qp*x-xq)**2.0_qp/eM &
&       +2.0_qp*xiP/eP-2.0_qp*xiM/eM &
&       -xiP**2.0_qp/eP**3.0_qp*(2.0_qp*x+xq)**2.0_qp &
&       +xiM**2.0_qp/eM**3.0_qp*(2.0_qp*x-xq)**2.0_qp

END SUBROUTINE derivpt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE oangpp
!Calcule la position des trois points anguleux opp(1),opp(2),opp(3) de la densité spectrale
!ainsi que la valeur de x, alias k, où min_x E(u=1) est atteint (valeur inscrite dans opp(4)) 
USE recettes
IMPLICIT NONE
REAL(QP) d2en0,x2
opp(3)=2.0_qp*sqrt(1.0_qp+(xq**2.0_qp/4.0_qp-x0)**2.0_qp)
if(xq**2.0_qp/4.0_qp>x0)then
 ptbranchmtpp=1
 opp(1)=opp(3)
 opp(2)=opp(3)
 opp(4)=-1e-31_qp
else
 opp(1)=2.0_qp !min_k,uE=2Delta
 d2en0=-64*x0*(1+x0**2)+48*(1+x0**2)*xq**2-12*x0*xq**4+xq**6    !derivee seconde en k=0
 if(d2en0>0.0_qp)then           !Si elle est positive, le min est en zero
  ptbranchmtpp=2
  opp(2)=opp(3)
  opp(4)=-1e-31_qp
 else
  ptbranchmtpp=3
  x2=rtsafeq(derivpp,(/bidon/),1.e-14_qp,sqrt(x0)+2.0_qp,1.e-19_qp)!Le min est entre 0 et rac(x0)
  opp(2)=Enpp(x2,1.0_qp)
  opp(4)=x2
 endif
endif
!write(6,*)"opp=",opp(1:3)
END SUBROUTINE oangpp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE oangpt
!Calcule la position du point anguleux opt(1) de la densité spectrale
!ainsi que la valeur de x, alias k, où min_x E(u=1) est atteint (valeur inscrite dans opt(2)) 
USE recettes
IMPLICIT NONE
if(xq**2.0_qp/4.0_qp>x0)then
 ptbranchmtpt=0
 opt(1)=0.0_qp
 opt(2)=0.0_qp
else
  ptbranchmtpt=1
  opt(2)=rtsafeq(derivpt,(/bidon/),1.e-17_qp,sqrt(x0-xq**2.0_qp/4.0_qp),1.e-19_qp)!Le max est entre 0 et rac(x0-q^2/4), pts où En=0
  opt(1)=abs(Enpt(opt(2),1.0_qp))
endif
END SUBROUTINE oangpt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calcxqjoin
!Calcule la valeur de xq où opp(2) et opp(3) se rejoingnent
IMPLICIT NONE
if(x0<0)stop "xqjoin n'existe que pour x0>0"
xqjoin=2*sqrt(x0+(sqrt(1+x0**2)-x0)**(1.0_qp/3.0_qp)-(sqrt(1+x0**2)+x0)**(1.0_qp/3.0_qp))
END SUBROUTINE calcxqjoin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE mat_pairfield_pttq(om0,e,det,M,fr)
USE recettes
!mat_pairfield computes the 2*2 fluctuation matrix in the phase-modulus basis in the limit xq<<1/x0 (in the BCS regime)
!input: same for mat

!output: 
!M(1:2,1:2) a 2*2 symmetric matrix in the basis 1:phase, 2:modulus
!det=M(1,1)*M(2,2)-M(1,2)**2, the determinant of the matrix
!fr=1/M the pair propagator

IMPLICIT NONE
COMPLEX(QPC), DIMENSION(1:2,1:2), INTENT(OUT) :: M,fr
COMPLEX(QPC), INTENT(OUT) :: det
REAL(QP), INTENT(IN)  :: om0,e

REAL(QP) ss,ll,nn1,nn2,mm,re_ellpi1,re_ellpi2,im_ellpi1,t0,zeta,rem11,rem22,imm11,imm22
COMPLEX(QPC) fp,fm

call oangpp
call oangpt

if(om0<opp(1)) stop "om n’est pas dans le continuum de brisure de paires"

ll=asinh(x0)
mm =-exp(2.0_qp*ll)

ss=acosh(om0/opp(1))
!write(6,*)"ss=",ss

nn1= exp(ll+ss)
nn2=-exp(ll-ss)
re_ellpi1=ellpi_q(PI/2.0_qp,nn1,mm)
re_ellpi2=ellpi_q(PI/2.0_qp,nn2,mm)

t0=asin(1.0_qp/sqrt(nn1))
im_ellpi1=-PI/((1.0_qp-mm*sin(t0)**2)**(0.5_qp)*(-nn1*sin(2.0_qp*t0)))
!write(6,*)"im_ellpi1=",im_ellpi1

fp=(sinh(ss)+sinh(ll))*(re_ellpi1+iiq*im_ellpi1-re_ellpi2)+cosh(ss)*ellf(PI/2.0_qp,iiq*exp(ll))
!write(6,*)"fp=",fp

ss=-ss

!write(6,*)"ss=",ss
nn1= exp(ll+ss)
nn2=-exp(ll-ss)
re_ellpi1=ellpi_q(PI/2.0_qp,nn1,mm)
re_ellpi2=ellpi_q(PI/2.0_qp,nn2,mm)

if(om0<opp(3))then
 t0=asin(1.0_qp/sqrt(nn1))
 im_ellpi1=PI/((1.0_qp-mm*sin(t0)**2)**(0.5_qp)*(-nn1*sin(2.0_qp*t0)))
! write(6,*)"im_ellpi1=",im_ellpi1
endif

fm=(sinh(ss)+sinh(ll))*(re_ellpi1+iiq*im_ellpi1-re_ellpi2)+cosh(ss)*ellf(PI/2.0_qp,iiq*exp(ll))
!write(6,*)"fm=",fm

ss=-ss

M(1,2)=-PI*sqrt(2*exp(ll))*(fp+fm)
M(1,2)=cmplx(real(-PI*sqrt(2*exp(ll))*(fp+fm)),-PI*rhopp(3.5_qp,om0,-1.0_qp),kind=qpc)
if((opp(1)<om0).AND.(om0<opp(2)))then
 zeta=(om0-opp(1))/xq**2/x0
 rem11=-PI**2*log((1+sqrt(1-zeta))/sqrt(zeta))/xq
 rem22= PI**2*(sqrt(1-zeta)-zeta*log((1+sqrt(1-zeta))/sqrt(zeta)))*xq*x0/2.0_qp
 rem22= rem22+2.0_qp*PI*x0**(1.5_qp)*xq**2*(zeta-1.0_qp/3.0_qp)
 imm11=-PI**3/2/xq
 imm22=-PI**3*zeta*xq*x0/4
 M(1,1)=cmplx(rem11,imm11,kind=qpc)
 M(2,2)=cmplx(rem22,imm22,kind=qpc)
else
 M(1,1)=-PI*sqrt(2*exp(ll))*(fp-fm)/tanh(ss)
 M(1,1)=cmplx(real(-PI*sqrt(2*exp(ll))*(fp-fm)/tanh(ss)),-PI*rhopp(1.5_qp,om0,-1.0_qp),kind=qpc)
 M(2,2)=M(1,1)*tanh(ss)**2
 M(2,2)=cmplx(real(M(1,1)*tanh(ss)**2),-PI*rhopp(2.5_qp,om0,-1.0_qp),kind=qpc)
endif

det=M(1,1)*M(2,2)-M(1,2)**2.0_qp

fr(1,1)=M(2,2)/det
fr(2,2)=M(1,1)/det
fr(1,2)=-M(1,2)/det
fr(2,1)=fr(1,2)

END SUBROUTINE mat_pairfield_pttq
END MODULE dspec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    if(opp(1)>2*x0)then
!     arg(1,1:5)=         (/nsu,    nsu,   nsu,             nsu,     nsu/)
!     choix(1:5)=         (/msqu,   msql,  msqu,            msql,    rinf/)
!     intpp=decoupe(inter,(/0.0_qp, 2*x0,  x0+opp(1)/2.0_qp,opp(1),  2*opp(1),bmax/),arg(1:1,1:5),choix(1:5),EPSpp,bla1)
!    else
!     if(om0<opp(1))then 
!      arg(1,1:5)=         (/nsu,    nsu,   nsu,             nsu, nsu/)
!      choix(1:5)=         (/msqu,   msql,  msqu,            msql,rinf/)
!      intpp=decoupe(inter,(/0.0_qp, opp(1),x0+opp(1)/2.0_qp,2*x0,4*x0,bmax/),arg(1:1,1:5),choix(1:5),EPSpp,bla1)
!     elseif(om0<2*x0)then
!      arg(1,1:5)=         (/nsu,    su,     su,    nsu,  nsu/)
!      choix(1:5)=         (/msqu,   msql,   msqu,  msql, rinf/)
!      intpp=decoupe(inter,(/0.0_qp, opp(1), om0,   2*x0, 4*x0,bmax/),arg(1:1,1:5),choix(1:5),EPSpp,bla1)
!
!      Icomp=-r0*log((2*x0-om0)/(om0-opp(1)))
!      call ecrit(bla1,'Icomp=',Icomp)
!     else
!      arg(1,1:6)=         (/nsu,    nsu,    nsu,                 su,  su,  nsu/)
!      choix(1:6)=         (/msqu,   msql,   msqu,                msql,rinf,rinf/)
!      intpp=decoupe(inter,(/0.0_qp, opp(1), (2*x0+opp(1))/2.0_qp,2*x0,om0, 2*om0-2*x0,bmax/),arg(1:1,1:6),choix(1:6),EPSpp,bla1)
!     endif
!    endif
!    if(opp(2)<2*x0)then
!     if(om0<opp(1))then
!      arg(1,1:6)=         (/nsu,    nsu,    nsu,    nsu,              nsu,   nsu/)
!      choix(1:6)=         (/mpnt,   msqu,   msql,   msqu,             msql,  rinf/)
!      intpp=decoupe(inter,(/0.0_qp, opp(1), opp(2), x0+opp(2)/2.0_qp, 2*x0,  4*x0,bmax/),arg(1:1,1:6),choix(1:6),EPSpp,bla1)
!     elseif(om0<opp(2))then
!      arg(1,1:7)=         (/nsu,    su,     su,    nsu,    nsu,              nsu,   nsu/)
!      choix(1:7)=         (/mpnt,   mpnt,   msqu,  msql,   msqu,             msql,  rinf/)
!      intpp=decoupe(inter,(/0.0_qp, opp(1), om0,   opp(2), x0+opp(2)/2.0_qp, 2*x0,  4*x0,bmax/),arg(1:1,1:7),choix(1:7),EPSpp,bla1)
!  
!      Icomp=-r0*log((opp(2)-om0)/(om0-opp(1)))
!      call ecrit(bla1,'Icomp=',Icomp)
!     elseif(om0<2*x0)then
!      arg(1,1:6)=         (/nsu,    nsu,     su,     su,   nsu,   nsu/)
!      choix(1:6)=         (/mpnt,   msqu,    msql,   msqu, msql,  rinf/)
!      intpp=decoupe(inter,(/0.0_qp, opp(1),  opp(2), om0,  2*x0,  4*x0,bmax/),arg(1:1,1:6),choix(1:6),EPSpp,bla1)
!
!      Icomp=-r0*log((2*x0-om0)/(om0-opp(2)))
!      call ecrit(bla1,'Icomp=',Icomp)
!     else
!      arg(1,1:7)=         (/nsu,    nsu,     nsu,    nsu,             su,    su,  nsu/)
!      choix(1:7)=         (/mpnt,   msqu,    msql,   msqu,            msql,  rinf,rinf/)
!      intpp=decoupe(inter,(/0.0_qp, opp(1),  opp(2), x0+opp(2)/2.0_qp,2*x0,  om0, 2*om0,bmax/),arg(1:1,1:7),choix(1:7),EPSpp,bla1)
!
!      Icomp=-r0*log((2*om0-om0)/(2*x0-om0))
!      call ecrit(bla1,'Icomp=',Icomp)
!     endif
!    else
!     if(om0<opp(1))then
!      arg(1,1:6)=         (/nsu,    nsu,    nsu,    nsu,              nsu,   nsu/)
!      choix(1:6)=         (/mpnt,   msqu,   msql,   msqu,             msql,  rinf/)
!      intpp=decoupe(inter,(/0.0_qp, opp(1), opp(2), x0+opp(2)/2.0_qp, 2*x0,  4*x0,bmax/),arg(1:1,1:6),choix(1:6),EPSpp,bla1)
!     elseif(om0<2*x0)then
!      arg(1,1:7)=         (/nsu,    su,     su,    nsu,    nsu,              nsu,   nsu/)
!      choix(1:7)=         (/mpnt,   mpnt,   msqu,  msql,   msqu,             msql,  rinf/)
!      intpp=decoupe(inter,(/0.0_qp, opp(1), om0,   opp(2), x0+opp(2)/2.0_qp, 2*x0,  4*x0,bmax/),arg(1:1,1:7),choix(1:7),EPSpp,bla1)
!  
!      Icomp=-r0*log((opp(2)-om0)/(om0-opp(1)))
!      call ecrit(bla1,'Icomp=',Icomp)
!     elseif(om0<opp(2))then
!      arg(1,1:6)=         (/nsu,    nsu,     su,     su,   nsu,   nsu/)
!      choix(1:6)=         (/mpnt,   msqu,    msql,   msqu, msql,  rinf/)
!      intpp=decoupe(inter,(/0.0_qp, opp(1),  opp(2), om0,  2*x0,  4*x0,bmax/),arg(1:1,1:6),choix(1:6),EPSpp,bla1)
!
!      Icomp=-r0*log((2*x0-om0)/(om0-opp(2)))
!      call ecrit(bla1,'Icomp=',Icomp)
!     else
!      arg(1,1:7)=         (/nsu,    nsu,     nsu,    nsu,             su,    su,  nsu/)
!      choix(1:7)=         (/mpnt,   msqu,    msql,   msqu,            msql,  rinf,rinf/)
!      intpp=decoupe(inter,(/0.0_qp, opp(1),  opp(2), x0+opp(2)/2.0_qp,2*x0,  om0, 2*om0,bmax/),arg(1:1,1:7),choix(1:7),EPSpp,bla1)
!
!      Icomp=-r0*log((2*om0-om0)/(2*x0-om0))
!      call ecrit(bla1,'Icomp=',Icomp)
!     endif
!    endif
!   elseif(ptbranchmtpp==3)then
!    if(om0<opp(1))then
!     arg(1,1:8)=         (/nsu,    nsu,    nsu,    nsu,                   nsu,   nsu,              nsu,   nsu/)
!     choix(1:8)=         (/mpnt,   mpnt,   msql,   msqu,                  msql,  msqu,             msql,  rinf/)
!     if(floor(num)==5) choix(2)=msql
!     bornes(1:9)=        (/0.0_qp, opp(1), opp(2), (opp(2)+opp(3))/2.0_qp,opp(3),x0+opp(3)/2.0_qp, 2*x0,  4*x0,bmax/)
!     intpp=decoupe(inter,bornes(1:9),arg(1:1,1:8),choix(1:8),EPSpp,bla1)
!    elseif(om0<opp(2))then
!     arg(1,1:9)  =       (/nsu,    su,     su,    nsu,    nsu,                   nsu,   nsu,             nsu,   nsu/)
!     choix(1:9)  =       (/mpnt,   mpnt,   msqu,  msql,   msqu,                  msql,  msqu,            msql,  rinf/)
!     if(floor(num)==5) choix(3)=rinf
!     bornes(1:10)=       (/0.0_qp, opp(1), om0,   opp(2), (opp(3)+opp(2))/2.0_qp,opp(3),x0+opp(3)/2.0_qp,2*x0,  4*x0,bmax/)
!     intpp=decoupe(inter,bornes(1:10),arg(1:1,1:9),choix(1:9),EPSpp,bla1)
! 
!     Icomp=-r0*log((opp(2)-om0)/(om0-opp(1)))
!     call ecrit(bla1,'Icomp=',Icomp)
!    elseif(om0<opp(3))then
!     arg(1,1:8)=         (/nsu,    nsu,     su,     su,   nsu,   nsu,             nsu,   nsu/)
!     choix(1:8)=         (/mpnt,   msqu,    msql,   msqu, msql,  msqu,            msql,  rinf/)
!     bornes(1:9)=        (/0.0_qp, opp(1),  opp(2), om0,  opp(3),x0+opp(3)/2.0_qp,2*x0,  4*x0,bmax/)
!     intpp=decoupe(inter,bornes(1:9),arg(1:1,1:8),choix(1:8),EPSpp,bla1)

!     Icomp=-r0*log((opp(3)-om0)/(om0-opp(2)))
!     call ecrit(bla1,'Icomp=',Icomp)
!    elseif(om0<2*x0)then
!     arg(1,1:8)=         (/nsu,    nsu,     nsu,    nsu,                   su,    su,   nsu,   nsu/)
!     choix(1:8)=         (/mpnt,   msqu,    msql,   msqu,                  msql,  msqu, msql,  rinf/)
!     bornes(1:9)=        (/0.0_qp, opp(1),  opp(2), (opp(3)+opp(2))/2.0_qp,opp(3),om0,  2*x0,  4*x0,bmax/)
!     intpp=decoupe(inter,bornes(1:9),arg(1:1,1:8),choix(1:8),EPSpp,bla1)

!     Icomp=-r0*log((2*x0-om0)/(om0-opp(3)))
!     call ecrit(bla1,'Icomp=',Icomp)
!    else
!     arg(1,1:9)=         (/nsu,    nsu,     su,     su,                    nsu,   nsu,             su,  su,  nsu/)
!     choix(1:9)=         (/mpnt,   msqu,    msql,   msqu,                  msql,  msqu,            msql,rinf,rinf/)
!     bornes(1:10)=       (/0.0_qp, opp(1),  opp(2), (opp(3)+opp(2))/2.0_qp,opp(3),x0+opp(3)/2.0_qp,2*x0,om0, 2*om0,bmax/)
!     intpp=decoupe(inter,bornes(1:10),arg(1:1,1:9),choix(1:9),EPSpp,bla1)
!     
!     Icomp=-r0*log((2*om0-om0)/(om0-2*x0))
!     call ecrit(bla1,'Icomp=',Icomp)
!    endif
!   endif
