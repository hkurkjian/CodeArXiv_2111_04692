MODULE eqdetat
 USE nrtype
 USE recettes, ONLY : elle,ellf
 IMPLICIT NONE
 REAL(QP), PARAMETER ::  bidoneqdet=1.0e200_qp
 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION x1(x0)
 REAL(QP), INTENT(IN) :: x0
 REAL(QP) x1
 x1=sqrt(0.5_qp*(sqrt(1+x0**2.0_qp)+x0))
 END FUNCTION x1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION kappa(x0)
 REAL(QP), INTENT(IN) :: x0
 REAL(QP) kappa
 kappa=sqrt(x1(x0)**2.0_qp/(1.0_qp + x0**2.0_qp)**(0.5_qp))
 END FUNCTION kappa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION I4(x0)
 REAL(QP), INTENT(IN) :: x0
 REAL(QP) I4
 I4=PI/2.0_qp*x1(x0)
 END FUNCTION I4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION I5(x0)
 REAL(QP), INTENT(IN) :: x0
 REAL(QP) I5
 I5=(1.0_qp+x0**2.0_qp)**(0.25_qp)*elle(PI/2.0_qp,kappa(x0))&
 -0.25_qp*ellf(PI/2.0_qp,kappa(x0))/(x1(x0)**2.0_qp*(1+x0**2.0_qp)**(0.25_qp))
 END FUNCTION I5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION I6(x0)
 REAL(QP), INTENT(IN) :: x0
 REAL(QP) I6
 I6=0.5_qp*ellf(PI/2.0_qp,kappa(x0))/(1.0_qp+x0**2.0_qp)**(0.25_qp)
 END FUNCTION I6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION unsurkfa0(x0)
 REAL(QP), INTENT(IN) :: x0
 REAL(QP) unsurkfa0
 unsurkfa0=-4.0_qp/PI*(x0*I6(x0)-I5(x0))/(x0*I5(x0)+I6(x0))**(1.0_qp/3.0_qp)
 END FUNCTION unsurkfa0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION eF(x0)!E_F/Delta at zero temperature
 REAL(QP), INTENT(IN) :: x0
 REAL(QP) eF
 eF=(x0*I5(x0)+I6(x0))**(2.0_qp/3.0_qp)
 END FUNCTION eF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION mc2sDelta(x0)
 REAL(QP), INTENT(IN) :: x0
 REAL(QP) mc2sDelta
 mc2sDelta=2.0_qp*x0/3.0_qp*(I6(x0)/x0/I5(x0) + 1)/(I6(x0)**2/I5(x0)**2 + 1)
 END FUNCTION mc2sDelta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION mu(unsa)
 USE recettes
 REAL(QP), INTENT(IN) :: unsa
 REAL(QP) mu
 mu=rtsafeq(zero,(/bidon/),-1000.0_qp,1000.0_qp,1.0e-15_qp)
 CONTAINS
 SUBROUTINE zero(mmu,arg,f,df)
  REAL(QP), INTENT(IN) :: mmu
  REAL(QP), DIMENSION(:), INTENT(IN) :: arg
  REAL(QP), INTENT(OUT) :: f,df
  REAL(QP) eps
  eps=1.e-5_qp
  f=unsurkFa0(mmu)-unsa
  df=(unsurkFa0(mmu+eps)-unsurkFa0(mmu-eps))/(2.0_qp*eps)
 END SUBROUTINE zero
 END FUNCTION mu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!À température non nulle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION eFsurmu(bmu,bd) !Returns eF/mu(T) given values of mu/T (bmu) and Delta/T (bd)
 USE modsim
 REAL(QP), INTENT(IN) :: bmu,bd
 REAL(QP) eFsurmu,bmax,EPS
  
 EPS=1.0e-9_qp
 bmax=1.0e10_qp
 if(bmu>0.0_qp)then
  eFsurmu=qromo(intxi,0.0_qp,1.0_qp,(/bidoneqdet/),midpntq,EPS)
  eFsurmu=eFsurmu+qromo(intx,sqrt(2.0_qp),bmax,(/bidoneqdet/),midinfq,EPS)
 else
  eFsurmu=qromo(intx,0.0_qp,1.0_qp,(/bidoneqdet/),midpntq,EPS)
  eFsurmu=eFsurmu+qromo(intx,1.0_qp,bmax,(/bidoneqdet/),midinfq,EPS)
 endif
 eFsurmu=(3.0_qp*eFsurmu/2.0_qp)**(2.0_qp/3.0_qp)
 CONTAINS 
  FUNCTION intxi(xi,arg)
   REAL(QP), INTENT(IN), DIMENSION (:) :: xi,arg
   REAL(QP), DIMENSION(size(xi))       ::  intxi,eps
   eps=sqrt(xi**2+(bd/bmu)**2)
   intxi=sqrt(1.0_qp+xi)/2.0_qp*(1.0_qp-xi*tanh(bmu*eps/2.0_qp)/eps)
   intxi=intxi+sqrt(1.0_qp-xi)/2.0_qp*(1.0_qp+xi*tanh(bmu*eps/2.0_qp)/eps)
  END FUNCTION intxi
  FUNCTION intx(x,arg)
   REAL(QP), INTENT(IN), DIMENSION (:) :: x,arg
   REAL(QP), DIMENSION(size(x))       ::  intx,xi,eps
   xi=x**2-abs(bmu)/bmu
   eps=sqrt(xi**2+(bd/bmu)**2)
   intx=x**2*(1.0_qp-xi*tanh(abs(bmu)*eps/2.0_qp)/eps)
  END FUNCTION intx
 END FUNCTION eFsurmu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION unsurkFa(bmu,bd) !Returns 1/kF*a given values of mu/T (bmu) and Delta/T (bd)
 USE modsim
 REAL(QP), INTENT(IN) :: bmu,bd
 REAL(QP) unsurkFa,bmax,EPS
 EPS=1.0e-9_qp
 bmax=1.0e10_qp
 if(bmu>0.0_qp)then
  unsurkFa=qromo(intxi,0.0_qp,1.0_qp,(/bidoneqdet/),midsquq,EPS)
  unsurkFa=unsurkFa+qromo(intx,sqrt(2.0_qp),bmax,(/bidoneqdet/),midinfq,EPS)
 else
  unsurkFa=qromo(intx,0.0_qp,1.0_qp,(/bidoneqdet/),midpntq,EPS)
  unsurkFa=unsurkFa+qromo(intx,1.0_qp,bmax,(/bidoneqdet/),midinfq,EPS)
 endif
 unsurkFa=unsurkFa*4.0_qp/(PI*sqrt(eFsurmu(bmu,bd)))
 CONTAINS 
  FUNCTION intxi(xi,arg)
   REAL(QP), INTENT(IN), DIMENSION (:) :: xi,arg
   REAL(QP), DIMENSION(size(xi))       ::  intxi,eps
   eps=sqrt(xi**2+(bd/bmu)**2)
   intxi=      sqrt(1.0_qp+xi)/2.0_qp*(1.0_qp/(2.0_qp*(1.0_qp+xi))-tanh(bmu*eps/2.0_qp)/(2.0_qp*eps))
   intxi=intxi+sqrt(1.0_qp-xi)/2.0_qp*(1.0_qp/(2.0_qp*(1.0_qp-xi))-tanh(bmu*eps/2.0_qp)/(2.0_qp*eps))
  END FUNCTION intxi
  FUNCTION intx(x,arg)
   REAL(QP), INTENT(IN), DIMENSION (:) :: x,arg
   REAL(QP), DIMENSION(size(x))       ::  intx,xi,eps
   xi=x**2-abs(bmu)/bmu
   eps=sqrt(xi**2+(bd/bmu)**2)
   intx=x**2*(1.0_qp/(2.0_qp*x**2)-tanh(abs(bmu)*eps/2.0_qp)/(2.0_qp*eps))
  END FUNCTION intx
 END FUNCTION unsurkFa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION Ftstc(bmu,bd) !Returns T/Tc given values of mu/T (bmu) and Delta/T (bd)
 USE modsim
 REAL(QP), INTENT(IN) :: bmu,bd
 REAL(QP) Ftstc
 REAL(QP) us,ef,efc,bmuc
  
 us=unsurkFa(bmu,bd)
 ef=eFsurmu(bmu,bd)
! write(6,*)"bmu,bd,ef,us=",bmu,bd,ef,us

 bmuc=bmuTc(us)
! write(6,*)"bmuc=",bmuc
 efc=eFsurmuTc(bmuc)
! write(6,*)"efc=",efc

 Ftstc=efc*bmuc/ef/bmu
 END FUNCTION Ftstc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Près de Tc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION eFsurmuTc(bmu) !Returns eF/mu(Tc) given a value of mu(Tc)/Tc
 USE modsim
 REAL(QP), INTENT(IN) :: bmu
 REAL(QP) eFsurmuTc,bmax,EPS
  
 EPS=1.0e-9_qp
 bmax=1.0e10_qp
 if(bmu>0.0_qp)then
  eFsurmuTc=qromo(intxi,0.0_qp,1.0_qp,(/bidoneqdet/),midpntq,EPS)
  eFsurmuTc=eFsurmuTc+qromo(intx,sqrt(2.0_qp),bmax,(/bidoneqdet/),midinfq,EPS)
 else
  eFsurmuTc=qromo(intx,0.0_qp,1.0_qp,(/bidoneqdet/),midpntq,EPS)
  eFsurmuTc=eFsurmuTc+qromo(intx,1.0_qp,bmax,(/bidoneqdet/),midinfq,EPS)
 endif
 eFsurmuTc=(3.0_qp*eFsurmuTc/2.0_qp)**(2.0_qp/3.0_qp)
 CONTAINS 
  FUNCTION intxi(xi,arg)
   REAL(QP), INTENT(IN), DIMENSION (:) :: xi,arg
   REAL(QP), DIMENSION(size(xi))       ::  intxi
   intxi=sqrt(1.0_qp+xi)/2.0_qp*(1.0_qp-tanh(bmu*xi/2.0_qp))
   intxi=intxi+sqrt(1.0_qp-xi)/2.0_qp*(1.0_qp+tanh(bmu*xi/2.0_qp))
  END FUNCTION intxi
  FUNCTION intx(x,arg)
   REAL(QP), INTENT(IN), DIMENSION (:) :: x,arg
   REAL(QP), DIMENSION(size(x))       ::  intx,xi
   xi=x**2-abs(bmu)/bmu
   intx=x**2*(1.0_qp-abs(xi)/xi*tanh(abs(bmu)*xi/2.0_qp))
  END FUNCTION intx
 END FUNCTION eFsurmuTc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION unsurkFaTc(bmu) !Returns 1/kF*a given a value of mu(Tc)/Tc
 USE modsim
 REAL(QP), INTENT(IN) :: bmu
 REAL(QP) unsurkFaTc,bmax,EPS
  
 EPS=1.0e-9_qp
 bmax=1.0e10_qp
 if(bmu>0.0_qp)then
  unsurkFaTc=qromo(intxi,0.0_qp,1.0_qp,(/bidoneqdet/),midsquq,EPS)
  unsurkFaTc=unsurkFaTc+qromo(intx,sqrt(2.0_qp),bmax,(/bidoneqdet/),midinfq,EPS)
 else
  unsurkFaTc=qromo(intx,0.0_qp,1.0_qp,(/bidoneqdet/),midpntq,EPS)
  unsurkFaTc=unsurkFaTc+qromo(intx,1.0_qp,bmax,(/bidoneqdet/),midinfq,EPS)
 endif
 unsurkFaTc=unsurkFaTc*4.0_qp/(PI*sqrt(eFsurmuTc(bmu)))
 CONTAINS 
  FUNCTION intxi(xi,arg)
   REAL(QP), INTENT(IN), DIMENSION (:) :: xi,arg
   REAL(QP), DIMENSION(size(xi))       ::  intxi
   intxi=      sqrt(1.0_qp+xi)/2.0_qp*(1.0_qp/(2.0_qp*(1.0_qp+xi))-tanh(bmu*xi/2.0_qp)/(2.0_qp*xi))
   intxi=intxi+sqrt(1.0_qp-xi)/2.0_qp*(1.0_qp/(2.0_qp*(1.0_qp-xi))-tanh(bmu*xi/2.0_qp)/(2.0_qp*xi))
  END FUNCTION intxi
  FUNCTION intx(x,arg)
   REAL(QP), INTENT(IN), DIMENSION (:) :: x,arg
   REAL(QP), DIMENSION(size(x))       ::  intx,xi
   xi=x**2-abs(bmu)/bmu
   intx=x**2*(1.0_qp/(2.0_qp*x**2)-tanh(abs(bmu)*xi/2.0_qp)/(2.0_qp*xi))
  END FUNCTION intx
 END FUNCTION unsurkFaTc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION bmuTc(unskfa) !Returns mu(Tc)/Tc given a value of 1/kF*a (inversion of the preceding routine)
 USE recettes
 REAL(QP), INTENT(IN) :: unskfa
 REAL(QP) bmuTc,EPS
  

 if(unskfa>3100.0_qp)STOP "unskfa trop grand"
 if(unskfa<-10.0_qp) STOP "unskfa trop petit"

 EPS=1.e-8_qp
 if(unskfa<-1.0_qp)then
  bmuTc=rtsafeq(zero,(/bidoneqdet/),10.0_qp,1.e5_qp,EPS) !
 elseif(unskfa>10.0_qp)then
  bmuTc=rtsafeq(zero,(/bidoneqdet/),-20.0_qp,-4.0_qp,EPS) !
 else
  bmuTc=rtsafeq(zero,(/bidoneqdet/),-5.0_qp,10.0_qp,EPS)
 endif
 CONTAINS 
  SUBROUTINE zero(bmu,arg,f,df)
  REAL(QP), INTENT(IN) :: bmu
  REAL(QP), DIMENSION(:), INTENT(IN) :: arg
  REAL(QP), INTENT(OUT) :: f,df
  REAL(QP) eps 
  eps=1.e-5_qp
  f=unsurkFaTc(bmu)-unskfa
  df=(unsurkFaTc(bmu+eps)-unsurkFaTc(bmu-eps))/(2.0_qp*eps)
  END SUBROUTINE zero
 END FUNCTION bmuTc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Conversion de unskfa, TsTc à bmu,bd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE bmubd(unskfa,TsTc,bmu,bd)
 USE recettes
 REAL(QP), INTENT(IN)  :: unskfa,TsTc
 REAL(QP), INTENT(INOUT) :: bmu,bd
 REAL(QP) xdep(1:2),bmuc,efc,TsTF
 LOGICAL check
  
! if(unskfa>3100.0_qp)STOP "unskfa trop grand"
! if(unskfa<-10.0_qp) STOP "unskfa trop petit"
 xdep=(/bmu,bd/)
 bmuc=bmuTc(unskfa)
 efc=eFsurmuTc(bmuc)
 TsTF=TsTc/efc/bmuc
 call newt_q(xdep,1.0e-9_qp,check,fnewt)
 bmu=xdep(1)
 bd =xdep(2)
 CONTAINS
  FUNCTION fnewt(x)
  REAL(QP), DIMENSION(:), INTENT(IN) :: x
  REAL(QP), DIMENSION(size(x))  :: fnewt
  fnewt(1)=unskfa-unsurkFa(x(1),x(2))
  fnewt(2)=TsTF-1/x(1)/eFsurmu(x(1),x(2))
!  write(6,*)"x,fnewt=",x,fnewt
  END FUNCTION fnewt
END SUBROUTINE bmubd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE bmubdTF(unskfa,TsTF,bmu,bd)
 USE recettes
 REAL(QP), INTENT(IN)  :: unskfa,TsTF
 REAL(QP), INTENT(INOUT) :: bmu,bd
 REAL(QP) xdep(1:2),bmuc,efc
 LOGICAL check
 xdep=(/bmu,bd/)
 call newt_q(xdep,1.0e-9_qp,check,fnewt)
 bmu=xdep(1)
 bd =xdep(2)
 CONTAINS
  FUNCTION fnewt(x)
  REAL(QP), DIMENSION(:), INTENT(IN) :: x
  REAL(QP), DIMENSION(size(x))  :: fnewt
  fnewt(1)=unskfa-unsurkFa(x(1),x(2))
  fnewt(2)=TsTF-1/x(1)/eFsurmu(x(1),x(2))
  END FUNCTION fnewt
END SUBROUTINE bmubdTF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION sommes(bmu,bd)
 USE recettes
 USE modsim
 REAL(QP), INTENT(IN)  :: bmu,bd
 REAL(QP) :: sommes(1:5)

 REAL(QP) x0,bmax,EPS,S1,S2,S3,S4,S5
  
 EPS=1.0e-9_qp
 bmax=1.0e9_qp
 x0=bmu/bd
 S1=   qromoq(integr,0.0_qp,sqrt(x0),(/1.5_qp/),midsquq,EPS)
 S1=S1+qromoq(integr,sqrt(x0),bmax,  (/1.5_qp/),  midinfq,EPS)
 S1=4*PI*S1
! write(6,*)"S1=",S1
 S2=   qromoq(integr,0.0_qp,sqrt(x0),(/2.5_qp/),midsquq,EPS)
 S2=S2+qromoq(integr,sqrt(x0),bmax,  (/2.5_qp/),  midinfq,EPS)
 S2=4*PI*S2
! write(6,*)"S2=",S2
 S3=   qromoq(integr,0.0_qp,sqrt(x0),(/3.5_qp/),midsquq,EPS)
 S3=S3+qromoq(integr,sqrt(x0),bmax,  (/3.5_qp/),  midexpq,EPS)
 S3=4*PI*bd*S3/2
! write(6,*)"S3=",S3
 S4=   qromoq(integr,0.0_qp,sqrt(x0),(/4.5_qp/),midsquq,EPS)
 S4=S4+qromoq(integr,sqrt(x0),bmax,  (/4.5_qp/),  midinfq,EPS)
 S4=4*PI*bd*S4/2
! write(6,*)"S4=",S4
 S5=   qromoq(integr,0.0_qp,sqrt(x0),(/5.5_qp/),midsquq,EPS)
 S5=S5+qromoq(integr,sqrt(x0),bmax,  (/5.5_qp/),  midinfq,EPS)
 S5=4*PI*bd*S5/2
! write(6,*)"S5=",S5
 sommes=(/S1,S2,S3,S4,S5/)
 CONTAINS
  FUNCTION integr(x,arg)
   REAL(QP), INTENT(IN), DIMENSION (:) :: x,arg
   REAL(QP), DIMENSION(size(x))       ::  integr,xi,eps,Xx

   INTEGER n,is
   xi=x**2-x0
   eps=sqrt(xi**2+1)
   n=floor(arg(1))
   integr(:)=0.0_qp
   if(n<3)then
    Xx=tanh(bd*eps/2)
   else
    Xx=1/cosh(bd*eps/2)**2
   endif
   if(n==1) integr=x**2*Xx*xi   /eps**3
   if(n==2) integr=x**2*Xx      /eps**3
   if(n==3) integr=x**2*Xx*xi   /eps**2
   if(n==4) integr=x**2*Xx*xi**2/eps**2
   if(n==5) integr=x**2*Xx      /eps**2
  END FUNCTION integr
END FUNCTION sommes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION dmudd(bmu,bd)
 USE recettes
 USE modsim
 REAL(QP), INTENT(IN)  :: bmu,bd
 REAL(QP) :: dmudd

 REAL(QP) s(1:5)
 LOGICAL check
  
 s=sommes(bmu,bd)
 dmudd=(s(3)-s(1))/(s(2)+s(4))
END FUNCTION dmudd
END MODULE eqdetat
