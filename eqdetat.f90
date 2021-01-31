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
 FUNCTION eF(x0)
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
!À température non nulle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION eFsurmu(bmu,bd)
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
 FUNCTION unsurkFa(bmu,bd)
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
! Près de Tc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION eFsurmuTc(bmu)
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
 FUNCTION unsurkFaTc(bmu)
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
 FUNCTION bmuTc(unskfa)
 USE recettes
 REAL(QP), INTENT(IN) :: unskfa
 REAL(QP) bmuTc,EPS
  

 if(unskfa>3100.0_qp)STOP "unskfa trop grand"
 if(unskfa<-10.0_qp) STOP "unskfa trop petit"

 EPS=1.e-7_qp
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
! SUBROUTINE bmubd(unskfa,TsTc,bmu,bd)
! USE recettes
! REAL(QP), INTENT(IN)  :: unskfa,TsTc
! REAL(QP), INTENT(OUT) :: bmu,bd
! REAL(QP) bmax,EPS
!  
! bmuc   = bmuTc(unskfa)
! eFsmuC = eFsmuTc(bmuc)
!
! if(unskfa>3100.0_qp)STOP "unskfa trop grand"
! if(unskfa<-10.0_qp) STOP "unskfa trop petit"
!
! if(unskfa<-1.0_qp)then
!  bmuTc=rtsafeq(zero,(/bidoneqdet/),10.0_qp,1.e5_qp,1.e-7_qp) !
! elseif(unskfa>10.0_qp)then
!  bmuTc=rtsafeq(zero,(/bidoneqdet/),-20.0_qp,-4.0_qp,1.e-7_qp) !
! else
!  bmuTc=rtsafeq(zero,(/bidoneqdet/),-5.0_qp,10.0_qp,1.e-7_qp)
! endif
! CONTAINS 
!  SUBROUTINE zero(bmu,arg,f,df)
!  REAL(QP), INTENT(IN) :: bmu
!  REAL(QP), DIMENSION(:), INTENT(IN) :: arg
!  REAL(QP), INTENT(OUT) :: f,df
!  REAL(QP) eps 
!  eps=1.e-5_qp
!  f=unsurkFaTc(bmu)-unskfa
!  df=(unsurkFaTc(bmu+eps)-unsurkFaTc(bmu-eps))/(2.0_qp*eps)
!  END SUBROUTINE zero
! END FUNCTION bmuTc

END MODULE eqdetat
