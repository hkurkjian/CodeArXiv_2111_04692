MODULE intldc !Contribution of the branch cuts ("ldc") to the fermionic self-energy
USE vars ! For variables share by all subroutines here and in the dspec module (in particular x0=mu/Delta and xq=q/q_Delta, opp(1:4), location of the angular points of the qp-qp branch cut)
USE dspec
IMPLICIT NONE
REAL(QP) EPSom,EPSu,EPSq
LOGICAL bla0
CHARACTER(len=15) donnees
CHARACTER(len=9) donneesq
CHARACTER(len=5) suffixe
CONTAINS
!FUNCTION grille(qmin,qmax)
!REAL(QP) qmin,qmax
!call mat_pairfield(ome,e,det,Mat,Gam)
!END FUNCTION grille
FUNCTION selfE(k,zk)
 REAL(QP), INTENT(IN) :: k,zk
 REAL(QP) selfE(1:3)
 REAL(QP) Iq1(1:3),Iq2(1:3)
 REAL(QP) argq(1:1),e
 LOGICAL err

! !Read input file
! open(10,file="intldc.inp")
!  read(10,*)suffixe
!  read(10,*)x0
!  read(10,*)k
!  read(10,*)zk
! close(10)

 temperaturenulle=.TRUE.

 e=0.0_qp
 bla0=.TRUE.


 call system("rm "// "intq"//suffixe//".dat")

 open(13,file="intq"//suffixe//".dat",POSITION="APPEND")
  write(13,*)"!Valeurs de q,intq pour x0,zk,k=",x0,zk,k
 close(13)

! qromovq(f,a,b,dim,arg,varchange,EPS): Romberg integration of the function f
! with values in R^dim from a to b with precision EPS. arg is a vector of 
! static parameters of the function,
! varchange is a subroutine that performs a change of variable for improper integrals:
! varchange=midpntvq when the function is integrated over a compact interval where it takes finite values
! varchange=midinfvq when b -> +oo (b can be as large as allowed by machine precision) and f decays at least as 1/x^2
! varchange=racinfvq when b -> +oo and f decays as 1/x^(3/2)
! varchange=midsquvq/midsqlvq f has a 1/sqrt(b-x) or 1/sqrt(x-a) (integrable) divergence at the upper/lower bound of the integration interval

 open(11,file="grillex0_10.dat")

 Iq1(:)=0.0_qp
 Iq2(:)=0.0_qp

 argq(1)=bidon !Nothing in argq here
 Iq1=qromovq(intq,0.0_qp ,        10.0_qp,       3,argq,midpntvq,EPSq,err)
 write(6,*)"Iq1=",Iq1
  if(err) write(6,*) "convergence non atteinte dans l’intégrale sur q"
! Iq2=qromovq(intq,10.0_qp,        1000.0_qp,     3,argq,midpntvq,EPSq,err)
  if(err) write(6,*) "convergence non atteinte dans l’intégrale sur q"

 selfE=2.0_qp*PI*(Iq1+Iq2) !Integration sur phi

 open(13,file="intq"//suffixe//".dat",POSITION="APPEND")
  write(13,*)
 close(13)

 close(11)

CONTAINS
  
 FUNCTION intq(q,argq,m) !Computes size(q) points of the function to be integrate over q
  INTEGER,  INTENT(IN) :: m !Always called with m=3: the 3 coefficients of the self-energy matrix
  REAL(QP), INTENT(IN), DIMENSION(:)  ::  q,argq
  REAL(QP)  intq(size(q),m)

  REAL(QP), DIMENSION(1:3) ::  I,Ia,Ib,Ic,Id,Ie
  REAL(QP) bmax,qs
  INTEGER is
 
  write(6,*)"q=",q
  intq(:,:)=0.0_qp
 
  do is=1,size(q)
   qs=q(is) !Current value of q, passed on to the intom function in its "arg" argument 
   xq=qs
   call oangpp

   write(6,*)"ptbranchmtpp, qs=",ptbranchmtpp,qs
   write(6,*)"opp=",opp
  
   Ia(:)=0.0_qp
   Ib(:)=0.0_qp
   Ic(:)=0.0_qp
   Id(:)=0.0_qp
   Ie(:)=0.0_qp
   I(:) =0.0_qp
   intq(is,:)=0.0_qp
   
   bmax=1.e6_qp
  
    if(ptbranchmtpp==1)then !BEC-like behavior: integrated from branch cut lower-edge opp(1) to infinity
      Ib=qromovq(intom,opp(1)         ,bmax                ,3,(/qs/),racinfvq,EPSom,err) !deals with the 1/om^(3/2) decay at large om
      call ecrit(bla0,'Ib=',Ib)
      if(err) write(6,*) "convergence non atteinte dans l’intégrale sur omega"
    elseif(ptbranchmtpp==2)then !One angular point opp(2) besides the lower-edge
      Ib=qromovq(intom,opp(1)         ,opp(2)              ,3,(/qs/),midpntvq,EPSom,err) !Integrate from the edge to the angular point
      call ecrit(bla0,'Ib=',Ib)
      if(err) write(6,*) "convergence non atteinte dans l’intégrale sur omega"
      Ic=qromovq(intom,opp(2)         ,2.0_qp*opp(2)       ,3,(/qs/),midpntvq,EPSom,err) !then from opp(2) to 2*opp(2), this circumscribes the numerical difficulty around opp(2)
      call ecrit(bla0,'Ic=',Ic)
      if(err) write(6,*) "convergence non atteinte dans l’intégrale sur omega"
      Id=qromovq(intom,2.0_qp*opp(2)  ,bmax                ,3,(/qs/),racinfvq,EPSom,err) !then from 2*opp(2) to infinity
      call ecrit(bla0,'Id=',Id)
      if(err) write(6,*) "convergence non atteinte dans l’intégrale sur omega"
    elseif(ptbranchmtpp==3)then !Two angular points opp(2) and opp(3) besides the lower-edge
      Ib=qromovq(intom,opp(1)         ,opp(2)              ,3,(/qs/),midpntvq,EPSom,err)
      call ecrit(bla0,'Ib=',Ib)
      if(err) write(6,*) "convergence non atteinte dans l’intégrale sur omega"
      Ic=qromovq(intom,opp(2)         ,opp(3)              ,3,(/qs/),midpntvq,EPSom,err)
      call ecrit(bla0,'Ic=',Ic)
      if(err) write(6,*) "convergence non atteinte dans l’intégrale sur omega"
      Id=qromovq(intom,opp(3)         ,2.0_qp*opp(3)       ,3,(/qs/),midpntvq,EPSom,err)
      call ecrit(bla0,'Id=',Id)
      if(err) write(6,*) "convergence non atteinte dans l’intégrale sur omega"
      Ie=qromovq(intom,2.0_qp*opp(3)         ,bmax         ,3,(/qs/),racinfvq,EPSom,err)
      call ecrit(bla0,'Ie=',Ie)
      if(err) write(6,*) "convergence non atteinte dans l’intégrale sur omega"
    else
     STOP "Erreur de ptbranchmntpp"
    endif
  
   I=Ib+Ic+Id+Ie !Combines the integration intervals
   write(6,*)"I=",I

   open(13,file="intq"//suffixe//".dat",POSITION="APPEND")
    write(13,*)qs,I
   close(13)

   intq(is,:)=I(:)*qs**2 !Jacobian of the q integration
  
  enddo
 END FUNCTION intq

 FUNCTION intom(om,arg,m) !Computes size(om) points of the function  to be integrate over om
  IMPLICIT NONE
  INTEGER,  INTENT(IN) :: m !m=3 here
  REAL(QP), INTENT(IN), DIMENSION(:)  ::  om,arg !arg(1) should be the value of q
  REAL(QP), DIMENSION(size(om),m)       ::  intom

  COMPLEX(QPC) Gam(1:2,1:2),Mat(1:2,1:2),Mat2(1:2,1:2),MatCat(1:2,1:2),det
  REAL(QP) reM11,reM22,reM12,reM21,imM11,imM22,imM12,imM21,omfi,xqfi
  REAL(QP) q,Iu(1:3),argintu(1:2),rho(1:2,1:2),ome
  REAL(QP) deb,fin
  INTEGER is

  q=arg(1) !value of q passed on to the intu function
  argintu(2)=q
  intom(:,:)=0.0_qp

  do is=1,size(om)
   ome=om(is)

   argintu(1)=ome-zk !value of the energy denominator passed on to the intu function

   Iu(:)=0.0_qp
!  call cpu_time(deb)
   Iu=qromovq2(intu,-1.0_qp,1.0_qp,3,argintu,midpntvq,EPSu,err) !computes int_-1^1 du (V^2,U^2,UV)/(ome-z+eps)
!   call cpu_time(fin)
!   write(6,*)"Iu, temps écoulé=",deb-fin
!   call cpu_time(deb)
   read(11,*)xqfi,omfi,reM11,reM12,reM21,reM22,imM11,imM12,imM21,imM22
   Mat(1,1)=cmplx(reM11,imM11,kind=qpc)
   Mat(2,2)=cmplx(reM22,imM22,kind=qpc)
   Mat(1,2)=cmplx(reM12,imM12,kind=qpc)
   Mat(2,1)=Mat(1,2)
   det=Mat(1,1)*Mat(2,2)-Mat(1,2)**2
!   call mat_pairfield(ome,e,det,Mat2,Gam)
   if(abs(xq -xqfi)>1.e-20_qp)stop "xq -xqfi"
   if(abs(ome-omfi)>1.e-20_qp)stop "ome-omfi"
   
!   write(6,*)xq,ome
!   write(6,*)xqfi,omfi
!   write(6,*)real(Mat)
!   write(6,*)real(Mat2)
!   write(6,*)imag(Mat)
!   write(6,*)imag(Mat2)

   MatCat(1,1)=(Mat(1,1)+Mat(2,2))/2.0_qp+Mat(1,2)
   MatCat(2,2)=(Mat(1,1)+Mat(2,2))/2.0_qp-Mat(1,2)
   MatCat(1,2)=(Mat(2,2)-Mat(1,1))/2.0_qp

   Gam(1,1)=MatCat(2,2)/det
   Gam(2,2)=MatCat(1,1)/det
   Gam(1,2)=-MatCat(1,2)/det

   rho=-imag(Gam)/PI
   intom(is,1)=rho(1,1)*Iu(1)
   intom(is,2)=rho(2,2)*Iu(2)
   intom(is,3)=rho(1,2)*Iu(3)
!   call cpu_time(fin)
!   write(6,*)"lecture, temps écoulé=",deb-fin
!   write(6,*)"ome,intom=",ome,intom(is,:)!*ome**(3.0_qp/2.0_qp)

!   open(13,file=donneesq,POSITION="APPEND")
!    write(13,*)ome,intom(is,:)
!   close(13)

!   write(6,*)"ome,rho =",ome,rho(1,1),rho(2,2),rho(1,2)
  enddo
 END FUNCTION intom

 FUNCTION intu(u,arg,m) !Computes size(u) points of the function to be integrate over u
  IMPLICIT NONE
  INTEGER,  INTENT(IN) :: m ! m=3 here
  REAL(QP), INTENT(IN), DIMENSION(:)  ::  u,arg !arg(1) should be the energy in the denominator, arg(2) should be the value of q
  REAL(QP), DIMENSION(size(u),m)       ::  intu

  REAL(QP) kmq2,q,en,xi,eps,U2,V2,UV,us
  INTEGER is

  en=arg(1)
  q=arg(2)
  intu(:,:)=0.0_qp

  do is=1,size(u)

   us   =u(is)
   kmq2 =k**2+q**2-2.0_qp*k*q*us
   xi   =kmq2-x0
   eps  =sqrt(xi**2+1.0_qp)
   U2   =0.5_qp*(1.0_qp+xi/eps)
   V2   =1.0_qp-U2
   UV   =0.5_qp/eps

   intu(is,1)=V2   /(en+eps)
   intu(is,2)=U2   /(en+eps)
   intu(is,3)=UV   /(en+eps)

  enddo
 END FUNCTION intu

END FUNCTION selfE
END MODULE intldc
