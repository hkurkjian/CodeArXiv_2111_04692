MODULE intldc !Contribution of the branch cuts ("ldc") to the fermionic self-energy
USE vars ! For variables share by all subroutines here and in the dspec module (in particular x0=mu/Delta and xq=q/q_Delta, opp(1:4), location of the angular points of the qp-qp branch cut)
USE dspec
USE modsim
IMPLICIT NONE
REAL(QP) EPSom,EPSu,EPSq
LOGICAL lecture,ecriture,bla0,st
INTEGER profondeur
CHARACTER(len=15) donnees
CHARACTER(len=30) fichierlec1,fichierlec2,fichierlec3,fichom2,fichom2p
CHARACTER(len=9) donneesq
CHARACTER(len=5) suffixe
REAL(QP) q1,q2,q3,q4
REAL(QP) xi0,xiP,xiM,epsP,epsM,xmin,xmax
REAL(QP) k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12
REAL(QP) l1,l2,l3,l4,l5,l6,l7,l8
INTEGER, PARAMETER :: al=1,bet=2,gam=3,delt=4,epsi=5,alti=6,betti=7,deltti=8,epsiti=9
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION selfEldc(k,zk)
 REAL(QP), INTENT(IN) :: k,zk
 COMPLEX(QPC) selfEldc(1:6)
 COMPLEX(QPC) Iq1(1:6),Iq2(1:6),Iq3(1:6)
 REAL(QP) Iqinf(1:6)
 REAL(QP) argq(1:1),e
 LOGICAL err
 CHARACTER(len=250) chainebidon

 temperaturenulle=.TRUE.

 e=0.0_qp


! qromovq(f,a,b,dim,arg,varchange,EPS,jmax): Romberg integration of the function f
! with values in R^dim from a to b with precision EPS. arg is a vector of 
! static parameters of the function,
! varchange is a subroutine that performs a change of variable for improper integrals:
! varchange=midpntvq when the function is integrated over a compact interval where it takes finite values
! varchange=midinfvq when b -> +oo (b can be as large as allowed by machine precision) and f decays at least as 1/x^2
! varchange=racinfvq when b -> +oo and f decays as 1/x^(3/2)
! varchange=midsquvq/midsqlvq f has a 1/sqrt(b-x) or 1/sqrt(x-a) (integrable) divergence at the upper/lower bound of the integration interval

 write(6,*)"fichierlec1,fichierlec2: ",fichierlec1," ",fichierlec2
 if(lecture)then
  open(11,file=trim(fichierlec1))
  read(11,*)chainebidon
  write(6,*)chainebidon
  open(12,file=trim(fichierlec2))
  read(12,*)chainebidon
  write(6,*)chainebidon
!  open(13,file=trim(fichierlec3))
!  read(13,*)
 endif

 if(ecriture)then
  open (16,file=trim(fichierlec1))
  write(16,*)"! grille de valeur de q,om et Mat de qmin=",q1,"à qmax=",q2," avec profondeur=",profondeur
  open (17,file=trim(fichierlec2))
  write(17,*)"! grille de valeur de q,om et Mat de qmin=",q2,"à qmax=",q3," avec profondeur=",profondeur
  open (18,file=trim(fichierlec3))
  write(18,*)"! grille de valeur de q,om et Mat de qmin=",q3,"à qmax=",q4," avec profondeur=",profondeur
 endif

 Iq1(:)=0.0_qp
 Iq2(:)=0.0_qp
 Iq3(:)=0.0_qp

 argq(1)=1.5_qp 
 Iq1=qromovfixed(intq,q1 ,   q2,   6,argq,midpntvcq,EPSq,profondeur,err)
 write(6,*)"re Iq1=",real(Iq1)
 write(6,*)"im Iq1=",imag(Iq1)
 if(err)  call erreur("q")

 argq(1)=2.5_qp 
 Iq2=qromovfixed(intq,q2,    q3,   6,argq,midpntvcq,EPSq,profondeur,err)
 write(6,*)"re Iq2=",real(Iq2)
 write(6,*)"im Iq2=",imag(Iq2)
 if(err)  call erreur("q")

! argq(1)=3.5_qp 
! Iq3=qromovfixed(intq,q3,    q4,   6,argq,midinfvcq,EPSq,profondeur,err)
! write(6,*)"re Iq3=",real(Iq3)
! write(6,*)"im Iq3=",imag(Iq3)
! if(err)  call erreur("q")

 Iqinf(:)=0.0_qp
 Iqinf(1)=1.0_qp/(2.0_qp*sqrt(3.0_qp)*PI**3*q3**4)
 Iqinf(5)=-Iqinf(1)
! For the other integrals, the large q contribution (vanishing at least of 1/q3**6) is neglected

 selfEldc=2.0_qp*PI*(Iq1+Iq2+Iq3+Iqinf) !Integration sur phi

 open(20,file="intq"//suffixe//".dat",POSITION="APPEND")
  write(20,*)
 close(20)

 close(11)
 close(12)
 close(13)

 close(16)
 close(17)
 close(18)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION intq(q,argq,m) !Computes size(q) points of the function to be integrate over q
  USE nrutil
  INTEGER,  INTENT(IN) :: m !m=6: the 3 coefficients of the 1<->3 self-energy matrix and the 3 coeff of the 4<->0 process
  REAL(QP), INTENT(IN), DIMENSION(:)  ::  q,argq
  COMPLEX(QPC)  intq(size(q),m)

  COMPLEX(QPC), DIMENSION(1:6) ::  I,Ia,Ib,Ic,Id,Ie,It
  REAL(QPC)   , DIMENSION(1:6) ::  Iinf
  REAL(QP) bmax,qs,fich,bmax2
  INTEGER is
 
  fich=argq(1)
  intq(:,:)=0.0_qp
 
  do is=1,size(q)
   qs=q(is) !Current value of q, passed on to the intom function in its "arg" argument 
   xq=qs
   call oangpp

   call ecrit(bla0,"qs=",qs)
   call ecrit(bla0,"ptbranchmtpp=",real(ptbranchmtpp,kind=qpc))
   call ecrit(bla0,"opp=",opp)
  
   xi0=k**2+qs**2-x0
   xiP=k**2+qs**2+2.0_qp*k*qs-x0
   xiM=k**2+qs**2-2.0_qp*k*qs-x0
   epsP=sqrt(xiP**2+1.0_qp)
   epsM=sqrt(xiM**2+1.0_qp)

   xmin=epsP-xiP
   xmax=epsM-xiM

   Ia(:)=0.0_qp
   Ib(:)=0.0_qp
   Ic(:)=0.0_qp
   Id(:)=0.0_qp
   Ie(:)=0.0_qp
   I(:) =0.0_qp
   intq(is,:)=0.0_qp
   
   bmax =1.e6_qp

   Iinf(1)=(xmax-xmin)            /(       sqrt(2.0_qp)*PI**3*k*qs*bmax**(1.0_qp/2.0_qp))
   Iinf(2)=(xmin**(-1)-xmax**(-1))/(9.0_qp*sqrt(2.0_qp)*PI**3*k*qs*bmax**(9.0_qp/2.0_qp))
   Iinf(3)=log(xmax/xmin)         /(5.0_qp*sqrt(2.0_qp)*PI**3*k*qs*bmax**(5.0_qp/2.0_qp))
   Iinf(4)=-Iinf(2)
   Iinf(5)=-Iinf(1)
   Iinf(6)= Iinf(3)
   if(bla0)then
    write(6,*)"Iinf=",Iinf
   endif
  
   if(ptbranchmtpp==1)then !BEC-like behavior: integrated from branch cut lower-edge opp(1) to infinity
     Ib=qromovfixed(intom,opp(1)         ,bmax                ,6,(/qs,fich/),racinfvcq,EPSom,profondeur,err) !deals with the 1/om^(3/2) decay at large om
     call ecrit(bla0,'Ib=',Ib)
     if(err)  call erreur("omega")
   elseif(ptbranchmtpp==2)then !One angular point opp(2) besides the lower-edge
    Ib=qromovfixed(intom,opp(1)         ,opp(2)              ,6,(/qs,fich/),midpntvcq,EPSom,profondeur,err) !Integrate from the edge to the angular point
    call ecrit(bla0,'Ib=',Ib)
    if(err)  call erreur("omega")
    Ic=qromovfixed(intom,opp(2)         ,2.0_qp*opp(2)       ,6,(/qs,fich/),midpntvcq,EPSom,profondeur,err) !then from opp(2) to 2*opp(2), this circumscribes the numerical difficulty around opp(2)
    call ecrit(bla0,'Ic=',Ic)
    if(err)  call erreur("omega")
    Id=qromovfixed(intom,2.0_qp*opp(2)  ,bmax                ,6,(/qs,fich/),racinfvcq,EPSom,profondeur,err) !then from 2*opp(2) to infinity
    call ecrit(bla0,'Id=',Id)
    if(err)  call erreur("omega")
   elseif(ptbranchmtpp==3)then !Two angular points opp(2) and opp(3) besides the lower-edge
    Ib=qromovfixed(intom,opp(1)         ,opp(2)              ,6,(/qs,fich/),midpntvcq,EPSom,profondeur,err)
    call ecrit(bla0,'Ib=',Ib)
    if(err)  call erreur("omega")
    Ic=qromovfixed(intom,opp(2)         ,opp(3)              ,6,(/qs,fich/),midpntvcq,EPSom,profondeur,err)
    call ecrit(bla0,'Ic=',Ic)
    if(err)  call erreur("omega")
    Id=qromovfixed(intom,opp(3)         ,2.0_qp*opp(3)       ,6,(/qs,fich/),midpntvcq,EPSom,profondeur,err)
    call ecrit(bla0,'Id=',Id)
    if(err)  call erreur("omega")
    Ie=qromovfixed(intom,2.0_qp*opp(3)  ,bmax                ,6,(/qs,fich/),racinfvcq,EPSom,profondeur,err)
    call ecrit(bla0,'Ie=',Ie)
    if(err)  call erreur("omega")
   endif
  
   I=Ib+Ic+Id+Ie+Iinf !Combine the integration intervals
   if(bla0)then
    write(6,*)"qs,I=",qs,real(I)
!    write(6,*)"qs,Iinf(1),I(1)=",qs,Iinf(1),real(I(1))
   endif

!   open(20,file="intq"//suffixe//".dat",POSITION="APPEND")
!    write(20,*)qs,real(I(1:3))
!   close(20)

   intq(is,:)=I(:)*qs**2 !Jacobian of the q integration
  
  enddo

 END FUNCTION intq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION intom(om,arg,m) !Computes size(om) points of the function  to be integrate over om
  IMPLICIT NONE
  INTEGER,  INTENT(IN) :: m !m=6 here
  REAL(QP), INTENT(IN), DIMENSION(:)  ::  om,arg !arg(1) should be the value of q
  COMPLEX(QPC), DIMENSION(size(om),m)       ::  intom

  COMPLEX(QPC) Gam(1:2,1:2),Mat(1:2,1:2),Mat2(1:2,1:2),MatCat(1:2,1:2),det
  REAL(QP) reM11,reM22,reM12,reM21,imM11,imM22,imM12,imM21,omfi,xqfi
  REAL(QP) q,argintu(1:2),rho(1:2,1:2),ome,enM,enP
  COMPLEX(QPC) IuP(1:3),IuM(1:3)
  REAL(QP) deb,fin,Iu(1:3)
  INTEGER is,fich

  q=arg(1) !value of q passed on to the intu function
  fich=floor(arg(2))
  argintu(2)=q
  intom(:,:)=0.0_qp

  do is=1,size(om)
   ome=om(is)

!value of the energy denominator passed on to the intu and Iuanaly functions
   enM= ome-zk
   enP= ome+zk

   IuP(:)=0.0_qp
   IuM(:)=0.0_qp
   argintu(1)=enM
!   Iu=qromovq(intu,-1.0_qp,1.0_qp,3,argintu,midpntvq,EPSu) !computes int_-1^1 du (V^2,U^2,UV)/(ome-z+eps)
   IuM=Iuanaly(enM,k,q)
   IuP=Iuanaly(enP,k,q)

   if(lecture)then
    read(10+fich,*)xqfi,omfi,reM11,reM12,reM21,reM22,imM11,imM12,imM21,imM22
    Mat(1,1)=cmplx(reM11,imM11,kind=qpc)
    Mat(2,2)=cmplx(reM22,imM22,kind=qpc)
    Mat(1,2)=cmplx(reM12,imM12,kind=qpc)
    Mat(2,1)=Mat(1,2)
    det=Mat(1,1)*Mat(2,2)-Mat(1,2)**2
    if(abs(xq -xqfi)>1.e-13_qp)then
     stop "xq -xqfi"
    endif
    if(abs(ome-omfi)>1.e-13_qp)then
     write(6,*)"ome,omfi=",ome,omfi
     stop "ome-omfi"
    endif
   else
    call mat_pairfield(ome,e,det,Mat,Gam)
    if(ecriture)then
     write(15+fich,*)q,ome,real(Mat),imag(Mat)
    endif
   endif
   
   MatCat(1,1)=(Mat(1,1)+Mat(2,2))/2.0_qp+Mat(1,2)
   MatCat(2,2)=(Mat(1,1)+Mat(2,2))/2.0_qp-Mat(1,2)
   MatCat(1,2)=(Mat(2,2)-Mat(1,1))/2.0_qp

   Gam(1,1)=MatCat(2,2)/det
   Gam(2,2)=MatCat(1,1)/det
   Gam(1,2)=-MatCat(1,2)/det

   rho=-imag(Gam)/PI

   intom(is,1)=-rho(1,1)*IuM(1)
   intom(is,2)=-rho(2,2)*IuM(2)
   intom(is,3)=-rho(1,2)*IuM(3)

   intom(is,4)= rho(2,2)*IuP(2)
   intom(is,5)= rho(1,1)*IuP(1)
   intom(is,6)=-rho(1,2)*IuP(3)

   if(bla0)then
    write(6,FMT="(A20,8G20.10)")"q,ome,real(intom)=",q,ome,real(intom(is,:))!*ome**(3.0_qp/2.0_qp)
   endif

  enddo
 END FUNCTION intom

! FUNCTION intu(u,arg,m) !Computes size(u) points of the function to be integrate over u
!  IMPLICIT NONE
!  INTEGER,  INTENT(IN) :: m ! m=3 here
!  REAL(QP), INTENT(IN), DIMENSION(:)  ::  u,arg !arg(1) should be the energy in the denominator, arg(2) should be the value of q
!  REAL(QP), DIMENSION(size(u),m)       ::  intu
!
!  REAL(QP) kmq2,q,en,xi,eps,U2,V2,UV,us
!  INTEGER is
!
!  en=arg(1)
!  q=arg(2)
!  intu(:,:)=0.0_qp
!
!  do is=1,size(u)
!
!   us   =u(is)
!   kmq2 =k**2+q**2-2.0_qp*k*q*us
!   xi   =kmq2-x0
!   eps  =sqrt(xi**2+1.0_qp)
!   U2   =0.5_qp*(1.0_qp+xi/eps)
!   V2   =1.0_qp-U2
!   UV   =0.5_qp/eps
!
!   intu(is,1)=V2   /(en+eps)
!   intu(is,2)=U2   /(en+eps)
!   intu(is,3)=UV   /(en+eps)
!
!  enddo
! END FUNCTION intu
END FUNCTION selfEldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION Iuanaly(en,k,q) !Computes Iu analytically
 USE recettes
 IMPLICIT NONE
 REAL(QP), INTENT(IN) :: en,k,q
 COMPLEX(QPC) Iuanaly(1:3)
 COMPLEX(QPC) x1,x2

 REAL(QP) I1,I2,I3,imI1,imI2,rac

 Iuanaly(:)=0.0_qp
 imI1=0.0_qp
 imI2=0.0_qp
 
 if(abs(en)>1.0_qp)then
  rac=sqrt(en**2-1.0_qp)
  x1=cmplx(-en+rac,0.0_qp,kind=qpc)
  x2=cmplx(-en-rac,0.0_qp,kind=qpc)
 else
  rac=sqrt(1.0_qp-en**2)
  x1=cmplx(-en,+rac,kind=qpc)
  x2=cmplx(-en,-rac,kind=qpc)
 endif
 I1=log(abs((xmax-x1)*(xmax-x2)/(xmin-x1)/(xmin-x2)))
 I3=log(xmax/xmin)
 if(abs(en)<1.0_qp)then
  I2=2.0_qp*(argum(xmax-x1)-argum(xmin-x1))/rac
 else
  I2=log(abs((xmax-x1)*(xmin-x2)/(xmax-x2)/(xmin-x1)))/rac
  if((xmax>real(x1)).AND.(xmin<real(x1)))then
   imI1=imI1+PI*sign(1.0_qp,en)
   imI2=imI2+PI*sign(1.0_qp,en)
  endif
  if((xmax>real(x2)).AND.(xmin<real(x2)))then
   imI1=imI1-PI*sign(1.0_qp,en)
   imI2=imI2+PI*sign(1.0_qp,en)
  endif
  imI2=imI2/rac
 endif
 Iuanaly(1)=( (I1+iiq*imI1)-en*(I2+iiq*imI2))          /(4*k*q) !I_V^2
 Iuanaly(2)=(-(I1+iiq*imI1)-en*(I2+iiq*imI2)+2.0_qp*I3)/(4*k*q) !I_U^2
 Iuanaly(3)=  (I2+iiq*imI2)                            /(4*k*q) !I_UV
END FUNCTION Iuanaly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE erreur(var)
CHARACTER(len=*), INTENT(IN) :: var
  if(st)then
   write(6,*) "convergence non atteinte dans l’intégrale sur "//var
   stop
  else
   write(6,*) "convergence non atteinte dans l’intégrale sur "//var
  endif
END SUBROUTINE erreur

! @@
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION intim(k,zk)
 REAL(QP), INTENT(IN) :: k,zk
 COMPLEX(QPC) intim(1:3)
 COMPLEX(QPC) Iq1(1:3),Iq2(1:3),Iq3(1:3)
 REAL(QP) Iqinf(1:3)
 REAL(QP) argq(1:1),e,bq(1:8)
 CHARACTER(len=250) chainebidon

 INTEGER taille
 INTEGER, DIMENSION(1:7) :: config

 CHARACTER(len=2) reg

 temperaturenulle=.TRUE.
 e=0.0_qp

 Iq1(:)=0.0_qp
 Iq2(:)=0.0_qp
 Iq3(:)=0.0_qp

 call lignesenergie(k)
 call region (k,zk-2.0_qp,reg,taille,config)
 if(bla0)then
  write(6,*)"reg=",reg,"  config=",ecritconfig(taille,config) 
 endif
 call bornesq(k,zk,bq)
 if(bla0)then
  write(6,"(A15,8G20.10)")"zk,bq=",zk,bq(1:taille)
  write(6,*)
 endif
!itai=1
!grecque=config(itai)
!Iq1=qromovq(intimq,bq(itai),

END FUNCTION intim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE bornesom(k,zk,q,grecque,res,bom)
USE recettes
REAL(QP), INTENT(IN) :: k,zk,q
INTEGER, INTENT(IN) :: grecque
INTEGER, INTENT(OUT) :: res(1:4)
REAL(QP), INTENT(OUT) :: bom(1:4)
REAL(QP) s

bom(1)=max(2.0_qp,2*epsBCS(q/2))

if(grecque<6)then 
 s=+1.0_qp
else
 s=-1.0_qp
endif

if((grecque==al).OR.(grecque==alti))then
 res(1)=1
 bom(2)=zk-epsBCS(k-s*q)
elseif((grecque==bet).OR.(grecque==betti))then
 res(1)=1
 bom(2)=zk-epsBCS(k-s*q)
 res(2)=2
 bom(3)=zk-1.0_qp
elseif(grecque==gam)then
 res(1)=2
 bom(2)=zk-1.0_qp
elseif((grecque==delt).OR.(grecque==deltti))then
 res(1)=0
 bom(2)=zk-epsBCS(k+s*q)
 res(2)=1
 bom(3)=zk-epsBCS(k-s*q)
elseif((grecque==epsi).OR.(grecque==epsiti))then
 res(1)=0
 bom(2)=zk-epsBCS(k+s*q)
 res(2)=1
 bom(3)=zk-epsBCS(k-s*q)
 res(3)=2
 bom(4)=zk-1.0_qp
endif

END SUBROUTINE bornesom
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE bornesq(k,zk,bq)
USE recettes
REAL(QP), INTENT(IN) :: k,zk
REAL(QP), INTENT(OUT) :: bq(1:8)

REAL(QP) zkt,km,kp,qm,q1m,q2m,q3m,q4m,q1p,q2p,q1mbis,q3mbis,q2pbis,qc,qd
REAL(QP) vecq(1:13)
REAL(QP) sC,dsC

zkt=   zk-2.0_qp
km=    1.0e50_qp  
kp=    1.0e50_qp  
q1m=   1.0e50_qp  
q2m=   1.0e50_qp  
q3m=   1.0e50_qp  
q4m=   1.0e50_qp  

q1mbis=1.0e50_qp  
q3mbis=1.0e50_qp  

q1p=   1.0e50_qp  
q2p=   1.0e50_qp 

qc=    1.0e50_qp
qm=    1.0e50_qp

if(k<k0)then
 km=k0-k
 kp=k+k0
 if((l5>zkt).AND.(zkt>l1)) q1m=k-sqrt(k0**2-sqrt(zkt**2-1))
 if(l5>zkt)                q2m=k+sqrt(k0**2-sqrt(zkt**2-1))
                           q3m=k+sqrt(k0**2+sqrt(zkt**2-1))


 q3mbis=rtsafe(soleC,(/-1.0_qp,1.0_qp/),k+k0,1.e18_qp,1.e-18_qp)
 q3m=min(q3m,q3mbis)

 if(l1>zkt) q1p=-k+sqrt(k0**2-sqrt(zkt**2-1))
 q2p=-k+sqrt(k0**2+sqrt(zkt**2-1))

 q2pbis=rtsafe(soleC,(/ 1.0_qp,1.0_qp/),k0-k,1.e18_qp,1.e-18_qp)
 q2p=min(q2p,q2pbis)

 if(zkt>l4) qc=sqrt(k0**2-k**2)

else
 if(zkt>l1) q1p=rtsafe(soleC,(/ 1.0_qp,1.0_qp/),0.0_qp,1.e18_qp,1.e-18_qp)
 if(zkt>l6)then
  q4m=rtsafe(soleC,(/-1.0_qp,1.0_qp/),k+k0  ,1.e18_qp,1.e-18_qp)
  kp=k0+k
 elseif(zkt>l7)then
  qm=sqrt(4*k0**2+2*sqrt((-1.0_qp + zkt)*(3.0_qp + zkt)))
 endif
 if(k<2*k0)then
  km=k-k0
  if(zkt<l1) q1m=k-sqrt(k0**2+sqrt(zkt**2-1))
  if(zkt<l5) q2m=k-sqrt(k0**2-sqrt(zkt**2-1))
  if((l6>zkt).AND.(zkt>l4))then
    qd   =rtsafe(soleC,(/-1.0_qp,-1.0_qp/),k+1.0e-17_qp   ,k+k0    ,1.e-18_qp)
    q3m   =rtsafe(soleC,(/-1.0_qp, 1.0_qp/),k   ,qd     ,1.e-18_qp)
    q3mbis=rtsafe(soleC,(/-1.0_qp, 1.0_qp/),qd ,k+k0    ,1.e-18_qp)
  endif
  if((l5>zkt).AND.(zkt>l6)) q3m=rtsafe(soleC,(/-1.0_qp,1.0_qp/),k ,k+k0    ,1.e-18_qp)
 elseif(k<3*k0)then
  km=k-k0
  if(zkt<l1) q1m=k-sqrt(k0**2+sqrt(zkt**2-1))
  if(zkt<l6) q2m=rtsafe(soleC,(/-1.0_qp,1.0_qp/),k-k0 ,k+k0    ,1.e-18_qp)
 else
  if((zkt>l7).AND.(l6>zkt))then
   q1m=rtsafe(soleC,(/-1.0_qp,1.0_qp/),0.0_qp,k-k0    ,1.e-18_qp)
   q2m=rtsafe(soleC,(/-1.0_qp,1.0_qp/),k-k0  ,k+k0    ,1.e-18_qp)
   km=k-k0
  elseif((zkt>l6).AND.(l1>zkt))then
   q1m=rtsafe(soleC,(/-1.0_qp,1.0_qp/),0.0_qp,k-k0    ,1.e-18_qp)
   km=k-k0
  elseif((zkt>l8).AND.(l7>zkt))then
   qd=rtsafe(soleC,(/-1.0_qp,-1.0_qp/),0.0_qp,k-k0    ,1.e-18_qp)
   write(6,*)"qd=",qd
   q1m   =rtsafe(soleC,(/-1.0_qp,1.0_qp/),0.0_qp,qd    ,1.e-18_qp)
   q1mbis=rtsafe(soleC,(/-1.0_qp,1.0_qp/),qd    ,k-k0    ,1.e-18_qp)
  else
   km=k-k0
  endif
 endif
endif


vecq=(/0.0_qp,q1m,q2m,q3m,q4m,q1p,q2p,kp,km,qm,qc,q1mbis,q3mbis/)
!write(6,*)"vecq=",real(vecq,sp)
call tri(vecq)
bq=vecq(1:8)


CONTAINS 
  SUBROUTINE soleC(q,arg,x,dx)
  REAL(QP), INTENT(IN) :: q
  REAL(QP), DIMENSION(:), INTENT(IN) :: arg
  REAL(QP), INTENT(OUT) :: x,dx

  REAL(QP) sC,dsC,ddsC
  REAL(QP) dec,ddec,s,derivee
  
  s=arg(1)
  derivee=arg(2)

  sC=zk-ec(q)-epsBCS(k+s*q)
!  write(6,*)"q,sC=",q,sC
  if(q<2*k0)then
   dec =0.0_qp
   ddec=0.0_qp
  else
   dec = deps(q/2)
   ddec=ddeps(q/2)/2
  endif
  dsC  =-dec -s   * deps(k+s*q)
  ddsC =-ddec-s**2*ddeps(k+s*q)
  if(derivee>0.0_qp)then
   x=sC
   dx=dsC
  else
!   write(6,*)"q,dsC=",q,dsC
   x=dsc
   dx=ddsc
  endif  
  END SUBROUTINE soleC
END SUBROUTINE bornesq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE bornesk

k0=sqrt(x0)
k1=k0/sqrt(2.0_qp)
k2=3*k0/5
k3=k0/2
k4=k0/sqrt(5.0_qp)
k5=k0/3
k6=(sqrt(2.0_qp)-1)*k0/2
k7=k0/5

k8 =(1+sqrt(2.0_qp))*k0/2
k9 =sqrt(2.0_qp)*k0
k10=-k0+sqrt(4*k0**2+2*sqrt(k0**4+2*sqrt(1+k0**4)-2))
k11=2*k0
k12=3*k0

END SUBROUTINE bornesk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE region(k,zkt,reg,taille,config)
REAL(QP), INTENT(IN) :: k,zkt
!INTEGER reg
CHARACTER(len=2), INTENT(OUT) :: reg
INTEGER, INTENT(OUT) :: taille
INTEGER, DIMENSION(1:7), INTENT(OUT) :: config

if(k<k0)then
 if(zkt<1)then
  reg="00"
  taille=0
 elseif(min(l1,l4)>zkt)then
  reg="A0"
  taille=6
  config(1:taille)=(/0,alti,betti,gam,bet,al/)
 elseif((zkt>l1).AND.(min(l2,l4)>zkt))then
  reg="B0"
  taille=6
  config(1:taille)=(/deltti,alti,betti,gam,bet,al/)
 elseif((zkt>l4).AND.(min(l1,l2)>zkt))then
  reg="B1"
  taille=7
  config(1:taille)=(/0,alti,betti,epsiti,epsi,bet,al/)
 elseif((zkt>l2).AND.(l4>zkt))then
  reg="C0"
  taille=6
  config(1:taille)=(/deltti,epsiti,betti,gam,bet,al/)
 elseif((l2>zkt).AND.(zkt>max(l1,l4)))then
  reg="C1"
  taille=7
  config(1:taille)=(/deltti,alti,betti,epsiti,epsi,bet,al/)
 elseif((min(l1,l3)>zkt).AND.(zkt>l2))then
  reg="C2"
  taille=7
  config(1:taille)=(/0,alti,deltti,epsiti,epsi,bet,al/)
 elseif((zkt>l4).AND.(l5>zkt).AND.(k>k1))then
  reg="D0"
  taille=7
  config(1:taille)=(/deltti,epsiti,epsi,bet,gam,bet,al/)
 elseif((zkt>max(l2,l4)).AND.(l5>zkt).AND.(k>k3))then
  reg="D1"
  taille=7
  config(1:taille)=(/deltti,epsiti,betti,epsiti,epsi,bet,al/)
 elseif((zkt>max(l1,l2)).AND.(min(l3,l5)>zkt).AND.(k>k7))then
  reg="D2"
  taille=7
  config(1:taille)=(/deltti,alti,deltti,epsiti,epsi,bet,al/)
 elseif((l1>zkt).AND.(zkt> l3))then
  reg="D3"
  taille=7
  config(1:taille)=(/0,alti,deltti,epsiti,epsi,delt,al/)
 elseif((l3>zkt).AND.(zkt>l5))then
  reg="E0"
  taille=5
  config(1:taille)=(/deltti,epsiti,epsi,bet,al/)
 elseif((l5>zkt).AND.(zkt>max(l1,l3)))then
  reg="E1"
  taille=7
  config(1:taille)=(/deltti,alti,deltti,epsiti,epsi,delt,al/)
 elseif(zkt>max(l3,l5))then
  reg="F0"
  taille=5
  config(1:taille)=(/deltti,epsiti,epsi,delt,al/)
 else
  stop "Erreur dans region k<k0"
 endif
else
! write(6,*)"top"
 if(l8>zkt)then
  reg="00"
  taille=0
 elseif((min(l1,l5)>zkt).AND.(zkt>l6))then
  reg="A0"
  taille=6
  config(1:taille)=(/0,al,bet,gam,bet,al/)
 elseif((min(l2,l5)>zkt).AND.(zkt>l1))then
  reg="B0"
  taille=6
  config(1:taille)=(/delt,al,bet,gam,bet,al/)
 elseif((zkt>max(l6,l5)).AND.(l1>zkt))then
  reg="B1"
  taille=4
  config(1:taille)=(/0,al,bet,al/)
 elseif((zkt>l2).AND.(l5>zkt))then
  reg="D0"
  taille=6
  config(1:taille)=(/delt,epsi,bet,gam,bet,al/)
 elseif((l2>zkt).AND.(zkt>max(l1,l5)))then
  reg="D1"
  taille=4
  config(1:taille)=(/delt,al,bet,al/)
 elseif((zkt>max(l2,l5)).AND.(l3>zkt))then
  reg="E0"
  taille=4
  config(1:taille)=(/delt,epsi,bet,al/)
 elseif(zkt>l3)then
  reg="F0"
  taille=4
  config(1:taille)=(/delt,epsi,delt,al/)
 elseif((k<k11).AND.((zkt<l4).OR.((zkt>l5).AND.(l6>zkt))))then
  reg="G0"
  taille=4
  config(1:taille)=(/0,al,bet,gam/)
 elseif((k>k11).AND.(k12>k).AND.(l6>zkt))then
  reg="G0"
  taille=4
  config(1:taille)=(/0,al,bet,gam/)
 elseif((k>k12).AND.(l6>zkt).AND.(zkt>l7))then
  reg="G0"
  taille=4
  config(1:taille)=(/0,al,bet,gam/)
 elseif((k>k12).AND.(l7>zkt).AND.(zkt>l8))then
  reg="H0"
  taille=2
  config(1:taille)=(/0,al/)
 elseif((k<k11).AND.(min(l5,l6)>zkt).AND.(zkt>l4))then
  reg="J0"
  taille=6
  config(1:taille)=(/0,al,bet,gam,bet,gam/)
 else 
  stop "Erreur dans region k>k0"
 endif
endif
END SUBROUTINE region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION ecritconfig(taille,config)
INTEGER, INTENT(IN) :: taille
INTEGER, DIMENSION(:), INTENT(IN) :: config
CHARACTER(len=90) :: ecritconfig

INTEGER itai
ecritconfig=trim(ecritc(config(1)))
do itai=2,taille
 ecritconfig=trim(ecritconfig)//"  "//trim(ecritc(config(itai)))
enddo
CONTAINS 
  FUNCTION ecritc(c)
  INTEGER, INTENT(IN) :: c
  CHARACTER(len=7) :: ecritc
  if(c==0)then
   ecritc="0"
  elseif(c==1)then
   ecritc="al"
  elseif(c==2)then
   ecritc="bet"
  elseif(c==3)then
   ecritc="gam"
  elseif(c==4)then
   ecritc="delt"
  elseif(c==5)then
   ecritc="epsi"
  elseif(c==6)then
   ecritc="al_ti"
  elseif(c==7)then
   ecritc="bet_ti"
  elseif(c==8)then
   ecritc="delt_ti"
  elseif(c==9)then
   ecritc="epsi_ti"
  endif
  END FUNCTION ecritc
END FUNCTION ecritconfig
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE lignesenergie(k)

REAL(QP), INTENT(IN) :: k

l1=epsBCS(k)
l5=epsBCS(0.0_qp)
if(k<k0)then
 l2=epsBCS(2*k-k0)
 l3=epsBCS(2*k+k0)
 l4=epsBCS(k+sqrt(k0**2-k**2))
elseif(k<3*k0)then
 l2=epsBCS(2*k-k0)
 l3=epsBCS(2*k+k0)+ec(k+k0)-2
 if(k<2*k0)then
  l4=solom2(k,fichom2)
 endif
 l6=ec(k+k0)-1
 l8=1.0_qp
else
 l2=epsBCS(2*k-k0)+ec(k-k0)-2
 l3=epsBCS(2*k+k0)+ec(k+k0)-2
 l6=ec(k+k0)-1
 l7=ec(k-k0)-1
 l8=solom2(k,fichom2p)
endif

END SUBROUTINE lignesenergie 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION ec(q)

REAL(QP), INTENT(IN) :: q
REAL(QP) ec

if(q<2*k0)then
 ec=2.0_qp
else
 ec=2*sqrt(1+(q**2/4-x0)**2)
endif

END FUNCTION ec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION xiBCS(k)
REAL(QP), INTENT(IN) :: k
REAL(QP) xiBCS

xiBCS=k**2-x0
END FUNCTION xiBCS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION deps(k)
REAL(QP), INTENT(IN) :: k
REAL(QP) deps

deps=2*k*xiBCS(k)/epsBCS(k)
END FUNCTION deps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION ddeps(k)
REAL(QP), INTENT(IN) :: k
REAL(QP) ddeps

ddeps=4*k**2/epsBCS(k)**3+4*k**2/epsBCS(k)
END FUNCTION ddeps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION epsBCS(k)
REAL(QP), INTENT(IN) :: k
REAL(QP) epsBCS

epsBCS=sqrt((k**2-x0)**2+1)
END FUNCTION epsBCS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION solom2(k,fichdep)
USE recettes
REAL(QP), INTENT(IN) :: k
CHARACTER(len=*), INTENT(IN) :: fichdep
REAL(QP) solom2

REAL(QP) qdep,zdep,klu,vec(1:2)
INTEGER ik

open(25,file=trim(fichdep))
 do ik=1,10000
  read(25,*)klu,qdep,zdep
  if(klu>k)exit
 enddo
close(25)
vec=(/qdep,zdep/)
call mnewt_q(20,vec,1.e-11_qp,1.e-9_qp,solP)
solom2=vec(2)-2
CONTAINS
 SUBROUTINE solP(v,P,JP)
 REAL(QP), DIMENSION(:), INTENT(IN)  :: v
 REAL(QP), DIMENSION(:), INTENT(OUT) :: P
 REAL(QP), DIMENSION(:,:), INTENT(OUT) :: JP

 REAL(QP) q,z

 q=v(1)
 z=v(2)
 P(1)=144 - 96*k**4 + 16*k**8 + 192*k**2*k0**2 - 64*k**6*k0**2 + 288*k0**4 - &
         32*k**4*k0**4 + 192*k**2*k0**6 + 144*k0**8 - 160*z**2 - 32*k**4*z**2 + 64*k**2*k0**2*z**2 - &
         160*k0**4*z**2 + 16*z**4 + (384*k**3 - 128*k**7 - 384*k*k0**2 + 384*k**5*k0**2 + 128*k**3*k0**4 - &
         384*k*k0**6 + 128*k**3*z**2 - 128*k*k0**2*z**2)*q + &
         (-576*k**2 + 448*k**6 - 896*k**4*k0**2 - 320*k**2*k0**4 - 192*k**2*z**2 + 128*k0**2*z**2)*q**2 + &
         (384*k - 896*k**5 + 1024*k**3*k0**2 + 384*k*k0**4 + 128*k*z**2)*q**3 + &
         (-72 + 1112*k**4 - 560*k**2*k0**2 - 72*k0**4 - 40*z**2)*q**4 + (-864*k**3 + 96*k*k0**2)*q**5 + &
         400*k**2*q**6 - 96*k*q**7 + 9*q**8

 P(2)=384*k**3 - 128*k**7 - 384*k*k0**2 + 384*k**5*k0**2 + 128*k**3*k0**4 - 384*k*k0**6 + &
        5*(-864*k**3 + 96*k*k0**2)*q**4 + 2400*k**2*q**5 - 672*k*q**6 + 72*q**7 + 128*k**3*z**2 - &
        128*k*k0**2*z**2 + 4*q**3*(-72 + 1112*k**4 - 560*k**2*k0**2 - 72*k0**4 - 40*z**2) + &
        3*q**2*(384*k - 896*k**5 + 1024*k**3*k0**2 + 384*k*k0**4 + 128*k*z**2) + &
        2*q*(-576*k**2 + 448*k**6 - 896*k**4*k0**2 - 320*k**2*k0**4 - 192*k**2*z**2 + 128*k0**2*z**2)

 JP(1,1)=P(2)

 JP(1,2)=-320*z - 64*k**4*z + 128*k**2*k0**2*z - 320*k0**4*z + 256*k*q**3*z - 80*q**4*z + 64*z**3 + &
        q**2*(-384*k**2*z + 256*k0**2*z) + q*(256*k**3*z - 256*k*k0**2*z)

 JP(2,1)=20*(-864*k**3 + 96*k*k0**2)*q**3 + 12000*k**2*q**4 - 4032*k*q**5 + 504*q**6 + &
        12*q**2*(-72 + 1112*k**4 - 560*k**2*k0**2 - 72*k0**4 - 40*z**2) + &
        6*q*(384*k - 896*k**5 + 1024*k**3*k0**2 + 384*k*k0**4 + 128*k*z**2) + &
        2*(-576*k**2 + 448*k**6 - 896*k**4*k0**2 - 320*k**2*k0**4 - 192*k**2*z**2 + 128*k0**2*z**2)

 JP(2,2)=256*k**3*z - 256*k*k0**2*z + 768*k*q**2*z - 320*q**3*z + 2*q*(-384*k**2*z + 256*k0**2*z)

 END SUBROUTINE solP
END FUNCTION solom2
END MODULE intldc
