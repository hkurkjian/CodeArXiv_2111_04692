MODULE intldc !Contribution of the branch cuts ("ldc") to the fermionic self-energy
USE vars ! For variables share by all subroutines here and in the dspec module (in particular x0=mu/Delta and xq=q/q_Delta, opp(1:4), location of the angular points of the qp-qp branch cut)
USE dspec
USE modsim
USE angularint
IMPLICIT NONE
REAL(QP) EPSom,EPSu,EPSq
LOGICAL lecture,ecriture,bla0,bla00,st
INTEGER profondeur
CHARACTER(len=30) fichierlec1,fichierlec2,fichierlec3,fichom2,fichom2p
CHARACTER(len=50) suffixe_intq1,prefixe,suffixe
REAL(QP) q1,q2,q3,q4
REAL(QP) xiP,xiM,epsP,epsM,xmin,xmax
REAL(QP) k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12
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

 if(bla0)then
   write(6,*)"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   write(6,*)
   write(6,*)"+++++++++++++++++++++++++++++ selfEldc ++++++++++++++++++++++++++++"
   write(6,*)
   write(6,*)"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   write(6,*)
   write(6,*)"k,zk=",k,zk
   write(6,*)"lecture,ecriture=",lecture,ecriture
   write(6,*)"fichierlec1,fichierlec2: ",fichierlec1," ",fichierlec2
   write(6,*)"q1,q2,q3=",q1,q2,q3
   write(6,*)
   write(6,*)"EPSq,EPSom,profondeur=",EPSq,EPSom,profondeur
   write(6,*)
   write(6,*)"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
 endif

! qromovq(f,a,b,dim,arg,varchange,EPS,jmax): Romberg integration of the function f
! with values in R^dim from a to b with precision EPS. arg is a vector of 
! static parameters of the function,
! varchange is a subroutine that performs a change of variable for improper integrals:
! varchange=midpntvq when the function is integrated over a compact interval where it takes finite values
! varchange=midinfvq when b -> +oo (b can be as large as allowed by machine precision) and f decays at least as 1/x^2
! varchange=racinfvq when b -> +oo and f decays as 1/x^(3/2)
! varchange=midsquvq/midsqlvq f has a 1/sqrt(b-x) or 1/sqrt(x-a) (integrable) divergence at the upper/lower bound of the integration interval

 if(lecture)then
  open(11,file=trim(fichierlec1))
  read(11,*)chainebidon
  write(6,*)chainebidon
  open(12,file=trim(fichierlec2))
  read(12,*)chainebidon
  write(6,*)chainebidon
 endif

 if(ecriture)then
  open (16,file=trim(fichierlec1))
  write(16,*)"! grille de valeur de q,om et Mat de qmin=",q1,"à qmax=",q2," avec profondeur=",profondeur
  open (17,file=trim(fichierlec2))
  write(17,*)"! grille de valeur de q,om et Mat de qmin=",q2,"à qmax=",q3," avec profondeur=",profondeur
 endif

 Iq1(:)=0.0_qp
 Iq2(:)=0.0_qp
 Iq3(:)=0.0_qp

 suffixe="sEfixed"
 argq(1)=1.5_qp 
 Iq1=qromovfixed(intq,q1 ,   q2,   6,argq,midpntvcq,EPSq,profondeur,err)
 write(6,*)"re Iq1=",real(Iq1)
 write(6,*)"im Iq1=",imag(Iq1)
 if(err)  call erreur("q")

! argq(1)=2.5_qp 
! Iq2=qromovfixed(intq,q2,    q3,   6,argq,midpntvcq,EPSq,profondeur,err)
! write(6,*)"re Iq2=",real(Iq2)
! write(6,*)"im Iq2=",imag(Iq2)
! if(err)  call erreur("q")

 Iqinf(:)=0.0_qp
! Iqinf(1)=1.0_qp/(2.0_qp*sqrt(3.0_qp)*PI**3*q3**4)
! Iqinf(5)=-Iqinf(1)
! For the other integrals, the large q contribution (vanishing at least of 1/q3**6) is neglected

 selfEldc=2.0_qp*PI*(Iq1+Iq2+Iq3+Iqinf) !Integration sur phi

 open(20,file="intq"//trim(prefixe)//trim(suffixe)//".dat",POSITION="APPEND")
  write(20,*)
 close(20)

 close(11)
 close(12)

 close(16)
 close(17)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION intq(q,argq,m) !Computes size(q) points of the function to be integrate over q
  USE nrutil
  INTEGER,  INTENT(IN) :: m !m=6: the 3 coefficients of the 1<->3 self-energy matrix and the 3 coeff of the 4<->0 process
  REAL(QP), INTENT(IN), DIMENSION(:)  ::  q,argq
  COMPLEX(QPC)  intq(size(q),m)

  COMPLEX(QPC), DIMENSION(1:6) ::  I,Ia,Ib,Ic,Id,Ie
  REAL(QPC)   , DIMENSION(1:6) ::  Iinf
  REAL(QP) bmax,qs,fich
  INTEGER is
 
  fich=argq(1)
  intq(:,:)=0.0_qp
 
  do is=1,size(q)
   qs=q(is) !Current value of q, passed on to the intom function in its "arg" argument 
   xq=qs
   call oangpp

   if(bla0)then
    write(6,*)"----------------------------------------"
    write(6,*)
    write(6,*)"qs=",qs
    write(6,*)"ptbranchmtpp=",ptbranchmtpp
    write(6,*)"opp=",opp
    write(6,*)
   endif
  
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
     write(6,FMT="(A6,6G20.10)")"Iinf=",Iinf
   endif
  
   if(ptbranchmtpp==1)then !BEC-like behavior: integrated from branch cut lower-edge opp(1) to infinity
     Ib=qromovfixed(intom,opp(1)         ,bmax                ,6,(/qs,fich/),racinfvcq,EPSom,profondeur,err) !deals with the 1/om^(3/2) decay at large om
     if(bla0)then
       write(6,FMT="(A10,6G20.10)")'real(Ib)=',real(Ib)
       write(6,FMT="(A10,6G20.10)")'imag(Ib)=',imag(Ib)
     endif
     if(err)  call erreur("omega")

   elseif(ptbranchmtpp==2)then !One angular point opp(2) besides the lower-edge
    Ib=qromovfixed(intom,opp(1)         ,opp(2)              ,6,(/qs,fich/),midpntvcq,EPSom,profondeur,err) !Integrate from the edge to the angular point
    if(bla0)then
      write(6,FMT="(A10,6G20.10)")'real(Ib)=',real(Ib)
      write(6,FMT="(A10,6G20.10)")'imag(Ib)=',imag(Ib)
    endif
    if(err)  call erreur("omega")

    Ic=qromovfixed(intom,opp(2)         ,2.0_qp*opp(2)       ,6,(/qs,fich/),midpntvcq,EPSom,profondeur,err) !then from opp(2) to 2*opp(2), this circumscribes the numerical difficulty around opp(2)
    if(bla0)then
      write(6,FMT="(A10,6G20.10)")'real(Ic)=',real(Ic)
      write(6,FMT="(A10,6G20.10)")'imag(Ic)=',imag(Ic)
    endif
    if(err)  call erreur("omega")

    Id=qromovfixed(intom,2.0_qp*opp(2)  ,bmax                ,6,(/qs,fich/),racinfvcq,EPSom,profondeur,err) !then from 2*opp(2) to infinity
    if(bla0)then
      write(6,FMT="(A10,6G20.10)")'real(Id)=',real(Id)
      write(6,FMT="(A10,6G20.10)")'imag(Id)=',imag(Id)
    endif
    if(err)  call erreur("omega")

   elseif(ptbranchmtpp==3)then !Two angular points opp(2) and opp(3) besides the lower-edge
    Ib=qromovfixed(intom,opp(1)         ,opp(2)              ,6,(/qs,fich/),midpntvcq,EPSom,profondeur,err)
    if(bla0)then
      write(6,FMT="(A10,6G20.10)")'real(Ib)=',real(Ib)
      write(6,FMT="(A10,6G20.10)")'imag(Ib)=',imag(Ib)
    endif
    if(err)  call erreur("omega")

    Ic=qromovfixed(intom,opp(2)         ,opp(3)              ,6,(/qs,fich/),midpntvcq,EPSom,profondeur,err)
    if(bla0)then
      write(6,FMT="(A10,6G20.10)")'real(Ic)=',real(Ic)
      write(6,FMT="(A10,6G20.10)")'imag(Ic)=',imag(Ic)
    endif
    if(err)  call erreur("omega")

    Id=qromovfixed(intom,opp(3)         ,2.0_qp*opp(3)       ,6,(/qs,fich/),midpntvcq,EPSom,profondeur,err)
    if(bla0)then
      write(6,FMT="(A10,6G20.10)")'real(Id)=',real(Id)
      write(6,FMT="(A10,6G20.10)")'imag(Id)=',imag(Id)
    endif
    if(err)  call erreur("omega")

    Ie=qromovfixed(intom,2.0_qp*opp(3)  ,bmax                ,6,(/qs,fich/),racinfvcq,EPSom,profondeur,err)
    if(bla0)then
      write(6,FMT="(A10,6G20.10)")'real(Ie)=',real(Ie)
      write(6,FMT="(A10,6G20.10)")'imag(Ie)=',imag(Ie)
    endif
    if(err)  call erreur("omega")
   endif
  
   I=Ib+Ic+Id+Ie+Iinf !Combine the integration intervals
   if(bla0)then
    write(6,*)
    write(6,FMT="(A12,7G20.10)")"qs,real(I)=",qs,real(I)
    write(6,FMT="(A12,7G20.10)")"qs,imag(I)=",qs,imag(I)
   endif


   intq(is,:)=I(:)*qs**2 !Jacobian of the q integration

   open(20,file="intq"//trim(prefixe)//trim(suffixe)//".dat",POSITION="APPEND")
    write(20,*)qs,real(intq(is,1:6)),imag(intq(is,1:3))
   close(20)
  
  enddo

 END FUNCTION intq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION intom(om,arg,m) !Computes size(om) points of the function  to be integrate over om
  IMPLICIT NONE
  INTEGER,  INTENT(IN) :: m !m=6 here
  REAL(QP), INTENT(IN), DIMENSION(:)  ::  om,arg !arg(1) should be the value of q
  COMPLEX(QPC), DIMENSION(size(om),m)       ::  intom

  COMPLEX(QPC) Gam(1:2,1:2),Mat(1:2,1:2),MatCat(1:2,1:2),det
  REAL(QP) reM11,reM22,reM12,reM21,imM11,imM22,imM12,imM21,omfi,xqfi
  REAL(QP) q,argintu(1:2),rho(1:2,1:2),ome,enM,enP
  COMPLEX(QPC) IuP(1:3),IuM(1:3)
  INTEGER is,fich

  q=arg(1) !value of q passed on to the intu function
  fich=floor(arg(2))
  argintu(2)=q
  intom(:,:)=0.0_qp

  if((bla00).AND.(size(om)==1))then
   write(6,*)
   write(6,*)"****************"
  endif

  do is=1,size(om)
   ome=om(is)

!value of the energy denominator passed on to the intu and Iuanaly functions
   enM= ome-zk
   enP= ome+zk

   IuP(:)=0.0_qp
   IuM(:)=0.0_qp
   argintu(1)=enM
!   Iu=qromovq(intu,-1.0_qp,1.0_qp,3,argintu,midpntvq,EPSu) !computes int_-1^1 du (V^2,U^2,UV)/(ome-z+eps)
   IuM=conjg(Iuanaly(enM,k,q,xmin,xmax)) !enM a une (petite) partie imaginaire négative venant de -zk
   IuP=      Iuanaly(enP,k,q,xmin,xmax)

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

   Gam(1,1)= MatCat(2,2)/det
   Gam(2,2)= MatCat(1,1)/det
   Gam(1,2)=-MatCat(1,2)/det

   rho=-imag(Gam)/PI

   intom(is,1)=-rho(1,1)*IuM(1)
   intom(is,2)=-rho(2,2)*IuM(2)
   intom(is,3)=-rho(1,2)*IuM(3)

   intom(is,4)= rho(2,2)*IuP(2)
   intom(is,5)= rho(1,1)*IuP(1)
   intom(is,6)=-rho(1,2)*IuP(3)

   if(bla00)then
    write(6,FMT="(A19,8G20.10)")"q,ome,real(intom),imag(intom)=",q,ome,&
        real(intom(is,1:6)),imag(intom(is,1:3))!*ome**(3.0_qp/2.0_qp)
   endif

  enddo

  if(bla00)then
   write(6,*)
  endif

 END FUNCTION intom
END FUNCTION selfEldc
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
!FUNCTION selfEldc2(k,zk,interpolation,fich)
!USE estM
!REAL(QP), INTENT(IN) :: k,zk
!LOGICAL, INTENT(IN) :: interpolation
!CHARACTER(len=90), INTENT(IN) :: fich
!COMPLEX(QPC) selfEldc2(1:6)
!
!CHARACTER(len=2) reg
!INTEGER tconf,config(1:7)
!
!REAL(QP) qmax,int6(1:6),le(1:8),bq(1:8)
!
!call bornesk
!call bornesq (k,zk-2.0_qp,le,reg,tconf,config,bq)
!temperaturenulle=.TRUE.
!
!suffixe="sEres"
!open(25,file="intq"//trim(prefixe)//trim(suffixe)//".dat")
!close(25)
!selfEldc2(1:3)=intres(k,zk,interpolation,config(1:tconf),bq(1:tconf+1))
!int6=intpasres(k,zk,interpolation,2,bq(tconf+1),qmax)
!selfEldc2(1:3)=selfEldc2(1:3)+int6(1:3)
!
!suffixe="sEpasres"
!open(25,file="intq"//trim(prefixe)//trim(suffixe)//".dat")
!close(25)
!int6=intpasres(k,zk,interpolation,3,0.0_qp,qmax)
!selfEldc2(4:6)=int6(4:6)
!
!if(interpolation) call unload_data
!
!END FUNCTION selfEldc2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION intres(k,zk,interpolation,fich)
 USE bestM
 USE recettes
 REAL(QP), INTENT(IN) :: k,zk
 LOGICAL,  INTENT(IN) :: interpolation
 CHARACTER(len=90), INTENT(IN), DIMENSION(1:2) :: fich
 COMPLEX(QPC) intres(1:6)

 INTEGER,  ALLOCATABLE, DIMENSION(:) :: config
 REAL(QP), ALLOCATABLE, DIMENSION(:) :: bq

 CHARACTER(len=2) reg
 INTEGER tconf,configbis(1:7)
 REAL(QP) e,qmax,le(1:8),bqbis(1:8)

 CHARACTER(len=250) chainebidon
 INTEGER grecque,igr

 call bornesk
 call bornesq (k,zk-2.0_qp,le,reg,tconf,configbis,bqbis)
 temperaturenulle=.TRUE.
 e=0.0_qp

 allocate(bq(1:tconf+1))
 allocate(config(1:tconf))
 bq    =bqbis(1:tconf+1)
 config=configbis(1:tconf)

 if(interpolation) call load_data(fich(1))
 if(interpolation) call loadom2(fich(2))
 
 if(bla0)then
   write(6,*)"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   write(6,*)
   write(6,*)"++++++++++++++++++++++++++++++ intres +++++++++++++++++++++++++++++"
   write(6,*)
   write(6,*)"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   write(6,*)
   write(6,*)
   write(6,*)"k,zk=",k,zk
   write(6,FMT="(A21,8G20.10)")"lignes d’énergie=",le
   write(6,*)                  "reg=",reg
   write(6,*)                  "config=",ecritconfig(tconf-1,config) 
   write(6,*)                  "bq="    ,bq(1:tconf+1)
 endif

 suffixe="res"
! open(25,file="intq"//trim(prefixe)//trim(suffixe)//".dat")
! close(25)

 intres(:)=cmplx(0.0_qp,0.0_qp,kind=qpc)

 bmax =1.e6_qp
 qmax= 16.0_qp

 if(tconf==0)then
  intres=qromovcq(intresq,0.0_qp,qmax,6,(/bidon/),midpntvcq,EPSq)
 else
  do igr=5,size(config)
   grecque=config(igr)
   if(bla0)then
    write(6,*)"---------------------------------"
    write(6,*)
    write(6,*)"igr=",igr
    write(6,*)"bq(igr),bq(igr+1)=",bq(igr),bq(igr+1)
    write(6,*)
   endif
   intres=intres+qromovcq(intresq,bq(igr),bq(igr+1),6,(/bidon/),midpntvcq,EPSq)
   if(bla0)then
    write(6,*)
    write(6,*)"---------------------------------"
   endif
  enddo
  intres=intres+qromovcq(intresq,bq(igr+1),qmax,6,(/bidon/),midpntvcq,EPSq)
 endif
 intres=2.0_qp*PI*intres
 if(interpolation) call unload_data

 CONTAINS
  FUNCTION intresq(q,argq,m)
  USE nrutil
  USE recettes
  USE modsim
  INTEGER,  INTENT(IN) :: m 
  REAL(QP), INTENT(IN), DIMENSION(:)  ::  q,argq
  COMPLEX(QPC)  intresq(size(q),m)

  REAL(QP) qs,ommil
  REAL(QP) bom(1:3),bom2(1:6)
  INTEGER is,ttot,tres,trout,res(1:2),pos_bom(1:6),p1,p2
  REAL(QP), ALLOCATABLE, DIMENSION(:,:) :: arg
  REAL(QP), ALLOCATABLE, DIMENSION(:)   :: bomf
  INTEGER,  ALLOCATABLE, DIMENSION(:)   :: vres,routint
  INTEGER itai
  COMPLEX(QPC)  intmax(m)

  do is=1,size(q)
  
   qs=q(is)
   xq=qs
   call oangpp
   call bornesom(k,zk,qs,grecque,res,bom,tres)
 
   xiP=k**2+qs**2+2.0_qp*k*qs-x0
   xiM=k**2+qs**2-2.0_qp*k*qs-x0
   epsP=sqrt(xiP**2+1.0_qp)
   epsM=sqrt(xiM**2+1.0_qp)
 
   xmin=epsP-xiP
   xmax=epsM-xiM

   if(bla0)then
    write(6,*)"---------------------------------"
    write(6,*)
    write(6,*)"qs=",qs
    write(6,*)"grecque=",ecritc(grecque)
    write(6,*)"ptbranchmtpp=",ptbranchmtpp
    write(6,*)"opp=",opp
    write(6,*)"res=",res
    write(6,FMT="(A5,3G20.10)")"bom=",bom
   endif

!Combine les points anguleux de opp avec ceux de bom
   bom2(:)=1.e100_qp
   if(grecque==0)then
    ttot=ptbranchmtpp !ttot=nbr tot de points anguleux.
    bom2(1:ptbranchmtpp)=opp(1:ptbranchmtpp)
   elseif((grecque==al).OR.(grecque==alti).OR.(grecque==bet).OR.(grecque==betti).OR.(grecque==gam))then
    ttot=tres+ptbranchmtpp !Dans ce cas bom(1)=opp(1): on évite le double comptage
    bom2(1:ptbranchmtpp)=opp(1:ptbranchmtpp)
    bom2(ptbranchmtpp+1:ptbranchmtpp+tres)=bom(2:tres+1)
   else
    ttot=tres+ptbranchmtpp+1 !Dans ce cas bom(1).NE.opp(1)
    bom2(1:ptbranchmtpp)=opp(1:ptbranchmtpp)
    bom2(ptbranchmtpp+1:ptbranchmtpp+tres+1)=bom(1:tres+1)
   endif
   call tri_pos(bom2,pos_bom)
 
!Découpe les intervalles par le milieu, assigne routint (pour le changement de variable) et vres (pour le nombre d’angles de résonnance)
   trout=2*ttot-1 !trout=nombre d’intervalle d’integration. Nombre de bornes (avec les milieux et bmax: 2*ttot)
   allocate(bomf(1:trout+1))
   allocate(routint(1:trout))
   allocate(vres(1:trout))
   vres(:)   =0
   routint(:)=mpnt
   do itai=1,ttot-1
    p1=pos_bom(itai)
    p2=pos_bom(itai+1)
    ommil=(bom2(itai)+bom2(itai+1))/2
    bomf(2*itai-1)=bom2(itai)
    bomf(2*itai)  =ommil
    if(p1>ptbranchmtpp) routint(2*itai-1)=msql
    if(p2>ptbranchmtpp) routint(2*itai)  =msqu
!    if((ptbranchmtpp==3).AND.(p1==2)) routint(2*itai-1)  =msql
    if(grecque.NE.0)then 
     call locate(bom(1:tres+1),ommil,p1)
     if((p1>0).AND.(p1<tres+1))then
       vres(2*itai-1:2*itai)=res(p1)
     endif
    endif
   enddo
   bomf(trout)  =bom2(ttot)
   bomf(trout+1)=bmax
   routint(trout)=rinf
   vres(trout)=0

   allocate(arg(1:2,1:trout))
   arg(1,:)=qs
   arg(2,:)=vres(:)+0.5_qp

!   bomf(1:6)=bom2
!   trout=ttot-1
!   do itai=ttot-1,1,-1
!    p1=pos_bom(itai)
!    p2=pos_bom(itai+1)
!    if((p1>ptbranchmtpp).AND.(p2.LE.ptbranchmtpp))then
!     routint(itai)=msql
!    elseif((p1.LE.ptbranchmtpp).AND.(p2>ptbranchmtpp))then
!     routint(itai)=msqu
!    elseif((p1.LE.ptbranchmtpp).AND.(p2.LE.ptbranchmtpp))then
!     routint(itai)=mpnt
!    elseif((p1>ptbranchmtpp).AND.(p2>ptbranchmtpp))then
!     trout=trout+1
!     do jtai=ttot+3,itai+1,-1
!      bomf(jtai+1)=bomf(jtai)
!      routint(jtai+1)=routint(jtai)
!     enddo
!     bomf(itai+1)=(bom2(itai)+bom2(itai+1))/2
!     routint(itai+1)=msqu
!     routint(itai)  =msql
!    endif
!   enddo
!   allocate(vres(1:trout))
!   vres(:)=100
!   if(grecque==0)then 
!    vres(:)=0
!   else
!    do itai=1,trout
!     ommil=(bomf(itai)+bomf(itai+1))/2.0_qp
!     call locate(bom(1:tres+1),ommil,p1)
!     if((p1==0).OR.(p1==tres+1))then
!       vres(itai)=0
!     else
!       vres(itai)=res(p1)
!     endif
!    enddo
!   endif


   if(bla0)then
    write(6,FMT="(A6,9G20.10)")"bomf=",bomf(1:trout+1)
    write(6,*)"routint=",ecritrout(trout,routint(1:trout))
    write(6,*)"vres="   ,vres(1:trout)
    write(6,*)
    write(6,*)"************************"
    write(6,*)
    write(6,*)"        decoupe         "
    write(6,*)
   endif
   intresq(is,:)=decoupevcq(intresom,bomf(1:trout+1),6,arg,routint(1:trout),EPSom,bla0)
   intresq(is,:)=intresq(is,:)*qs**2
   open(25,file="intq"//trim(prefixe)//trim(suffixe)//".dat",POSITION="APPEND")
    write(25,*)qs,real(intresq(is,:)),imag(intresq(is,1:3))
   close(25)
   deallocate(bomf)
   deallocate(routint)
   deallocate(vres)
   deallocate(arg)
  enddo
  END FUNCTION intresq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION intresom(om,arg,m) !Computes size(om) points of the function  to be integrate over om
   IMPLICIT NONE
   INTEGER,  INTENT(IN) :: m !m=6 here
   REAL(QP), INTENT(IN), DIMENSION(:)  ::  om,arg !arg(1) should be the value of q
   COMPLEX(QPC), DIMENSION(size(om),m)       ::  intresom

   COMPLEX(QPC) Gam(1:2,1:2),Matt(1:2,1:2),MatCat(1:2,1:2),det
   REAL(QP) q,rho(1:2,1:2),ome,enP,enM
   COMPLEX(QPC) IuP(1:3),IuM(1:3)
   INTEGER is,r

   q=arg(1) 
   r=floor(arg(2))!Number of resonance angles
   intresom(:,:)=0.0_qp
!   if(r==0) return

   if(bla00.AND.(size(om)==1))then
    write(6,*)"************************"
   endif

   do is=1,size(om)
    ome=om(is)

    enM= ome-zk
    enP= ome+zk
  
    IuM(:)=0.0_qp
    IuP(:)=0.0_qp

    IuM=conjg(Iuanaly(enM,k,q,xmin,xmax)) !enM a une (petite) partie imaginaire négative venant de -zk
    IuP=      Iuanaly(enP,k,q,xmin,xmax)
  
    if(interpolation)then
     call estmat_pairfield(ome,e,det,Matt,Gam)
    else
     call mat_pairfield(ome,e,det,Matt,Gam)
    endif
    
    MatCat(1,1)=(Matt(1,1)+Matt(2,2))/2.0_qp+Matt(1,2)
    MatCat(2,2)=(Matt(1,1)+Matt(2,2))/2.0_qp-Matt(1,2)
    MatCat(1,2)=(Matt(2,2)-Matt(1,1))/2.0_qp
  
    Gam(1,1)= MatCat(2,2)/det
    Gam(2,2)= MatCat(1,1)/det
    Gam(1,2)=-MatCat(1,2)/det
  
    rho=-imag(Gam)/PI
  
    intresom(is,1)=-rho(1,1)*IuM(1)
    intresom(is,2)=-rho(2,2)*IuM(2)
    intresom(is,3)=-rho(1,2)*IuM(3)
  
    intresom(is,4)= rho(2,2)*IuP(2)
    intresom(is,5)= rho(1,1)*IuP(1)
    intresom(is,6)=-rho(1,2)*IuP(3)

    if(bla00)then
     write(6,FMT="(A21,8G20.10)")"q,ome,real(intresom)=",q,ome,real(intresom(is,1)),imag(intresom(is,1))
    endif

  enddo
  if(bla00)then
   write(6,*)
  endif
  END FUNCTION intresom
END FUNCTION intres
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
SUBROUTINE bornesq(k,zkt,le,reg,tconf,config,bq)
USE recettes
REAL(QP), INTENT(IN) :: k,zkt
CHARACTER(len=2), INTENT(OUT) :: reg
INTEGER, INTENT(OUT) :: tconf
INTEGER, DIMENSION(1:7), INTENT(OUT) :: config
REAL(QP), INTENT(OUT) :: bq(1:8)
REAL(QP), INTENT(OUT) :: le(1:8)

REAL(QP) km,kp,qm,q1m,q2m,q3m,q4m,q1p,q2p,q1mbis,q3mbis,q2pbis,qc,qd
REAL(QP) l1,l2,l3,l4,l5,l6,l7,l8
REAL(QP) vecq(1:13)

l1=1.0e50_qp; l2=1.0e50_qp; l3=1.0e50_qp; l4=1.0e50_qp; 
l5=1.0e50_qp; l6=1.0e50_qp; l7=1.0e50_qp; l8=1.0e50_qp;

!lignes d’énergie
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
le=(/l1,l2,l3,l4,l5,l6,l7,l8/)

!reg et configuration
if(k<k0)then
 if(zkt<1)then
  reg="00"
  tconf=0
 elseif(min(l1,l4)>zkt)then
  reg="A0"
  tconf=6
  config(1:tconf)=(/0,alti,betti,gam,bet,al/)
 elseif((zkt>l1).AND.(min(l2,l4)>zkt))then
  reg="B0"
  tconf=6
  config(1:tconf)=(/deltti,alti,betti,gam,bet,al/)
 elseif((zkt>l4).AND.(min(l1,l2)>zkt))then
  reg="B1"
  tconf=7
  config(1:tconf)=(/0,alti,betti,epsiti,epsi,bet,al/)
 elseif((zkt>l2).AND.(l4>zkt))then
  reg="C0"
  tconf=6
  config(1:tconf)=(/deltti,epsiti,betti,gam,bet,al/)
 elseif((l2>zkt).AND.(zkt>max(l1,l4)))then
  reg="C1"
  tconf=7
  config(1:tconf)=(/deltti,alti,betti,epsiti,epsi,bet,al/)
 elseif((min(l1,l3)>zkt).AND.(zkt>l2))then
  reg="C2"
  tconf=7
  config(1:tconf)=(/0,alti,deltti,epsiti,epsi,bet,al/)
 elseif((zkt>l4).AND.(l5>zkt).AND.(k>k1))then
  reg="D0"
  tconf=7
  config(1:tconf)=(/deltti,epsiti,epsi,bet,gam,bet,al/)
 elseif((zkt>max(l2,l4)).AND.(l5>zkt).AND.(k>k3))then
  reg="D1"
  tconf=7
  config(1:tconf)=(/deltti,epsiti,betti,epsiti,epsi,bet,al/)
 elseif((zkt>max(l1,l2)).AND.(min(l3,l5)>zkt).AND.(k>k7))then
  reg="D2"
  tconf=7
  config(1:tconf)=(/deltti,alti,deltti,epsiti,epsi,bet,al/)
 elseif((l1>zkt).AND.(zkt> l3))then
  reg="D3"
  tconf=7
  config(1:tconf)=(/0,alti,deltti,epsiti,epsi,delt,al/)
 elseif((l3>zkt).AND.(zkt>l5))then
  reg="E0"
  tconf=5
  config(1:tconf)=(/deltti,epsiti,epsi,bet,al/)
 elseif((l5>zkt).AND.(zkt>max(l1,l3)))then
  reg="E1"
  tconf=7
  config(1:tconf)=(/deltti,alti,deltti,epsiti,epsi,delt,al/)
 elseif(zkt>max(l3,l5))then
  reg="F0"
  tconf=5
  config(1:tconf)=(/deltti,epsiti,epsi,delt,al/)
 else
  stop "Erreur dans region k<k0"
 endif
else
 if(l8>zkt)then
  reg="00"
  tconf=0
 elseif((min(l1,l5)>zkt).AND.(zkt>l6))then
  reg="A0"
  tconf=6
  config(1:tconf)=(/0,al,bet,gam,bet,al/)
 elseif((min(l2,l5)>zkt).AND.(zkt>l1))then
  reg="B0"
  tconf=6
  config(1:tconf)=(/delt,al,bet,gam,bet,al/)
 elseif((zkt>max(l6,l5)).AND.(l1>zkt))then
  reg="B1"
  tconf=4
  config(1:tconf)=(/0,al,bet,al/)
 elseif((zkt>l2).AND.(l5>zkt))then
  reg="D0"
  tconf=6
  config(1:tconf)=(/delt,epsi,bet,gam,bet,al/)
 elseif((l2>zkt).AND.(zkt>max(l1,l5)))then
  reg="D1"
  tconf=4
  config(1:tconf)=(/delt,al,bet,al/)
 elseif((zkt>max(l2,l5)).AND.(l3>zkt))then
  reg="E0"
  tconf=4
  config(1:tconf)=(/delt,epsi,bet,al/)
 elseif(zkt>l3)then
  reg="F0"
  tconf=4
  config(1:tconf)=(/delt,epsi,delt,al/)
 elseif((k<k11).AND.((zkt<l4).OR.((zkt>l5).AND.(l6>zkt))))then
  reg="G0"
  tconf=4
  config(1:tconf)=(/0,al,bet,gam/)
 elseif((k>k11).AND.(k12>k).AND.(l6>zkt))then
  reg="G0"
  tconf=4
  config(1:tconf)=(/0,al,bet,gam/)
 elseif((k>k12).AND.(l6>zkt).AND.(zkt>l7))then
  reg="G0"
  tconf=4
  config(1:tconf)=(/0,al,bet,gam/)
 elseif((k>k12).AND.(l7>zkt).AND.(zkt>l8))then
  reg="H0"
  tconf=2
  config(1:tconf)=(/0,al/)
 elseif((k<k11).AND.(min(l5,l6)>zkt).AND.(zkt>l4))then
  reg="J0"
  tconf=6
  config(1:tconf)=(/0,al,bet,gam,bet,gam/)
 else 
  stop "Erreur dans region k>k0"
 endif
endif

!bornesq
km=    1.0e50_qp; kp=    1.0e50_qp; q1m=   1.0e50_qp; q2m=   1.0e50_qp; q3m=   1.0e50_qp; q4m=   1.0e50_qp; 
q1mbis=1.0e50_qp; q3mbis=1.0e50_qp; 
q1p=   1.0e50_qp; q2p=   1.0e50_qp; 
qc=    1.0e50_qp; qm=    1.0e50_qp;

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
  if(zkt<l1) write(6,*)"top"
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

  sC=zkt+2-ec(q)-epsBCS(k+s*q)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE bornesom(k,zk,q,grecque,res,bom,tres)
USE recettes
REAL(QP), INTENT(IN) :: k,zk,q
INTEGER, INTENT(IN) :: grecque
INTEGER, INTENT(OUT) :: res(1:2),tres
REAL(QP), INTENT(OUT) :: bom(1:3)
REAL(QP) s

bom(:)=1.e100_qp

if(grecque<6)then 
 s=+1.0_qp
else
 s=-1.0_qp
endif

if((grecque==al).OR.(grecque==alti))then
 tres=1
 bom(1)=ec(q)
 res(1)=1
 bom(2)=zk-epsBCS(k-s*q)
elseif((grecque==bet).OR.(grecque==betti))then
 tres=2
 bom(1)=ec(q)
 res(1)=1
 bom(2)=zk-epsBCS(k-s*q)
 res(2)=2
 bom(3)=zk-1.0_qp
elseif(grecque==gam)then
 tres=1
 bom(1)=ec(q)
 res(1)=2
 bom(2)=zk-1.0_qp
elseif((grecque==delt).OR.(grecque==deltti))then
 tres=1
 bom(1)=zk-epsBCS(k+s*q)
 res(1)=1
 bom(2)=zk-epsBCS(k-s*q)
elseif((grecque==epsi).OR.(grecque==epsiti))then
 tres=2
 bom(1)=zk-epsBCS(k+s*q)
 res(1)=1
 bom(2)=zk-epsBCS(k-s*q)
 res(2)=2
 bom(3)=zk-1.0_qp
elseif(grecque==0)then
 tres=0
endif

END SUBROUTINE bornesom
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION ecritconfig(tconf,config)
INTEGER, INTENT(IN) :: tconf
INTEGER, DIMENSION(:), INTENT(IN) :: config
CHARACTER(len=90) :: ecritconfig

INTEGER itai
ecritconfig=trim(ecritc(config(1)))
do itai=2,tconf
 ecritconfig=trim(ecritconfig)//"  "//trim(ecritc(config(itai)))
enddo
END FUNCTION ecritconfig
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
END MODULE intldc
