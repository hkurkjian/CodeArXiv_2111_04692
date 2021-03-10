PROGRAM selfcons
USE modsim
USE vars
USE selftot
USE recettes, ONLY : rtsafe,mnewt,rtsafeq
IMPLICIT NONE
LOGICAL blaSC

REAL(QP)  mu,k,kmin,kmax,dk,zkdep(1:1)
REAL(QP)  zk
REAL(QP)  EPS(1:3)
CHARACTER(len=90) fichiers(1:5),suffixe,suffeintq !1:3 fichldc, 3 fichlec, 4 fichpol
INTEGER nk,ik
LOGICAL nvofich,res,newton

REAL(QP) contK(1:2), xik, e0, Uk, Vk, OS
REAL(QP) detZero, detCont, cont, borneinf, bornesup
COMPLEX(QPC) sigbidon(1:2,1:6)

open(10,file='selfcons.inp')
 read(10,*)mu
 read(10,*)kmin
 read(10,*)kmax
 read(10,*)nk
 read(10,*)zkdep(1)
 read(10,*)fichiers(1) !pour charger bestM/donnees
 read(10,*)fichiers(2) !pour charger bestM/donnees_sup
 read(10,*)fichiers(3) !pour charger bestM/donnees_sup2
 read(10,*)fichiers(4) !pour intldc/intpasres
 read(10,*)fichiers(5) !pour intpole
 read(10,*)EPS(1)      !intldc/EPSq
 read(10,*)EPS(2)      !intldc/EPSom
 read(10,*)EPS(3)      !intpole/EPSq
 read(10,*)res         !TRUE to use detGres, FALSE for detG
 read(10,*)newton      !TRUE to use mnewt,   FALSE for rtsafe
 read(10,*)suffixe     !terminaison of output files
 read(10,*)nvofich     !TRUE to overwrite output files
close(10)

if(nvofich)then
 open(20,file="DONNEES/autocor"//trim(suffixe)//".dat")
  write(20,*)"!Valeurs de k,z_s(k) et bornes des continua pour x0=",x0
 close(20)
endif

if(nk==0)then
 dk=0.0
else
 dk=(kmax-kmin)/nk
endif

blaSC=.TRUE.
suffeintq="bidon"

call initialisation(mu,0,fichiers,0)
do ik=0,nk
 k=kmin+dk*ik
 if(blaSC) write(6,*)"-------------------------------------"
 if(blaSC) write(6,*)
 if(blaSC) write(6,*)"k=",k
 if(blaSC) write(6,*)


 ! Calculate continuum thresholds
 contK=thresholds(mu,k) 

 ! Mean-field functions
 xik=k**2-mu
 e0=sqrt(xik**2+1.0_qp)
 Uk=sqrt((1.0_qp+xik/e0)/2.0_qp)
 Vk=sqrt((1.0_qp-xik/e0)/2.0_qp)
    
 ! Look for a solution below the continuum
 if(newton)then
  if(blaSC) write(6,*)"mnewt, initial guess: zkdep=",zkdep
  call mnewt(20,zkdep,0.001_qp*EPS(1),0.1_qp*EPS(1),detGnewt)
 else
  OS=1.e-6_qp
  cont=contK(2)
  if(blaSC) write(6,*)"rtsafe: trying to bracket the root"
  if(res)then
   detZero=detGres(k,     OS,EPS,sigbidon,suffeintq)
  else
   detZero=detG   (k,     OS,EPS,sigbidon,suffeintq)
  endif
  if(blaSC) write(6,*)"zk,det=",OS         ,detZero
 
  if(res)then
   detCont=detGres(k,cont-OS,EPS,sigbidon,suffeintq)
  else
   detCont=detG   (k,cont-OS,EPS,sigbidon,suffeintq)
  endif
  if(blaSC) write(6,*)"zk,det=",cont-OS,detCont
 
 
  if(detZero*detCont<0.0_qp)then
    if(blaSC) write(6,*)"Looking for solution below continuum"
    borneinf=zkdep(1)
    borneinf=OS
    bornesup=cont-OS
    zkdep(1)=rtsafeq(detGrtsafe,(/bidon/),borneinf,bornesup,1.e-6_qp)
  else
    if(blaSC) write(6,*)"Solution is inside the continuum"
  endif
 endif

 open(20,file="DONNEES/solautocor"//trim(suffixe)//".dat",POSITION="APPEND")
  write(20,*)k,zkdep,contK
 close(20)

enddo

call desinit
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE detGrtsafe(zkIn,arg,det,ddet)
    REAL(QP), INTENT(IN) :: zkIn
    REAL(QP), DIMENSION(:), INTENT(IN) :: arg
    REAL(QP), INTENT(OUT) :: det,ddet

    REAL(QP) detP,detM,h

    h=zkIn*0.1_qp*OS

    if(res)then
     det =real(detGres(k,zkIn  ,EPS,sigbidon,suffeintq))
    else
     det =real(detG   (k,zkIn  ,EPS,sigbidon,suffeintq))
    endif
    if(blaSC) write(6,*)"zk,det=",zkIn  ,det

    if(res)then
     detP=real(detGres(k,zkIn+h,EPS,sigbidon,suffeintq))
    else
     detP=real(detG   (k,zkIn+h,EPS,sigbidon,suffeintq))
    endif
    if(blaSC) write(6,*)"zk,det=",zkIn+h,detP

    if(res)then
     detM=real(detGres(k,zkIn-h,EPS,sigbidon,suffeintq))
    else
     detM=real(detG   (k,zkIn-h,EPS,sigbidon,suffeintq))
    endif
    if(blaSC) write(6,*)"zk,det=",zkIn-h,detM

    ddet=(detP-detM)/(2.0_qp*h)

    if(blaSC)then
      write(6,*)"ddet=",ddet
      write(6,*)
    endif

  END SUBROUTINE detGrtsafe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE detGnewt(zk,det,ddet)
    REAL(QP), INTENT(IN),  DIMENSION(:) :: zk
    REAL(QP), INTENT(OUT), DIMENSION(:)   :: det
    REAL(QP), INTENT(OUT), DIMENSION(:,:)   :: ddet

    REAL(QP) zkIn
    REAL(QP) detP,detM,h

    ! Self-energy
    zkIn=zk(1)
    h=zkIn*0.0001_qp

    if(res)then
     det(1)=real(detGres(k,zkIn  ,EPS,sigbidon,suffeintq))
    else
     det(1)=real(detG   (k,zkIn  ,EPS,sigbidon,suffeintq))
    endif
    if(blaSC) write(6,*)"zk,det=",zkIn  ,det

    if(res)then
     detP  =real(detGres(k,zkIn+h,EPS,sigbidon,suffeintq))
    else
     detP  =real(detG   (k,zkIn+h,EPS,sigbidon,suffeintq))
    endif
    if(blaSC) write(6,*)"zk,det=",zkIn+h,detP

    if(res)then
     detM  =real(detGres(k,zkIn-h,EPS,sigbidon,suffeintq))
    else
     detM  =real(detG   (k,zkIn-h,EPS,sigbidon,suffeintq))
    endif
    if(blaSC) write(6,*)"zk,det=",zkIn-h,detM

    ddet(1,1)=(detP-detM)/(2.0_qp*h)

    if(blaSC)then
      write(6,*)"ddet=",ddet
      write(6,*)
    endif

  END SUBROUTINE detGnewt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END PROGRAM selfcons
