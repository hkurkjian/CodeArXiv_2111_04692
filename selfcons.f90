PROGRAM selfcons
USE modsim
USE vars
USE selftot
USE recettes, ONLY : rtsafe,mnewt
IMPLICIT NONE
LOGICAL blaSC

REAL(QP)  mu,k,kmin,kmax,dk,zkdep
REAL(QP)  zk
REAL(QP)  EPS(1:3)
CHARACTER(len=90) fichiers(1:6),suffixe !1:2 fichldc, 3 fichlec, 4:5 fichom2, 6 fichpol
INTEGER nk,ik
LOGICAL nvofich

REAL(QP) contK(1:2), xik, e0, Uk, Vk, OS
REAL(QP) detZero, detCont
COMPLEX(QPC) sigbidon(1:2,1:6)

open(10,file='selfcons.inp')
 read(10,*)mu
 read(10,*)kmin
 read(10,*)kmax
 read(10,*)nk
 read(10,*)zkdep
 read(10,*)fichiers(1) !pour charger bestM/donnees
 read(10,*)fichiers(2) !pour charger bestM/donnees_sup
 read(10,*)fichiers(3) !pour intldc/intpasres
 read(10,*)fichiers(4) !pour intldc/lignesenergie
 read(10,*)fichiers(5) !pour intldc/lignesenergie
 read(10,*)fichiers(6) !pour intpole
 read(10,*)EPS(1)      !intldc/EPSq
 read(10,*)EPS(2)      !intldc/EPSom
 read(10,*)EPS(3)      !intpole/EPSq
 read(10,*)suffixe !terminaison of output files
 read(10,*)nvofich !TRUE to overwrite output files
close(10)

if(nvofich)then
 open(20,file="DONNEES/solautocor"//trim(suffixe)//".dat")
  write(20,*)"!Valeurs de k,z_s(k) et bornes des continua pour x0=",x0
 close(20)
endif

if(nk==0)then
 dk=0.0
else
 dk=(kmax-kmin)/nk
endif

blaSC=.TRUE.

do ik=0,nk
 k=kmin+dk*ik
 if(blaSC) write(6,*)"k=",k

 ! Calculate continuum thresholds
 contK=thresholds(mu,k,(/fichiers(4:6)/)) 

 ! Mean-field functions
 xik=k**2-mu
 e0=sqrt(xik**2+1.0_qp)
 Uk=sqrt((1.0_qp+xik/e0)/2.0_qp)
 Vk=sqrt((1.0_qp-xik/e0)/2.0_qp)
    
 call initialisation(mu,0,.FALSE.,(/fichiers(1:2)/))
 ! Look for a solution below the continuum
 OS=1.e-6_qp
 detZero=detG(k,         OS,(/fichiers(3),fichiers(6)/),EPS,sigbidon)
 if(blaSC) write(6,*)"zk,det=",OS         ,detZero
 detCont=detG(k,contK(1)-OS,(/fichiers(3),fichiers(6)/),EPS,sigbidon)
 if(blaSC) write(6,*)"zk,det=",contK(1)-OS,detCont


 if(detZero*detCont<0.0_qp)then
   if(blaSC) write(6,*)"Looking for solution below continuum"
   zk=rtsafe(rootFunSC,(/bidon/),OS,MINVAL(contK)-OS,1.e-5_qp)
 else
   if(blaSC) write(6,*)"Solution is inside the continuum"
 endif

 open(20,file="DONNEES/solautocor"//trim(suffixe)//".dat",POSITION="APPEND")
  write(20,*)k,zk,contK
 close(20)

enddo

CONTAINS
  SUBROUTINE rootFunSC(zkIn,arg,det,ddet)
    IMPLICIT NONE
    REAL(QP), INTENT(IN) :: zkIn
    REAL(QP), DIMENSION(:), INTENT(IN) :: arg
    REAL(QP), INTENT(OUT) :: det,ddet

    REAL(QP) detP,detM,h

    ! Self-energy
    h=zkIn*0.0001_qp
    det =real(detG(k,zkIn  ,(/fichiers(3),fichiers(6)/),EPS,sigbidon))
    if(blaSC) write(6,*)"zk,det=",zkIn  ,det
    detP=real(detG(k,zkIn+h,(/fichiers(3),fichiers(6)/),EPS,sigbidon))
    if(blaSC) write(6,*)"zk,det=",zkIn+h,det
    detM=real(detG(k,zkIn-h,(/fichiers(3),fichiers(6)/),EPS,sigbidon))
    if(blaSC) write(6,*)"zk,det=",zkIn-h,det

    ddet=(detP-detM)/(2.0_qp*h)

    if(blaSC)then
      write(6,*)"zk=",zkIn
      write(6,*)"det,ddet=",det,ddet
    endif

  END SUBROUTINE rootFunSC
END PROGRAM selfcons
