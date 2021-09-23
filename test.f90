program test
USE modsim
USE recettes
USE dspec
USE vars
USE estM
USE eqdetat
USE OMP_LIB
IMPLICIT NONE

REAL(QP) om,dom,M(1:2,1:2),dM(1:2,1:2),A(1:6),dMm(1:3),q,ommax,ommin
COMPLEX(QPC) Mm2(1:6)
REAL(QP) Mmbid(1:6,1:3)
REAL(QP) bk(0:12),le(1:8)
REAL(QP) EPSq,EPSom,EPSpole,EPS(1:3)
REAL(QP) :: k,zk,testt,mu,th(1:2),rho(1:2,1:2)
COMPLEX(QPC) sigma(1:2,1:6),sigcomb(1:3)
CHARACTER(len=90) fich,fichiers(1:5),fichpole
COMPLEX(QPC) Gamm(1:2,1:2),Matt(1:2,1:2),MatCat(1:2,1:2),det

REAL(QP) ptq,ptom,ptM(1:3),ptdM(1:3)
INTEGER nivobla,itest

LOGICAL interpol,testest,testdspec,testselftot,lecture,ecriture

testest=.FALSE.
testest=.TRUE.

!Paramètres de dspec
temperaturenulle=.TRUE.
EPSpp=1.0e-8_qp
x0crit=1.0_qp
bla1=.TRUE.
bla1=.FALSE.
bla2=.TRUE.
bla2=.FALSE.

!Paramètres de estM
qpetit=0.1_qp/x0
qpetit=0.015_qp
omgrand=3000.0_qp
!Fichiers de données
fich="DONNEES/LU/BCS_LU"

x0=0.8604366861256786_qp

if(testest)then
 blaM=.TRUE.
 blaM=.FALSE.
 blaerr=.FALSE.
 call load_data(fich)
 
 call calcxqjoin
 xq=168.521030460055679625590380867127081_qp
 call oangpp
 write(6,*)"opp=",opp(1:3)
 ommin=28300
 ommax=28500
 dom=(ommax-ommin)/1000
 do itest=1,1000
  om=ommin+itest*dom
 
  call estmat_pairfield(om,0.0_qp,det,Matt,Gamm)
  write(6,*)"om,xq,real(Matt(1,1))=",om,xq,real(Matt(1,1))
  write(6,*)"om,xq,real(Matt(2,2))=",om,xq,real(Matt(2,2))
  write(6,*)"om,xq,real(Matt(1,2))=",om,xq,real(Matt(1,2))
  write(6,*)"om,xq,real(Matt(1,1))=",om,xq,imag(Matt(1,1))
  write(6,*)"om,xq,real(Matt(2,2))=",om,xq,imag(Matt(2,2))
  write(6,*)"om,xq,real(Matt(1,2))=",om,xq,imag(Matt(1,2))
  write(6,*)
  write(6,*)"det=",det
  write(6,*)
  
  call mat_pairfield(om,0.0_qp,det,Matt,Gamm)
  write(6,*)"om,xq,real(Matt(1,1))=",om,xq,real(Matt(1,1))
  write(6,*)"om,xq,real(Matt(2,2))=",om,xq,real(Matt(2,2))
  write(6,*)"om,xq,real(Matt(1,2))=",om,xq,real(Matt(1,2))
  write(6,*)"om,xq,real(Matt(1,1))=",om,xq,imag(Matt(1,1))
  write(6,*)"om,xq,real(Matt(2,2))=",om,xq,imag(Matt(2,2))
  write(6,*)"om,xq,real(Matt(1,2))=",om,xq,imag(Matt(1,2))
  write(6,*)
  write(6,*)"det=",det
  write(6,*)

 enddo
 
 call unload_data
 stop
endif

end program
