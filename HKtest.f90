program test
USE modsim
USE recettes
USE dspec
USE vars
USE Zerom
USE intpole
USE intldc
USE bestM
USE eqdetat
USE selftot
USE OMP_LIB
IMPLICIT NONE

REAL(QP) om,dom,M(1:2,1:2),dM(1:2,1:2),A(1:6),dMm(1:3),q
COMPLEX(QPC) Mm2(1:6)
REAL(QP) Mmbid(1:6,1:3)
REAL(QP) bk(0:12),le(1:8)
REAL(QP) EPSq,EPSom,EPSpole,EPS(1:3)
REAL(QP) :: k,zk,testt,mu,th(1:2),rho(1:2,1:2)
COMPLEX(QPC) sigma(1:2,1:6),sigcomb(1:3)
CHARACTER(len=90) fichom2(1:2),fichlec,fichgri(1:3),suffintq,suffsE,fichiers(1:5),fichpole
COMPLEX(QPC) Gamm(1:2,1:2),Matt(1:2,1:2),MatCat(1:2,1:2),det

REAL(QP) ptq,ptom,ptM(1:3),ptdM(1:3)
INTEGER nivobla

LOGICAL interpol,testpt,testdspec,testselftot,lecture,ecriture

blaM=.TRUE.
fichgri(1)="DONNEES/BCS_4_sup2"
fichgri(2)="DONNEES/BCS_4_sup3"

!!call combineom2(fichgri(1),fichgri(2))
!!stop
!
!Mmbid(1,:)=(/5.0_qp, 1.0_qp,1.0_qp/)
!Mmbid(2,:)=(/13.0_qp,2.0_qp,7.0_qp/)
!Mmbid(3,:)=(/3.0_qp, 3.0_qp,3.0_qp/)
!Mmbid(4,:)=(/6.0_qp, 9.0_qp,4.0_qp/)
!Mmbid(5,:)=(/15.0_qp,5.0_qp,5.0_qp/)
!Mmbid(6,:)=(/2.0_qp, 6.0_qp,6.0_qp/)
!call tricol_q(Mmbid,4,2)
!write(6,*)
!write(6,*)
!write(6,*) Mmbid(1,:)
!write(6,*) Mmbid(2,:)
!write(6,*) Mmbid(3,:)
!write(6,*) Mmbid(4,:)
!write(6,*) Mmbid(5,:)
!write(6,*) Mmbid(6,:)
!stop
!Mm2=cmplx(Mmbid,0.0_qp,kind=qpc)
!write(6,*)Mm2

testpt=.TRUE.
testpt=.FALSE.
testdspec=.TRUE.
testdspec=.FALSE.
testselftot=.FALSE.
testselftot=.TRUE.

!Paramètres physiques
x0=4.0_qp
k=2.8_qp
zk=sqrt((k**2-x0)**2+1)-0.0001_qp
k=3.8_qp
zk=3.4_qp
zk=22.0_qp
k=6.14124999999999999999999999999999982_qp         
zk=300.00039999999999999999999999999999998_qp
!k 0 to 3
!zk 1 to 6

!Paramètres de dspec
temperaturenulle=.TRUE.
EPSpp=1.0e-7_qp
x0crit=0.0_qp
bla1=.TRUE.
bla1=.FALSE.
bla2=.TRUE.
bla2=.FALSE.

!Paramètres de estM
blaM=.TRUE.
blaM=.FALSE.
blaerr=.TRUE.
blaerr=.FALSE.
qpetit=0.1_qp/x0
qpetit=0.03_qp
!Fichiers de données
fichgri(1)="DONNEES/BCS_4_2"
fichgri(2)="DONNEES/BCS_4_sup_comb"
fichgri(3)="DONNEES/BCS_4_sup4"


!Paramètres de intldc
EPSom=1.0e-5_qp
EPSq =1.0e-3_qp
EPSpole=1.0e-5_qp
EPS=(/EPSq,EPSom,EPSpole/)
bla0=.TRUE.
bla00=.FALSE.
bla00=.TRUE.
lecture=.FALSE.
lecture=.TRUE.
ecriture=.TRUE.
ecriture=.FALSE.
bq=(/0.0_qp,10.0_qp,100.0_qp/)
profondeur=5
suffintq="nvobestM"
interpol=.TRUE.
fichlec ="DONNEES/grillebis_x04"
fichpole="DONNEES/BCS_4_pole"


if(testselftot)then
nivobla=2
ecrintq=1
fichiers=(/fichgri,fichlec,fichpole/)
call initialisation(x0,nivobla,fichiers(1:5),ecrintq)
endif


if(testpt)then
 blaM=.TRUE.
 blaerr=.TRUE.
 call load_data(fichgri(1))
 call loadom2(fichgri(2))
 call loadom3(fichgri(3))
 
 call calcxqjoin
 xq=4.82291968861206930767704624488931714         
 om=4.14927348663679825437208934202557343
 call oangpp
 write(6,*)"opp=",opp(1:3)
 
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
 
 call unload_data
 stop
endif

if(testdspec)then
 bla1=.TRUE.
 bla2=.TRUE.
 bla2=.FALSE.
 
 
 call calcxqjoin
 x0=100.0_qp
 x0=4.0_qp
 xq=3.15277857142857142857142857142857182_qp
 om=3.63054415685928261130364771224973059_qp
 om=760915.550879716444940200860623913995_qp
 xq=4.29779714285714285714285714285714286_qp
 om=9995.0594221204099389844000076800279_qp
 xq=4.29578914285714285714285714285714324_qp
 xq=4.166666667_qp
 om=93367.19425_qp
 om=80.9166935138540908220669656417116605_qp
 xq=13.3333333333333333333333333333333333_qp
 xq=13.3333333333333333333333333333333323_qp
 om=80.9166935138540908220669656417116605_qp
 xq=60.85759778         
 om=2840000.1080184
 xq=14.12789396_qp         
 om= 369.8180720_qp 
 call oangpp
 write(6,*)"opp=",opp(1:3)
 
 call mat_pairfield(om,0.0_qp,det,Matt,Gamm)
 MatCat(1,1)=(Matt(1,1)+Matt(2,2))/2.0_qp+Matt(1,2)
 MatCat(2,2)=(Matt(1,1)+Matt(2,2))/2.0_qp-Matt(1,2)
 MatCat(1,2)=(Matt(2,2)-Matt(1,1))/2.0_qp

 Gamm(1,1)= MatCat(2,2)/det
 Gamm(2,2)= MatCat(1,1)/det
 Gamm(1,2)=-MatCat(1,2)/det

 rho=-imag(Gamm)/PI

 write(6,*)"om,xq,real(Matt(1,1))=",om,xq,real(Matt(1,1))
 write(6,*)"om,xq,real(Matt(2,2))=",om,xq,real(Matt(2,2))
 write(6,*)"om,xq,real(Matt(1,2))=",om,xq,real(Matt(1,2))
 write(6,*)"om,xq,real(Matt(1,1))=",om,xq,imag(Matt(1,1))
 write(6,*)"om,xq,real(Matt(2,2))=",om,xq,imag(Matt(2,2))
 write(6,*)"om,xq,real(Matt(1,2))=",om,xq,imag(Matt(1,2))
 write(6,*)
 write(6,*)"om,xq,real(rho(1,1))=",om,xq,rho(1,1)
 write(6,*)"om,xq,real(rho(2,2))=",om,xq,rho(2,2)
 write(6,*)"om,xq,real(rho(1,2))=",om,xq,rho(1,2)
 write(6,*)
 write(6,*)"det=",det
 write(6,*)
 call mat_pairfield_gom0(om,0.0_qp,det,Matt,Gamm)
 MatCat(1,1)=(Matt(1,1)+Matt(2,2))/2.0_qp+Matt(1,2)
 MatCat(2,2)=(Matt(1,1)+Matt(2,2))/2.0_qp-Matt(1,2)
 MatCat(1,2)=(Matt(2,2)-Matt(1,1))/2.0_qp

 Gamm(1,1)= MatCat(2,2)/det
 Gamm(2,2)= MatCat(1,1)/det
 Gamm(1,2)=-MatCat(1,2)/det

 rho=-imag(Gamm)/PI

 write(6,*)"om,xq,real(Matt(1,1))=",om,xq,real(Matt(1,1))
 write(6,*)"om,xq,real(Matt(2,2))=",om,xq,real(Matt(2,2))
 write(6,*)"om,xq,real(Matt(1,2))=",om,xq,real(Matt(1,2))
 write(6,*)"om,xq,real(Matt(1,1))=",om,xq,imag(Matt(1,1))
 write(6,*)"om,xq,real(Matt(2,2))=",om,xq,imag(Matt(2,2))
 write(6,*)"om,xq,real(Matt(1,2))=",om,xq,imag(Matt(1,2))
 write(6,*)
 write(6,*)"om,xq,real(rho(1,1))=",om,xq,rho(1,1)
 write(6,*)"om,xq,real(rho(2,2))=",om,xq,rho(2,2)
 write(6,*)"om,xq,real(rho(1,2))=",om,xq,rho(1,2)
 write(6,*)
 write(6,*)"det=",det
 write(6,*)
 stop

endif

mu=x0
th=thresholds(mu,k)
write(6,*)"th=",th
if(zk<th(2))then
 det=detG   (k,zk,EPS,sigma,suffintq)
else
 det=detGres(k,zk,EPS,sigma,suffintq)
endif

suffsE="hktest"
call bornesk(bk)
write(6,*)"bk=",bk
call lignesenergie(k,le)
write(6,*)"le=",le
Mm2=intres(k,zk,interpol,EPS(1:2),bk,le,suffintq)
write(6,*)"Mm2=",Mm2
open(17,file="selfE"//trim(suffsE)//".dat")
 write(17,*)real(Mm2),imag(Mm2)
close(17)

!Mm2=intpasres(k,zk,lecture,ecriture,profondeur,EPS,bq,fichlec,suffintq)
open(17,file="selfE"//trim(suffsE)//".dat")
 write(17,*)real(Mm2),imag(Mm2)
close(17)

stop
end program
