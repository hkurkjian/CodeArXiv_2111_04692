MODULE selfE
FUNCTION selfEtotk(mu,k,zk,
                   interpolation,EPSldc,fichom2,fichldc,
                   fichpol,EPSpole,
                   nivobla)
USE dspec
USE intldc
USE intpole

!Paramètres de dspec
temperaturenulle=.TRUE.
EPSpp=1.0e-8_qp
x0crit=0.0_qp
bla1=.TRUE.
bla1=.FALSE.
bla2=.TRUE.
bla2=.FALSE.

if(nivobla==0)then
 bla0 =.FALSE.
 bla00=.FALSE.
 blaM=.FALSE.
 blaerr=.FALSE.
 blaPole=.FALSE.
elseif(nivobla==1)then
 bla0 =.TRUE.
 bla00=.FALSE.
 blaM=.FALSE.
 blaerr=.FALSE.
 blaPole=.TRUE.
elseif(nivobla==2)then
 bla0 =.TRUE.
 bla00=.FALSE.
 blaM=.FALSE.
 blaerr=.TRUE.
 blaPole=.TRUE.
elseif(nivobla==3)then
 bla0 =.TRUE.
 bla00=.TRUE.
 blaM=.TRUE.
 blaerr=.TRUE.
 blaPole=.TRUE.
endif


if(nivobla>0) write(6,*)"-----------------------------------"
if(nivobla>0) write(6,*)
if(nivobla>0) write(6,*)"        Calcul de selfEldc"
if(nivobla>0) write(6,*)
call bornesk(bk)
if(nivobla>0) write(6,*)"bk=",bk
call lignesenergie(k,fichom2,le)
if(nivobla>0) write(6,*)"le=",le


END FUNCTION selfEtot
END MODULE selfE

PROGRAM encorr
USE dspec
IMPLICIT NONE

REAL(QP) :: k,zk
COMPLEX(QPC) sE(1:6)
REAL(QP) dk,dzk
REAL(QP) kmin,kmax,zkmin,zkmax
INTEGER nk,nzk,ik,izk
CHARACTER(len=30) grille
LOGICAL nvofich

!k 0 to 3
!zk 1 to 6
open(10,file='encorr.inp')
 read(10,*)kmin
 read(10,*)kmax
 read(10,*)nk
 read(10,*)zkmin
 read(10,*)zkmax
 read(10,*)nzk
 read(10,*)suffixe
 read(10,*)bla0
 read(10,*)bla00
 read(10,*)grille
 read(10,*)nvofich
close(10)


write(6,*)'--------------------'
write(6,*)
write(6,*)'Programme encorr'
write(6,*)
write(6,*)'suffixe=',suffixe
write(6,*)'kmin=',kmin
write(6,*)'kmax=',kmax
write(6,*)'nk=',nk
write(6,*)'zkmin=',zkmin
write(6,*)'zkmax=',zkmax
write(6,*)'nzk=',nzk



if(nk==0)then
 dk=0.0
else
 dk=(kmax-kmin)/nk
endif

if(nzk==0)then
 dzk=0.0
else
 dzk=(zkmax-zkmin)/nzk
endif

!Paramètres de estM
blaM=.TRUE.
blaM=.FALSE.
blaerr=.TRUE.
blaerr=.FALSE.
qpetit=0.1_qp/x0
qpetit=0.03_qp
!Fichiers de données
fichgri(1)="BCS_4_2"
fichgri(2)="BCS_4_sup"


!Paramètres de intpole
fichier="BCS_4_pole"

!Paramètres de intldc
EPSom=1.0e-5_qp
EPSq =1.0e-3_qp
EPSpole =1.0e-6_qp
EPS=(/EPSq,EPSom/)
bla0=.TRUE.
bla00=.FALSE.
!Paramètres de intldc/intpasres
lecture=.TRUE.
ecriture=.FALSE.
!Paramètres de intldc/intres
interpol=.TRUE.
fichom2(1) ="DONNEES/Tom1.dat"
fichom2(2) ="DONNEES/Tom1p.dat"

write(6,*)trim(grille)
open(15,file=trim(grille))
 read(15,*)x0,bq(1),bq(2),bq(3),fichlec(1),fichlec(2),profondeur
 write(6,*)
 write(6,*)"Pour les points non résonnants"
 write(6,*)
 write(6,*)"x0=",x0
 write(6,*)"bq(1),bq(2),bq(3)=",bq(1),bq(2),bq(3)
 write(6,*)"fichiers=",fichlec(1),fichlec(2)
 write(6,*)"profondeur=",profondeur
 write(6,*)
close(15)

if(nvofich)then
 open(14,file="selfE"//trim(suffixe)//".dat")
  write(14,*)"!Valeurs de k,zk et sE (l’autoénergie) pour x0=",x0
 close(14)
endif

do ik=0,nk
 k=kmin+dk*ik
 write(6,*)"k=",k
 call bornesk(bk)
 write(6,*)"bk=",bk
 do izk=0,nzk
  zk=zkmin+dzk*izk
  call lignesenergie(k,zk-2.0_qp,fichom2,le)
  write(prefixe,FMT="(A1,I2,A1,I2)")"_",ik,"_",izk
  write(6,*)trim(prefixe)
  write(6,*)"zk=",zk
  write(6,FMT="(A3,8G20.10)")"le=",le
  if((zk-2.0_qp)<min(le))then
   sE=intpasres(k,zk,lecture,ecriture,profondeur,EPS,bq,fichlec,suffintq)
  else
   sE=intres   (k,zk,interpol,EPS,fichgri,bk,le,prefixe)
  endif
 
  sEpole=selfEpole(k,zk,EPSpole)

  open(14,file="selfE"//trim(suffixe)//".dat",POSITION="APPEND")
   write(14,*)k,zk,real(sE),imag(sE)
  close(14)

 enddo
enddo


END PROGRAM encorr
