PROGRAM encorr
USE vars
USE dspec
USE bestM
USE intldc
USE intpole
USE recettes
IMPLICIT NONE

REAL(QP) :: mu,k,zk
COMPLEX(QPC) selfEtot(1:6),selfEldc(1:6),selfEpol(1:6)
REAL(QP) dk,dzk
REAL(QP) kmin,kmax,zkmin,zkmax,EPS(1:3)
REAL(QP) bk(0:12),le(1:8),bk2(0:12),le2(1:8),bqbidon(1:3)
INTEGER nk,nzk,ik,izk,nivobla,profondeurbidon
CHARACTER(len=90) fichldc(1:2),fichlec,fichom2(1:2),fichierpole,suffixe,suffintq
CHARACTER(len=5) cik,cizk
LOGICAL nvofich

!k 0 to 3
!zk 1 to 6
open(10,file='encorr.inp')
 read(10,*)mu
 read(10,*)kmin
 read(10,*)kmax
 read(10,*)nk
 read(10,*)zkmin
 read(10,*)zkmax
 read(10,*)nzk
 read(10,*)fichldc(1)
 read(10,*)fichldc(2)
 read(10,*)fichlec
 read(10,*)fichom2(1)
 read(10,*)fichom2(2)
 read(10,*)fichierpole
 read(10,*)EPS(1)!EPSq  pour intldc
 read(10,*)EPS(2)!EPSom pour intldc
 read(10,*)EPS(3)!EPSq  pour intpole
 read(10,*)nivobla
 read(10,*)suffixe
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
write(6,*)'fichiers:',trim(fichldc(1))," ",trim(fichldc(2))
write(6,*)'fichiers:',trim(fichom2(1))," ",trim(fichom2(2))
write(6,*)'fichiers:',trim(fichlec)," ",trim(fichierpole)
write(6,*)'precision:',EPS

!Niveaux de blabla (blaPole: intpole, bla0 et bla00: intldc, blaM et blaerr: estM)
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

!Initialisation de bestM
qpetit=0.03_qp
call load_data(fichldc(1))
call loadom2  (fichldc(2))

!Initialisation de dspec
x0=mu
temperaturenulle=.TRUE.
EPSpp=1.0e-7_qp
x0crit=0.0_qp

!Initialisation de intpole
fichpol=fichierpole

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

if(nvofich)then
 open(20,file="DONNEES/selfEldc"//trim(suffixe)//".dat")
  write(20,*)"!Valeurs de k,zk et selfEldc pour x0=",x0
 close(20)
 open(21,file="DONNEES/selfEpol"//trim(suffixe)//".dat")
  write(21,*)"!Valeurs de k,zk et selfEpol pour x0=",x0
 close(21)
 open(22,file="DONNEES/selfEtot"//trim(suffixe)//".dat")
  write(22,*)"!Valeurs de k,zk et selfEtot pour x0=",x0
 close(22)
endif

do ik=0,nk
 k=kmin+dk*ik
 write(6,*)"k=",k
 do izk=0,nzk
  zk=zkmin+dzk*izk
  write(6,*)"zk=",zk
  write(cik, FMT="(I2)")ik
  write(cizk,FMT="(I2)")izk
  cik=adjustl(cik)
  cizk=adjustl(cizk)
  write(suffintq,*)"_",trim(cik),"_",trim(cizk)
  suffintq=adjustl(suffintq)
  write(6,*)"suffintq:",trim(suffintq)

  if(nivobla>0) write(6,*)"-----------------------------------"
  if(nivobla>0) write(6,*)
  if(nivobla>0) write(6,*)"        Calcul de selfEldc"
  if(nivobla>0) write(6,*)
  
  call bornesk(bk)
  call lignesenergie(k,fichom2,le)
  bk2=bk; le2=le
  call tri(bk2)
  call tri(le2)
  if(nivobla>0) write(6,FMT="(A3,12G20.10)")"bk=",bk2
  if(nivobla>0) write(6,FMT="(A3,8G20.10)")"le=",le2
  
  profondeurbidon=intbidon
  bqbidon(:)=bidon
  if((zk-2.0_qp)<MINVAL(le))then
   selfEldc=intres   (k,zk,.TRUE.,EPS(1:2),bk,le,suffintq)
!   selfEldc   =cmplx(intpasres(k,zk,.TRUE.,.FALSE.,profondeurbidon,EPS(1:2),bqbidon,fichlec,suffintq),0.0_qp,kind=qpc)
  else
   selfEldc=intres   (k,zk,.TRUE.,EPS(1:2),bk,le,suffintq)
  endif
  
  open(20,file="DONNEES/selfEldc"//trim(suffixe)//".dat",POSITION="APPEND")
   write(20,*)k,zk,real(selfEldc),imag(selfEldc)
  close(20)

  if(nivobla>0) write(6,*)"-----------------------------------"
  if(nivobla>0) write(6,*)
  if(nivobla>0) write(6,*)"        Calcul de selfEpole"
  if(nivobla>0) write(6,*)
  
  selfEpol=selfEpole(k,zk,EPS(3))

  open(21,file="DONNEES/selfEpol"//trim(suffixe)//".dat",POSITION="APPEND")
   write(21,*)k,zk,real(selfEpol),imag(selfEpol)
  close(21)

  selfEtot=selfEpol+selfEldc

  open(22,file="DONNEES/selfEtot"//trim(suffixe)//".dat",POSITION="APPEND")
   write(22,*)k,zk,real(selfEtot),imag(selfEtot)
  close(22)

 enddo
enddo

call unload_data
END PROGRAM encorr

!CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!FUNCTION selfEtot(mu,k,zk, !paramètres physique
!                   fichom2,fichldc,suffintq   !pour intldc/intres
!                   fichlec,   !pour intldc/intpasres
!                   fichierpol !pour intpol
!                   EPS,nivobla)
!COMPLEX(QPC) selfEtot(1:6)
!LOGICAL, INTENT(IN) :: interpol,nivobla
!CHARACTER(len=*), INTENT(IN) fichierpol,fichlec,fichldc(1:2),fichom2(1:2)
!REAL(QP, INTENT(IN)  EPS(1:3)
!REAL(QP) bk(0:12),le(1:8)
!
!!Paramètres de dspec
!x0=mu
!temperaturenulle=.TRUE.
!EPSpp=1.0e-7_qp
!x0crit=0.0_qp
!bla1=.TRUE.
!bla1=.FALSE.
!bla2=.TRUE.
!bla2=.FALSE.
!
!!pour intpol
!fichpol=fichierpol
!EPSpole=1.0e-6_qp
!
!!Paramètres de estM
!qpetit=0.01_qp
!
!if(nivobla==0)then
! bla0 =.FALSE.
! bla00=.FALSE.
! blaM=.FALSE.
! blaerr=.FALSE.
! blaPole=.FALSE.
!elseif(nivobla==1)then
! bla0 =.TRUE.
! bla00=.FALSE.
! blaM=.FALSE.
! blaerr=.FALSE.
! blaPole=.TRUE.
!elseif(nivobla==2)then
! bla0 =.TRUE.
! bla00=.FALSE.
! blaM=.FALSE.
! blaerr=.TRUE.
! blaPole=.TRUE.
!elseif(nivobla==3)then
! bla0 =.TRUE.
! bla00=.TRUE.
! blaM=.TRUE.
! blaerr=.TRUE.
! blaPole=.TRUE.
!endif
!
!selfEtot(:)=0.0_qp
!
!if(nivobla>0) write(6,*)"-----------------------------------"
!if(nivobla>0) write(6,*)
!if(nivobla>0) write(6,*)"        Calcul de selfEldc"
!if(nivobla>0) write(6,*)
!
!call bornesk(bk)
!if(nivobla>0) write(6,*)"bk=",bk
!call lignesenergie(k,fichom2,le)
!if(nivobla>0) write(6,*)"le=",le
!
!if((zk-2.0_qp)<min(le))then
! selfEtot%re   =intpasres(k,zk,.TRUE.,.FALSE.,intbidon,EPS(1:2),(/bidon,bidon,bidon/),fichlec,suffintq)
! selfEtot%im(:)=0.0_qp
!else
! selfEtot=intres   (k,zk,.TRUE.,EPS(3),fichgri,bk,le,suffintq)
!endif
!
!if(nivobla>0) write(6,*)"-----------------------------------"
!if(nivobla>0) write(6,*)
!if(nivobla>0) write(6,*)"        Calcul de selfEpole"
!if(nivobla>0) write(6,*)
!
!selfEtot=selfEtot+selfEpole(k,zk,EPSpole)
!
!END FUNCTION selfEtot
!END MODULE selfE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
