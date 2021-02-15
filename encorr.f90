PROGRAM encorr
USE vars
USE recettes
USE selftot
IMPLICIT NONE

REAL(QP) :: mu,k,zk
COMPLEX(QPC) sigma(1:2,1:6),sigcomb(1:3),det
REAL(QP) dk,dzk
REAL(QP) kmin,kmax,zkmin,zkmax,EPS(1:3)
REAL(QP) bqbidon(1:3),th(1:2)
INTEGER nk,nzk,ik,izk,nivobla,profondeurbidon
CHARACTER(len=90) fichldc(1:2),fichiers(1:6),suffixe,suffintq
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
 read(10,*)fichiers(1) !pour charger bestM/donnees
 read(10,*)fichiers(2) !pour charger bestM/donnees_sup
 read(10,*)fichiers(3) !pour intldc/intpasres
 read(10,*)fichiers(4) !pour intldc/lignesenergie
 read(10,*)fichiers(5) !pour intldc/lignesenergie
 read(10,*)fichiers(6) !pour intpole
 read(10,*)EPS(1)      !intldc/EPSq
 read(10,*)EPS(2)      !intldc/EPSom
 read(10,*)EPS(3)      !intpole/EPSq
 read(10,*)nivobla !Verbose level, from 0 to 3
 read(10,*)suffixe !terminaison of output files
 read(10,*)nvofich !TRUE to overwrite output files
 read(10,*)ecrintq !0 to avoid writing intq files, 1 to write, 2 to overwrite
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
write(6,*)'fichiers ldc:',trim(fichiers(1))," ",trim(fichiers(2))
write(6,*)'fichiers om2:',trim(fichiers(4))," ",trim(fichiers(5))
write(6,*)'fichier  lec:',trim(fichiers(3))
write(6,*)'fichier pole:',trim(fichiers(6))
write(6,*)'precisions:',EPS

call initialisation(mu,nivobla,.TRUE.,fichiers(1:2))

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

  th=thresholds(mu,k,(/fichiers(4:6)/))
  if(zk<th(2))then
   det=detG(k,zk,(/fichiers(3),fichiers(6)/),EPS,sigma)
!   det=detGres(k,zk,(/fichiers(3),fichiers(6)/),EPS,sigcomb)
  else
   det=detGres(k,zk,(/fichiers(4:6)/),EPS,sigma)
  endif
  
  sigcomb(1:3)=sigma(1,1:3)+sigma(1,4:6)+sigma(2,1:3)+sigma(2,4:6)

  open(20,file="DONNEES/selfEldc"//trim(suffixe)//".dat",POSITION="APPEND")
  open(21,file="DONNEES/selfEpol"//trim(suffixe)//".dat",POSITION="APPEND")
  open(22,file="DONNEES/selfEtot"//trim(suffixe)//".dat",POSITION="APPEND")
  open(23,file="DONNEES/continuum"//trim(suffixe)//".dat",POSITION="APPEND")
   write(6,*)"k,zk,bords des continua=",k,zk,th
   write(6,*)"re(selfEldc)=",real(sigma(1,1:6))
   write(6,*)"im(selfEldc)=",imag(sigma(1,1:6))
   write(6,*)"re(selfEpol)=",real(sigma(2,1:6))
   write(6,*)"im(selfEpol)=",imag(sigma(2,1:6))
   write(6,*)"re(selfEtot)=",real(sigcomb(1:3))
   write(6,*)"im(selfEtot)=",imag(sigcomb(1:3))
   write(20,*)k,zk,real(sigma(1,1:6))
   write(20,*)k,zk,imag(sigma(1,1:6))
   write(21,*)k,zk,real(sigma(2,1:6))
   write(21,*)k,zk,imag(sigma(2,1:6))
   write(22,*)k,zk,real(sigcomb(1:3))
   write(22,*)k,zk,imag(sigcomb(1:3))
   write(23,*)k,th
  close(20)
  close(21)
  close(22)
  close(23)

 enddo
enddo

call desinit
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
