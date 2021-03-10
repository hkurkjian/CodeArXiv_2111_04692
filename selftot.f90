MODULE selftot
USE vars
USE dspec
USE intpole
USE intldc
USE bestM
IMPLICIT NONE
CONTAINS
SUBROUTINE initialisation(mu,nivobla,fichiers,eintq)
INTEGER, INTENT(IN) :: nivobla,eintq
CHARACTER(len=90), INTENT(IN) :: fichiers(1:5)
REAL(QP), INTENT(IN) :: mu

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
 bla00=.TRUE.
 blaM=.FALSE.
 blaerr=.FALSE.
 blaPole=.TRUE.
elseif(nivobla==3)then
 bla0 =.TRUE.
 bla00=.TRUE.
 blaM=.FALSE.
 blaerr=.TRUE.
 blaPole=.TRUE.
elseif(nivobla==4)then
 bla0 =.TRUE.
 bla00=.TRUE.
 blaM=.TRUE.
 blaerr=.TRUE.
 blaPole=.TRUE.
endif

!Initialisation de dspec
x0=mu
temperaturenulle=.TRUE.
EPSpp=1.0e-6_qp
x0crit=0.0_qp
!Initialisation de intpole
if(trim(fichiers(5)).NE."bidon")then
  call rdInfo(fichiers(5))
endif
!Initialisation de intldc
ecrintq=eintq
if(trim(fichiers(4)).NE."bidon")then
 call ini_intpasres(.TRUE.,.FALSE.,fichiers(4))
endif
!Initialisation de bestM
qpetit=0.03_qp

!chargement de donnees et donnees_sup ?
if(fichiers(1).NE."bidon")then
     call load_data(fichiers(1))
     call loadom2  (fichiers(2))
     call loadom3  (fichiers(3))
endif
END SUBROUTINE initialisation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE desinit
call unload_data
END SUBROUTINE desinit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION detG(k,zk,EPS,sigcomb,suff)
  REAL(QP), INTENT(IN) :: k,zk, EPS(1:3)
  CHARACTER(len=90), INTENT(IN) :: suff
  COMPLEX(QPC), INTENT(OUT) :: sigcomb(1:2,1:6)
  COMPLEX(QPC) detG

  COMPLEX(QPC) sig(1:3)
  REAL(QP) xik

  xik=k**2-x0

  sigcomb(1,:)=cmplx(intpasres(k,zk,.TRUE.,.FALSE.,EPS(1:2),suff),0.0_qp)
  sigcomb(2,:)=0.0_qp
  sigcomb(2,:)=selfEpole(k,zk,EPS(3))
!  sigcomb(1,:)=0.0_qp

  sig=sigcomb(1,1:3)+sigcomb(1,4:6)+sigcomb(2,1:3)+sigcomb(2,4:6)
  detG=(-zk+xik-sig(1))*(-zk-xik-sig(2))-(1.0_qp-sig(3))**2
END FUNCTION detG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION detGres(k,zk,EPS,sigcomb,suff)
  REAL(QP), INTENT(IN) :: k,zk, EPS(1:3)
  CHARACTER(len=90), INTENT(IN) :: suff
  COMPLEX(QPC), INTENT(OUT) :: sigcomb(1:2,1:6)
  COMPLEX(QPC) detGres

  REAL(QP) bk(0:12),le(1:8)
  COMPLEX(QPC) sig(1:3)
  REAL(QP) xik

  xik=k**2-x0

  call bornesk(bk)
  call lignesenergie(k,le)

  sigcomb(2,:)=0.0_qp
  sigcomb(2,:)=selfEpole(k,zk,EPS(3))
  sigcomb(1,:)=intres(k,zk,.TRUE.,EPS(1:2),bk,le,suff)

  sig=sigcomb(1,1:3)+sigcomb(1,4:6)+sigcomb(2,1:3)+sigcomb(2,4:6)
  detGres=(-zk+xik-sig(1))*(-zk-xik-sig(2))-(1.0_qp-sig(3))**2
END FUNCTION detGres
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION thresholds(mu,k)
  REAL(QP), INTENT(IN) :: mu,k
  REAL(QP) :: thresholds(1:2)
  thresholds(1)=contPole(k)
  if(k<3*sqrt(x0))then
   thresholds(2)=3.0_qp
  else
   thresholds(2)=3*epsBCS(k/3.0_qp)
  endif

END FUNCTION thresholds
END MODULE selftot
