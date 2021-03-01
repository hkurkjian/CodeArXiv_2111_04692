MODULE selftot
USE vars
USE dspec
USE intpole
USE intldc
USE bestM
IMPLICIT NONE
LOGICAL chargbestM

CONTAINS
SUBROUTINE initialisation(mu,nivobla,chargbestM,fichldc)
INTEGER, INTENT(IN) :: nivobla
LOGICAL, INTENT(IN) :: chargbestM
CHARACTER(len=90), INTENT(IN) :: fichldc(1:2)
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
EPSpp=1.0e-7_qp
x0crit=0.0_qp
!Initialisation de intldc
ecrintq=0
!Initialisation de bestM
qpetit=0.03_qp

!chargement de donnees et donnees_sup ?
if(chargbestM)then
     call load_data(fichldc(1))
     call loadom2  (fichldc(2))
endif
END SUBROUTINE initialisation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE desinit
call unload_data
END SUBROUTINE desinit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION detG(k,zk,fichiers,EPS,sigcomb)
  REAL(QP), INTENT(IN) :: k,zk, EPS(1:3)
  CHARACTER(len=90), INTENT(IN) :: fichiers(1:2)
  COMPLEX(QPC), INTENT(OUT) :: sigcomb(1:2,1:6)
  COMPLEX(QPC) detG

  COMPLEX(QPC) sig(1:3)
  REAL(QP) bqbidon(1:3)
  REAL(QP) xik
  INTEGER profondeurbidon

  xik=k**2-x0

  profondeurbidon=intbidon
  bqbidon(:)=bidon

  sigcomb(1,:)=cmplx(intpasres(k,zk,.TRUE.,.FALSE.,profondeurbidon,EPS(1:2),bqbidon,fichiers(1),""),0.0_qp)
  sigcomb(2,:)=0.0_qp
!  sigcomb(2,:)=selfEpole(k,zk,EPS(3),fichiers(2))
!  sigcomb(1,:)=0.0_qp

  sig=sigcomb(1,1:3)+sigcomb(1,4:6)+sigcomb(2,1:3)+sigcomb(2,4:6)
  detG=(-zk+xik-sig(1))*(-zk-xik-sig(2))-(1.0_qp-sig(3))**2
END FUNCTION detG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION detGres(k,zk,fichpol,EPS,sigcomb)
  REAL(QP), INTENT(IN) :: k,zk, EPS(1:3)
  CHARACTER(len=90), INTENT(IN) :: fichpol
  COMPLEX(QPC), INTENT(OUT) :: sigcomb(1:2,1:6)
  COMPLEX(QPC) detGres

  REAL(QP) bk(0:12),le(1:8)
  COMPLEX(QPC) sig(1:3)
  REAL(QP) xik

  xik=k**2-x0

  call bornesk(bk)
  call lignesenergie(k,le)

  sigcomb(2,:)=0.0_qp
!  sigcomb(2,:)=selfEpole(k,zk,EPS(3),fichpol)
  sigcomb(1,:)=intres(k,zk,.TRUE.,EPS(1:2),bk,le,"")

  sig=sigcomb(1,1:3)+sigcomb(1,4:6)+sigcomb(2,1:3)+sigcomb(2,4:6)
  detGres=(-zk+xik-sig(1))*(-zk-xik-sig(2))-(1.0_qp-sig(3))**2
END FUNCTION detGres
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION thresholds(mu,k,fichpol)
  REAL(QP), INTENT(IN) :: mu,k
  CHARACTER(len=90), INTENT(IN) :: fichpol
  REAL(QP) :: thresholds(1:2)
  thresholds(1)=contPole(k,fichpol)
  thresholds(2)=3*epsBCS(k/3.0_qp)

END FUNCTION thresholds
END MODULE selftot
