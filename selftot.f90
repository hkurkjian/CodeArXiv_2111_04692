MODULE selftot
USE modsim
USE vars
USE dspec
USE intpole
USE intldc
USE bestM
IMPLICIT NONE
LOGICAL usebestM

CONTAINS
SUBROUTINE initialisation(chargbestM,fichldc)
!Initialisation de dspec
x0=mu
temperaturenulle=.TRUE.
EPSpp=1.0e-7_qp
x0crit=0.0_qp
bla1=.FALSE.;bla2=.FALSE.
!Initialisation de intldc
ecrintq=0
bla0=.FALSE.;bla00=.FALSE.
!Initialisation de intpole
blaPole=.FALSE.
!Initialisation de bestM
qpetit=0.03_qp
blaM=.FALSE.;blaerr=.FALSE.;
if(chargbestM)then
     call load_data(fichldc(1))
     call loadom2  (fichldc(2))
endif
END SUBROUTINE(initialisation)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION detG(k,zk,fichiers,EPS,sigcomb)
  REAL, INTENT(IN) :: k,zk, EPS(1:3)
  CHARACTER(len=90), INTENT(IN) :: fichiers(1:2)
  COMPLEX(QPC), INTENT(OUT) :: sigcomb(1:6)
  COMPLEX(QPC) detG

  COMPLEX(QPC) sig(1:3)

  call initialisation(.FALSE.,(/"",""/))
  profondeurbidon=intbidon
  bqbidon(:)=bidon

  sigldc=cmplx(intpasres(k,zk,.TRUE.,.FALSE.,profondeurbidon,EPS(1:2),bqbidon,fichiers(1),""),0.0_qp)
  sigpol=selfEpole(k,zk,EPS(3),fichier(2))

  sig=sigpol(1:3)+sigpol(4:6)+sigldc(1:3)+sigldc(4:6)
  detG=(-zk+xik-sig(1))*(-zk-xik-sig(2))-(1.0_qp-sig(3))**2
END FUNCTION detG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION detGres(k,zk,fichiers,EPS,sigcomb,chargbestM)
  REAL, INTENT(IN) :: k,zk, EPS(1:3)
  CHARACTER(len=90), INTENT(IN) :: fichiers(1:5)
  COMPLEX(QPC), INTENT(OUT) :: sigcomb(1:6)
  LOGICAL, INTENT(IN) :: chargbestM
  COMPLEX(QPC) detGres

  REAL(QP) bk(0:12),le(1:8)
  COMPLEX(QPC) sig(1:3)

  call initialisation(chargbestM,fichiers(1:2))

  call bornesk(bk)
  call lignesenergie(k,fichiers(3:4),le)

  sigcomb(1,:)=intres(k,zk,.TRUE.,EPS(1:2),bk,le,"")
  sigcomb(2,:)=selfEpole(k,zk,EPS(3),fichier(3))

  sig=sigcomb(1,1:3)+sigcomb(1,4:6)+sigcomb(2,1:3)+sigcomb(2,4:6)
  detGres=(-zk+xik-sig(1))*(-zk-xik-sig(2))-(1.0_qp-sig(3))**2
END FUNCTION detGres
END MODULE selftot
