program test2
USE modsim
USE recettes
USE dspec
USE vars
USE Zerom
USE intpole
USE eqdetat
IMPLICIT NONE

COMPLEX(QPC) selfEpol(1:6)
REAL(QP) k,zk

k= 2.1_qp
zk=3.4_qp
blaPole=.TRUE.
fichpol="BCS_4_pole"

selfEpol=selfEpole(k,zk,1.0e-6_qp)
write(6,*)"selfEpol=",selfEpol
END program test2
