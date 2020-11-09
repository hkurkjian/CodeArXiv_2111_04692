program test
USE nrtype
USE modsim
USE recettes
USE dspec
USE vars
IMPLICIT NONE

REAL(QP) om
COMPLEX(QPC) Gam(1:2,1:2),Matt(1:2,1:2),det
COMPLEX(QPC) Gam3(1:3,1:3),Mat3(1:3,1:3)

om=7.12_qp
xq=0.5_qp
x0=10._qp
beta=0.0_qp

bla1=.FALSE.
bla1=.TRUE.

temperaturenulle=.TRUE.
!temperaturenulle=.FALSE.

EPSpp=1.0e-8_qp
EPSrpp=1.0e-10_qp
EPSpt=1.0e-6_qp
EPSrpt=1.0e-7_qp

xq=12.11111111111111111111111111111111154_qp
om=150.095123261642178399749904036442505_qp

x0crit=9.0_qp
call mat(om,0.0_qp,det,Mat3,Gam3)
write(6,*)"quasip-quasip    angular points=",opp(1:3)
write(6,*)"quasip-quasihole angular point =",opt(1)
write(6,*)"re Mat Delta  =",real(Mat3(1,1)),real(Mat3(1,2)),real(Mat3(2,2))
write(6,*)"re Mat density=",real(Mat3(3,3)),real(Mat3(1,3)),real(Mat3(2,3))
write(6,*)"im Mat Delta  =",imag(Mat3(1,1)),imag(Mat3(1,2)),imag(Mat3(2,2))
write(6,*)"im Mat density=",imag(Mat3(3,3)),imag(Mat3(1,3)),imag(Mat3(2,3))

write(6,*)
write(6,*)
x0crit=19.0_qp
call mat(om,0.0_qp,det,Mat3,Gam3)
write(6,*)"re Mat Delta  =",real(Mat3(1,1)),real(Mat3(1,2)),real(Mat3(2,2))
write(6,*)"re Mat density=",real(Mat3(3,3)),real(Mat3(1,3)),real(Mat3(2,3))
write(6,*)"im Mat Delta  =",imag(Mat3(1,1)),imag(Mat3(1,2)),imag(Mat3(2,2))
write(6,*)"im Mat density=",imag(Mat3(3,3)),imag(Mat3(1,3)),imag(Mat3(2,3))
end program
