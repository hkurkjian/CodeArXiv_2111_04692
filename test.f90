program test
USE modsim
USE recettes
USE dspec
USE vars
USE Zerom
USE intpole
IMPLICIT NONE

REAL(QP) om,dom,M(1:2,1:2),dM(1:2,1:2),A(1:6),Mm(1:3),dMm(1:3)
REAL(QP) :: xik,epsk,Uk2,Vk2,k
COMPLEX(QPC) Gam(1:2,1:2),Matt(1:2,1:2),det,Gamp(1:2,1:2),Matp(1:2,1:2),Gamm(1:2,1:2),Matm(1:2,1:2)
COMPLEX(QPC) Gam3(1:3,1:3),Mat3(1:3,1:3)

x0=0.860436686125678599999999999999999961_qp
xq=0.5_qp
om=0.5_qp
bla2=.TRUE.

EPSpp=1.0e-8_qp
EPSrpp=1.0e-10_qp
EPSpt=1.0e-6_qp
EPSrpt=1.0e-7_qp

call Zero(om,Mm,dMm)

write(6,*) "om=",om
stop

fichier="lu"
A=selfEpole(sqrt(x0),1.0_qp)

k=sqrt(x0)
xik=k**2.0_qp-x0
epsk=sqrt(xik**2+1)
Uk2=(1+xik/epsk)/2
Vk2=(1-xik/epsk)/2
write(6,*)"deltaeps=",A(3)/epsk-Uk2*A(1)-Vk2*A(2)
write(6,*)"deltaeps=",A(6)/epsk-Uk2*A(4)-Vk2*A(5)

stop















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

x0crit=19.0_qp
xq=1.0_qp
om=1.995123261642178399749904036442505_qp
dom=1.0e-4_qp
bla1=.FALSE.
bla2=.TRUE.

call Zero(om,M,dM)
call Zero(om,M,dM)
write(6,*)"om,M,dM=",om,M,dM
stop

call der_mat(om,0.0_qp,M,dM)
write(6,*)"quasip-quasip    angular points=",opp(1:3)
write(6,*)"quasip-quasihole angular point =",opt(1)
write(6,*)"dM   =",dM(1,1),dM(2,2),dM(1,2)

stop
call mat_pairfield(om    ,0.0_qp,det,Matt,Gam )
call mat_pairfield(om+dom,0.0_qp,det,Matp,Gamp)
call mat_pairfield(om-dom,0.0_qp,det,Matm,Gamm)

write(6,*)"dM1   =",(Matp(1,1)-Matm(1,1))/dom/2,(Matp(2,2)-Matm(2,2))/dom/2,(Matp(1,2)-Matm(1,2))/dom/2
write(6,*)"dM2   =",dM(1,1),dM(2,2),dM(1,2)

stop

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
