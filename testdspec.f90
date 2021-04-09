program test
USE dspec
USE vars
USE OMP_LIB
IMPLICIT NONE

REAL(QP) om,dom,M(1:2,1:2),dM(1:2,1:2),A(1:6),dMm(1:3),q
COMPLEX(QPC) Mm2(1:6)
REAL(QP) Mmbid(1:6,1:3),rho(1:2,1:2)
COMPLEX(QPC) Gamm(1:2,1:2),Matt(1:2,1:2),MatCat(1:2,1:2),det

LOGICAL testdspec

testdspec=.FALSE.
testdspec=.TRUE.

!Paramètres de dspec
temperaturenulle=.TRUE.
EPSpp=1.0e-9_qp
x0crit=0.0_qp
x0crit=1.0_qp
bla1=.FALSE.
bla1=.TRUE.
bla2=.TRUE.
bla2=.FALSE.

 
call calcxqjoin
x0=4.0_qp
x0=100.0_qp
x0=0.8604366861256786_qp
xq=0.270614857
om=2.00008678567156939262684483498551258_qp
om=2.05751104503364324991542327941859321_qp
xq=1.57347200000000000000000000000000021_qp
xq=2.03009440145508893240000000000000020
om=2.263373185916696019076155463109632531
om=2.17678420571626069824921576610125391
xq=1.31290057142857142857142857142857159
om=2.00042250152115481492663972570040443_qp
xq=1.85523657142857142857142857142857159_qp
xq=0.921174857142857142857142857142857258
om=9.03277526534274179047385969699288689
om=2.25922195715940769053137863442953025_qp
xq=1.15768685714285714285714285714285729_qp
call oangpp
write(6,*)"opp=",opp(1:3)
 
!!$OMP PARALLEL DO
! do itest=1,4
! xq=0.1_qp+itest
! call oangpp
! write(6,*)"opp=",opp(1:3)
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
! call mat_pairfield_gom0(om,0.0_qp,det,Matt,Gamm)
! MatCat(1,1)=(Matt(1,1)+Matt(2,2))/2.0_qp+Matt(1,2)
! MatCat(2,2)=(Matt(1,1)+Matt(2,2))/2.0_qp-Matt(1,2)
! MatCat(1,2)=(Matt(2,2)-Matt(1,1))/2.0_qp
!
! Gamm(1,1)= MatCat(2,2)/det
! Gamm(2,2)= MatCat(1,1)/det
! Gamm(1,2)=-MatCat(1,2)/det
!
! rho=-imag(Gamm)/PI
!
! write(6,*)"om,xq,real(Matt(1,1))=",om,xq,real(Matt(1,1))
! write(6,*)"om,xq,real(Matt(2,2))=",om,xq,real(Matt(2,2))
! write(6,*)"om,xq,real(Matt(1,2))=",om,xq,real(Matt(1,2))
! write(6,*)"om,xq,real(Matt(1,1))=",om,xq,imag(Matt(1,1))
! write(6,*)"om,xq,real(Matt(2,2))=",om,xq,imag(Matt(2,2))
! write(6,*)"om,xq,real(Matt(1,2))=",om,xq,imag(Matt(1,2))
! write(6,*)
! write(6,*)"om,xq,real(rho(1,1))=",om,xq,rho(1,1)
! write(6,*)"om,xq,real(rho(2,2))=",om,xq,rho(2,2)
! write(6,*)"om,xq,real(rho(1,2))=",om,xq,rho(1,2)
! write(6,*)
! write(6,*)"det=",det
! write(6,*)
! enddo
!!$OMP END PARALLEL DO
end program
