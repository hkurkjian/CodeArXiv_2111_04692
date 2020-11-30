!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Zerom
USE dspec
USE vars
IMPLICIT NONE
LOGICAL bla2
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Zero(om,M,dM)
USE dspec
USE vars
USE recettes, ONLY : mnewt
IMPLICIT NONE
REAL(QP), INTENT(INOUT) :: om
REAL(QP), INTENT(OUT) :: M(1:3),dM(1:3)

REAL(QP) :: Mb(1:2,1:2),dMb(1:2,1:2)
REAL(QP), DIMENSION(1:1) :: omdep
REAL(QP) tolx

temperaturenulle=.TRUE.


tolx=4.0_qp*xq**2.0_qp*1.e-11_qp
if(bla2)then
 write(6,*)'--------------------'
 write(6,*)
 write(6,*)'Zerodspec'
 write(6,*)
 
 write(6,*)'--------------------'
 write(6,*)'x0=',x0
 write(6,*)'xq=',xq
 write(6,*)'om=',om
endif

omdep=(/om/)
call mnewt(20,omdep,tolx,1.e-9_qp,derivee)

om=omdep(1)

call der_mat(om,0.0_qp,Mb,dMb)
M =(/Mb(1,1),Mb(2,2),Mb(1,2)/)
dM=(/dMb(1,1),dMb(2,2),dMb(1,2)/)


END SUBROUTINE Zero
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE derivee(oom,det,Jdet)
USE dspec
IMPLICIT NONE
REAL(QP), DIMENSION(:), INTENT(IN) :: oom
REAL(QP), DIMENSION(:), INTENT(OUT) :: det
REAL(QP), DIMENSION(:,:), INTENT(OUT) :: Jdet
REAL(QP) Mr(1:2,1:2),dMbis(1:2,1:2)

if(bla2)then
 write(6,*)'----------------------'
 write(6,*)'oom=',oom
endif

call der_mat(oom(1),0.0_qp,Mr,dMbis)

det=Mr(1,1)*Mr(2,2)-Mr(1,2)**2
Jdet=dMbis(1,1)*Mr(2,2)+Mr(1,1)*dMbis(2,2)-2.0_qp*Mr(1,2)*dMbis(1,2)
if(bla2)then
 write(6,*)'om,det,Jdet=',oom(1),det,Jdet
 write(6,*)'Jdet=',Jdet(1,1)
endif
END SUBROUTINE derivee
END MODULE Zerom
