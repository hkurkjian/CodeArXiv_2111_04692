!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM Zero 
USE nrtype
USE dspec
USE recettes, ONLY : mnewt
IMPLICIT NONE

INTERFACE 
 SUBROUTINE derivee(zvec,detvec,Jdet)
 USE nrtype
 IMPLICIT NONE
 REAL(QP), DIMENSION(:), INTENT(IN) :: zvec
 REAL(QP), DIMENSION(:), INTENT(OUT) :: detvec
 REAL(QP), DIMENSION(:,:), INTENT(OUT) :: Jdet
 END SUBROUTINE derivee
END INTERFACE

REAL(QP), DIMENSION(1:1) :: det
REAL(QP), DIMENSION(1:1,1:1) :: Jdet
REAL(QP), DIMENSION(1:1) :: omdep
CHARACTER(len=90) fichier

REAL(QP) xqmin,xqmax,dxq,tolx,coup
REAL(QP) x0min,x0max,dx0
INTEGER ixq,nxq
INTEGER ix0,nx0

open(10,file='Zerodspec.inp')
 read(10,*)x0min
 read(10,*)x0max
 read(10,*)nx0
 read(10,*)xqmin
 read(10,*)xqmax
 read(10,*)nxq
 read(10,*)omdep(1)
 read(10,*)coup
 read(10,*)fichier
close(10)

call system('rm ' // fichier)

open(14,file=fichier)
  write(14,*)'# Résultats pour z (col 3 et 4) et alpha (col 5 et 6)'
  write(14,*)'# fct de x0, xq (col 1 et 2)'
close(14)

bla1=.FALSE.
temperaturenulle=.TRUE.

write(6,*)'--------------------'
write(6,*)
write(6,*)'Programme Zerodspec'
write(6,*)
write(6,*)'fichier=',fichier
write(6,*)'x0max=',x0max
write(6,*)'x0min=',x0min
write(6,*)'nx0=',nx0
write(6,*)'xqmin=',xqmin
write(6,*)'xqmax=',xqmax
write(6,*)'nxq=',nxq
write(6,*)'coup=',coup

if(nx0==0)then
 dx0=0.0
else
 dx0=(x0max-x0min)/nx0
endif

if(nxq==0)then
 dxq=0.0
else
 dxq=(xqmax-xqmin)/nxq
endif

do ix0=0,nx0
x0=x0min+dx0*ix0
do ixq=0,nxq
 xq=xqmin+dxq*ixq
 tolx=4.0_qp*xq**2.0_qp*1.e-11_qp
 write(6,*)'--------------------'
 write(6,*)'x0=',x0
 write(6,*)'xq=',xq
 write(6,*)'omdep=',omdep

 call mnewt(20,omdep,tolx,1.e-9_qp,derivee)

 open(14,file=fichier,position='APPEND')
  write(14,*)x0,xq,omdep
 close(14)

enddo
enddo
END PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE derivee(oom,det,Jdet)
USE dspec
IMPLICIT NONE
REAL(QP), DIMENSION(:), INTENT(IN) :: oom
REAL(QP), DIMENSION(:), INTENT(OUT) :: det
REAL(QP), DIMENSION(:,:), INTENT(OUT) :: Jdet
REAL(QP) omP,omM,detP,detM
REAL(QP) ddet1,ddet2 
!REAL(QP) domega,bidon(1:2,1:2)
REAL(QP) domega
COMPLEX(QPC) detc,Mbidon(1:3,1:3),invMbid(1:3,1:3)

write(6,*)'----------------------'
write(6,*)'oom=',oom

domega=(oom(1)-2.0_qp)*1.e-3_qp
domega=(oom(1))*1.e-4_qp

omP=oom(1)+domega
omM=oom(1)-domega


call mat_pairfield(oom(1),0.0_qp,detc,Mbidon,invMbid)
det(1)=real(detc)
write(6,*)'om,det=',oom(1),det
call mat_pairfield(omP,0.0_qp,detc,Mbidon,invMbid)
detP=real(detc)
write(6,*)'om,det=',omP,detP
call mat_pairfield(omM,0.0_qp,detc,Mbidon,invMbid)
detM=real(detc)
write(6,*)'om,det=',omM,detM

Jdet(1,1)=(detP-detM)/(2.0_qp*domega)

write(6,*)'det=',det(1)
write(6,*)'Jdet=',Jdet(1,1)
END SUBROUTINE
