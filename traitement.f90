PROGRAM traitement
USE nrtype
IMPLICIT NONE

REAL(QP) :: mu,k,zk,xik,sigr(1:6),kvieux,redet,imdet
COMPLEX(QPC) sig(1:3),det
REAL(QP) bqbidon(1:3),th(1:2)
INTEGER ik,err

mu=4.0_qp
!open(10,file='DONNEES/selfEtotabove3.dat')
!open(11,file='DONNEES/cartebove3.dat')
!open(12,file='DONNEES/selfEtotabove4.dat')
open(13,file='DONNEES/detbelow4.dat')
open(14,file='DONNEES/detbelow6.dat')
 do ik=1,1000000
  kvieux=k
!  read(10,*,end=11,iostat=err)k,zk,sigr
  read(13,*,end=11,iostat=err)k,zk,redet,imdet
  if(err.NE.0) cycle
  write(6,*)"k,kvieux,abs(kvieux-k)=",k,kvieux,abs(kvieux-k)
  if(abs(kvieux-k)>0.0001_qp) write(14,*)
!  if(abs(kvieux-k)>0.01_qp) write(12,*)
!  if(abs(kvieux-k)>0.01_qp) write(11,*)
!  write(12,*)k,zk,sigr
!  sig=cmplx(sigr(1:3),sigr(4:6),kind=qpc)
!  xik=k**2-mu
!  det=(-zk+xik-sig(1))*(-zk-xik-sig(2))-(1.0_qp-sig(3))**2
!  write(11,*)k,zk,real(det),imag(det)
  write(14,*)k,zk,redet,imdet
 enddo
!close(10)
!close(11)
!close(12)
close(13)
11 CONTINUE
write(6,*)"Fini!"
END PROGRAM traitement
