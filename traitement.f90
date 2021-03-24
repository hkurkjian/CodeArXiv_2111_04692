PROGRAM traitement
USE nrtype
IMPLICIT NONE

REAL(QP) :: mu,k,zk,xik,sigr(1:6),kvieux,redet,imdet
COMPLEX(QPC) sig(1:3),det,fr(1:3)
REAL(QP) bqbidon(1:3),th(1:2)
INTEGER ik,err

mu=4.0_qp
open(13,file='DONNEES/detpk2.dat')
open(14,file='DONNEES/detpk4.dat')
 do ik=1,1000000
  kvieux=k
  read(13,*,end=11,iostat=err)k,zk,redet,imdet
  if(err.NE.0) cycle
  write(6,*)"k,kvieux,abs(kvieux-k)=",k,kvieux,abs(kvieux-k)
  if(abs(kvieux-k)>0.0001_qp) write(14,*)
  write(14,*)k,zk,redet,imdet
 enddo
11 CONTINUE
close(13)
close(14)
write(6,*)"Fini1!"
open(23,file='DONNEES/dett3gk2.dat')
open(24,file='DONNEES/dett3gk4.dat')
 do ik=1,1000000
  kvieux=k
  read(23,*,end=12,iostat=err)k,zk,redet,imdet
  if(err.NE.0) cycle
  write(6,*)"k,kvieux,abs(kvieux-k)=",k,kvieux,abs(kvieux-k)
  if(abs(kvieux-k)>0.0001_qp) write(24,*)
  write(24,*)k,zk,redet,imdet
 enddo
12 CONTINUE
close(23)
close(24)
write(6,*)"Fini2!"
open(33,file='DONNEES/selfEtott3gk2.dat')
open(34,file='DONNEES/frt3gk.dat')
 do ik=1,1000000
  kvieux=k
  read(33,*,end=13,iostat=err)k,zk,sigr
  if(abs(kvieux-k)>0.01_qp) write(34,*)
  sig=cmplx(sigr(1:3),sigr(4:6),kind=qpc)
  xik=k**2-mu
  det=(-zk+xik-sig(1))*(-zk-xik-sig(2))-(1.0_qp-sig(3))**2
  fr(1)=sig(1)/det
  fr(2)=sig(2)/det
  fr(3)=-sig(3)/det
  write(6,*)"k,kvieux,abs(kvieux-k)=",k,kvieux,abs(kvieux-k)
  write(34,*)k,zk,real(fr),imag(fr)
 enddo
13 CONTINUE
close(12)
write(6,*)"Fini3!"
END PROGRAM traitement
