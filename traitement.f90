PROGRAM traitement
USE nrtype
USE intpole
IMPLICIT NONE

REAL(QP) :: mu,k,zk,xik,siggr(1:6),kvieux,redet,imdet
COMPLEX(QPC) sigg(1:3),det,fr(1:3)
REAL(QP) bqbidon(1:3),th(1:2)
INTEGER ik,err,iq
CHARACTER(len=90) fich

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
open(13,file='DONNEES/detgk2.dat')
open(14,file='DONNEES/detgk4.dat')
 do ik=1,1000000
  kvieux=k
  read(13,*,end=15,iostat=err)k,zk,redet,imdet
  if(err.NE.0) cycle
  write(6,*)"k,kvieux,abs(kvieux-k)=",k,kvieux,abs(kvieux-k)
  if(abs(kvieux-k)>0.0001_qp) write(14,*)
  write(14,*)k,zk,redet,imdet
 enddo
15 CONTINUE
close(13)
close(14)
write(6,*)"Fini1bis!"
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
  read(33,*,end=13,iostat=err)k,zk,siggr
  if(abs(kvieux-k)>0.01_qp) write(34,*)
  sigg=cmplx(siggr(1:3),siggr(4:6),kind=qpc)
  xik=k**2-mu
  det=(-zk+xik-sigg(1))*(-zk-xik-sigg(2))-(1.0_qp-sigg(3))**2
  fr(1)=sigg(1)/det
  fr(2)=sigg(2)/det
  fr(3)=-sigg(3)/det
  write(6,*)"k,kvieux,abs(kvieux-k)=",k,kvieux,abs(kvieux-k)
  write(34,*)k,zk,real(fr),imag(fr)
 enddo
13 CONTINUE
close(12)
write(6,*)"Fini3!"
fich="DONNEES/BCS_4_pole"
call rdInfo(fich)
open(55,file='DONNEES/dispLU.dat')
do iq=1,nqeff
 write(55,*)donpol(1:2,iq)
enddo
close(55)
END PROGRAM traitement
