PROGRAM encorr
USE vars
USE recettes
USE selftot
USE OMP_LIB
IMPLICIT NONE

REAL(QP) :: mu,k,zk
COMPLEX(QPC) sigma(1:2,1:6),sigcomb(1:3),det
REAL(QP) dk,dzk
REAL(QP) kmin,kmax,zkmin,zkmax,EPS(1:3)
REAL(QP) bqbidon(1:3),th(1:4)
INTEGER nk,nzk,nkdeb,ik,izk,nzkdeb,ic,nivobla,eintq,profondeurbidon
CHARACTER(len=90) fichiers(1:5),suffixe,suffintq
CHARACTER(len=5) cik,cizk
LOGICAL nvofich

!k 0 to 3
!zk 1 to 6
open(10,file='encorr.inp')
 read(10,*)mu          !interaction regime =x0 in dspec
 read(10,*)kmin        !grid of points in k
 read(10,*)kmax
 read(10,*)nk,nkdeb
 read(10,*)zkmin       !grid of points in zk
 read(10,*)zkmax
 read(10,*)nzk,nzkdeb
 read(10,*)fichiers(1) !pour charger estM/donnees
 read(10,*)fichiers(2) !pour intldc/intpasres
 read(10,*)fichiers(3) !pour intpole
 read(10,*)EPS(1)      !intldc/EPSq
 read(10,*)EPS(2)      !intldc/EPSom
 read(10,*)EPS(3)      !intpole/EPSqpole
 read(10,*)nivobla !Verbose level, from 0 to 3
 read(10,*)suffixe !terminaison of output files
 read(10,*)nvofich !TRUE to overwrite output files
 read(10,*)eintq !0 to avoid writing intq files, 1 to write, 2 to overwrite
close(10)


write(6,*)'--------------------'
write(6,*)
write(6,*)'Programme encorr'
write(6,*)
write(6,*)'suffixe=',suffixe
write(6,*)'kmin=',kmin
write(6,*)'kmax=',kmax
write(6,*)'nk,nkdeb=',nk,nkdeb
write(6,*)'zkmin=',zkmin
write(6,*)'zkmax=',zkmax
write(6,*)'nzk,nzkdeb=',nzk,nzkdeb
write(6,*)'fichiers ldc:',trim(fichiers(1))
write(6,*)'fichier  lec:',trim(fichiers(2))
write(6,*)'fichier pole:',trim(fichiers(3))
write(6,*)'precisions:',EPS

call initialisation(mu,nivobla,fichiers(1:3),eintq)

if(nk==0)then
 dk=0.0
else
 dk=(kmax-kmin)/nk
endif

if(nzk==0)then
 dzk=0.0
else
 dzk=(zkmax-zkmin)/nzk
endif

if(nvofich)then
 open(20,file="DONNEES/selfEldc"//trim(suffixe)//".dat")
  write(20,*)"!Valeurs de k,zk et selfEldc pour x0=",x0
 close(20)
 open(21,file="DONNEES/selfEpol"//trim(suffixe)//".dat")
  write(21,*)"!Valeurs de k,zk et selfEpol pour x0=",x0
 close(21)
 open(22,file="DONNEES/selfEtot"//trim(suffixe)//".dat")
  write(22,*)"!Valeurs de k,zk et selfEtot pour x0=",x0
 close(22)
 open(23,file="DONNEES/det"//trim(suffixe)//".dat")
  write(23,*)"!Valeurs de k,zk et detG pour x0=",x0
 close(23)
 open(24,file="DONNEES/continuum"//trim(suffixe)//".dat")
  write(24,*)"!Valeurs de k et des bords de continuum pour x0=",x0
 close(24)
endif

izk=-1+nzkdeb
ik=nkdeb
!$OMP  PARALLEL DO &
!$OMP& PRIVATE(k,zk,cik,cizk,th,suffintq,det,sigcomb,sigma) &
!$OMP& SHARED(ik,izk,fichiers,EPS,mu,kmin,zkmin,dk,dzk) SCHEDULE(DYNAMIC)
do ic=1,(nk-nkdeb+1)*(nzk-nzkdeb+1)
!$OMP CRITICAL
 izk=izk+1
 if(izk==nzk+1)then
  izk=0
  ik=ik+1
 endif
 write(6,*)"ik,izk,fil=",ik,izk,OMP_GET_THREAD_NUM()
 k=kmin+dk*ik
 zk=zkmin+dzk*izk
 write(6,*)"k,zk=",k,zk

 write(cik, FMT="(I2)")ik
 write(cizk,FMT="(I2)")izk
 cik=adjustl(cik)
 cizk=adjustl(cizk)
 write(suffintq,*)"_",trim(cik),"_",trim(cizk)
 suffintq=adjustl(suffintq)
 write(6,*)"fil,suffintq:",omp_get_thread_num(),trim(suffintq)
 !$OMP END CRITICAL

 th=thresholds(mu,k)
 write(6,*)"th=",th(1),th(4)
 if(zk<th(4))then
  det=detG   (k,zk,EPS,sigma,suffintq)
 else
  det=detGres(k,zk,EPS,sigma,suffintq)
 endif
 
 sigcomb(1:3)=sigma(1,1:3)+sigma(1,4:6)+sigma(2,1:3)+sigma(2,4:6)

!$OMP CRITICAL
 open(20,file="DONNEES/selfEldc"//trim(suffixe)//".dat",POSITION="APPEND")
 open(21,file="DONNEES/selfEpol"//trim(suffixe)//".dat",POSITION="APPEND")
 open(22,file="DONNEES/selfEtot"//trim(suffixe)//".dat",POSITION="APPEND")
 open(23,file="DONNEES/det"     //trim(suffixe)//".dat",POSITION="APPEND")
  write(6,*)"k,zk,bords des continua=",k,zk,th(1),th(4)
  write(6,*)"re(selfEldc)=",real(sigma(1,1:6))
  write(6,*)"im(selfEldc)=",imag(sigma(1,1:6))
  write(6,*)"re(selfEpol)=",real(sigma(2,1:6))
  write(6,*)"im(selfEpol)=",imag(sigma(2,1:6))
  write(6,*)"re(selfEtot)=",real(sigcomb(1:3))
  write(6,*)"im(selfEtot)=",imag(sigcomb(1:3))
  write(6,*)"det=",real(det),imag(det)
  write(6,*)
  write(20,*)k,zk,real(sigma(1,1:6)),imag(sigma(1,1:6))
  write(21,*)k,zk,real(sigma(2,1:6)),imag(sigma(2,1:6))
  write(22,*)k,zk,real(sigcomb(1:3)),imag(sigcomb(1:3))
  write(23,*)k,zk,real(det),         imag(det)
  if(izk==nzk)then
   open(24,file="DONNEES/continuum"//trim(suffixe)//".dat",POSITION="APPEND")
    write(24,*)k,th(1),th(4),th(2:3)
   close(24)
  endif
 close(20)
 close(21)
 close(22)
 close(23)
 !$OMP END CRITICAL
enddo
!$OMP END PARALLEL DO

call desinit
call desini(.TRUE.,.FALSE.)
END PROGRAM encorr
