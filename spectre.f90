PROGRAM spectre
 USE nrtype
 USE modsim
 USE vars
 USE dspec
 USE Zerom
 USE eqdetat
 IMPLICIT NONE
 REAL(QP) :: om,M(1:3),dM(1:3),x
 INTEGER nn,n,taille
 CHARACTER(len=90) fichier
 REAL(QP) xqmin,xqmax,dxq,tolx,omdep
 REAL(QP) ccheck,Xx,Theta
 INTEGER ixq,ixqdep,nxq
 LOGICAL nouveauf

 open(10,file='spectre.inp')
  read(10,*)x0
  read(10,*)xqmin
  read(10,*)xqmax
  read(10,*)nxq
  read(10,*)om
  read(10,*)fichier
 close(10)

 bla1=.TRUE.
 bla1=.FALSE.
 bla3=.FALSE.
 bla3=.TRUE.

 EPSpp=1.0e-8_qp
 EPSrpp=1.0e-10_qp
 x0crit=10.0_qp
 beta=bidon

 if(nxq==0)then
  dxq=0.0
 else
  dxq=(xqmax-xqmin)/nxq
 endif

 Theta =4.0_qp*PI*I5(x0)
 Xx    =4.0_qp*PI*I6(x0)
 ccheck=mc2sDelta(x0)

 xq=xqmin
 write(6,*)'xq=',xq
 write(6,*)'om=',om
 write(6,*)'fichier=',trim(fichier)
 inquire(iolength=nn)xq,om,M,dM
 write(6,*)"Taille d un enregistrement en octets:",nn

 nouveauf=.FALSE.
 nouveauf=.TRUE.
 if(nouveauf)then
  call system("rm "//trim(fichier)//".dat")
  call system("rm "//trim(fichier)//".info")
  open(10,file=trim(fichier)//".info")
   write(10,*)"!x0,Theta,Xx,ccheck,xqmin,xqmax,nxq,taille enregistrement"
   write(10,*)x0,Theta,Xx,ccheck,xqmin,xqmax,nxq,nn
  close(10)
 endif

 ixqdep=0
 do ixq=ixqdep,nxq
  xq=xqmin+dxq*ixq
  tolx=4.0_qp*xq**2.0_qp*1.e-11_qp
  write(6,*)'--------------------'
  write(6,*)
  write(6,*)"ixq=",ixq
  write(6,*)
  call Zero(om,M,dM)
  write(6,*)
  write(6,*)"om,M=",om,M
  write(6,*)"dM=",dM
 
  open(11,file=trim(fichier)//".dat",ACTION="WRITE",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   write(11,REC=ixq+1)xq,om,M,dM
  close(11)

 enddo
END PROGRAM spectre
