PROGRAM pointsM
USE nrtype
USE vars
USE dspec
IMPLICIT NONE
COMPLEX(QPC) Gg(1:2,1:2),Mm(1:2,1:2),det
CHARACTER(len=90) fichier
REAL(QP) xqmin,xq1,xq2,xqmax,dxq,qsep(1:4),xqdep,xqfin
REAL(QP) om,ommin,ommax,dom,dom2,omsep(1:5)
REAL(QP) y,ymin,ymax,dy
REAL(QP) bmax,Mmv(1:6),donnees(1:7,1:10000)
INTEGER ixq,nq,nq1,nq2,nq3,ifenq,nqfen(1:3),nfenetres
INTEGER iom,nom,nom1,nom2,nom3,nom4,ifen,nomfen(1:4),nominf
INTEGER compteur,compteurq,nn
LOGICAL nouveauf

open(10,file='pointsM.inp')
 read(10,*)x0
 read(10,*)xqmin
 read(10,*)nq1
 read(10,*)nq2
 read(10,*)nq3
 read(10,*)nom1
 read(10,*)nom2
 read(10,*)nom3
 read(10,*)nom4
 read(10,*)nominf
 read(10,*)fichier
close(10)

bla1=.TRUE.
bla1=.FALSE.

EPSpp=1.0e-8_qp
x0crit=4.0_qp
beta=bidon
temperaturenulle=.TRUE.

call calcxqjoin
xq1  =xqjoin
xq2  =2*sqrt(x0)
xqmax=4*xq2
qsep=(/xqmin,xq1,xq2,xqmax/)
nomfen=(/nom1,nom2,nom3,nom4/)
nqfen =(/nq1,nq2,nq3/)

write(6,*)'xqmin=' ,xqmin
write(6,*)'xq1='   ,xq1
write(6,*)'xq2='   ,xq2
write(6,*)'xqmax=' ,xqmax
write(6,*)'nom1='  ,nom1
write(6,*)'nom2='  ,nom2
write(6,*)'nom3='  ,nom3
write(6,*)'nom4='  ,nom4
write(6,*)'nominf=',nominf
write(6,*)'fichier=',trim(fichier)

inquire(iolength=nn)donnees(1:7,1:nom1+nom2+nom3+nom4+nominf)
write(6,*)"Taille de l’enregistrement en octets:",nn

nouveauf=.FALSE.
nouveauf=.TRUE.
bmax=1.0e6_qp

if(nouveauf)then
 call system("rm "//trim(fichier)//".dat")
 call system("rm "//trim(fichier)//"grilleq.dat")
 call system("rm "//trim(fichier)//".info")
 open(12,file=trim(fichier)//".info")
  write(12,*)"!x0,xqmin,xq1,xq2,xqmax,nq1,nq2,nq3,nom1,nom2,nom3,nom4,nominf,nn,bmax"
  write(12,*)  x0,xqmin,xq1,xq2,xqmax,nq1,nq2,nq3,nom1,nom2,nom3,nom4,nominf,nn,bmax
 close(12)
 open(13,file=trim(fichier)//"grilleq.dat")
  write(13,*)"!ixq,xq,compteur,taille de l’enregistrement"
 close(13)
endif



compteurq=0
do ifenq=1,3

 if(ifenq==1)then
  nfenetres=4
 elseif(ifenq==2)then
  nfenetres=2
 else
  nfenetres=1
 endif

 xqdep=qsep(ifenq)
 xqfin=qsep(ifenq+1)
 nq=nqfen(ifenq)
 dxq=(xqfin-xqdep)/nq
 write(6,*)'--------------------'
 write(6,*)
 write(6,*)"xqdep,xqfin,dxq=",xqdep,xqfin,dxq
 do ixq=1,nq
  xq=xqdep+dxq*(ixq-0.5_qp)
  call oangpp
  write(6,*)'--------------------'
  write(6,*)
  write(6,*)"xq,ixq+compteurq=",xq,ixq+compteurq
  write(6,*)"opp=",opp(1:3)
  write(6,*) 

  compteur=0
  donnees(:,:)=0.0_qp

  dom2=min(opp(2)-opp(1),(opp(3)-opp(2))/4)
  if(ifenq==1)then
   omsep=(/opp(1),opp(2),opp(2)+dom2,opp(3),4*opp(3)/)
  elseif(ifenq==2)then
   omsep(1:3)=(/opp(1),opp(3),4*opp(3)/)
  else
   omsep(1:2)=(/opp(3),4*opp(3)/)
  endif

  do ifen=1,nfenetres
   ommin=omsep(ifen)
   ommax=omsep(ifen+1)
!Nbr de points
   if(ifenq==1) nom=nomfen(ifen)  
   if(ifenq==2)then
     if(ifen==1) nom=nom1 !quand om2=om3 la premiere fenetre est de type 1 pas 2
     if(ifen==2) nom=nom3 
   endif
   if(ifenq==3)  nom=nom3

   dom=(ommax-ommin)/nom
   write(6,*)"ommim,ommax,dom=",ommin,ommax,dom
   do iom=1,nom
    om=ommin+dom*(iom-0.5_qp)
    call mat_pairfield(om,0.0_qp,det,Mm,Gg)
    Mmv=(/real(Mm(1,1)),real(Mm(2,2)),real(Mm(1,2)),imag(Mm(1,1)),imag(Mm(2,2)),imag(Mm(1,2))/)
    write(6,*)"iom+compteur,om,Mmv(1)=",iom+compteur,om,Mmv(1)
    donnees(1,iom+compteur)=om
    donnees(2:7,iom+compteur)=Mmv
   enddo
   compteur=compteur+nom
  enddo

  ommin=4*opp(3)
  ommax=bmax
  ymin=1/ommax**(1.5_qp)
  ymax=1/ommin**(1.5_qp)
  dy=(ymax-ymin)/nominf
  write(6,*)"ommim,ommax,ymin,ymax,dy=",ommin,ommax,ymin,ymax,dy
  do iom=1,nominf
   y=ymin+dy*(iom-0.5_qp)
   om=1.0_qp/y**(2.0_qp/3.0_qp)
   call mat_pairfield(om,0.0_qp,det,Mm,Gg)
   Mmv=(/real(Mm(1,1)),real(Mm(2,2)),real(Mm(1,2)),imag(Mm(1,1)),imag(Mm(2,2)),imag(Mm(1,2))/)
   write(6,*)"iom+compteur,om,Mmv(1)=",iom+compteur,om,Mmv(1)
   donnees(1,iom+compteur)=om
   donnees(2:7,iom+compteur)=Mmv
  enddo
  compteur=compteur+nominf

  open(11,file=trim(fichier)//".dat",ACTION="WRITE",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   write(11,REC=ixq+compteurq)donnees(1:7,1:nom1+nom2+nom3+nom4+nominf)
!   write(11,REC=ixq+compteurq)ixq,om
 close(11)
  open(13,file=trim(fichier)//"grilleq.dat",POSITION="APPEND")
   write(13,*)ixq+compteurq,xq,compteur,nn
  close(13)
  write(6,*)
 enddo
 compteurq=compteurq+nq
enddo
END PROGRAM pointsM
