PROGRAM pointsM
USE nrtype
USE vars
USE dspec
IMPLICIT NONE
COMPLEX(QPC) Gg(1:2,1:2),Mm(1:2,1:2),det
CHARACTER(len=90) fichier
REAL(QP) xqmin,xq1a,xq1b,xq2,xqmax,dxq,qsep(1:5),xqdep,xqfin
REAL(QP) om,ommin,ommax,dom,dom2,dom2p,dom3,dom3p,omsep(1:8)
REAL(QP) y,ymin,ymax,dy
REAL(QP) bmax,Mmv(1:6),donnees(1:7,1:10000)
INTEGER ixq,nq,nq1a,nq1b,nq2,nq3,ifenq,nqfen(1:4),nfenetres
INTEGER iom,nom,nom1a,nom1b,nom2a,nom2b,nom2c,nom3a,nom3b,ifen,nomfen(1:8),nominf
INTEGER compteur,compteurq,nn

open(10,file='pointsM.inp')
 read(10,*)x0
 read(10,*)xqmin
 read(10,*)nq1a
 read(10,*)nq1b
 read(10,*)nq2
 read(10,*)nq3
 read(10,*)nom1a
 read(10,*)nom1b
 read(10,*)nom2a
 read(10,*)nom2b
 read(10,*)nom2c
 read(10,*)nom3a
 read(10,*)nom3b
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
xq1a  =xqjoin/4.0_qp
xq1b  =xqjoin
xq2   =2*sqrt(x0)
xqmax =4*xq2
qsep=(/xqmin,xq1a,xq1b,xq2,xqmax/)
nomfen=(/nom1a,nom1b,nom2a,nom2b,nom2c,nom3a,nom3b,nominf/)
nqfen =(/nq1a,nq1b,nq2,nq3/)

write(6,*)'xqmin='  ,xqmin
write(6,*)'xq1a='   ,xq1a
write(6,*)'xq1b='   ,xq1b
write(6,*)'xq2='    ,xq2
write(6,*)'xqmax='  ,xqmax
write(6,*)'nom1a='  ,nom1a
write(6,*)'nom1b='  ,nom1b
write(6,*)'nom2a='  ,nom2a
write(6,*)'nom2b='  ,nom2b
write(6,*)'nom2c='  ,nom2c
write(6,*)'nom3a='  ,nom3a
write(6,*)'nom3b='  ,nom3b
write(6,*)'nominf=' ,nominf
write(6,*)'fichier=',trim(fichier)

inquire(iolength=nn)donnees(1:7,1:sum(nomfen))
write(6,*)"Taille de l’enregistrement en octets:",nn

bmax=1.0e6_qp

call system("rm "//trim(fichier)//".dat")
call system("rm "//trim(fichier)//"grilleq.dat")
call system("rm "//trim(fichier)//".info")
open(12,file=trim(fichier)//".info")
 write(12,*)"!x0,xqmin,xq1a,xq1b,xq2,xqmax,nq1a,nq1b,nq2,nq3,nom1a,nom1b,nom2a,nom2b,nom2c,nom3a,nom3a,nominf,nn,bmax"
 write(12,*)  x0,xqmin,xq1a,xq1b,xq2,xqmax,nqfen,nomfen,nominf,nn,bmax
close(12)
open(13,file=trim(fichier)//"grilleq.dat")
 write(13,*)"!ixq,xq,compteur,taille de l’enregistrement"
close(13)

compteurq=0
do ifenq=1,4

 if(ifenq.LE.2)then
  nfenetres=7
 elseif(ifenq==3)then
  nfenetres=4
 else
  nfenetres=2
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

  dom2  =min((opp(2)-opp(1))/2,(opp(3)-opp(2))/4)
  dom2p =min((opp(2)-opp(1))/8,5.0e-2_qp)
  if(ifenq.LE.2)then
   dom3  =min((opp(3)-opp(2))/8,5.0e-1_qp)
   dom3p =min((opp(3)-opp(2))/8,5.0e-2_qp)
   omsep=(/opp(1),opp(2)-dom2p,opp(2),opp(2)+dom2,opp(3)-dom3p,opp(3),opp(3)+dom3,4*opp(3)/)
  elseif(ifenq==3)then
   dom3  =min((opp(3)-opp(1))/8,5.0e-1_qp)
   dom3p =min((opp(3)-opp(1))/8,5.0e-2_qp)
   omsep(1:5)=(/opp(1),opp(3)-dom3p,opp(3),opp(3)+dom3,4*opp(3)/)
  else
   dom3  =5.0e-1_qp
   omsep(1:3)=(/opp(3),opp(3)+dom3,4*opp(3)/)
  endif

  do ifen=1,nfenetres
   ommin=omsep(ifen)
   ommax=omsep(ifen+1)
!Nbr de points
   if(ifenq.LE.2) nom=nomfen(ifen)  
   if(ifenq==3)then
     if(ifen==1) nom=nom1a !quand om2=om3 la premiere fenetre est de type 1 pas 2
     if(ifen==2) nom=nom1b 
     if(ifen==3) nom=nom3a 
     if(ifen==4) nom=nom3b 
   endif
   if(ifenq==4)  nom=nomfen(ifen+5)

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

!  ymin=1/ommax**(1.5_qp)
!  ymax=1/ommin**(1.5_qp)

  ymin=1/ommax**(0.5_qp)
  ymax=1/ommin**(0.5_qp)

  dy=(ymax-ymin)/nominf
  write(6,*)"ommim,ommax,ymin,ymax,dy=",ommin,ommax,ymin,ymax,dy
  do iom=1,nominf
   y=ymin+dy*(iom-0.5_qp)
!   om=1.0_qp/y**(2.0_qp/3.0_qp)
   om=1.0_qp/y**(2.0_qp)
   call mat_pairfield(om,0.0_qp,det,Mm,Gg)
   Mmv=(/real(Mm(1,1)),real(Mm(2,2)),real(Mm(1,2)),imag(Mm(1,1)),imag(Mm(2,2)),imag(Mm(1,2))/)
   write(6,*)"iom+compteur,om,Mmv(1)=",iom+compteur,om,Mmv(1)
   donnees(1,iom+compteur)=om
   donnees(2:7,iom+compteur)=Mmv
  enddo
  compteur=compteur+nominf

  open(11,file=trim(fichier)//".dat",ACTION="WRITE",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   write(11,REC=ixq+compteurq)donnees(1:7,1:sum(nomfen))
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
