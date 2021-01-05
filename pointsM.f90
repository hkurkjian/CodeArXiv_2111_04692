PROGRAM pointsM
USE nrtype
USE vars
USE dspec
IMPLICIT NONE
COMPLEX(QPC) Gg(1:2,1:2),Mm(1:2,1:2),det
CHARACTER(len=90) fichier
REAL(QP) xqmin,xqmax,dxq,om,ommin,ommax,dom,y,ymin,ymax,dy,bmax,Mmv(1:6),donnees(1:7,1:10000)
INTEGER ixq,ixqdep,nxq,iom,nom1,nom2,nom3,nominf,compteur,nn
LOGICAL nouveauf

open(10,file='pointsM.inp')
 read(10,*)x0
 read(10,*)xqmin
 read(10,*)xqmax
 read(10,*)nxq
 read(10,*)nom1
 read(10,*)nom2
 read(10,*)nom3
 read(10,*)nominf
 read(10,*)fichier
close(10)

bla1=.TRUE.
bla1=.FALSE.

EPSpp=1.0e-8_qp
x0crit=4.0_qp
beta=bidon
temperaturenulle=.TRUE.

if(nxq==0)then
 dxq=0.0
else
 dxq=(xqmax-xqmin)/nxq
endif

write(6,*)'xqmin=',xqmin
write(6,*)'xqmax=',xqmax
write(6,*)'nom1='  ,nom1
write(6,*)'nom2='  ,nom2
write(6,*)'nom3='  ,nom3
write(6,*)'nominf=',nominf
write(6,*)'fichier=',trim(fichier)

nouveauf=.FALSE.
nouveauf=.TRUE.
if(nouveauf)then
 call system("rm "//trim(fichier)//".dat")
 call system("rm "//trim(fichier)//"grilleq.dat")
 call system("rm "//trim(fichier)//".info")
 open(12,file=trim(fichier)//".info")
  write(12,*)"!x0,xqmin,xqmax,nxq,nom1,nom2,nom3,nominf"
  write(12,*)  x0,xqmin,xqmax,nxq,nom1,nom2,nom3,nominf
 close(12)
 open(12,file=trim(fichier)//"grilleq.dat")
 close(12)
endif

bmax=1.0e6_qp
ixqdep=0
do ixq=ixqdep,nxq
 xq=xqmin+dxq*ixq
 call oangpp
 write(6,*)'--------------------'
 write(6,*)
 write(6,*)"xq,ixq=",xq,ixq
 write(6,*)"opp=",opp(1:3)
 write(6,*) 

 compteur=0

 if(ptbranchmtpp>1)then
  ommin=opp(1)
  ommax=opp(2)
  dom=(ommax-ommin)/(nom1+1)
  write(6,*)"ommim,ommax,dom=",ommin,ommax,dom
  do iom=1,nom1
   om=ommin+dom*iom
   call mat_pairfield(om,0.0_qp,det,Mm,Gg)
   Mmv=(/real(Mm(1,1)),real(Mm(2,2)),real(Mm(1,2)),imag(Mm(1,1)),imag(Mm(2,2)),imag(Mm(1,2))/)
   write(6,*)"iom+compteur,om,Mmv(1)=",iom+compteur,om,Mmv(1)
   donnees(1,iom+compteur)=om
   donnees(2:7,iom+compteur)=Mmv
  enddo
 compteur=compteur+nom1
 endif

 if(ptbranchmtpp>2)then
  ommin=opp(2)
  ommax=opp(3)
  dom=(ommax-ommin)/(nom1+1)
  write(6,*)"ommim,ommax,dom=",ommin,ommax,dom
  do iom=1,nom2
   om=ommin+dom*iom
   call mat_pairfield(om,0.0_qp,det,Mm,Gg)
   Mmv=(/real(Mm(1,1)),real(Mm(2,2)),real(Mm(1,2)),imag(Mm(1,1)),imag(Mm(2,2)),imag(Mm(1,2))/)
   write(6,*)"iom+compteur,om,Mmv(1)=",iom+compteur,om,Mmv(1)
   donnees(1,iom+compteur)=om
   donnees(2:7,iom+compteur)=Mmv
  enddo
 compteur=compteur+nom2
 endif

 ommin=opp(3)
 ommax=4*opp(3)
 dom=(ommax-ommin)/(nom1+1)
 write(6,*)"ommim,ommax,dom=",ommin,ommax,dom
 do iom=1,nom3
  om=ommin+dom*iom
  call mat_pairfield(om,0.0_qp,det,Mm,Gg)
  Mmv=(/real(Mm(1,1)),real(Mm(2,2)),real(Mm(1,2)),imag(Mm(1,1)),imag(Mm(2,2)),imag(Mm(1,2))/)
  write(6,*)"iom+compteur,om,Mmv(1)=",iom+compteur,om,Mmv(1)
  donnees(1,iom+compteur)=om
  donnees(2:7,iom+compteur)=Mmv
 enddo
 compteur=compteur+nom3

 ommin=4*opp(3)
 ommax=bmax
 ymin=1/ommax**(1.5_qp)
 ymax=1/ommin**(1.5_qp)
 dy=(ymax-ymin)/(nominf+1)
 write(6,*)"ommim,ommax,ymin,ymax,dy=",ommin,ommax,ymin,ymax,dy
 do iom=1,nominf
  y=ymin+dy*iom
  om=1.0_qp/y**(2.0_qp/3.0_qp)
  call mat_pairfield(om,0.0_qp,det,Mm,Gg)
  Mmv=(/real(Mm(1,1)),real(Mm(2,2)),real(Mm(1,2)),imag(Mm(1,1)),imag(Mm(2,2)),imag(Mm(1,2))/)
  write(6,*)"iom+compteur,om,Mmv(1)=",iom+compteur,om,Mmv(1)
  donnees(1,iom+compteur)=om
  donnees(2:7,iom+compteur)=Mmv
 enddo
 compteur=compteur+nominf

  inquire(iolength=nn)donnees(1:7,1:compteur)
  write(6,*)"Taille de l’enregistrement en octets:",nn
  open(11,file=trim(fichier)//".dat",ACTION="WRITE",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   write(11,REC=ixq+1)donnees(1:7,1:compteur)
  close(11)
  open(12,file=trim(fichier)//"grilleq.dat",POSITION="APPEND")
   write(12,*)ixq,xq,compteur,nn
  close(12)
  write(6,*)

enddo
END PROGRAM pointsM
