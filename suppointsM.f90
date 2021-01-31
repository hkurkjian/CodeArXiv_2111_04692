PROGRAM suppointsM
USE nrtype
USE vars
USE dspec
IMPLICIT NONE
COMPLEX(QPC) Gg(1:2,1:2),Mm(1:2,1:2),det
CHARACTER(len=90) fichier
REAL(QP) xqmin,xq2,dxq
REAL(QP) om,ommin,ommax,dom,dom2,dom2p,y,ymax,dy,xq1
REAL(QP) Mmv(1:6)
REAL(QP),ALLOCATABLE, DIMENSION(:,:) :: donnees
INTEGER  ixq,nq
INTEGER  iom,nom
INTEGER  nn

open(10,file='suppointsM.inp')
 read(10,*)x0
 read(10,*)xqmin
 read(10,*)nq
 read(10,*)nom
 read(10,*)fichier
close(10)

bla1=.TRUE.
bla1=.FALSE.

EPSpp=1.0e-8_qp
x0crit=4.0_qp
beta=bidon
temperaturenulle=.TRUE.

xq2   =2*sqrt(x0)
call calcxqjoin
xq1 = xqjoin

write(6,*)'xq2='    ,xq2
write(6,*)'nom='    ,nom
write(6,*)'nq ='    ,nq
write(6,*)'xq1 ='   ,xq1
write(6,*)'fichier=',trim(fichier)

allocate(donnees(1:7,1:2*nom))
inquire(iolength=nn)donnees
write(6,*)"Taille de l’enregistrement en octets:",nn

call system("rm "//trim(fichier)//".dat")
call system("rm "//trim(fichier)//"grilleq.dat")
call system("rm "//trim(fichier)//".info")

open(12,file=trim(fichier)//".info")
 write(12,*)"!x0,xq2,nq,nom,nn"
 write(12,*)  x0,xq2,nq,nom,nn
close(12)

open(13,file=trim(fichier)//"grilleq.dat")
 write(13,*)"!ixq,xq"
close(13)

dxq=(xq2-xqmin)/nq
write(6,*)'--------------------'
write(6,*)
write(6,*)"dxq=",dxq
do ixq=1,nq
 xq=xqmin+dxq*(ixq-0.5_qp)
 call oangpp
 write(6,*)'--------------------'
 write(6,*)
 write(6,*)"xq,ixq=",xq,ixq
 write(6,*)"opp=",opp(1:3)
 write(6,*) 

 donnees(:,:)=0.0_qp

 dom2 =(opp(2)-opp(1))/4
 if(xq<xq1) dom2p=min((opp(2)-opp(1))/4,(opp(3)-opp(2))/4)
 if(xq>xq1) dom2p=(opp(2)-opp(1))/4
 dom=dom2p/nom
 ommin=opp(2)
 ymax=sqrt(dom2)
 dy=ymax/nom
 write(6,*)"ommim,opp(2),ommax=",opp(2)-dom2,opp(2),opp(2)+dom2p
!$OMP PARALLEL DO PRIVATE(om,det,Mm,Gg,Mmv) SCHEDULE(DYNAMIC)
 do iom=1,2*nom
  if(iom.LE.nom)then
   y=dy*(iom-0.5_qp)
   om=opp(2)-y**2
  else
   om=ommin+dom*(iom-nom-0.5_qp)
  endif
  call mat_pairfield(om,0.0_qp,det,Mm,Gg)
  Mmv=(/real(Mm(1,1)),real(Mm(2,2)),real(Mm(1,2)),imag(Mm(1,1)),imag(Mm(2,2)),imag(Mm(1,2))/)
  write(6,*)"iom,om-opp(2),Mmv(1)=",iom,om-opp(2),Mmv(1)
  donnees(1  ,iom)=om
  donnees(2:7,iom)=Mmv
 enddo
!$OMP END PARALLEL DO

 open(11,file=trim(fichier)//".dat",ACTION="WRITE",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
  write(11,REC=ixq)donnees
 close(11)
 open(13,file=trim(fichier)//"grilleq.dat",POSITION="APPEND")
  write(13,*)ixq,xq
 close(13)
 write(6,*)
enddo
END PROGRAM suppointsM
