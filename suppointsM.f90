PROGRAM suppointsM
USE nrtype
USE vars
USE dspec
USE OMP_LIB
IMPLICIT NONE
COMPLEX(QPC) Gg(1:2,1:2),Mm(1:2,1:2),det
CHARACTER(len=90) fichier
REAL(QP) xqmin,xqmax,xq2,dxq
REAL(QP) om,ommin,ommax,dom,dom1,dom2,dom3,dom2p,y,ymax,dy,xq1
REAL(QP) Mmv(1:6)
REAL(QP),ALLOCATABLE, DIMENSION(:,:,:) :: donnees
INTEGER, ALLOCATABLE, DIMENSION(:)   :: fini
INTEGER  ixq,ixqdep,ixqfin,nq
INTEGER  iom,nom,icb,nmax
INTEGER  nn

open(10,file='suppointsM.inp')
 read(10,*)x0
 read(10,*)xqmin
 read(10,*)xqmax
 read(10,*)ixqdep
 read(10,*)ixqfin
 read(10,*)nq
 read(10,*)nom
 read(10,*)fichier
close(10)

bla1=.TRUE.
bla1=.FALSE.


EPSpp=1.0e-8_qp
x0crit=0.0_qp
beta=bidon
temperaturenulle=.TRUE.

xq2   =2*sqrt(x0)
call calcxqjoin
xq1 = xqjoin

write(6,*)'xq2='    ,xq2
write(6,*)'nom='    ,nom
write(6,*)'nq ='    ,nq
write(6,*)'ixqdep=' ,ixqdep
write(6,*)'ixqfin=' ,ixqfin
write(6,*)'xq1 ='   ,xq1
write(6,*)'fichier=',trim(fichier)

allocate(donnees(1:7,1:4*nom,1:nq))
allocate(fini(1:nq))
inquire(iolength=nn)donnees(:,:,1)
write(6,*)"Taille de l’enregistrement en octets:",nn

!call system("rm "//trim(fichier)//".dat")
!call system("rm "//trim(fichier)//"grilleq.dat")
!call system("rm "//trim(fichier)//".info")

!open(12,file=trim(fichier)//".info")
! write(12,*)"!x0,xq2,nq,nom,nn"
! write(12,*)  x0,xq2,nq,nom,nn
!close(12)
!
!open(13,file=trim(fichier)//"grilleq.dat")
! write(13,*)"!ixq,xq"
!close(13)
!
dxq=(xqmax-xqmin)/nq
write(6,*)'--------------------'
write(6,*)
write(6,*)"dxq=",dxq
donnees(:,:,:)=0.0_qp
fini(:)=0
!$OMP PARALLEL DO PRIVATE(om,det,Mm,Gg,Mmv,ixq,iom,dom,dom1,dom2,dom3,nmax,y,dy) SHARED(fini,donnees) SCHEDULE(DYNAMIC)
do icb=(ixqdep-1)*4*nom+1,ixqfin*4*nom
 ixq=icb/(4*nom)+1
 iom=MODULO(icb,4*nom)
 if(iom==0)then 
  iom=4*nom
  ixq=ixq-1
 endif
 xq=xqmin+dxq*(ixq-0.5_qp)
 call oangpp
 if(xq>xq2)then
  nmax=nom
 elseif(xq>xq1)then
  nmax=3*nom
 else
  nmax=4*nom
 endif


 dom1 =(opp(2)-opp(1))/2
 dom2 =min(0.3_qp,(opp(3)-opp(2)))
 dom3 =0.05_qp*opp(3)
 dy=sqrt(dom1)/nom
 if(iom==1)then
  write(6,*)'--------------------'
  write(6,*)
  write(6,*)"xq,ixq=",xq,ixq
  write(6,*)"opp=",opp(1:3)
  write(6,*) 
  write(6,*)"départ,1er arrêt,2e arrêt,fin,2e dépar,fint=",opp(1),opp(2)-dom1,opp(2),opp(2)+dom2,opp(3),opp(3)+dom3
 endif
 if(iom.LE.nom)then
  if(xq>xq2) cycle
  dom=dom1/nom
  om=opp(1)+dom*(iom-0.5_qp)
 elseif(iom.LE.2*nom)then
  if(xq>xq2) cycle
  y=dy*(iom-nom-0.5_qp)
  om=opp(2)-y**2
 elseif(iom.LE.3*nom)then
  if(xq>xq1) cycle
  dom=dom2/nom
  om=opp(2)+dom*(iom-2*nom-0.5_qp)
 else
  dom=dom3/nom
  om=opp(3)+dom*(iom-3*nom-0.5_qp)
 endif
 call mat_pairfield(om,0.0_qp,det,Mm,Gg)
 Mmv=(/real(Mm(1,1)),real(Mm(2,2)),real(Mm(1,2)),imag(Mm(1,1)),imag(Mm(2,2)),imag(Mm(1,2))/)
 donnees(1  ,iom,ixq)=om
 donnees(2:7,iom,ixq)=Mmv
 !$OMP CRITICAL
 fini(ixq)=fini(ixq)+1
 !$OMP END CRITICAL
 write(6,*)"ixq,iom,om,Mmv(1),pt=",ixq,iom,om,Mmv(1),"  ",fini(ixq)," sur ",nmax
 if(fini(ixq)==nmax)then
  write(6,*)ixq,"est fini"
  open(11,file=trim(fichier)//".dat",ACTION="WRITE",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   write(11,REC=ixq)donnees(:,:,ixq)
  close(11)
  open(13,file=trim(fichier)//"grilleq.dat",POSITION="APPEND")
   write(13,*)ixq,xq
   write(6,*)ixq,xq
  close(13)
  write(6,*)
 endif
enddo
!$OMP END PARALLEL DO

END PROGRAM suppointsM
