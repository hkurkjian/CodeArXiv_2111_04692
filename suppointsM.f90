PROGRAM suppointsM
USE nrtype
USE vars
USE dspec
USE modpol
USE OMP_LIB
IMPLICIT NONE
COMPLEX(QPC) Gg(1:2,1:2),Mm(1:2,1:2),det
CHARACTER(len=90) fichier
REAL(QP) xqmin,xqmax,xq1,xq2,dxq
REAL(QP) om,dom,dom1,dom2,dom3,dom2p,dom3p,y,dymax,dy,bmax
REAL(QP) Mmv(1:6)
REAL(QP),ALLOCATABLE, DIMENSION(:,:,:) :: donnees
REAL(QP), ALLOCATABLE, DIMENSION(:) :: xqfen,bornes 
INTEGER, ALLOCATABLE, DIMENSION(:)   :: tfen,fenom,passe,fini
INTEGER  ixq,ixqdep,ixqfin,nq
INTEGER  iom,nom1,nom2,icb,nmax
INTEGER  ifen,nfen,nfen1,nfen2,nn,npoints
LOGICAL  nvofich

open(10,file='suppointsM.inp')
 read(10,*)x0
 read(10,*)xqmin
 read(10,*)xqmax
 read(10,*)ixqdep
 read(10,*)ixqfin
 read(10,*)nq
 read(10,*)nfen1
 read(10,*)nfen2
 read(10,*)nom1
 read(10,*)nom2
 read(10,*)bmax
 read(10,*)nvofich
 read(10,*)fichier
close(10)

bla1=.TRUE.
bla1=.FALSE.

EPSpp=1.0e-8_qp
x0crit=1.0_qp
beta=bidon
temperaturenulle=.TRUE.

call calcxqjoin
xq1 = xqjoin
xq2   =2*sqrt(x0)

npoints=nfen1*nom1+nfen2*nom2
nfen=nfen1+nfen2

write(6,*)'nom1,nom2,nfen1,nfen2,npoints='    ,nom1,nom2,nfen1,nfen2,npoints
write(6,*)'ixqdep,ixqfin;nq=' ,ixqdep,ixqfin,"; ",nq
write(6,*)'xq1,xq2 ='   ,xq1,xq2
write(6,*)'bmax='   ,bmax
write(6,*)'nvofich=',nvofich
write(6,*)'fichier=',trim(fichier)

allocate(bornes(1:nfen+1),fenom(1:nfen+1),passe(1:nfen+1))
allocate(xqfen(1:nfen),tfen(1:nfen))
allocate(donnees(1:7,1:npoints,1:nq))
allocate(fini(1:nq))

tfen=(/1,3,2,1,3,2,1,3,4/)
xqfen=(/xq2,xq2,xq2,xq1,xq1,xq1,1.0e100_qp,1.0e100_qp,1.0e100_qp/)
fenom=(/1,nom1,nom2,nom1,nom1,nom2,nom1,nom1,nom2,nom2-1/)
do ifen=1,nfen+1
 passe(ifen)=sum(fenom(1:ifen))
enddo

inquire(iolength=nn)donnees(:,:,1)
write(6,*)"Taille de l’enregistrement en octets:",nn

if(nvofich)then
   open(12,file=trim(fichier)//".info")
    write(12,*)"!x0,nq,nom1,nom2,nn"
    write(12,*)  x0,nq,nom1,nom2,nn
   close(12)
   
   open(13,file=trim(fichier)//"grilleq.dat")
    write(13,*)"!ixq,xq,opp(1:3)"
   close(13)

   call system("rm "//trim(fichier)//".dat")
endif

dxq=(xqmax-xqmin)/nq
write(6,*)'--------------------'
write(6,*)
write(6,*)"dxq=",dxq
donnees(:,:,:)=0.0_qp
fini(:)=0
!$OMP  PARALLEL DO &
!$OMP& PRIVATE(om,det,Mm,Gg,Mmv,icb,ixq,iom,ifen,dom,dom1,dom2,dom3,bornes,nmax,y,dy,dymax) &
!$OMP& SHARED(fini,donnees) SCHEDULE(DYNAMIC)
do icb=(ixqdep-1)*npoints+1,ixqfin*npoints
 ixq=icb/npoints+1
 iom=MODULO(icb,npoints)
 if(iom==0)then 
  iom=npoints
  ixq=ixq-1
 endif
 xq=xqmin+dxq*(ixq-0.5_qp)
 call oangpp

 if(xq>xq2)then
  nmax=npoints-2*(2*nom1+nom2)
 elseif(xq>xq1)then
  nmax=npoints-(2*nom1+nom2)
 else
  nmax=npoints
 endif


 dom1  =min((opp(2)-opp(1))/10,0.1_qp)
 if((opp(3)-opp(1)).LE.0.1_qp) dom1=(opp(2)-opp(1))/4
 dom2  =min(0.1_qp,(opp(3)-opp(2))/4)
 dom3  =min(0.05_qp*opp(3),(opp(3)-opp(2))/4)
 dom3p =0.05_qp*opp(3)

 bornes=(/opp(1),opp(1)+dom1,opp(2)-dom1,opp(2),opp(2)+dom2,opp(3)-dom3,opp(3),opp(3)+dom3p,4*opp(3),bmax/)

 if(iom==1)then
  write(6,*)'--------------------'
  write(6,*)
  write(6,*)"xq,ixq=",xq,ixq
  write(6,*)"opp=",opp(1:3)
  write(6,*) 
  write(6,*)"découpe intervalle II:" ,bornes(1:4)
  write(6,*)"découpe intervalle III:",bornes(4:7)
  write(6,*)"découpe intervalle IV:" ,bornes(7:10)
 endif

 call locatei(passe,iom,ifen)
 if(xq>xqfen(ifen)) cycle

 dy=sqrt(bornes(ifen+1)-bornes(ifen))/nom1
 dymax=(1/sqrt(bornes(ifen))-1/sqrt(bornes(ifen+1)))/nom2
 dom=(bornes(ifen+1)-bornes(ifen))/nom2

 y=dy*(iom-passe(ifen)+1-0.5_qp)
 if(tfen(ifen)==4) y=1/sqrt(bornes(ifen+1))+dymax*(iom-passe(ifen)+1-0.5_qp)

 if(tfen(ifen)==1) om=bornes(ifen)  +y**2
 if(tfen(ifen)==2) om=bornes(ifen+1)-y**2
 if(tfen(ifen)==3) om=bornes(ifen)  +dom*(iom-passe(ifen)+1-0.5_qp)
 if(tfen(ifen)==4) om=1.0_qp/y**2
 
 if(om>200000.0_qp)then
  call mat_pairfield_gom0(om,0.0_qp,det,Mm,Gg)
 else
  call mat_pairfield     (om,0.0_qp,det,Mm,Gg)
 endif
 Mmv=(/real(Mm(1,1)),real(Mm(2,2)),real(Mm(1,2)),imag(Mm(1,1)),imag(Mm(2,2)),imag(Mm(1,2))/)
 donnees(1  ,iom,ixq)=om
 donnees(2:7,iom,ixq)=Mmv
 !$OMP CRITICAL
 fini(ixq)=fini(ixq)+1
 if(fini(ixq)==nmax)then
  write(6,*)ixq,"est fini"
  open(11,file=trim(fichier)//".dat",ACTION="WRITE",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   write(11,REC=ixq)donnees(:,:,ixq)
  close(11)
  open(13,file=trim(fichier)//"grilleq.dat",POSITION="APPEND")
   write(13,*)ixq,xq,opp(1:3)
   write(6,*)ixq,xq
  close(13)
  write(6,*)
 endif
 !$OMP END CRITICAL
 if(iom==passe(ifen)) write(6,*)"bornes:",bornes(ifen),bornes(ifen+1)
 write(6,*)"ixq,xq,ifen,iom,om,Mmv(1),pt=",ixq,xq,ifen,iom,om,Mmv(1),"  ",fini(ixq)," sur ",nmax
enddo
!$OMP END PARALLEL DO

END PROGRAM suppointsM
