MODULE estM
USE nrtype
USE vars
USE dspec
USE modpol
IMPLICIT NONE
REAL(QP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: donnees,donneesgq
REAL(QP), DIMENSION(:,:),   ALLOCATABLE, TARGET :: vecq,vecqgq
REAL(QP) bmax,qpetit,omgrand,gq
INTEGER, TARGET :: nom1,nom2,nq
INTEGER, TARGET :: nom1gq,nom2gq,nqgq
LOGICAL blaM,blaMM,blaerr,blablaerr
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE estmat_pairfield(om0,e,det,M,fr)
COMPLEX(QPC), DIMENSION(1:2,1:2), INTENT(OUT) :: M,fr
COMPLEX(QPC), INTENT(OUT) :: det
REAL(QP), INTENT(IN)  :: om0,e

REAL(QP),  DIMENSION(1:6) :: Mmv
REAL(QP) :: q

q=xq
Mmv=interpolM_recerr(q,om0)
M(1,1)=cmplx(Mmv(1),Mmv(4),kind=qpc)
M(2,2)=cmplx(Mmv(2),Mmv(5),kind=qpc)
M(1,2)=cmplx(Mmv(3),Mmv(6),kind=qpc)
M(2,1)=M(1,2)

det=M(1,1)*M(2,2)-M(1,2)**2.0_qp

fr(1,1)=M(2,2)/det
fr(2,2)=M(1,1)/det
fr(1,2)=-M(1,2)/det
fr(2,1)=fr(1,2)
END SUBROUTINE estmat_pairfield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION interpolM_recerr(q,om)
USE modsim
REAL(QP), INTENT(IN) :: om,q
REAL(QP), DIMENSION(1:6) :: interpolM_recerr

LOGICAL :: err
REAL(QP),  DIMENSION(1:6) :: dinterpolM
COMPLEX(QPC) Gg(1:2,1:2),Mm(1:2,1:2),det
INTEGER iMm

interpolM_recerr=interpoloM(q,om,dinterpolM,err)

if(blaMM.OR.blablaerr)then
 do iMm=1,6
  write(6,*)"om,q,interpolM,dinterpolM=",om,q,interpolM_recerr(iMm),dinterpolM(iMm)
 enddo
 write(6,*)
endif

if(err)then
!if(.TRUE.)then
   if(blaerr) write(6,*)"Erreur sèche,q,om=",q,om
   if(q<qpetit)then
    if(blaerr) write(6,*)"Appel de mat_pairfield_pttq"
    call mat_pairfield_pttq(om,0.0_qp,det,Mm,Gg)
   elseif(om>omgrand)then
    if(blaerr) write(6,*)"Appel de mat_pairfield_gom0"
    call mat_pairfield_gom0(om,0.0_qp,det,Mm,Gg)
   else
    if(blaerr) write(6,*)"Appel de mat_pairfield"
    call mat_pairfield(om,0.0_qp,det,Mm,Gg)
   endif
   interpolM_recerr(1:6)=(/real(Mm(1,1)),real(Mm(2,2)),real(Mm(1,2)),imag(Mm(1,1)),imag(Mm(2,2)),imag(Mm(1,2))/)

   if(blaerr)then
    do iMm=1,6
     write(6,*)"om,q,interpol,dinterpol(récupéré)=",om,xq,interpolM_recerr(iMm),dinterpolM(iMm)
    enddo
    write(6,*)
   endif
endif

END FUNCTION interpolM_recerr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION interpoloM(q,om,dinterpolom,err)
USE modsim
REAL(QP), INTENT(IN) :: om,q
LOGICAL,  INTENT(INOUT) :: err
REAL(QP), DIMENSION(1:6), INTENT(OUT) :: dinterpolom
REAL(QP), DIMENSION(1:6) :: interpolom

INTEGER posq,posom,fen,nom,npasse,posc
INTEGER iposq,iMm,no1,no2,nqq
REAL(QP) ptsq(1:4),ptsom(1:4,1:4),ptsy(1:4,1:4),ptsM(1:3,1:4,1:4)
REAL(QP) y,oppref(1:4,1:3),oppr(1:3)
REAL(QP), DIMENSION(:,:), ALLOCATABLE :: sgrille,sgrilley
REAL(QP), DIMENSION(:,:,:), POINTER :: don
REAL(QP), DIMENSION(:,:), POINTER   :: vq
LOGICAL :: err2

err   =.FALSE.
blablaerr=.FALSE.

ptsM (:,:,:)=0.0_qp
ptsom(:,:)  =0.0_qp
ptsq (:)    =0.0_qp

if(q<gq)then
 don => donnees
 vq  => vecq
 no1=nom1
 no2=nom2
 nqq=nq
else
 don => donneesgq
 vq  => vecqgq
 no1=nom1gq
 no2=nom2gq
 nqq=nqgq
endif

xq=q
call oangpp
call locate(vq(:,1),q,posq)
call locate(opp(1:3),om,fen)

if(fen==1)then
  nom=2*no1+1*no2
  npasse=0
elseif(fen==2)then 
  nom=2*no1+1*no2
  npasse=2*no1+1*no2
elseif(fen==3)then
  nom=1*no1+2*no2
  npasse=4*no1+2*no2
endif
allocate(sgrille (1:4,1:nom),sgrilley(1:4,1:nom))

if((posq<1).OR.(posq.GE.nqq))then
  err=.TRUE.
  return
endif

if((posq<2).OR.(posq.GE.nqq-1))then
  err2=.TRUE.
  if(blaerr) blablaerr=.TRUE.
  ptsq  (2:3)    =vq(posq:posq+1,1)
  oppref(2:3,1:3)=vq(posq:posq+1,2:4)
else
  ptsq  (1:4)    =vq(posq-1:posq+2,1)
  oppref(1:4,1:3)=vq(posq-1:posq+2,2:4)
endif

do iposq=posq-1,posq+2
 posc=iposq-posq+2
 oppr=oppref(posc,:)
 if((iposq<1).OR.(iposq>nqq)) cycle

 sgrille=don(1:4,npasse+1:npasse+nom,iposq)

 if(om<(opp(1)+opp(2))/2)then 
  if((opp(2)-opp(1))<0.01_qp)then
   y=(om-opp(1))/(opp(2)-opp(1))
   sgrilley(1,:)=(sgrille(1,:)-oppr(1))/(oppr(2)-oppr(1))
  else
   y=om-opp(1)
   sgrilley(1,:)=sgrille(1,:)-oppr(1)
  endif
 elseif(om<opp(2))then 
  y=sqrt(opp(2)-om)
  sgrilley(1,:)  =sqrt(oppr(2)-sgrille(1,:))
 elseif(om<opp(3))then
  y=om-opp(2)
  sgrilley(1,:)=sgrille(1,:)-oppr(2)
 elseif(om<4*opp(3))then
  y=om-opp(3)
  sgrilley(1,:)=sgrille(1,:)-oppr(3)
 else
  y=1.0_qp/om**(0.5_qp)
  sgrilley(1,:)  =1.0_qp/sgrille(1,:)**(0.5_qp)
 endif
 call locate(sgrilley(1,:),y,posom)

 if((posom<1).OR.(posom.GE.nom))then
  err2=.TRUE.
  if(blaerr) write(6,*)"*********************************************"
  if(blaerr) write(6,*)
  if(blaerr) write(6,*)"om hors grille à posq,iposq=",posq,iposq
  if(blaerr) write(6,*)"om,posom(nom),pts extrêmes=",om,posom,"(",nom,")",sgrille(1,1),sgrille(1,size(sgrille,dim=2))
  if(blaerr) write(6,*)
  if(blaerr) write(6,*)"*********************************************"
  if(blaerr) blablaerr=.TRUE.
  if((iposq==posq).OR.(iposq==posq+1)) err=.TRUE.
  cycle
 endif

 if((posom<2).OR.(posom.GE.nom-1))then
  err2=.TRUE.
  ptsM(1:3,2:3,posc)=sgrille (2:4,posom:posom+1)
  ptsom   (2:3,posc)=sgrille (1,  posom:posom+1)
  ptsy    (2:3,posc)=sgrilley(1,  posom:posom+1)
 else
  ptsM(1:3,1:4,posc)=sgrille (2:4,posom-1:posom+2)
  ptsom   (1:4,posc)=sgrille (1,  posom-1:posom+2)
  ptsy    (1:4,posc)=sgrilley(1,  posom-1:posom+2)
 endif
enddo

if(err.AND.(blaerr)) blablaerr=.TRUE.

if(err)then
 if(blaerr) write(6,*)"Erreur de type 0"
 if(blaerr) write(6,*)"*********************************************"
endif


if(blaMM.OR.blablaerr)then
 write(6,*)
 write(6,*)"------------- interpolom -------------------"
 write(6,*)
 write(6,*)"q,om,y=",q,om,y
 write(6,*)"opp=",opp(1:3)
 write(6,*)
 write(6,*)"posq,  vecq(posq)  =",posq,  vq(posq,1)
 write(6,*)"posq+1,vecq(posq+1)=",posq+1,vq(posq+1,1)
 write(6,*)
 write(6,*)"--------------- Points retenus ---------------"
 write(6,*)
 write(6,*)"ptsq=",ptsq
 write(6,*)
 write(6,*)"ptsom(1,:)=",ptsom(1,:)
 write(6,*)"ptsom(2,:)=",ptsom(2,:)
 write(6,*)"ptsom(3,:)=",ptsom(3,:)
 write(6,*)"ptsom(4,:)=",ptsom(4,:)
 write(6,*)
 write(6,*)"ptsM(1,1,:)=",ptsM(1,1,:)
 write(6,*)"ptsM(1,2,:)=",ptsM(1,2,:)
 write(6,*)"ptsM(1,3,:)=",ptsM(1,3,:)
 write(6,*)"ptsM(1,4,:)=",ptsM(1,4,:)
 write(6,*)
endif

if(err) return

do iMm=1,3
  if(err2)then
   call polin2(ptsq(2:3),ptsy(2:3,2:3),ptsM(iMm,2:3,2:3),q,y,interpolom(iMm),dinterpolom(iMm))
  else
   call polin2(ptsq,ptsy,ptsM(iMm,:,:),q,y,interpolom(iMm),dinterpolom(iMm))
  endif
  if(abs(dinterpolom(iMm)/interpolom(iMm))>0.01_qp)then
    if(blaerr) blablaerr=.TRUE.
    if(blaerr) write(6,*)"grosse erreur d’interpolation iMm=",iMm
  endif
  if(abs(dinterpolom(iMm)/interpolom(iMm))>0.2_qp)then
    err=.TRUE.
    if(blaerr) write(6,*)"trop grosse erreur d’interpolation iMm=",iMm
  endif
enddo

interpolom(4)=-PI*rhopp(1.5_qp,om,-1.0_qp)
interpolom(5)=-PI*rhopp(2.5_qp,om,-1.0_qp)
interpolom(6)=-PI*rhopp(3.5_qp,om,-1.0_qp)

dinterpolom(4:6)=EPSpp

deallocate(sgrille)
deallocate(sgrilley)

END FUNCTION interpoloM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE load_data(fich)
INTEGER ixq,ixqbis,nn,nngq,npoints
CHARACTER(len=90), INTENT(IN) :: fich
REAL(QP) xqlec,qbidon
open(212,file=trim(fich)//".info")
 read(212,*)
 read(212,*) x0,nq,qbidon,qbidon,nom1,nom2,nn
 if(blaM)then
  write(6,*)"----------------------------------------"
  write(6,*)
  write(6,*)"Chargement du tableau de données primaire"
  write(6,*)
  write(6,*)"x0=",x0
  write(6,*)"nq=",nq
  write(6,*)"nn=",nn
  write(6,*)
 endif
close(212)

npoints=5*nom1+4*nom2
allocate(donnees(1:7,1:npoints,1:nq))
allocate(vecq(1:nq,1:4))

donnees(1,:,:)=1.0e50_qp
vecq(:,:)=1.0e50_qp

open(213,file=trim(fich)//"grilleq.dat")
read(213,*)
do ixq=1,nq
  read(213,*)ixqbis,vecq(ixq,:)
  open(211,file=trim(fich)//".dat",ACTION="READ",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   read(211,REC=ixq)donnees(1:7,1:npoints,ixq)
   donnees(:,1*nom1+1*nom2+1:2*nom1+1*nom2,ixq)=donnees(:,2*nom1+1*nom2:1*nom1+1*nom2+1:-1,ixq)
   donnees(:,3*nom1+2*nom2+1:4*nom1+2*nom2,ixq)=donnees(:,4*nom1+2*nom2:3*nom1+2*nom2+1:-1,ixq)
   donnees(:,5*nom1+3*nom2+1:5*nom1+4*nom2,ixq)=donnees(:,5*nom1+4*nom2:5*nom1+3*nom2+1:-1,ixq)
   if(blaMM)then 
    write(6,*)"ixq,donnees(1:3,75,ixq)=",ixq,vecq(ixq,1),vecq(ixq,3),donnees(1:3,3450,ixq)
   endif
  close(211)
enddo
close(213)

if(blaM)then
 write(6,*)"Chargement terminé"
 write(6,*)
 write(6,*)"----------------------------------------"
endif
! ************************************************************** !
open(212,file=trim(fich)//"gq.info")
 read(212,*)
 read(212,*) x0,nqgq,gq,qbidon,nom1gq,nom2gq,nngq
 if(blaM)then
  write(6,*)"----------------------------------------"
  write(6,*)
  write(6,*)"Chargement du tableau de données grands q"
  write(6,*)
  write(6,*)"x0=",x0
  write(6,*)"nqgq=",nqgq
  write(6,*)"nngq=",nngq
  write(6,*)
 endif
close(212)

npoints=5*nom1gq+4*nom2gq
allocate(donneesgq(1:7,1:npoints,1:nqgq))
allocate(vecqgq(1:nqgq,1:4))

donneesgq(1,:,:)=1.0e50_qp
vecqgq(:,:)=1.0e50_qp

open(213,file=trim(fich)//"gqgrilleq.dat")
read(213,*)
do ixq=1,nqgq
  read(213,*)ixqbis,vecqgq(ixq,:)
  open(211,file=trim(fich)//"gq.dat",ACTION="READ",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nngq)
   read(211,REC=ixq)donneesgq(1:7,1:npoints,ixq)
   donneesgq(:,5*nom1gq+3*nom2gq+1:5*nom1gq+4*nom2gq,ixq)=donneesgq(:,5*nom1gq+4*nom2gq:5*nom1gq+3*nom2gq+1:-1,ixq)
   if(blaMM)then 
    write(6,*)"ixq,q,opp(2),donneesgq(1:3,75,ixq)=",ixq,vecqgq(ixq,1),vecqgq(ixq,3),donneesgq(1:3,2000,ixq)
   endif
  close(211)
enddo
close(213)



if(blaM)then
 write(6,*)"Chargement gq terminé"
 write(6,*)
 write(6,*)"----------------------------------------"
endif
END SUBROUTINE load_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE unload_data
deallocate(donnees)
deallocate(vecq)
END SUBROUTINE unload_data
END MODULE estM
