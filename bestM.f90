MODULE bestM
USE nrtype
USE vars
USE dspec
USE modpol
IMPLICIT NONE
REAL(QP), DIMENSION(:,:,:), ALLOCATABLE :: donnees,donnees_sup,donnees_sup2
REAL(QP), DIMENSION(:,:),   ALLOCATABLE :: vecq,vecq_sup,vecq_sup2
REAL(QP) qsep(1:5),bmax,qpetit,xq2
INTEGER nqfen(1:4),nomfen(1:8),nbidon,nqsup,nomsup,nqsup2,nomsup2
LOGICAL blaM,blaerr,blablaerr
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

LOGICAL :: err(0:2),recom2
REAL(QP),  DIMENSION(1:6) :: dinterpolM
INTEGER iMm
COMPLEX(QPC) Gg(1:2,1:2),Mm(1:2,1:2),det
REAL(QP) dom2,dom2p

interpolM_recerr=interpolM(q,om,dinterpolM,err)

if(blaM.OR.blablaerr)then
 do iMm=1,6
  write(6,*)"om,q,interpolM,dinterpolM=",om,q,interpolM_recerr(iMm),dinterpolM(iMm)
 enddo
 write(6,*)
endif

if(err(0).OR.err(1))then
 if(blaerr) write(6,*)"Tentative de récupération par interpolom2"
 recom2=.FALSE.

 interpolM_recerr=interpolom2(q,om,dinterpolM,err)
 if(.NOT.(err(0).OR.err(1)))then 
  if(blaerr) write(6,*)"Erreur récupérée par intom2"
 else
  interpolM_recerr=interpolom3(q,om,dinterpolM,err)
  if(.NOT.(err(0).OR.err(1)))then 
   if(blaerr) write(6,*)"Erreur récupérée par intom3"
  else
   if(blaerr) write(6,*)"Erreur sèche,q,om=",q,om
   if(q<qpetit)then
    call mat_pairfield_pttq(om,0.0_qp,det,Mm,Gg)
   else
    call mat_pairfield(om,0.0_qp,det,Mm,Gg)
   endif
   interpolM_recerr(1:6)=(/real(Mm(1,1)),real(Mm(2,2)),real(Mm(1,2)),imag(Mm(1,1)),imag(Mm(2,2)),imag(Mm(1,2))/)
  endif
 endif

 if(blaerr)then
  do iMm=1,6
   write(6,*)"om,q,interpol,dinterpol(récupéré)=",om,xq,interpolM_recerr(iMm),dinterpolM(iMm)
  enddo
  write(6,*)
 endif
endif

END FUNCTION interpolM_recerr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION interpolM(q,om,dinterpolM,err)
USE modsim
REAL(QP), INTENT(IN) :: om,q
LOGICAL,  INTENT(OUT) :: err(0:2)
REAL(QP), DIMENSION(1:6), INTENT(OUT) :: dinterpolM
REAL(QP), DIMENSION(1:6) :: interpolM

INTEGER posq,posom,posy
INTEGER nom,decalageom,fenom
INTEGER nqred(1:3),nomred(1:4)
INTEGER iposq,iMm
REAL(QP) y,qsepred(1:4),omsepred(1:5)
REAL(QP) ptsq(1:4),ptsy(1:4,1:4),ptsM(1:3,1:4,1:4)
REAL(QP) opp1,opp2,opp3,oppref(1:4)
REAL(QP), DIMENSION(:,:),   ALLOCATABLE :: sousgrille,sousgrilley

err(:)   =.FALSE.
blablaerr=.FALSE.

ptsM (:,:,:)=0.0_qp
ptsy (:,:)  =0.0_qp
ptsq (:)    =0.0_qp

nqred =(/nqfen(1)+nqfen(2),nqfen(3),nqfen(4)/)
nomred=(/sum(nomfen(1:2)),sum(nomfen(3:5)),sum(nomfen(6:7)),nomfen(8)/)

xq=q
call oangpp

qsepred =(/qsep(1),qsep(3),qsep(4),qsep(5)/)
omsepred=(/opp(1),opp(2),opp(3),4*opp(3),bmax/)

call locate(omsepred,om,fenom)

if(fenom==0)then
 if(blaerr) write(6,*)"om n'est pas dans le continuum de brisure de paires"
 if(blaerr) blablaerr=.TRUE.
 return
elseif(fenom==5)then
 if(blaerr) write(6,*)"om dépasse bmax"
 if(blaerr) blablaerr=.TRUE.
 return
endif

nom=nomred(fenom)

allocate(sousgrille (1:4,1:nom))
allocate(sousgrilley(1:4,1:nom))

call locate(vecq(:,1),q,posq)

if((posq<1).OR.(posq.GE.sum(nqred)))then
  err(0:1)=.TRUE.
  if(blaerr) write(6,*)"Erreur de type 0 venant de q"
  return
endif

if((posq<2).OR.(posq.GE.sum(nqred)-1))then
  err(2)=.TRUE.
  if(blaerr) blablaerr=.TRUE.
  ptsq(2:3)=vecq(posq:posq+1,1)
else
  ptsq(1:4)=vecq(posq-1:posq+2,1)
endif

do iposq=posq-1,posq+2
 if((iposq<1).OR.(iposq>sum(nqred))) cycle
 opp1=vecq(iposq,2)
 opp2=vecq(iposq,3)
 opp3=vecq(iposq,4)

 decalageom=0
 if(ptsq(iposq-posq+2)<qsepred(2))then
  if(fenom.GE.2) decalageom=sum(nomred(1:fenom-1))
 elseif(ptsq(iposq-posq+2)<qsepred(3))then
  if(fenom==3) decalageom=nomred(1)
  if(fenom==4) decalageom=nomred(1)+nomred(3)
  if(fenom==2) err(2)=.TRUE.
  if(fenom==2.AND.(iposq==posq)) err(0)=.TRUE.
  if(fenom==2.AND.(iposq==posq+1)) err(1)=.TRUE.
  if(fenom==2) cycle
 else
  if(fenom==4) decalageom=nomred(3)
  if(fenom==1) err(2)=.TRUE.
  if(fenom==1.AND.(iposq==posq)) err(0)=.TRUE.
  if(fenom==1.AND.(iposq==posq+1)) err(1)=.TRUE.
  if(fenom==1) cycle
 endif

 sousgrille=donnees(1:4,decalageom+1:decalageom+nom,iposq)
 sousgrilley(2:4,:)=sousgrille(2:4,:)

 if(om<(opp(1)+opp(2))/2.0_qp)then
  y=om-opp(1)
  oppref(iposq-posq+2)=opp1
  sousgrilley(1,:)  =sousgrille(1,:)-opp1
 elseif(om<(opp(2)+opp(3))/2.0_qp)then
  y=om-opp(2)
  oppref(iposq-posq+2)=opp2
  sousgrilley(1,:)  =sousgrille(1,:)-opp2
 elseif(om<4*opp(3))then
  y=om-opp(3)
  oppref(iposq-posq+2)=opp3
  sousgrilley(1,:)  =sousgrille(1,:)-opp3
 else
  y=1.0_qp/om**(0.5_qp)
  sousgrilley(1,:)  =1.0_qp/sousgrille(1,:)**(0.5_qp)
 endif
 call locate(sousgrilley(1,1:nom),y,posy)

 if(blaM) write(6,*) "iposq,posom=",iposq,posom
 if((posy<1).OR.(posy.GE.nom))then
  if(blaerr) write(6,*)"*********************************************"
  if(blaerr) write(6,*)
  if(blaerr) write(6,*)"om hors grille à posq,iposq=",posq,iposq
  if(blaerr) write(6,*)"om,posom(nom),pts extrêmes=",om,posy,"(",nom,")",sousgrille(1,1),sousgrille(1,nom)
  if(blaerr) write(6,*)
  if(blaerr) write(6,*)"*********************************************"
  if(blaerr) blablaerr=.TRUE.
  if(iposq==posq)then
   err(0)=.TRUE.
  elseif(iposq==posq+1)then
   err(1)=.TRUE.
  else
   err(2)=.TRUE.
  endif
  cycle
 endif

 if((posy<2).OR.(posy.GE.nom-1))then
  err(2)=.TRUE.
  ptsM(1:3,2:3,iposq-posq+2)=sousgrilley(2:4,posy:posy+1)
  ptsy    (2:3,iposq-posq+2)=sousgrilley(1,  posy:posy+1)

 else
  ptsM(1:3,1:4,iposq-posq+2)=sousgrilley(2:4,posy-1:posy+2)
  ptsy    (1:4,iposq-posq+2)=sousgrilley(1,  posy-1:posy+2)
 endif
enddo

if(any(err).AND.(blaerr)) blablaerr=.TRUE.

if(err(0).AND.err(1))then
 if(blaerr) write(6,*)"Erreur de type 0"
 if(blaerr) write(6,*)"*********************************************"
elseif(err(2))then
 if(blaerr) write(6,*)"Erreur de type 2 (le carré extérieur de 12 est incomplet)"
 if(blaerr) write(6,*)"*********************************************"
endif


if(blaM.OR.blablaerr)then
 write(6,*)
 write(6,*)"------------- interpolM -------------------"
 write(6,*)
 write(6,*)"q,om=",q,om
 write(6,*)"opp=",opp
 write(6,*)
 write(6,*)"qsepred=",qsepred
 write(6,*)"omsepred=",omsepred
 write(6,*)
 write(6,*)"fenom=",fenom
 write(6,*)
 write(6,*)"posq,  vecq(posq)  =",posq,  vecq(posq,:)
 write(6,*)"posq+1,vecq(posq+1)=",posq+1,vecq(posq+1,:)
 write(6,*)
 write(6,*)"--------------- Points retenus ---------------"
 write(6,*)
 write(6,*)"ptsq=",ptsq
 write(6,*)
 if(fenom==4)then
  write(6,*)"ptsom(1,:)=",1.0_qp/ptsy(1,:)**(2.0_qp)
  write(6,*)"ptsom(2,:)=",1.0_qp/ptsy(2,:)**(2.0_qp)
  write(6,*)"ptsom(3,:)=",1.0_qp/ptsy(3,:)**(2.0_qp)
  write(6,*)"ptsom(4,:)=",1.0_qp/ptsy(4,:)**(2.0_qp)
 else
  write(6,*)"ptsom(1,:)=",oppref(:)+ptsy(1,:)
  write(6,*)"ptsom(2,:)=",oppref(:)+ptsy(2,:)
  write(6,*)"ptsom(3,:)=",oppref(:)+ptsy(3,:)
  write(6,*)"ptsom(4,:)=",oppref(:)+ptsy(4,:)
 endif

 write(6,*)
 write(6,*)"ptsM(1,1,:)=",ptsM(1,1,:)
 write(6,*)"ptsM(1,2,:)=",ptsM(1,2,:)
 write(6,*)"ptsM(1,3,:)=",ptsM(1,3,:)
 write(6,*)"ptsM(1,4,:)=",ptsM(1,4,:)
 write(6,*)
endif

if(err(0).OR.err(1)) return

do iMm=1,3
 if(err(2))then
  call polin2(ptsq(2:3),ptsy(2:3,2:3),ptsM(iMm,2:3,2:3),q,y,interpolM(iMm),dinterpolM(iMm))
 else
  call polin2(ptsq,ptsy,ptsM(iMm,:,:),q,y,interpolM(iMm),dinterpolM(iMm))
 endif
 if(abs(dinterpolM(iMm)/interpolM(iMm))>0.01_qp)then
   if(blaerr) blablaerr=.TRUE.
   if(blaerr) write(6,*)"grosse erreur d’interpolation"
 endif
 if(abs(dinterpolM(iMm)/interpolM(iMm))>0.2_qp)then
   err(0:1)=.TRUE.
   if(blaerr) write(6,*)"trop grosse erreur d’interpolation"
 endif
enddo

interpolM(4)=-PI*rhopp(1.5_qp,om,-1.0_qp)
interpolM(5)=-PI*rhopp(2.5_qp,om,-1.0_qp)
interpolM(6)=-PI*rhopp(3.5_qp,om,-1.0_qp)

dinterpolM(4:6)=EPSpp

deallocate(sousgrille)
deallocate(sousgrilley)

END FUNCTION interpolM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION interpolom2(q,om,dinterpolom2,err)
USE modsim
REAL(QP), INTENT(IN) :: om,q
LOGICAL,  INTENT(INOUT) :: err(0:2)
REAL(QP), DIMENSION(1:6), INTENT(OUT) :: dinterpolom2
REAL(QP), DIMENSION(1:6) :: interpolom2

INTEGER posq,posom
INTEGER iposq,iMm
REAL(QP) ptsq(1:4),ptsom(1:4,1:4),ptsM(1:3,1:4,1:4)
REAL(QP) y,opp1,opp2,opp3,oppref(1:4)
REAL(QP), DIMENSION(:,:), ALLOCATABLE :: sousgrille

err(:)   =.FALSE.
blablaerr=.FALSE.

ptsM (:,:,:)=0.0_qp
ptsom(:,:)  =0.0_qp
ptsq (:)    =0.0_qp

allocate(sousgrille (1:4,1:nomsup))

xq=q
call oangpp

call locate(vecq_sup(:,1),q,posq)

if(((posq)<1).OR.((posq).GE.nqsup))then
  err(0:1)=.TRUE.
  if(blaerr)then 
   write(6,*)"*********************************************"
   write(6,*)
   write(6,*)"Dans interpolom2: Erreur de type 0 venant de q"
   write(6,*)"q,posq=",q,posq
   write(6,*)
   write(6,*)"*********************************************"
  endif
  return
endif

if(((posq)<2).OR.((posq).GE.nqsup-1))then
  err(2)=.TRUE.
  if(blaerr) blablaerr=.TRUE.
  ptsq(2:3)=vecq_sup(posq:posq+1,1)
else
  ptsq(1:4)=vecq_sup(posq-1:posq+2,1)
endif



do iposq=posq-1,posq+2
 if((iposq<1).OR.(iposq>nqsup)) cycle
 opp1=vecq_sup(iposq,2)
 opp2=vecq_sup(iposq,3)
 opp3=vecq_sup(iposq,4)
 if(om<(opp(1)+opp(2))/2)then 
  oppref(iposq-posq+2)=opp1
  sousgrille=donnees_sup(1:4,0*nomsup+1:1*nomsup,iposq)

  if((opp(2)-opp(1))<0.01_qp)then
   y=(om-opp(1))/(opp(2)-opp(1))
   sousgrille(1,:)=(sousgrille(1,:)-opp1)/(opp2-opp1)
  else
   y=om-opp(1)
   sousgrille(1,:)=sousgrille(1,:)-opp1
  endif
 elseif(om<opp(2))then 
  sousgrille=donnees_sup(1:4,1*nomsup+1:2*nomsup,iposq)

  y=sqrt(opp(2)-om)
  sousgrille(1,:)  =sqrt(opp2-sousgrille(1,:))
  sousgrille(2:4,:)=sousgrille(2:4,:)
 elseif(om<opp(3))then
  oppref(iposq-posq+2)=opp2
  sousgrille=donnees_sup(1:4,2*nomsup+1:3*nomsup,iposq)

  y=om-opp(2)
  sousgrille(1,:)=sousgrille(1,:)-opp2
 else
  oppref(iposq-posq+2)=opp3
  sousgrille=donnees_sup(1:4,3*nomsup+1:4*nomsup,iposq)

  y=om-opp(3)
  sousgrille(1,:)=sousgrille(1,:)-opp3
 endif
 call locate(sousgrille(1,:),y,posom)

 if(blaM) write(6,*) "iposq,posom=",iposq,posom
 if((posom<1).OR.(posom.GE.nomsup))then
  if(blaerr) write(6,*)"*********************************************"
  if(blaerr) write(6,*)
  if(blaerr) write(6,*)"om hors grille à posq,iposq=",posq,iposq
  if(blaerr) write(6,*)"om,posom(nomsup),pts extrêmes=",om,posom,"(",nomsup,")",donnees_sup(1,1,iposq)&
  ,donnees_sup(1,2*nomsup,iposq)
  if(blaerr) write(6,*)
  if(blaerr) write(6,*)"*********************************************"
  if(blaerr) blablaerr=.TRUE.
  if(iposq==posq)then
   err(0)=.TRUE.
  elseif(iposq==posq+1)then
   err(1)=.TRUE.
  else
   err(2)=.TRUE.
  endif
  cycle
 endif

 if((posom<2).OR.(posom.GE.nomsup-1))then
  err(2)=.TRUE.
  ptsM(1:3,2:3,iposq-posq+2)=sousgrille(2:4,posom:posom+1)
  ptsom   (2:3,iposq-posq+2)=sousgrille(1,  posom:posom+1)

 else
  ptsM(1:3,1:4,iposq-posq+2)=sousgrille(2:4,posom-1:posom+2)
  ptsom   (1:4,iposq-posq+2)=sousgrille(1,  posom-1:posom+2)
 endif
enddo

if(any(err).AND.(blaerr)) blablaerr=.TRUE.

if(err(0).AND.err(1))then
 if(blaerr) write(6,*)"Erreur de type 0"
 if(blaerr) write(6,*)"*********************************************"
elseif(err(2))then
 if(blaerr) write(6,*)"Erreur de type 2 (le carré extérieur de 12 est incomplet)"
 if(blaerr) write(6,*)"*********************************************"
endif


if(blaM.OR.blablaerr)then
 write(6,*)
 write(6,*)"------------- interpolom2 -------------------"
 write(6,*)
 write(6,*)"q,om,y=",q,om,y
 write(6,*)"opp=",opp(1:3)
 write(6,*)
 write(6,*)"posq,  vecq_sup(posq)  =",posq,  vecq_sup(posq,1)
 write(6,*)"posq+1,vecq_sup(posq+1)=",posq+1,vecq_sup(posq+1,1)
 write(6,*)
 write(6,*)"--------------- Points retenus ---------------"
 write(6,*)
 write(6,*)"ptsq=",ptsq
 write(6,*)
 if((om>(opp(1)+opp(2))/2.0_qp).AND.(om<opp(2)))then
  write(6,*)"ptsom(1,:)=",opp2-ptsom(1,:)**2
  write(6,*)"ptsom(2,:)=",opp2-ptsom(2,:)**2
  write(6,*)"ptsom(3,:)=",opp2-ptsom(3,:)**2
  write(6,*)"ptsom(4,:)=",opp2-ptsom(4,:)**2
 else
  if((opp(2)-opp(1))<0.01_qp)then
   write(6,*)"ptsom(1,:)=",ptsom(1,:)*(vecq_sup(posq-1:posq+2,3)-vecq_sup(posq-1:posq+2,2))+vecq_sup(posq-1:posq+2,2)
   write(6,*)"ptsom(2,:)=",ptsom(2,:)*(vecq_sup(posq-1:posq+2,3)-vecq_sup(posq-1:posq+2,2))+vecq_sup(posq-1:posq+2,2)
   write(6,*)"ptsom(3,:)=",ptsom(3,:)*(vecq_sup(posq-1:posq+2,3)-vecq_sup(posq-1:posq+2,2))+vecq_sup(posq-1:posq+2,2)
   write(6,*)"ptsom(4,:)=",ptsom(4,:)*(vecq_sup(posq-1:posq+2,3)-vecq_sup(posq-1:posq+2,2))+vecq_sup(posq-1:posq+2,2)
  else
   write(6,*)"ptsom(1,:)=",ptsom(1,:)+oppref(:)
   write(6,*)"ptsom(2,:)=",ptsom(2,:)+oppref(:)
   write(6,*)"ptsom(3,:)=",ptsom(3,:)+oppref(:)
   write(6,*)"ptsom(4,:)=",ptsom(4,:)+oppref(:)
  endif
 endif

 write(6,*)
 write(6,*)"ptsM(1,1,:)=",ptsM(1,1,:)
 write(6,*)"ptsM(1,2,:)=",ptsM(1,2,:)
 write(6,*)"ptsM(1,3,:)=",ptsM(1,3,:)
 write(6,*)"ptsM(1,4,:)=",ptsM(1,4,:)
 write(6,*)
endif

if(err(0).OR.err(1)) return

do iMm=1,3
  if(err(2))then
   call polin2(ptsq(2:3),ptsom(2:3,2:3),ptsM(iMm,2:3,2:3),q,y,interpolom2(iMm),dinterpolom2(iMm))
  else
   call polin2(ptsq,ptsom,ptsM(iMm,:,:),q,y,interpolom2(iMm),dinterpolom2(iMm))
  endif
  if(abs(dinterpolom2(iMm)/interpolom2(iMm))>0.01_qp)then
    if(blaerr) blablaerr=.TRUE.
    if(blaerr) write(6,*)"grosse erreur d’interpolation iMm=",iMm
  endif
  if(abs(dinterpolom2(iMm)/interpolom2(iMm))>0.2_qp)then
    err(0:1)=.TRUE.
    if(blaerr) write(6,*)"trop grosse erreur d’interpolation iMm=",iMm
  endif
enddo

interpolom2(4)=-PI*rhopp(1.5_qp,om,-1.0_qp)
interpolom2(5)=-PI*rhopp(2.5_qp,om,-1.0_qp)
interpolom2(6)=-PI*rhopp(3.5_qp,om,-1.0_qp)

dinterpolom2(4:6)=EPSpp

deallocate(sousgrille)

END FUNCTION interpolom2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION interpolom3(q,om,dinterpolom3,err)
USE modsim
REAL(QP), INTENT(IN) :: om,q
LOGICAL,  INTENT(INOUT) :: err(0:2)
REAL(QP), DIMENSION(1:6), INTENT(OUT) :: dinterpolom3
REAL(QP), DIMENSION(1:6) :: interpolom3

INTEGER posq,posom
INTEGER iposq,iMm
REAL(QP) ptsq(1:4),ptsom(1:4,1:4),ptsM(1:3,1:4,1:4)
REAL(QP) y,opp1,opp2,opp3,oppref(1:4)
REAL(QP), DIMENSION(:,:), ALLOCATABLE :: sousgrille

err(:)   =.FALSE.
blablaerr=.FALSE.

ptsM (:,:,:)=0.0_qp
ptsom(:,:)  =0.0_qp
ptsq (:)    =0.0_qp

allocate(sousgrille (1:4,1:nomsup2))

xq=q
call oangpp
call locate(vecq_sup2(:,1),q,posq)

if(((posq)<1).OR.((posq).GE.nqsup2))then
  err(0:1)=.TRUE.
  if(blaerr)then 
   write(6,*)"*********************************************"
   write(6,*)
   write(6,*)"Dans interpolom3: Erreur de type 0 venant de q"
   write(6,*)"q,posq=",q,posq
   write(6,*)
   write(6,*)"*********************************************"
  endif
  return
endif

if(((posq)<2).OR.((posq).GE.nqsup2-1))then
  err(2)=.TRUE.
  if(blaerr) blablaerr=.TRUE.
  ptsq(2:3)=vecq_sup2(posq:posq+1,1)
else
  ptsq(1:4)=vecq_sup2(posq-1:posq+2,1)
endif



do iposq=posq-1,posq+2
 if((iposq<1).OR.(iposq>nqsup2)) cycle
 opp1=vecq_sup2(iposq,2)
 opp2=vecq_sup2(iposq,3)
 opp3=vecq_sup2(iposq,4)
 if(om<(opp(1)+opp(2))/2)then 
  oppref(iposq-posq+2)=opp1
  sousgrille=donnees_sup2(1:4,0*nomsup2+1:1*nomsup2,iposq)

  if((opp(2)-opp(1))<0.01_qp)then
   y=(om-opp(1))/(opp(2)-opp(1))
   sousgrille(1,:)=(sousgrille(1,:)-opp1)/(opp2-opp1)
  else
   y=om-opp(1)
   sousgrille(1,:)=sousgrille(1,:)-opp1
  endif
 elseif(om<opp(2))then 
  sousgrille=donnees_sup2(1:4,1*nomsup2+1:2*nomsup2,iposq)

  y=sqrt(opp(2)-om)
  sousgrille(1,:)  =sqrt(opp2-sousgrille(1,:))
  sousgrille(2:4,:)=sousgrille(2:4,:)
 elseif(om<(opp(3)+opp(2))/2)then
  oppref(iposq-posq+2)=opp2
  sousgrille=donnees_sup2(1:4,2*nomsup2+1:3*nomsup2,iposq)

  y=om-opp(2)
  sousgrille(1,:)=sousgrille(1,:)-opp2
 elseif(om<opp(3))then
  oppref(iposq-posq+2)=opp3
  sousgrille=donnees_sup2(1:4,3*nomsup2+1:4*nomsup2,iposq)

  y=opp(3)-om
  sousgrille(1,:)=opp3-sousgrille(1,:)
 else
  oppref(iposq-posq+2)=opp3
  sousgrille=donnees_sup2(1:4,4*nomsup2+1:5*nomsup2,iposq)

  y=om-opp(3)
  sousgrille(1,:)=sousgrille(1,:)-opp3
 endif
 call locate(sousgrille(1,:),y,posom)

 if(blaM.OR.blaerr) write(6,*)"iposq,posom=",iposq,posom

 if((posom<1).OR.(posom.GE.nomsup2))then
  if(blaerr) write(6,*)"*********************************************"
  if(blaerr) write(6,*)
  if(blaerr) write(6,*)"om hors grille à posq,iposq=",posq,iposq
  if(blaerr) write(6,*)"om,posom(nomsup2),pts extrêmes=",om,posom,"(",nomsup2,")",donnees_sup2(1,1,iposq)&
  ,donnees_sup2(1,2*nomsup2,iposq)
  if(blaerr) write(6,*)
  if(blaerr) write(6,*)"*********************************************"
  if(blaerr) blablaerr=.TRUE.
  if(iposq==posq)then
   err(0)=.TRUE.
  elseif(iposq==posq+1)then
   err(1)=.TRUE.
  else
   err(2)=.TRUE.
  endif
  cycle
 endif

 if((posom<2).OR.(posom.GE.nomsup2-1))then
  err(2)=.TRUE.
  ptsM(1:3,2:3,iposq-posq+2)=sousgrille(2:4,posom:posom+1)
  ptsom   (2:3,iposq-posq+2)=sousgrille(1,  posom:posom+1)

 else
  ptsM(1:3,1:4,iposq-posq+2)=sousgrille(2:4,posom-1:posom+2)
  ptsom   (1:4,iposq-posq+2)=sousgrille(1,  posom-1:posom+2)
 endif
enddo

if(any(err).AND.(blaerr)) blablaerr=.TRUE.

if(err(0).AND.err(1))then
 if(blaerr) write(6,*)"Erreur de type 0"
 if(blaerr) write(6,*)"*********************************************"
elseif(err(2))then
 if(blaerr) write(6,*)"Erreur de type 2 (le carré extérieur de 12 est incomplet)"
 if(blaerr) write(6,*)"*********************************************"
endif


if(blaM.OR.blablaerr)then
 write(6,*)
 write(6,*)"------------- interpolom3 -------------------"
 write(6,*)
 write(6,*)"q,om,y=",q,om,y
 write(6,*)"opp=",opp(1:3)
 write(6,*)
 write(6,*)"posq,  vecq_sup2(posq)  =",posq,  vecq_sup2(posq,1)
 write(6,*)"posq+1,vecq_sup2(posq+1)=",posq+1,vecq_sup2(posq+1,1)
 write(6,*)
 write(6,*)"-------------- Points retenus ---------------"
 write(6,*)
 write(6,*)"ptsq=",ptsq
 write(6,*)
 if((om>(opp(1)+opp(2))/2.0_qp).AND.(om<opp(2)))then
  write(6,*)"ptsom(1,:)=",opp2-ptsom(1,:)**2
  write(6,*)"ptsom(2,:)=",opp2-ptsom(2,:)**2
  write(6,*)"ptsom(3,:)=",opp2-ptsom(3,:)**2
  write(6,*)"ptsom(4,:)=",opp2-ptsom(4,:)**2
 else
  if((opp(2)-opp(1))<0.01_qp)then
   write(6,*)"ptsom(1,:)=",ptsom(1,:)*(vecq_sup2(posq-1:posq+2,3)-vecq_sup2(posq-1:posq+2,2))+vecq_sup2(posq-1:posq+2,2)
   write(6,*)"ptsom(2,:)=",ptsom(2,:)*(vecq_sup2(posq-1:posq+2,3)-vecq_sup2(posq-1:posq+2,2))+vecq_sup2(posq-1:posq+2,2)
   write(6,*)"ptsom(3,:)=",ptsom(3,:)*(vecq_sup2(posq-1:posq+2,3)-vecq_sup2(posq-1:posq+2,2))+vecq_sup2(posq-1:posq+2,2)
   write(6,*)"ptsom(4,:)=",ptsom(4,:)*(vecq_sup2(posq-1:posq+2,3)-vecq_sup2(posq-1:posq+2,2))+vecq_sup2(posq-1:posq+2,2)
  else
   write(6,*)"ptsom(1,:)=",ptsom(1,:)+oppref(:)
   write(6,*)"ptsom(2,:)=",ptsom(2,:)+oppref(:)
   write(6,*)"ptsom(3,:)=",ptsom(3,:)+oppref(:)
   write(6,*)"ptsom(4,:)=",ptsom(4,:)+oppref(:)
  endif
 endif

 write(6,*)
 write(6,*)"ptsM(1,1,:)=",ptsM(1,1,:)
 write(6,*)"ptsM(1,2,:)=",ptsM(1,2,:)
 write(6,*)"ptsM(1,3,:)=",ptsM(1,3,:)
 write(6,*)"ptsM(1,4,:)=",ptsM(1,4,:)
 write(6,*)
endif

if(err(0).OR.err(1)) return

do iMm=1,3
  if(err(2))then
   call polin2(ptsq(2:3),ptsom(2:3,2:3),ptsM(iMm,2:3,2:3),q,y,interpolom3(iMm),dinterpolom3(iMm))
  else
   call polin2(ptsq,ptsom,ptsM(iMm,:,:),q,y,interpolom3(iMm),dinterpolom3(iMm))
  endif
  if(abs(dinterpolom3(iMm)/interpolom3(iMm))>0.01_qp)then
    if(blaerr) blablaerr=.TRUE.
    if(blaerr) write(6,*)"grosse erreur d’interpolation iMm=",iMm
  endif
  if(abs(dinterpolom3(iMm)/interpolom3(iMm))>0.2_qp)then
    err(0:1)=.TRUE.
    if(blaerr) write(6,*)"trop grosse erreur d’interpolation iMm=",iMm
  endif
enddo

interpolom3(4)=-PI*rhopp(1.5_qp,om,-1.0_qp)
interpolom3(5)=-PI*rhopp(2.5_qp,om,-1.0_qp)
interpolom3(6)=-PI*rhopp(3.5_qp,om,-1.0_qp)

dinterpolom3(4:6)=EPSpp

deallocate(sousgrille)

END FUNCTION interpolom3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE load_data(fich)
INTEGER ixq,ixqbis,compteur,nn
CHARACTER(len=90), INTENT(IN) :: fich
REAL(QP) xqlec
open(212,file=trim(fich)//".info")
 read(212,*)
 read(212,*) x0,qsep,nqfen,nomfen,nn,bmax
 if(blaM)then
  write(6,*)"----------------------------------------"
  write(6,*)
  write(6,*)"Chargement du tableau de données"
  write(6,*)
  write(6,*)"x0,qsep=",x0,qsep
  write(6,*)"nq1a,nq1b,nq2,nq3=",nqfen
  write(6,*)"nom1a,nom1b,nom2a,nom2b,nom2c,nom3a,nom3b=",nomfen
  write(6,*)"nn,bmax=",nn,bmax
  write(6,*)
 endif
close(212)

allocate(donnees(1:7,1:sum(nomfen),1:sum(nqfen)))
allocate(vecq(1:sum(nqfen),1:4))

donnees(1,:,:)=1.0e50_qp
vecq(:,:)=1.0e50_qp

open(213,file=trim(fich)//"grilleq.dat")
read(213,*)
do ixq=1,sum(nqfen)
!do ixq=1,700
  read(213,*)ixqbis,xqlec,compteur
  xq=xqlec
  call oangpp
  vecq(ixq,:)=(/xqlec,opp(1:3)/)
  open(211,file=trim(fich)//".dat",ACTION="READ",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   read(211,REC=ixq)donnees(1:7,1:sum(nomfen),ixq)
   if(blaM)then 
    write(6,*)"ixq,compteur,donnees(1:3,370,ixq)=",ixq,compteur,donnees(1:3,370,ixq)
   endif
   donnees(:,compteur+1:sum(nomfen),ixq)=1.0e50_qp
  close(211)
enddo
close(213)

if(blaM)then
 write(6,*)"Chargement terminé"
 write(6,*)
 write(6,*)"----------------------------------------"
endif

END SUBROUTINE load_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE loadom2(fich)
INTEGER ixq,ixqbis,nn
CHARACTER(len=90), INTENT(IN) :: fich
REAL(QP) xqlec
open(212,file=trim(fich)//".info")
 read(212,*)
 read(212,*) x0,xq2,nqsup,nomsup,nn
 if(blaM)then
  write(6,*)"----------------------------------------"
  write(6,*)
  write(6,*)"Chargement du tableau de données supplémentaires"
  write(6,*)
  write(6,*)"x0,xq2=",x0,xq2
  write(6,*)"nqsup =",nqsup
  write(6,*)"nomsup=",nomsup
  write(6,*)"nn=",nn
  write(6,*)
 endif
close(212)

allocate(donnees_sup(1:7,1:4*nomsup,1:nqsup))
allocate(vecq_sup(1:nqsup,1:4))

donnees_sup(1,:,:)=1.0e50_qp
vecq_sup(:,:)=1.0e50_qp

open(213,file=trim(fich)//"grilleq.dat")
read(213,*)
do ixq=1,nqsup
!do ixq=1,7990
  read(213,*)ixqbis,xqlec
  xq=xqlec
  call oangpp
  vecq_sup(ixq,:)=(/xqlec,opp(1:3)/)
  open(211,file=trim(fich)//".dat",ACTION="READ",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   read(211,REC=ixq)donnees_sup(1:7,1:4*nomsup,ixq)
   if(blaM)then 
    write(6,*)"ixq,donnees_sup(1:3,75,ixq)=",ixq,xqlec,opp(2),donnees_sup(1:3,75,ixq)
   endif
  close(211)
enddo
close(213)

if(blaM)then
 write(6,*)"Chargement terminé"
 write(6,*)
 write(6,*)"----------------------------------------"
endif

END SUBROUTINE loadom2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE loadom3(fich)
INTEGER ixq,ixqbis,nn
CHARACTER(len=90), INTENT(IN) :: fich
REAL(QP) xqlec
open(212,file=trim(fich)//".info")
 read(212,*)
 read(212,*) x0,xq2,nqsup2,nomsup2,nn
 if(blaM)then
  write(6,*)"----------------------------------------"
  write(6,*)
  write(6,*)"Chargement du tableau de données supplémentaires"
  write(6,*)
  write(6,*)"x0,xq2=",x0,xq2
  write(6,*)"nqsup2 =",nqsup2
  write(6,*)"nomsup2=",nomsup2
  write(6,*)"nn=",nn
  write(6,*)
 endif
close(212)

allocate(donnees_sup2(1:7,1:5*nomsup2,1:nqsup2))
allocate(vecq_sup2(1:nqsup2,1:4))

donnees_sup2(1,:,:)=1.0e50_qp
vecq_sup2(:,:)=1.0e50_qp

open(213,file=trim(fich)//"grilleq.dat")
read(213,*)
do ixq=1,nqsup2
  read(213,*)ixqbis,xqlec
  xq=xqlec
  call oangpp
  vecq_sup2(ixq,:)=(/xqlec,opp(1:3)/)
  open(211,file=trim(fich)//".dat",ACTION="READ",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   read(211,REC=ixq)donnees_sup2(1:7,1:5*nomsup2,ixq)
   if(blaM)then 
    write(6,*)"ixq,donnees_sup2(1:3,75,ixq)=",ixq,xqlec,opp(2),donnees_sup2(1:3,75,ixq)
   endif
  close(211)
enddo
close(213)

if(blaM)then
 write(6,*)"Chargement terminé"
 write(6,*)
 write(6,*)"----------------------------------------"
endif

END SUBROUTINE loadom3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE unload_data
deallocate(donnees)
deallocate(donnees_sup)
deallocate(donnees_sup2)
deallocate(vecq)
deallocate(vecq_sup)
deallocate(vecq_sup2)
END SUBROUTINE unload_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE combineom2(fich,fich2)
INTEGER ixq,ixqbis,compteur,nn
CHARACTER(len=90), INTENT(IN) :: fich,fich2
REAL(QP) xqlec
REAL(QP), DIMENSION(:,:), ALLOCATABLE :: donnees_temp
INTEGER nmin1,nmin2,nmax1,nmax2

open(216,file=trim(fich)//".info")
 read(216,*)
 read(216,*) x0,xq2,nqsup,nomsup,nn
close(216)

allocate(donnees_temp(1:7,1:4*nomsup))
open(213,file=trim(fich) //"grilleq.dat",ACTION="READ")
 read(213,*)
 read(213,*)nmin1,xqlec
 write(6,*)nmin1
 do 
  read(213, *, End = 1 ) nmax1,xqlec
!  write(6,*)nmax1
 enddo
 write(6,*)nmax1
1 Continue
close(213)
open(214,file=trim(fich2)//"grilleq.dat",ACTION="READ")
 read(214,*)
 read(214,*)nmin2,xqlec
 write(6,*)nmin2,xqlec
 do 
  read(214, *, End = 2 ) nmax2,xqlec
!  write(6,*)nmax2
 enddo
 write(6,*)nmax2
2 Continue
close(214)
write(6,*)nmin1,nmin2,nmax1,nmax2
do ixq=nmin1,nmax1-1

  open(211,file=trim(fich)//".dat",ACTION="READ",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   read(211,REC=ixq)donnees_temp(1:7,1:4*nomsup)
   write(6,*)"ixq,donnees_temp(1,1050)=",ixq,donnees_temp(1,1050)
  close(211)

  open(215,file=trim(fich)//"_comb.dat",ACTION="WRITE",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   write(215,REC=ixq)donnees_temp(1:7,1:4*nomsup)
  close(215)

enddo
write(6,*)
do ixq=nmin2,nmax2

  open(212,file=trim(fich2)//".dat",ACTION="READ",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   read(212,REC=ixq)donnees_temp(1:7,1:4*nomsup)
   write(6,*)"ixq,donnees_temp(1,1050)=",ixq,donnees_temp(1,1050)
  close(212)

  open(215,file=trim(fich)//"_comb.dat",ACTION="WRITE",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   write(215,REC=ixq)donnees_temp(1:7,1:4*nomsup)
  close(215)

enddo
END SUBROUTINE combineom2
END MODULE bestM
