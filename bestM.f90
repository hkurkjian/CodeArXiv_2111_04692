MODULE bestM
USE nrtype
USE vars
USE dspec
USE modpol
IMPLICIT NONE
REAL(QP), DIMENSION(:,:,:), ALLOCATABLE :: donnees,donnees_sup
REAL(QP), DIMENSION(:,:),   ALLOCATABLE :: vecq,vecq_sup
REAL(QP) qsep(1:5),bmax,qpetit,xq2
INTEGER nqfen(1:4),nomfen(1:8),nbidon,nqsup,nomsup
LOGICAL blaM,blaerr,blablaerr
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE load_data(fich)
INTEGER ixq,ixqbis,compteur,nn
CHARACTER(len=90), INTENT(IN) :: fich
REAL(QP) xqlec
open(12,file=trim(fich)//".info")
 read(12,*)
 read(12,*) x0,qsep,nqfen,nomfen,nn,bmax
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
close(12)

allocate(donnees(1:7,1:sum(nomfen),1:sum(nqfen)))
allocate(vecq(1:sum(nqfen),1:4))

donnees(1,:,:)=1.0e50_qp
vecq(:,:)=1.0e50_qp

open(13,file=trim(fich)//"grilleq.dat")
read(13,*)
!do ixq=1,sum(nqfen)
do ixq=1,700
  read(13,*)ixqbis,xqlec,compteur
  xq=xqlec
  call oangpp
  vecq(ixq,:)=(/xqlec,opp(1:3)/)
  open(11,file=trim(fich)//".dat",ACTION="READ",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   read(11,REC=ixq)donnees(1:7,1:sum(nomfen),ixq)
   if(blaM)then 
    write(6,*)"ixq,compteur,donnees(1:3,370,ixq)=",ixq,compteur,donnees(1:3,370,ixq)
   endif
   donnees(:,compteur+1:sum(nomfen),ixq)=0.0_qp
  close(11)
enddo
close(13)

if(blaM)then
 write(6,*)"Chargement terminé"
 write(6,*)
 write(6,*)"----------------------------------------"
endif

call system("rm pts_horsgrille.dat")
END SUBROUTINE load_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE loadom2(fich)
INTEGER ixq,ixqbis,compteur,nn
CHARACTER(len=90), INTENT(IN) :: fich
REAL(QP) xqlec
open(12,file=trim(fich)//".info")
 read(12,*)
 read(12,*) x0,xq2,nqsup,nomsup,nn
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
close(12)

allocate(donnees_sup(1:7,1:2*nomsup,1:nqsup))
allocate(vecq_sup(1:nqsup,1:4))

donnees_sup(1,:,:)=1.0e50_qp
vecq_sup(:,:)=1.0e50_qp

open(13,file=trim(fich)//"grilleq.dat")
read(13,*)
!do ixq=1,nqsup
do ixq=1,7000
  read(13,*)ixqbis,xqlec
  xq=xqlec
  call oangpp
  vecq_sup(ixq,:)=(/xqlec,opp(1:3)/)
  open(11,file=trim(fich)//".dat",ACTION="READ",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   read(11,REC=ixq)donnees_sup(1:7,1:2*nomsup,ixq)
   if(blaM)then 
    write(6,*)"ixq,donnees_sup(1:3,75,ixq)=",ixq,xqlec,opp(2),donnees_sup(1:3,75,ixq)
   endif
  close(11)
enddo
close(13)

if(blaM)then
 write(6,*)"Chargement terminé"
 write(6,*)
 write(6,*)"----------------------------------------"
endif

END SUBROUTINE loadom2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE unload_data
deallocate(donnees)
deallocate(vecq)
deallocate(donnees_sup)
deallocate(vecq_sup)
END SUBROUTINE unload_data
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

LOGICAL :: err(0:2)
REAL(QP),  DIMENSION(1:6) :: dinterpolM
INTEGER iMm
COMPLEX(QPC) Gg(1:2,1:2),Mm(1:2,1:2),det

interpolM_recerr=interpolM(q,om,dinterpolM,err)

if(err(0).OR.err(1))then
 if(blaerr) write(6,*)"Récupération d’erreur de type 0"
 if((om>(opp(1)+opp(2))/2).AND.(om<(opp(2)+opp(3))/2))then
!  if(blaerr) write(6,*)"Tentative de récupération par interpolom2"
  interpolM_recerr=interpolom2(q,om,dinterpolM,err)
 endif
 if(err(0).OR.err(1))then
!  write(6,*)"Erreur sèche,q,om=",q,om
  open(15,file="pts_horsgrille.dat",POSITION="APPEND")
   write(15,*)q,opp(2),om
  close(15)
  if(q<qpetit)then
   call mat_pairfield_pttq(om,0.0_qp,det,Mm,Gg)
  else
   call mat_pairfield(om,0.0_qp,det,Mm,Gg)
  endif
  interpolM_recerr(1:6)=(/real(Mm(1,1)),real(Mm(2,2)),real(Mm(1,2)),imag(Mm(1,1)),imag(Mm(2,2)),imag(Mm(1,2))/)
  if(blaerr)then
  do iMm=1,3
   write(6,*)"om,q,interpolM(1),dinterpolM(1)=",om,xq,interpolM_recerr(iMm),dinterpolM(iMm)
  enddo
  endif
 else
!  write(6,*)"Erreur récupérée par intom2"
 endif
endif

if(blaM.OR.blablaerr)then
 do iMm=4,6
  write(6,*)"om,q,im(interpol)                  =",om,q,interpolM_recerr(iMm)
 enddo
 write(6,*)
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
INTEGER nq,decalageq,fenq,nom,decalageom,fenom
INTEGER nqred(1:3),nomred(1:4)
INTEGER iposq,iMm
REAL(QP) y,qsepred(1:4),omsepred(1:5)
REAL(QP) ptsq(1:4),ptsom(1:4,1:4),ptsy(1:4,1:4),ptsM(1:3,1:4,1:4)
REAL(QP), DIMENSION(:,:),   ALLOCATABLE :: sousgrille,sousgrilley
REAL(QP), DIMENSION(:),     ALLOCATABLE :: sousgrilleq

err(:)   =.FALSE.
blablaerr=.FALSE.

ptsM (:,:,:)=0.0_qp
ptsom(:,:)  =0.0_qp
ptsy (:,:)  =0.0_qp
ptsq (:)    =0.0_qp

nqred =(/nqfen(1)+nqfen(2),nqfen(3),nqfen(4)/)
nomred=(/sum(nomfen(1:2)),sum(nomfen(3:5)),sum(nomfen(6:7)),nomfen(8)/)

xq=q
call oangpp

qsepred =(/qsep(1),qsep(3),qsep(4),qsep(5)/)
omsepred=(/opp(1),opp(2),opp(3),4*opp(3),bmax/)

call locate(qsepred ,q ,fenq)
call locate(omsepred,om,fenom)

if(fenq==0)then
 if(blaerr) write(6,*) "q inférieur à qmin"
 if(blaerr) blablaerr=.TRUE.
 return
elseif(fenq==4)then
 if(blaerr) write(6,*) "q dépasse qmax"
 if(blaerr) blablaerr=.TRUE.
 return
endif

if(fenom==0)then
 if(blaerr) write(6,*)"om n'est pas dans le continuum de brisure de paires"
 if(blaerr) blablaerr=.TRUE.
 return
elseif(fenom==5)then
 if(blaerr) write(6,*)"om dépasse bmax"
 if(blaerr) blablaerr=.TRUE.
 return
endif

nq =nqred(fenq)
nom=nomred(fenom)

if(fenom==4) y=1.0_qp/om**(0.5_qp)

decalageom=0
if(fenq==1)then
 decalageq=0
 if(fenom.GE.2) decalageom=sum(nomred(1:fenom-1))
elseif(fenq==2)then
 decalageq=nqfen(1)
 if(fenom==3) decalageom=nomred(1)
 if(fenom==4) decalageom=nomred(1)+nomred(3)
elseif(fenq==3)then
 decalageq=nqfen(1)+nqfen(2)
 if(fenom==4) decalageom=nomred(3)
endif

allocate(sousgrilleq(1:nq))
allocate(sousgrille (1:4,1:nom))
allocate(sousgrilley(1:4,1:nom))

sousgrilleq=vecq(decalageq+1:decalageq+nq,1)
call locate(sousgrilleq,q,posq)

if(((posq)<1).OR.((posq).GE.nq))then
  err(0:1)=.TRUE.
  if(blaerr)then 
   write(6,*)"*********************************************"
   write(6,*)
   write(6,*)"Erreur de type 0 venant de q"
   write(6,*)"q,posq,fenetre=",q,posq,qsepred(fenq),qsepred(fenq+1)
   write(6,*)
   write(6,*)"*********************************************"
  endif
  return
endif

if(((posq)<2).OR.((posq).GE.nq-1))then
  err(2)=.TRUE.
  if(blaerr) blablaerr=.TRUE.
  ptsq(2:3)=sousgrilleq(posq:posq+1)
else
  ptsq(1:4)=sousgrilleq(posq-1:posq+2)
endif



do iposq=posq-1,posq+2
 if((iposq<1).OR.(iposq>nq)) cycle
 sousgrille=donnees(1:4,decalageom+1:decalageom+nom,iposq)
 if(fenom==4)then
  sousgrilley(1,:)  =1.0_qp/sousgrille(1,:)**(0.5_qp)
  sousgrilley(2:4,:)=sousgrille(2:4,:)
  call locate(sousgrilley(1,1:nom),y,posy)
  posom=posy
 else
  call locate(sousgrille(1,1:nom),om,posom)
 endif

 if((posom<1).OR.(posom.GE.nom))then
  if(blaerr) write(6,*)"*********************************************"
  if(blaerr) write(6,*)
  if(blaerr) write(6,*)"om hors grille à posq,iposq=",posq,iposq
  if(blaerr) write(6,*)"om,posom(nom)=",om,posom,"(",nom,")"
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

 if(fenom==4)then
  if((posy<2).OR.(posy.GE.nom-1))then
   err(2)=.TRUE.
   ptsM(1:3,2:3,iposq-posq+2)=sousgrilley(2:4,posy:posy+1)
   ptsy    (2:3,iposq-posq+2)=sousgrilley(1,  posy:posy+1)
 
  else
   ptsM(1:3,1:4,iposq-posq+2)=sousgrilley(2:4,posy-1:posy+2)
   ptsy    (1:4,iposq-posq+2)=sousgrilley(1,  posy-1:posy+2)
  endif
 else
  if((posom<2).OR.(posom.GE.nom-1))then
   err(2)=.TRUE.
   ptsM(1:3,2:3,iposq-posq+2)=sousgrille(2:4,posom:posom+1)
   ptsom   (2:3,iposq-posq+2)=sousgrille(1,  posom:posom+1)
 
  else
   ptsM(1:3,1:4,iposq-posq+2)=sousgrille(2:4,posom-1:posom+2)
   ptsom   (1:4,iposq-posq+2)=sousgrille(1,  posom-1:posom+2)
  endif
 endif
enddo

if(any(err).AND.(blaerr)) blablaerr=.TRUE.

if(err(0).AND.err(1))then
 if(blaerr) write(6,*)"Erreur de type 0"
 if(blaerr) write(6,*)"*********************************************"
 return
elseif(err(0))then
 if(blaerr) write(6,*)"Erreur de type 1 (pts(2:3) manquant à posq)"
 if(blaerr) write(6,*)"*********************************************"
 return
elseif(err(1))then
 if(blaerr) write(6,*)"Erreur de type 1 (pts(2:3) manquant à posq+1)"
 if(blaerr) write(6,*)"*********************************************"
 return
elseif(err(2))then
 if(blaerr) write(6,*)"Erreur de type 2 (le carré extérieur de 12 est incomplet)"
 if(blaerr) write(6,*)"*********************************************"
endif



do iMm=1,3
  if(fenom==4)then
!     if(err(0))then
!      call polint(ptsy(2:3,3),ptsM(iMm,2:3,3),y,interpolM(iMm),dinterpolM(iMm))
!     elseif(err(1))then
!      call polint(ptsy(2:3,2),ptsM(iMm,2:3,2),y,interpolM(iMm),dinterpolM(iMm))
     if(err(2))then
      call polin2(ptsq(2:3),ptsy(2:3,2:3),ptsM(iMm,2:3,2:3),q,y,interpolM(iMm),dinterpolM(iMm))
     else
      call polin2(ptsq,ptsy,ptsM(iMm,:,:),q,y,interpolM(iMm),dinterpolM(iMm))
     endif
  else
!     if(err(0))then
!      call polint(ptsom(2:3,3),ptsM(iMm,2:3,3),om,interpolM(iMm),dinterpolM(iMm))
!     elseif(err(1))then
!      call polint(ptsom(2:3,2),ptsM(iMm,2:3,2),om,interpolM(iMm),dinterpolM(iMm))
     if(err(2))then
      call polin2(ptsq(2:3),ptsom(2:3,2:3),ptsM(iMm,2:3,2:3),q,om,interpolM(iMm),dinterpolM(iMm))
     else
      call polin2(ptsq,ptsom,ptsM(iMm,:,:),q,om,interpolM(iMm),dinterpolM(iMm))
     endif
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

if(blaM.OR.blablaerr)then
 write(6,*)
 write(6,*)"------------- interpolM -------------------"
 write(6,*)
 write(6,*)"q,om=",q,om
 write(6,*)"opp=",xq,opp
 write(6,*)
 write(6,*)"qsepred=",qsepred
 write(6,*)"omsepred=",omsepred
 write(6,*)
 write(6,*)"fenq =",fenq
 write(6,*)"fenom=",fenom
 write(6,*)
 write(6,*)"posq,  sousgrilleq(posq)  =",posq,  sousgrilleq(posq)
 write(6,*)"posq+1,sousgrilleq(posq+1)=",posq+1,sousgrilleq(posq+1)
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
  write(6,*)"ptsom(1,:)=",ptsom(1,:)
  write(6,*)"ptsom(2,:)=",ptsom(2,:)
  write(6,*)"ptsom(3,:)=",ptsom(3,:)
  write(6,*)"ptsom(4,:)=",ptsom(4,:)
 endif

 write(6,*)
 write(6,*)"ptsM(1,1,:)=",ptsM(1,1,:)
 write(6,*)"ptsM(1,2,:)=",ptsM(1,2,:)
 write(6,*)"ptsM(1,3,:)=",ptsM(1,3,:)
 write(6,*)"ptsM(1,4,:)=",ptsM(1,4,:)
 write(6,*)
 do iMm=1,3
  write(6,*)"om,q,interpolM(1),dinterpolM(1)=",om,q,interpolM(iMm),dinterpolM(iMm)
 enddo
endif

interpolM(4)=-PI*rhopp(1.5_qp,om,-1.0_qp)
interpolM(5)=-PI*rhopp(2.5_qp,om,-1.0_qp)
interpolM(6)=-PI*rhopp(3.5_qp,om,-1.0_qp)

dinterpolM(4:6)=EPSpp

deallocate(sousgrille)
deallocate(sousgrilleq)
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
REAL(QP) y,opp2
REAL(QP), DIMENSION(:,:), ALLOCATABLE :: sousgrille
REAL(QP), DIMENSION(:,:),   ALLOCATABLE :: sousgrilleq

err(:)   =.FALSE.
blablaerr=.FALSE.

ptsM (:,:,:)=0.0_qp
ptsom(:,:)  =0.0_qp
ptsq (:)    =0.0_qp

allocate(sousgrilleq(1:nqsup,1:3))
allocate(sousgrille (1:4,1:nomsup))

sousgrilleq=vecq_sup(1:nqsup,1:3)
call locate(sousgrilleq(:,1),q,posq)

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
  ptsq(2:3)=sousgrilleq(posq:posq+1,1)
else
  ptsq(1:4)=sousgrilleq(posq-1:posq+2,1)
endif



do iposq=posq-1,posq+2
 opp2=sousgrilleq(iposq,3)
! write(6,*)"m=",om
! write(6,*)"opp2=",opp2
 if((iposq<1).OR.(iposq>nqsup)) cycle
 if(om<opp(2))then 
  sousgrille=donnees_sup(1:4,1:nomsup      ,iposq)
  y=sqrt(opp(2)-om)
!  write(6,*)"xq,y=",xq,y

  sousgrille(1,:)  =sqrt(opp2-sousgrille(1,:))
  sousgrille(2:4,:)=sousgrille(2:4,:)

  call locate(sousgrille(1,:),y,posom)
 else
  sousgrille=donnees_sup(1:4,nomsup+1:2*nomsup,iposq)
!  do iMm=1,nomsup
!   write(6,*)sousgrille(1,iMm)
!  enddo
  call locate(sousgrille(1,:),om,posom)
 endif

 if((posom<1).OR.(posom.GE.nomsup))then
  if(blaerr) write(6,*)"*********************************************"
  if(blaerr) write(6,*)
  if(blaerr) write(6,*)"om hors grille à posq,iposq=",posq,iposq
  if(blaerr) write(6,*)"om,posom(nomsup)=",om,posom,"(",nomsup,")"
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
 return
elseif(err(0))then
 if(blaerr) write(6,*)"Erreur de type 1 (pts(2:3) manquant à posq)"
 if(blaerr) write(6,*)"*********************************************"
 return
elseif(err(1))then
 if(blaerr) write(6,*)"Erreur de type 1 (pts(2:3) manquant à posq+1)"
 if(blaerr) write(6,*)"*********************************************"
 return
elseif(err(2))then
 if(blaerr) write(6,*)"Erreur de type 2 (le carré extérieur de 12 est incomplet)"
 if(blaerr) write(6,*)"*********************************************"
endif



do iMm=1,3
 if(om<opp(2))then 
  if(err(2))then
   call polin2(ptsq(2:3),ptsom(2:3,2:3),ptsM(iMm,2:3,2:3),q,y,interpolom2(iMm),dinterpolom2(iMm))
  else
   call polin2(ptsq,ptsom,ptsM(iMm,:,:),q,y,interpolom2(iMm),dinterpolom2(iMm))
  endif
 else
  if(err(2))then
   call polin2(ptsq(2:3),ptsom(2:3,2:3),ptsM(iMm,2:3,2:3),q,om,interpolom2(iMm),dinterpolom2(iMm))
  else
   call polin2(ptsq,ptsom,ptsM(iMm,:,:),q,om,interpolom2(iMm),dinterpolom2(iMm))
  endif
 endif
  if(abs(dinterpolom2(iMm)/interpolom2(iMm))>0.01_qp)then
    if(blaerr) blablaerr=.TRUE.
    if(blaerr) write(6,*)"grosse erreur d’interpolation"
  endif
  if(abs(dinterpolom2(iMm)/interpolom2(iMm))>0.2_qp)then
    err(0:1)=.TRUE.
    if(blaerr) write(6,*)"trop grosse erreur d’interpolation"
  endif
enddo

if(blaM.OR.blablaerr)then
 write(6,*)
 write(6,*)"------------- interpolom2 -------------------"
 write(6,*)
 write(6,*)"q,om=",q,om
 write(6,*)"opp=",xq,opp
 write(6,*)
 write(6,*)"posq,  sousgrilleq(posq)  =",posq,  sousgrilleq(posq,1)
 write(6,*)"posq+1,sousgrilleq(posq+1)=",posq+1,sousgrilleq(posq+1,1)
 write(6,*)
 write(6,*)"--------------- Points retenus ---------------"
 write(6,*)
 write(6,*)"ptsq=",ptsq
 write(6,*)
 if(om<opp(2))then
  write(6,*)"ptsom(1,:)=",opp(2)-ptsom(1,:)**2
  write(6,*)"ptsom(2,:)=",opp(2)-ptsom(2,:)**2
  write(6,*)"ptsom(3,:)=",opp(2)-ptsom(3,:)**2
  write(6,*)"ptsom(4,:)=",opp(2)-ptsom(4,:)**2
 else
  write(6,*)"ptsom(1,:)=",ptsom(1,:)
  write(6,*)"ptsom(2,:)=",ptsom(2,:)
  write(6,*)"ptsom(3,:)=",ptsom(3,:)
  write(6,*)"ptsom(4,:)=",ptsom(4,:)
 endif

 write(6,*)
 write(6,*)"ptsM(1,1,:)=",ptsM(1,1,:)
 write(6,*)"ptsM(1,2,:)=",ptsM(1,2,:)
 write(6,*)"ptsM(1,3,:)=",ptsM(1,3,:)
 write(6,*)"ptsM(1,4,:)=",ptsM(1,4,:)
 write(6,*)
 do iMm=1,3
  write(6,*)"om,q,interpolom2(1),dinterpolom2(1)=",om,q,interpolom2(iMm),dinterpolom2(iMm)
 enddo
endif

interpolom2(4)=-PI*rhopp(1.5_qp,om,-1.0_qp)
interpolom2(5)=-PI*rhopp(2.5_qp,om,-1.0_qp)
interpolom2(6)=-PI*rhopp(3.5_qp,om,-1.0_qp)

dinterpolom2(4:6)=EPSpp

deallocate(sousgrille)
deallocate(sousgrilleq)

END FUNCTION interpolom2
END MODULE bestM
