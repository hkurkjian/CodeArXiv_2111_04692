MODULE estM
USE nrtype
USE vars
USE dspec
USE modpol
IMPLICIT NONE
REAL(QP), DIMENSION(:,:,:), ALLOCATABLE :: donnees
REAL(QP), DIMENSION(:),     ALLOCATABLE :: vecq
REAL(QP) xqmin,xq1,xq2,xqmax,bmax
INTEGER nq1,nq2,nq3,nom1,nom2,nom3,nominf
LOGICAL blaM
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE load_data(fich)
INTEGER ixq,ixqbis,compteur,nn
CHARACTER(len=90), INTENT(IN) :: fich
open(12,file=trim(fich)//".info")
 read(12,*)
 read(12,*) x0,xqmin,xq1,xq2,xqmax,nq1,nq2,nq3,nom1,nom2,nom3,nominf,nn,bmax
 if(blaM)then
  write(6,*)
  write(6,*)"Chargement du tableau de données"
  write(6,*)
  write(6,*)"x0,xqmin,xq1,xq2,xqmax=",x0,xqmin,xq1,xq2,xqmax
  write(6,*)"nq1,nq2,nq3,nom1,nom2,nom3,nominf=",nq1,nq2,nq3,nom1,nom2,nom3,nominf
  write(6,*)"nn,bmax=",nn,bmax
  write(6,*)
 endif
close(12)
allocate(donnees(1:7,1:nom1+nom2+nom3+nominf,1:nq1+nq2+nq3))
allocate(vecq(1:nq1+nq2+nq3))
donnees(1,:,:)=1.0e50_qp
open(13,file=trim(fich)//"grilleq.dat")
read(13,*)
do ixq=1,nq1+nq2+nq3
  read(13,*)ixqbis,xq,compteur
  vecq(ixq)=xq
!  write(6,*)"ixq,xq,compteur=",ixq,xq,compteur
  open(11,file=trim(fich)//".dat",ACTION="READ",ACCESS="DIRECT",FORM="UNFORMATTED",RECL=nn)
   read(11,REC=ixq)donnees(1:7,1:nom1+nom2+nom3+nominf,ixq)
   donnees(:,compteur+1:nom1+nom2+nom3+nominf,ixq)=0.0_qp
   if(blaM)then 
    write(6,*)"ixq,donnees(1:3,15,ixq)=",ixq,donnees(1:3,15,ixq) 
   endif
  close(11)
enddo
close(13)
if(blaM)then
 write(6,*)"Chargement terminé"
 write(6,*)
endif
!write(6,*) donnees(1:2,15,15)
!write(6,*) vecq
END SUBROUTINE load_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE unload_data
deallocate(donnees)
deallocate(vecq)
END SUBROUTINE unload_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE estmat_pairfield(om0,e,det,M,fr)
COMPLEX(QPC), DIMENSION(1:2,1:2), INTENT(OUT) :: M,fr
COMPLEX(QPC), INTENT(OUT) :: det
REAL(QP), INTENT(IN)  :: om0,e

LOGICAL :: errtype1,errtype2
REAL(QP),  DIMENSION(1:6) :: Mmv
REAL(QP) :: q

q=xq
Mmv=interpolM_recerr(q,om0)
M(1,1)=cmplx(Mmv(1),Mmv(4),kind=qpc)
M(2,2)=cmplx(Mmv(2),Mmv(5),kind=qpc)
M(1,2)=cmplx(Mmv(3),Mmv(5),kind=qpc)
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

LOGICAL :: errtype1,errtype2
REAL(QP),  DIMENSION(1:6) :: dinterpolM

COMPLEX(QPC) Gg(1:2,1:2),Mm(1:2,1:2),det

interpolM_recerr=interpolM(q,om,dinterpolM,errtype1,errtype2)

if(errtype1)then
 xq=q
 if(blaM)then
  write(6,*)"Récupération d’erreur de type 1"
 endif
 call mat_pairfield(om,0.0_qp,det,Mm,Gg)
 interpolM_recerr(1:6)=(/real(Mm(1,1)),real(Mm(2,2)),real(Mm(1,2)),imag(Mm(1,1)),imag(Mm(2,2)),imag(Mm(1,2))/)
 dinterpolM(1:6)=EPSpp
endif

if(blaM)then
 write(6,*)"interpolM_recerr(1),dinterpolM(1)=",interpolM_recerr(1),dinterpolM(1)
endif
END FUNCTION interpolM_recerr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION interpolM(q,om,dinterpolM,errtype1,errtype2)
USE modsim
REAL(QP), INTENT(IN) :: om,q
LOGICAL, INTENT(OUT) :: errtype1,errtype2
REAL(QP),  DIMENSION(1:6), INTENT(OUT) :: dinterpolM
REAL(QP), DIMENSION(1:6) :: interpolM

REAL(QP), DIMENSION(1:6,1:4,1:4) :: points
INTEGER posq,posom,posy,ixqbis,ixq,nq,decalageq,nom,decalageom,fen,iposq,npts,iMm
REAL(QP) ommin,ommax,dom,ymin,ymax,dy,y,dxq,xqdep,xqfin
REAL(QP) ptsq(1:4),ptsom(1:4,1:4),ptsy(1:4,1:4),ptsM(1:6,1:4,1:4)
REAL(QP), DIMENSION(:),     ALLOCATABLE :: grilleom
REAL(QP), DIMENSION(:,:),   ALLOCATABLE :: sousgrille,sousgrilley

errtype1=.FALSE.
errtype2=.FALSE.

xq=q
call oangpp
if(om<opp(1))then
 stop "interpolM: om n'est pas dans le continuum de brisure de paires"
elseif(om<opp(2))then
 fen=1
elseif(om<opp(3))then
 fen=2
elseif(om<4*opp(3))then
 fen=3
else
 fen=4
 y=1.0_qp/om**(1.5_qp)
endif

if(q<xq1)then
 xqdep=xqmin
 xqfin=xq1
 nq=nq1
 decalageq=0
 if(fen==1)then
  decalageom=0
  nom=nom1
 elseif(fen==2)then
  decalageom=nom1
  nom=nom2
 elseif(fen==3)then
  decalageom=nom1+nom2
  nom=nom3
 else
  decalageom=nom1+nom2+nom3
  nom=nominf
 endif
elseif(q<xq2)then
 xqdep=xq1
 xqfin=xq2
 nq=nq2
 decalageq=nq1
 if(fen==1)then
  decalageom=0
  nom=nom1
 elseif(fen==3)then
  decalageom=nom1
  nom=nom3
 else
  decalageom=nom1+nom3
  nom=nominf
 endif
else
 xqdep=xq2
 xqfin=4*xq2
 nq=nq3
 decalageq=nq1+nq2
 if(fen==3)then
  decalageom=0
  nom=nom3
 else
  decalageom=nom3
  nom=nominf
 endif
endif


dxq=(xqfin-xqdep)/nq
posq=floor((q-xqdep+0.5_qp*dxq)/dxq)+decalageq

if(((posq-decalageq)<1).OR.((posq-decalageq).GE.nq))then
  errtype1=.TRUE.
  errtype2=.TRUE.
  if(blaM) write(6,*)"Erreur de type 1"
  return
endif

if(((posq-decalageq)<2).OR.((posq-decalageq).GE.nq-1))then
  errtype2=.TRUE.
  ptsq(2:3)=vecq(posq:posq+1)
else
  ptsq(1:4)=vecq(posq-1:posq+2)
endif


if((posq>0).AND.(posq.LE.(nq1+nq2+nq3)).AND.blaM)then
 write(6,*)
 write(6,*)"------------- interpolM -------------------"
 write(6,*)
 write(6,*)"q,om=",q,om
 write(6,*)
 write(6,*)"posq,vecq(posq)=",posq,vecq(posq)
endif

allocate(sousgrille (1:7,1:nom))
allocate(sousgrilley(1:7,1:nom))

do iposq=posq-1,posq+2
 if(((iposq-decalageq)<1).OR.((iposq-decalageq)>nq)) cycle
 sousgrille=donnees(1:7,decalageom+1:decalageom+nom,iposq)
 if(fen==4)then
  sousgrilley(1,:)  =1.0_qp/sousgrille(1,:)**(1.5_qp)
  sousgrilley(2:7,:)=sousgrille(2:7,:)
  dy=sousgrilley(1,2)-sousgrilley(1,1)
  ymin=sousgrilley(1,1)-0.5_qp*dy
  posy=(y-ymin+0.5_qp*dy)/dy
  posom=posy
 else
  dom=sousgrille(1,2)-sousgrille(1,1)
  ommin=sousgrille(1,1)-0.5_qp*dom
  posom=(om-ommin+0.5_qp*dom)/dom
 endif

 if(((posom<1).OR.(posom.GE.nom)).AND.((iposq==posq).OR.(iposq==posq+1)))then
  errtype1=.TRUE.
  errtype2=.TRUE.
  if(blaM)then
   write(6,*)"Erreur de type 1"
   write(6,*)"posq,iposq=",posq,iposq
  endif
  return
 endif

 if(fen==4)then
  if((posy<2).OR.(posy.GE.nom-1))then
   errtype2=.TRUE.
   ptsM(1:6,2:3,iposq-posq+2)=sousgrilley(2:7,posy:posy+1)
   ptsy   (2:3,iposq-posq+2) =sousgrilley(1,  posy:posy+1)
 
  else
   ptsM(1:6,1:4,iposq-posq+2)=sousgrilley(2:7,posy-1:posy+2)
   ptsy   (1:4,iposq-posq+2) =sousgrilley(1,  posy-1:posy+2)
  endif
 else
  if((posom<2).OR.(posom.GE.nom-1))then
   errtype2=.TRUE.
   ptsM(1:6,2:3,iposq-posq+2)=sousgrille(2:7,posom:posom+1)
   ptsom   (2:3,iposq-posq+2)=sousgrille(1,  posom:posom+1)
 
  else
   ptsM(1:6,1:4,iposq-posq+2)=sousgrille(2:7,posom-1:posom+2)
   ptsom   (1:4,iposq-posq+2)=sousgrille(1,  posom-1:posom+2)
  endif
 endif
enddo


if(blaM)then
 write(6,*)"ptsq=",ptsq
 write(6,*)
 if(fen==4)then
  write(6,*)"ptsom(1,:)=",1.0_qp/ptsy(1,:)**(2.0_qp/3.0_qp)
  write(6,*)"ptsom(2,:)=",1.0_qp/ptsy(2,:)**(2.0_qp/3.0_qp)
  write(6,*)"ptsom(3,:)=",1.0_qp/ptsy(3,:)**(2.0_qp/3.0_qp)
  write(6,*)"ptsom(4,:)=",1.0_qp/ptsy(4,:)**(2.0_qp/3.0_qp)
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
endif

do iMm=1,6
  if(fen==4)then
     if(errtype2)then
       call polin2(ptsq(2:3),ptsy(2:3,2:3),ptsM(iMm,2:3,2:3),q,y,interpolM(iMm),dinterpolM(iMm))
     else
      call polin2(ptsq,ptsy,ptsM(iMm,:,:),q,y,interpolM(iMm),dinterpolM(iMm))
     endif
  else
     if(errtype2)then
      call polin2(ptsq(2:3),ptsom(2:3,2:3),ptsM(iMm,2:3,2:3),q,om,interpolM(iMm),dinterpolM(iMm))
     else
      call polin2(ptsq,ptsom,ptsM(iMm,:,:),q,om,interpolM(iMm),dinterpolM(iMm))
     endif
  endif
enddo
if(blaM)then
 do iMm=1,6
  write(6,*)"interpolM(1),dinterpolM(1)=",interpolM(iMm),dinterpolM(iMm)
 enddo
 write(6,*)
endif

END FUNCTION interpolM
END MODULE estM
