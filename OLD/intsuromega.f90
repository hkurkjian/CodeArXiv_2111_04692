!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE modpro
 USE nrtype
 USE nrutil, ONLY : ecrit
 REAL(QP) xq,x0,om0,e,r0(1:3),o(1:4)
 REAL(QP), PARAMETER :: bidon=-1.0e300_qp
 LOGICAL bla1,bla2,axereel
 COMPLEX(QPC) z,zbascule,xibasc(1:2),phibasc(1:2),r3basc(1:2),al,omil 
 INTEGER ptbranchmt
 INTERFACE ff
   MODULE PROCEDURE ffr, ffc
 END INTERFACE ff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION  intsuromega()
USE nrtype
USE recettes
USE modsim

IMPLICIT NONE

REAL(QP) bmax,num,s,om,dom
REAL(QP) s1,s2,s3,sc,vsec(1:3),singu,nonsingu
PARAMETER (s1=1.5_qp,s2=2.5_qp,s3=3.5_qp,singu=1.0_qp,nonsingu=-1.0_qp)
REAL(QP)     Ia(1:3),Ib(1:3),Ic(1:3),Id(1:3),Ie(1:3),Ig(1:3),Iinf(1:3)
COMPLEX(QPC) intsuromega(1:3),Ica(1:3),Icb(1:3),Icc(1:3),Icd(1:3),Ice(1:3),Icg(1:3),Icinf(1:3),Icorr(1:3),a
REAL(QP) r(1:3)
REAL(QP) EPS

INTEGER isec

EPS=1.0e-9_qp !Précision des intégrales

om0=real(z)
e=imag(z)

!write(6,*)'bla1=',bla1

call oang
call ecrit(bla1,'xq,x0,om0,e=',(/xq,x0,om0,e/))
call ecrit(bla1,'o=',o(1:3))
r0=dspec(om0)
call ecrit(bla1,'r0=',r0(1:3))
call ecrit(bla1,'')

bmax=9500000
Ia(:)=0.0_qp
Ib(:)=0.0_qp
Ic(:)=0.0_qp
Id(:)=0.0_qp
Ie(:)=0.0_qp
Ig(:)=0.0_qp
Iinf(:)=0.0_qp

Ica(:)  =cmplx(0.0_qp,0.0_qp,QPC)
Icb(:)  =cmplx(0.0_qp,0.0_qp,QPC)
Icc(:)  =cmplx(0.0_qp,0.0_qp,QPC)
Icd(:)  =cmplx(0.0_qp,0.0_qp,QPC)
Ice(:)  =cmplx(0.0_qp,0.0_qp,QPC)
Icg(:)  =cmplx(0.0_qp,0.0_qp,QPC)
Icinf(:)=cmplx(0.0_qp,0.0_qp,QPC)

if(axereel)then
 Iinf(1)= PI*(x0+xq**2.0_qp/4.0_qp)*sqrt(2.0_qp/bmax)
 Iinf(2)= Iinf(1)
 Iinf(3)=-PI*om0*sqrt(2.0_qp/bmax)*(1.0_qp+2.0_qp*(x0-xq**2.0_qp/4.0_qp)/(3.0_qp*bmax))
 call ecrit(bla1,'Iinf=',Iinf)
 
 vsec=(/s1,s2,s3/)
 do isec=1,3
  sc=vsec(isec)
  if(om0<o(1))then
   Ia(isec)=qromo(inter,0.0_qp         ,o(1)           ,(/sc,nonsingu/),midpntq,EPS)
   call ecrit(bla1,'Ia=',Ia(isec))
   if(ptbranchmt>1)then
    Ib(isec)=qromo(inter,o(1)           ,o(2)           ,(/sc,nonsingu/),midpntq,EPS)
    call ecrit(bla1,'Ib=',Ib(isec))
   endif
   if(ptbranchmt>2)then
    Ic(isec)=qromo(inter,o(2)           ,o(3)           ,(/sc,nonsingu/),midsqlq,EPS)
    call ecrit(bla1,'Ic=',Ic(isec))
   endif
   Id(isec)=qromo(inter,o(3)           ,bmax           ,(/sc,nonsingu/),racinfq,EPS)
   call ecrit(bla1,'Id=',Id(isec))
  elseif(om0<(o(1)+o(2))/2.0_qp)then
   Ia(isec)=qromo(inter,0.0_qp         ,o(1)           ,(/sc,nonsingu/),midpntq,EPS)
   call ecrit(bla1,'Ia=',Ia(isec))
   Ib(isec)=qromo(inter,o(1)           ,om0            ,(/sc,singu   /),midpntq,EPS)
   call ecrit(bla1,'Ib=',Ib(isec))
   Ic(isec)=qromo(inter,om0            ,2.0_qp*om0-o(1),(/sc,singu   /),midpntq,EPS)
   call ecrit(bla1,'Ic=',Ic(isec))
   Id(isec)=qromo(inter,2.0_qp*om0-o(1),o(2)           ,(/sc,nonsingu/),midpntq,EPS)
   call ecrit(bla1,'Id=',Id(isec))
   if(ptbranchmt>2)then
    Ie(isec)=qromo(inter,o(2)           ,o(3)           ,(/sc,nonsingu/),midsqlq,EPS)
    call ecrit(bla1,'Ie=',Ie(isec))
   endif
   Ig(isec)=qromo(inter,o(3)           ,bmax           ,(/sc,nonsingu/),racinfq,EPS)
   call ecrit(bla1,'Ig=',Ig(isec))
   call ecrit(bla1,'')
  elseif(om0<o(2))then
   Ia(isec)=qromo(inter,0.0_qp         ,o(1)           ,(/sc,nonsingu/),midpntq,EPS)
   call ecrit(bla1,'Ia=',Ia(isec))
   Ib(isec)=qromo(inter,o(1)           ,2.0_qp*om0-o(2),(/sc,nonsingu/),midpntq,EPS)
   call ecrit(bla1,'Ib=',Ib(isec))
   Ic(isec)=qromo(inter,2.0_qp*om0-o(2),om0            ,(/sc,singu   /),midpntq,EPS)
   call ecrit(bla1,'Ic=',Ic(isec))
   Id(isec)=qromo(inter,om0            ,o(2)           ,(/sc,singu   /),midpntq,EPS)
   call ecrit(bla1,'Id=',Id(isec))
   if(ptbranchmt>2)then
    Ie(isec)=qromo(inter,o(2)           ,o(3)           ,(/sc,nonsingu/),midsqlq,EPS)
    call ecrit(bla1,'Ie=',Ie(isec))
   endif
   Ig(isec)=qromo(inter,o(3)           ,bmax           ,(/sc,nonsingu/),racinfq,EPS)
   call ecrit(bla1,'Ig=',Ig(isec))
   call ecrit(bla1,'')
  elseif(om0<(o(2)+o(3))/2.0_qp)then
   Ia(isec)=qromo(inter,0.0_qp         ,o(1)           ,(/sc,nonsingu/),midpntq,EPS)
   call ecrit(bla1,'Ia=',Ia(isec))
   Ib(isec)=qromo(inter,o(1)           ,o(2)           ,(/sc,nonsingu/),midpntq,EPS)
   call ecrit(bla1,'Ib=',Ib(isec))
   Ic(isec)=qromo(inter,o(2)           ,om0            ,(/sc,singu   /),midpntq,EPS)
   call ecrit(bla1,'Ic=',Ic(isec))
   Id(isec)=qromo(inter,om0            ,2.0_qp*om0-o(2),(/sc,singu   /),midpntq,EPS)
   call ecrit(bla1,'Id=',Id(isec))
   Ie(isec)=qromo(inter,2.0_qp*om0-o(2),o(3)           ,(/sc,nonsingu/),midsqlq,EPS)
   call ecrit(bla1,'Ie=',Ie(isec))
   Ig(isec)=qromo(inter,o(3)           ,bmax           ,(/sc,nonsingu/),racinfq,EPS)
   call ecrit(bla1,'Ig=',Ig(isec))
   call ecrit(bla1,'')
  elseif(om0<o(3))then
   Ia(isec)=qromo(inter,0.0_qp         ,o(1)           ,(/sc,nonsingu/),midpntq,EPS)
   call ecrit(bla1,'Ia=',Ia(isec))
   Ib(isec)=qromo(inter,o(1)           ,o(2)           ,(/sc,nonsingu/),midpntq,EPS)
   call ecrit(bla1,'Ib=',Ib(isec))
   Ic(isec)=qromo(inter,o(2)           ,2.0_qp*om0-o(3),(/sc,nonsingu/),midpntq,EPS)
   call ecrit(bla1,'Ic=',Ic(isec))
   Id(isec)=qromo(inter,2.0_qp*om0-o(3),om0            ,(/sc,singu   /),midpntq,EPS)
   call ecrit(bla1,'Id=',Id(isec))
   Ie(isec)=qromo(inter,om0            ,o(3)           ,(/sc,singu   /),midpntq,EPS)
   call ecrit(bla1,'Ie=',Ie(isec))
   Ig(isec)=qromo(inter,o(3)           ,bmax           ,(/sc,nonsingu/),racinfq,EPS)
   call ecrit(bla1,'Ig=',Ig(isec))
   call ecrit(bla1,'')
  else
   Ia(isec)=qromo(inter,0.0_qp         ,o(1)           ,(/sc,nonsingu/),midpntq,EPS)
   call ecrit(bla1,'Ia=',Ia(isec))
   if(ptbranchmt>1)then
    Ib(isec)=qromo(inter,o(1)           ,o(2)           ,(/sc,nonsingu/),midpntq,EPS)
    call ecrit(bla1,'Ib=',Ib(isec))
   endif
   if(ptbranchmt>2)then
    Ic(isec)=qromo(inter,o(2)           ,o(3)           ,(/sc,nonsingu/),midpntq,EPS)
    call ecrit(bla1,'Ic=',Ic(isec))
   endif
   Id(isec)=qromo(inter,o(3)           ,om0            ,(/sc,singu   /),midsqlq,EPS)
   call ecrit(bla1,'Id=',Id(isec))
   Ie(isec)=qromo(inter,om0            ,2.0_qp*om0-o(3),(/sc,singu   /),midpntq,EPS)
   call ecrit(bla1,'Ie=',Ie(isec))
   Ig(isec)=qromo(inter,2.0_qp*om0-o(3),bmax          ,(/sc,nonsingu/),racinfq,EPS)
   call ecrit(bla1,'Ig=',Ig(isec))
   call ecrit(bla1,'')
  endif
 enddo
 intsuromega=cmplx(Ia+Ib+Ic+Id+Ie+Ig+Iinf,0.0_qp,QPC)
 call ecrit(bla1,'Ia+Ib+Ic+Id+Ie+Ig+Iinf=',Ia+Ib+Ic+Id+Ie+Ig+Iinf)
! stop

else
 Icinf(1)= PI*(x0+xq**2.0_qp/4.0_qp)*sqrt(2.0_qp/bmax)
 call ecrit(bla1,'Icinf(1)=',Icinf(1))
 Icinf(2)= Icinf(1)
 call ecrit(bla1,'Icinf(2)=',Icinf(2))
 Icinf(3)=-PI*z*sqrt(2.0_qp/bmax)*(1.0_qp+2.0_qp*(x0-xq**2.0_qp/4.0_qp)/(3.0_qp*bmax))
 call ecrit(bla1,'Icinf(3)=',Icinf(3))
 
 vsec=(/s1,s2,s3/)
 do isec=1,3
  sc=vsec(isec)
  if(om0<o(1))then
   a=0.0_qp
   Ica(isec)=qromo(intec,0.0_qp         ,o(1)           ,(/sc,nonsingu/),midpntcq,EPS)
   call ecrit(bla1,'Ica=',Ica(isec))
   if(ptbranchmt>1)then
    Icb(isec)=qromo(intec,o(1)           ,o(2)           ,(/sc,nonsingu/),midpntcq,EPS)
    call ecrit(bla1,'Icb=',Icb(isec))
   endif
   if(ptbranchmt>2)then
    Icc(isec)=qromo(intec,o(2)           ,o(3)           ,(/sc,nonsingu/),midsqlcq,EPS)
    call ecrit(bla1,'Icc=',Icc(isec))
   endif
   Icd(isec)=qromo(intec,o(3)           ,bmax           ,(/sc,nonsingu/),racinfcq,EPS)
   call ecrit(bla1,'Icd=',Icd(isec))
  elseif(om0<(o(1)+o(2))/2.0_qp)then
   a=om0-o(1)
   Ica(isec)=qromo(intec,0.0_qp         ,o(1)           ,(/sc,nonsingu/),midpntcq,EPS)
   call ecrit(bla1,'Ica=',Ica(isec))
   Icb(isec)=qromo(intec,o(1)           ,om0            ,(/sc,singu   /),midpntcq,EPS)
   call ecrit(bla1,'Icb=',Icb(isec))
   Icc(isec)=qromo(intec,om0            ,2.0_qp*om0-o(1),(/sc,singu   /),midpntcq,EPS)
   call ecrit(bla1,'Icc=',Icc(isec))
   Icd(isec)=qromo(intec,2.0_qp*om0-o(1),o(2)           ,(/sc,nonsingu/),midpntcq,EPS)
   call ecrit(bla1,'Icd=',Icd(isec))
   if(ptbranchmt>2)then
    Ice(isec)=qromo(intec,o(2)           ,o(3)           ,(/sc,nonsingu/),midsqlcq,EPS)
    call ecrit(bla1,'Ice=',Ice(isec))
   endif
   Icg(isec)=qromo(intec,o(3)           ,bmax           ,(/sc,nonsingu/),racinfcq,EPS)
   call ecrit(bla1,'Icg=',Icg(isec))
   call ecrit(bla1,'')
  elseif(om0<o(2))then
   a=o(2)-om0
   Ica(isec)=qromo(intec,0.0_qp         ,o(1)           ,(/sc,nonsingu/),midpntcq,EPS)
   call ecrit(bla1,'Ica=',Ica(isec))
   Icb(isec)=qromo(intec,o(1)           ,2.0_qp*om0-o(2),(/sc,nonsingu/),midpntcq,EPS)
   call ecrit(bla1,'Icb=',Icb(isec))
   Icc(isec)=qromo(intec,2.0_qp*om0-o(2),om0            ,(/sc,singu   /),midpntcq,EPS)
   call ecrit(bla1,'Icc=',Icc(isec))
   Icd(isec)=qromo(intec,om0            ,o(2)           ,(/sc,singu   /),midpntcq,EPS)
   call ecrit(bla1,'Icd=',Icd(isec))
   if(ptbranchmt>2)then
    Ice(isec)=qromo(intec,o(2)           ,o(3)           ,(/sc,nonsingu/),midsqlcq,EPS)
    call ecrit(bla1,'Ice=',Ice(isec))
   endif
   Icg(isec)=qromo(intec,o(3)           ,bmax           ,(/sc,nonsingu/),racinfcq,EPS)
   call ecrit(bla1,'Icg=',Icg(isec))
   call ecrit(bla1,'')
  elseif(om0<(o(2)+o(3))/2.0_qp)then
   a=om0-o(2)
   Ica(isec)=qromo(intec,0.0_qp         ,o(1)           ,(/sc,nonsingu/),midpntcq,EPS)
   call ecrit(bla1,'Ica=',Ica(isec))
   Icb(isec)=qromo(intec,o(1)           ,o(2)           ,(/sc,nonsingu/),midpntcq,EPS)
   call ecrit(bla1,'Icb=',Icb(isec))
   Icc(isec)=qromo(intec,o(2)           ,om0            ,(/sc,singu   /),midpntcq,EPS)
   call ecrit(bla1,'Icc=',Icc(isec))
   Icd(isec)=qromo(intec,om0            ,2.0_qp*om0-o(2),(/sc,singu   /),midpntcq,EPS)
   call ecrit(bla1,'Icd=',Icd(isec))
   Ice(isec)=qromo(intec,2.0_qp*om0-o(2),o(3)           ,(/sc,nonsingu/),midsqlcq,EPS)
   call ecrit(bla1,'Ice=',Ice(isec))
   Icg(isec)=qromo(intec,o(3)           ,bmax           ,(/sc,nonsingu/),racinfcq,EPS)
   call ecrit(bla1,'Icg=',Icg(isec))
   call ecrit(bla1,'')
  elseif(om0<o(3))then
   a=o(3)-om0
   Ica(isec)=qromo(intec,0.0_qp         ,o(1)           ,(/sc,nonsingu/),midpntcq,EPS)
   call ecrit(bla1,'Ica=',Ica(isec))
   Icb(isec)=qromo(intec,o(1)           ,o(2)           ,(/sc,nonsingu/),midpntcq,EPS)
   call ecrit(bla1,'Icb=',Icb(isec))
   Icc(isec)=qromo(intec,o(2)           ,2.0_qp*om0-o(3),(/sc,nonsingu/),midpntcq,EPS)
   call ecrit(bla1,'Icc=',Icc(isec))
   Icd(isec)=qromo(intec,2.0_qp*om0-o(3),om0            ,(/sc,singu   /),midpntcq,EPS)
   call ecrit(bla1,'Icd=',Icd(isec))
   Ice(isec)=qromo(intec,om0            ,o(3)           ,(/sc,singu   /),midpntcq,EPS)
   call ecrit(bla1,'Ice=',Ice(isec))
   Icg(isec)=qromo(intec,o(3)           ,bmax           ,(/sc,nonsingu/),racinfcq,EPS)
   call ecrit(bla1,'Icg=',Icg(isec))
   call ecrit(bla1,'')
  else
   a=om0-o(3)
   Ica(isec)=qromo(intec,0.0_qp         ,o(1)           ,(/sc,nonsingu/),midpntcq,EPS)
   call ecrit(bla1,'Ica=',Ica(isec))
   if(ptbranchmt>1)then
    Icb(isec)=qromo(intec,o(1)           ,o(2)           ,(/sc,nonsingu/),midpntcq,EPS)
    call ecrit(bla1,'Icb=',Icb(isec))
   endif
   if(ptbranchmt>2)then
    Icc(isec)=qromo(intec,o(2)           ,o(3)           ,(/sc,nonsingu/),midpntcq,EPS)
    call ecrit(bla1,'Icc=',Icc(isec))
   endif
   Icd(isec)=qromo(intec,o(3)           ,om0            ,(/sc,singu   /),midsqlcq,EPS)
   call ecrit(bla1,'Icd=',Icd(isec))
   Ice(isec)=qromo(intec,om0            ,2.0_qp*om0-o(3),(/sc,singu   /),midpntcq,EPS)
   call ecrit(bla1,'Ice=',Ice(isec))
   Icg(isec)=qromo(intec,2.0_qp*om0-o(3),bmax          ,(/sc,nonsingu/),racinfcq,EPS)
   call ecrit(bla1,'Icg=',Icg(isec))
   call ecrit(bla1,'')
  endif
  Icorr(isec)=r0(isec)*log((iiq*e+a)/(iiq*e-a))
!  Icorr=cmplx(0.0_qp,0.0_qp,QPC)
 enddo
 intsuromega=Ica+Icb+Icc+Icd+Ice+Icg+Icinf+Icorr
! call ecrit(bla1,'real(Ica+Icb+Icc+Icd+Ice+Icg+Icinf)=',real(Ica+Icb+Icc+Icd+Ice+Icg+Icinf+Icorr))
! call ecrit(bla1,'imag(Ica+Icb+Icc+Icd+Ice+Icg+Icinf)=',imag(Ica+Icb+Icc+Icd+Ice+Icg+Icinf+Icorr))
endif
END FUNCTION intsuromega
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION inter(om,arg)
USE nrtype
IMPLICIT NONE

REAL(QP), INTENT(IN), DIMENSION(:)  ::  om,arg
REAL(QP), DIMENSION(size(om))       ::  inter
REAL(QP) r(1:3),sec,s,ome
INTEGER is

sec=arg(1)
s=arg(2)

do is=1,size(om)
 ome=om(is)
 r =dspec(ome)
 if(s>0.0_qp)then
  if(sec<2.0_qp)then
   inter(is)=(r(1)-r0(1))/(om0-ome)-r(1)/(om0+ome)+PI*sqrt(ome)/sqrt(8.0_qp*(1.0_qp+(ome/2.0_qp-x0)**2.0_qp))
  elseif(sec<3.0_qp)then
   inter(is)=(r(2)-r0(2))/(om0-ome)-r(2)/(om0+ome)+PI*sqrt(ome)/sqrt(8.0_qp*(1.0_qp+(ome/2.0_qp-x0)**2.0_qp))
  else
   inter(is)=(r(3)-r0(3))/(om0-ome)+r(3)/(om0+ome)
  endif
 else
  if(sec<2.0_qp)then
   inter(is)=r(1)*(1.0_qp/(om0-ome)-1.0_qp/(om0+ome))+PI*sqrt(ome)/sqrt(8.0_qp*(1.0_qp+(ome/2.0_qp-x0)**2.0_qp))
  elseif(sec<3.0_qp)then
   inter(is)=r(2)*(1.0_qp/(om0-ome)-1.0_qp/(om0+ome))+PI*sqrt(ome)/sqrt(8.0_qp*(1.0_qp+(ome/2.0_qp-x0)**2.0_qp))
  else
   inter(is)=r(3)/(om0-ome)+r(3)/(om0+ome)
  endif
 endif
 if(isnan(inter(is)))STOP 'isnan(inter(is))'
enddo
!write(6,*)'om,inter=',om,inter
END FUNCTION inter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION intec(om,arg)
USE nrtype
IMPLICIT NONE

REAL(QP),     INTENT(IN), DIMENSION(:)  ::  om,arg
COMPLEX(QPC), DIMENSION(size(om))       ::  intec
REAL(QP) r(1:3),sec,s,ome
INTEGER is

sec=arg(1)
s=arg(2)

do is=1,size(om)
 ome=om(is)
 r =dspec(ome)
 if(s>0.0_qp)then
  if(sec<2.0_qp)then
   intec(is)=(r(1)-r0(1))/(z-ome)-r(1)/(z+ome)+PI*sqrt(ome)/sqrt(8.0_qp*(1.0_qp+(ome/2.0_qp-x0)**2.0_qp))
  elseif(sec<3.0_qp)then
   intec(is)=(r(2)-r0(2))/(z-ome)-r(2)/(z+ome)+PI*sqrt(ome)/sqrt(8.0_qp*(1.0_qp+(ome/2.0_qp-x0)**2.0_qp))
  else
   intec(is)=(r(3)-r0(3))/(z-ome)+r(3)/(z+ome)
  endif
 else
  if(sec<2.0_qp)then
   intec(is)=r(1)*(1.0_qp/(z-ome)-1.0_qp/(z+ome))+PI*sqrt(ome)/sqrt(8.0_qp*(1.0_qp+(ome/2.0_qp-x0)**2.0_qp))
!   write(6,*)'om,intec(is)',om(is),abs(intec(is))
  elseif(sec<3.0_qp)then
   intec(is)=r(2)*(1.0_qp/(z-ome)-1.0_qp/(z+ome))+PI*sqrt(ome)/sqrt(8.0_qp*(1.0_qp+(ome/2.0_qp-x0)**2.0_qp))
  else
   intec(is)=r(3)/(z-ome)+r(3)/(z+ome)
  endif
 endif
 if(isnan(real(intec(is))).OR.isnan(imag(intec(is))))STOP 'isnan(intec(is))'
enddo
END FUNCTION intec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION dspec(om)
USE nrtype
USE recettes
IMPLICIT NONE
REAL(QP), INTENT(IN)  :: om
REAL(QP) dspec(1:3)

REAL(QP) r1,r2,r3,rr3,xi(1:2),x
REAL(QP) sinht,cosht,tanht

call oang
cosht=om/2.0_qp
sinht=sqrt(cosht**2.0_qp-1.0_qp)
tanht=sinht/cosht

r1=0.0_qp
r2=0.0_qp
r3=0.0_qp

if(om<o(1))then
  call ecrit(bla2,'Secteur 0')
  dspec(:)=0.0_qp
elseif(om<o(2))then
  call ecrit(bla2,'Secteur 1')
  r1=2.0_qp*ellf(PI/2.0_qp,tanht)/cosht/4.0_qp
  r2=2.0_qp*elle(PI/2.0_qp,tanht)*cosht/2.0_qp
  dspec(1)=PI*r2/xq
  dspec(2)=PI*(r2-2.0_qp*r1)/xq
  dspec(3)=0.0_qp
elseif(om<o(3))then
  call ecrit(bla2,'Secteur 2')
  call bornesbrutal(om,2,xi)
 ! Partie décroissante
  x=-xi(1)/sinht
  rr3=-sinht*sqrt((1.0_qp-x**2.0_qp)/(1.0_qp-x**2.0_qp*tanht**2.0_qp))/2.0_qp
  r3=rr3
  r1=(ellf(PI/2.0_qp,tanht)-ellf(asin(x),tanht))/cosht/4.0_qp
  r2=(elle(PI/2.0_qp,tanht)-elle(asin(x),tanht))*cosht/2.0_qp-x*rr3*tanht
 ! Partie croissante
  x=xi(2)/sinht
  rr3=sinht*sqrt((1.0_qp-x**2.0_qp)/(1.0_qp-x**2.0_qp*tanht**2.0_qp))/2.0_qp
  r3=r3+rr3
  r1=r1+(ellf(PI/2.0_qp,tanht)-ellf(asin(x),tanht))/cosht/4.0_qp
  r2=r2+(elle(PI/2.0_qp,tanht)-elle(asin(x),tanht))*cosht/2.0_qp+x*rr3*tanht
  dspec(1)=PI*r2/xq
  dspec(2)=PI*(r2-2.0_qp*r1)/xq
  dspec(3)=PI*r3/xq
else
  call ecrit(bla2,'Secteur 3')
  call ecrit(bla2,'om')
  call bornesbrutal(om,3,xi)
  x=xi(1)/sinht
  if(isnan(asin(x)))then
   write(6,*)'isnan(asin(x))'
   write(6,*)'xi(1)=',xi(1)
   write(6,*)'sinht=',sinht
   write(6,*)'om=',om
   r1=0.0_qp
   r3=xq/2.0_qp*sqrt(om/2.0_qp)+xq*(x0-xq**2.0_qp/4.0_qp)/sqrt(8.0_qp*om)
   r2=r3*tanht*x
   write(6,*)'r3=',r3
  else
   r3=sinht*sqrt((1.0_qp-x**2.0_qp)/(1.0_qp-x**2.0_qp*tanht**2.0_qp))/2.0_qp
   r1=(ellf(PI/2.0_qp,tanht)-ellf(asin(x),tanht))/cosht/4.0_qp
   r2=(elle(PI/2.0_qp,tanht)-elle(asin(x),tanht))*cosht/2.0_qp+x*r3*tanht
  endif
  if(isnan(r3))then
   write(6,*)'x=',x
   write(6,*)'(1.0_qp-x**2.0_qp*tanht**2.0_qp)=',(1.0_qp-x**2.0_qp*tanht**2.0_qp)
   write(6,*)'(1.0_qp-x**2.0_qp)=',(1.0_qp-x**2.0_qp)
   STOP 'isnan(r3))'
 endif
 dspec(1)=PI*r2/xq
 dspec(2)=PI*(r2-2.0_qp*r1)/xq
 dspec(3)=PI*r3/xq
endif
END FUNCTION dspec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE bornesbrutal(om,sect,xi)
!Trouve les bornes d'intégration sur k ap integrale
!angulaire dans les secteurs 2 et 3
USE nrtype
USE recettes
IMPLICIT NONE
REAL(QP), INTENT(IN)  :: om
INTEGER, INTENT(IN) :: sect
REAL(QP), INTENT(OUT) :: xi(1:2)
REAL(QP) rac

xi=(/-1.e100_qp,-1.e100_qp/)

if(sect==2)then
 rac=rtsafeq(emax,(/om/),1.e-20_qp,o(4),1.e-18_qp*om)
 xi(1)=rac**2.0_qp-x0+xq**2.0_qp/4.0_qp
 rac=rtsafeq(emax,(/om/),o(4),1.e9_qp,1.e-18_qp*om)
 xi(2)=rac**2.0_qp-x0+xq**2.0_qp/4.0_qp

elseif(sect==3)then
 rac=rtsafe(emax,(/om/),1.e-20_qp,1.e13_qp,1.e-18_qp*om)
 xi(1)=rac**2.0_qp-x0+xq**2.0_qp/4.0_qp
 xi(2)=-1.e200_qp

else
 STOP 'Mauvais secteur'
endif
END SUBROUTINE bornesbrutal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION dspecpro(sect)
USE nrtype
USE recettes, ONLY : cacae,cacaf
IMPLICIT NONE

INTEGER, INTENT(IN)     :: sect
COMPLEX(QPC) dspecpro(1:3)

COMPLEX(QPC) r1(1:2),r2(1:2),r3(1:2),xi(1:2)
COMPLEX(QPC) tanht,sinht,cosht,bidon,bidon1,bidon2,phi(1:2),ellipcF,ellipcE
COMPLEX(QPC) zdep,zfin
LOGICAL axer
REAL(QP) d2en0

call ecrit(bla1,'secteur',real(sect,DP))
call oang
call ecrit(bla1,'o=',o(1:3))

cosht=z/2.0_qp
if(real(z).GE.0.0_qp) sinht= sqrt(cosht**2.0_qp-1.0_qp)
if(real(z)<0.0_qp) sinht=-sqrt(cosht**2.0_qp-1.0_qp)
tanht=sinht/cosht

if(sect==0)then
 dspecpro(1)=cmplx(0.0_qp,0.0_qp,QPC)
 dspecpro(2)=cmplx(0.0_qp,0.0_qp,QPC)
 dspecpro(3)=cmplx(0.0_qp,0.0_qp,QPC)
elseif(sect==1)then
 if(xq**2.0_dp/4.0_dp>x0) STOP 'Le secteur 1 est réduit à un point'
 r1(1)=2.0_qp*cacaf(pis2c,tanht)/cosht/4.0_qp
 r2(1)=2.0_qp*cacae(pis2c,tanht)*cosht/2.0_qp
 dspecpro(1)=PI*r2(1)/xq
 dspecpro(2)=PI*(r2(1)-2.0_qp*r1(1))/xq
 dspecpro(3)=cmplx(0.0_qp,0.0_qp,QPC)
elseif(sect==2)then
 d2en0=-64*x0*(1+x0**2)+48*(1+x0**2)*xq**2-12*x0*xq**4+xq**6    !derivee seconde en k=0
 if(d2en0>0.0_qp) STOP 'Le secteur 2 est réduit à un point'

 if((real(al)<real(omil)).AND.(imag(z)<imag(al)).AND.(real(z)<real(al)).AND.(.FALSE.))then
   write(6,*)
   write(6,*)"Bascule"
   write(6,*)
   zdep=zbascule
   xi=xibasc
   phi=phibasc
   r3=r3basc
   zfin=z
   axer=.FALSE.
   call suivi(zdep,zfin,axer,2,xi,phi,r3,.FALSE.)
  else
   zdep=omil
   zfin=z
   axer=.TRUE.
   call suivi(zdep,zfin,axer,2,xi,phi,r3,.FALSE.)
  endif
  call ecrit(bla1,'xi=',xi(1:2))

 ellipcE=cacae(pis2c,tanht)
 ellipcF=cacaf(pis2c,tanht)
! Partie décroissante
 r1(1)=(ellipcF-cacaf(phi(1),tanht))/cosht/4.0_qp
 r2(1)=(ellipcE-cacae(phi(1),tanht))*cosht/2.0_qp+xi(1)*r3(1)/cosht

! Partie croissante
 r1(2)=(ellipcF-cacaf(phi(2),tanht))/cosht/4.0_qp
 r2(2)=(ellipcE-cacae(phi(2),tanht))*cosht/2.0_qp+xi(2)*r3(2)/cosht

 open(16,file='test2.dat',POSITION='APPEND')
  write(16,*)real(z),&
  real(xi(1)),imag(xi(1)),real(phi(1)),imag(phi(1)),&
  real(r1(1)),imag(r1(1)),real(r2(1)),imag(r2(1)),real(r3(1)),imag(r3(1)),&
  real(bidon),imag(bidon),real(bidon1),imag(bidon1),real(bidon2),imag(bidon2)
 close(16)

 dspecpro(1)=PI*sum(r2)/xq
 dspecpro(2)=PI*sum(r2-2.0_qp*r1)/xq
 dspecpro(3)=PI*sum(r3)/xq

else !secteur 3

 if(imag(z)<0.0_qp)then
  if((real(al)<real(omil)).AND.(imag(z)<imag(al)).AND.(real(z)<real(al)).AND.(.FALSE.))then
   write(6,*)
   write(6,*)"Bascule"
   write(6,*)
   zdep=zbascule
   xi=xibasc
   phi=phibasc
   r3=r3basc
   zfin=z
   axer=.FALSE.
   call suivi(zdep,zfin,axer,3,xi,phi,r3,.FALSE.)
  else
   axer=.TRUE.
   zdep=cmplx(2.0_qp*o(3),0.0_qp,kind=qpc)
   zfin=z
   call suivi(zdep,zfin,axer,3,xi,phi,r3,.FALSE.)
  endif
  call ecrit(bla1,'xi=',xi(1:2))

  ellipcE=cacae(pis2c,tanht)
  ellipcF=cacaf(pis2c,tanht)

  r1(1)=(ellipcF-cacaf(phi(1),tanht))/cosht/4.0_qp
  r2(1)=(ellipcE-cacae(phi(1),tanht))*cosht/2.0_qp+xi(1)*r3(1)/cosht


!  write(6,*)"----------------------------"
!  write(6,*)"tanht=",tanht,"phi(1)=",phi(1)
! write(6,*)
!
!  bidon=cacae2(pis2c,tanht)
!  bidon1=cacaf2(pis2c,tanht)
!  write(6,*)"cacae2(pis2c,tanht)=",bidon
!  write(6,*)"cacaf2(pis2c,tanht)=",bidon1
!
!  bidon=cacae(pis2c,tanht)
!  bidon1=cacaf(pis2c,tanht)
!  write(6,*)"cacae(pis2c,tanht)=",bidon
!  write(6,*)"cacaf(pis2c,tanht)=",bidon1
!  
!  bidon=cacae2(phi(1),tanht)
!  bidon1=cacaf2(phi(1),tanht)
!  write(6,*)"cacae2(phi(1),tanht)=",bidon
!  write(6,*)"cacaf2(phi(1),tanht)=",bidon1
!  write(6,*)
!  
!  bidon=cacae(phi(1),tanht)
!  bidon1=cacaf(phi(1),tanht)
!  write(6,*)"cacae(phi(1),tanht)=",bidon
!  write(6,*)"cacaf(phi(1),tanht)=",bidon1
!  write(6,*)"----------------------------"



! bidon=xi(1)*r3(1)/cosht
! bidon1=(ellipcE-cacae2(phi(1),tanht))*cosht/2.0_qp
! bidon2=(ellipcE-cacae2(phi(1),tanht))*cosht/2.0_qp+xi(1)*r3(1)/cosht
 open(17,file='test3.dat',POSITION='APPEND')
  write(17,*)real(z),imag(z),&
  real(xi(1)),imag(xi(1)),real(phi(1)),imag(phi(1)),&
!  real(r1(1)),imag(r1(1)),&
  real(r2(1)),imag(r2(1)),&
!  real(r3(1)),imag(r3(1)),&
  real(bidon),imag(bidon),&
  real(bidon1),imag(bidon1),&
  real(bidon2),imag(bidon2)
 close(17)

  dspecpro(1)=PI*r2(1)/xq
  dspecpro(2)=PI*(r2(1)-2.0_qp*r1(1))/xq
  dspecpro(3)=PI*r3(1)/xq
 endif

endif
END FUNCTION dspecpro
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE suivi(zdep,zfin,axer,sect,xi,phi,r3,ecrire)
USE recettes
IMPLICIT NONE
INTEGER, INTENT(IN) :: sect
LOGICAL, INTENT(IN) :: axer,ecrire
COMPLEX(QPC), INTENT(IN)  :: zdep,zfin
COMPLEX(QPC), INTENT(INOUT) :: xi(1:2),phi(1:2),r3(1:2)

INTEGER it,nt,p
REAL(QP) t,dt,xir(1:3),rephi,imphi
COMPLEX(QPC) :: zc,xic(1:3),xia(1:3),r1(1:2),r2(1:2)
COMPLEX(QPC) tanht,sinht,cosht
COMPLEX(QPC) x(1:2),phic(1:2),r3c(1:2),bidon
COMPLEX(QPC) ellipcE,ellipcF

nt=10000


if(ecrire)then
 nt=10000
endif

dt=1.0_qp/real(nt)

if(sect==2)then

 if(axer)then
!On part de l'axe réel, il faut initialiser xi, r et phi
  if(abs(imag(zdep))>1.e-15_qp) stop "Le point de départ de suivi n'est pas sur l axe réel"
  call bornesbrutal(real(zdep),2,xir)
  xi(1:2)=cmplx(xir(1:2),0.0_qp,QPC)
 
  cosht=real(zdep)/2.0_qp
  sinht= sqrt(cosht**2.0_qp-1.0_qp)
  tanht=sinht/cosht
! Partie décroissante
  x(1)=-xi(1)/sinht
  phi(1)=asin(x(1))
  r3(1)= -sinht*sqrt((1.0_qp-x(1)**2.0_qp)/(1.0_qp-x(1)**2.0_qp*tanht**2.0_qp))/2.0_qp
! Partie croissante
  x(2)=xi(2)/sinht
  phi(2)=asin(x(2))
  r3(2)=  sinht*sqrt((1.0_qp-x(2)**2.0_qp)/(1.0_qp-x(2)**2.0_qp*tanht**2.0_qp))/2.0_qp
 endif

 do it=0,nt
  t=it*dt

! suivi de xi
  zc=(1.0_qp-t)*zdep+t*zfin
  call racinescubc(zc,xic,bidon)
  xic=xic-x0+xq**2.0_qp/4.0_qp
  xi(1)=xic(minloc(abs(xic-xi(1)),DIM=1))
  xi(2)=xic(minloc(abs(xic-xi(2)),DIM=1))

! suivi de r et phi
  cosht=zc/2.0_qp
! Compensation du saut de signe di sinht lorsque real(zc)<0
  if(real(zc).GE.0.0_qp) sinht= sqrt(cosht**2.0_qp-1.0_qp)
  if(real(zc)<0.0_qp) sinht=-sqrt(cosht**2.0_qp-1.0_qp)
  tanht=sinht/cosht
! Partie décroissante
  x(1)=-xi(1)/sinht
  p=minloc(abs((/phi(1)-asin(x(1)),phi(1)-(-PI-asin(x(1))),phi(1)-(PI-asin(x(1)))/)),DIM=1)
  if(p==1)then
   phi(1)=asin(x(1))
  elseif(p==2)then
   phi(1)=-PI-asin(x(1))
  elseif(p==3)then
   phi(1)=PI-asin(x(1))
  endif
  r3c(1)= -sinht*sqrt((1.0_qp-x(1)**2.0_qp)/(1.0_qp-x(1)**2.0_qp*tanht**2.0_qp))/2.0_qp
! Partie croissante
  x(2)=xi(2)/sinht
  p=minloc(abs((/phi(2)-asin(x(2)),phi(2)-(-PI-asin(x(2))),phi(2)-(PI-asin(x(2)))/)),DIM=1)
  if(p==1)then
   phi(2)=asin(x(2))
  elseif(p==2)then
   phi(2)=-PI-asin(x(2))
  elseif(p==3)then
   phi(2)=PI-asin(x(2))
  endif
  r3c(2)=  sinht*sqrt((1.0_qp-x(2)**2.0_qp)/(1.0_qp-x(2)**2.0_qp*tanht**2.0_qp))/2.0_qp

! Compensation des sauts de signe de r3
  r3(1)=r3c(1)*(-1.0_qp)**(minloc(abs((/r3(1)+r3c(1),r3(1)-r3c(1)/)),DIM=1))
  r3(2)=r3c(2)*(-1.0_qp)**(minloc(abs((/r3(2)+r3c(2),r3(2)-r3c(2)/)),DIM=1))

  if(ecrire)then

   ellipcE=cacae(pis2c,tanht)
   ellipcF=cacaf(pis2c,tanht)
! Partie décroissante
   r1(1)=(ellipcF-cacaf(phi(1),tanht))/cosht/4.0_qp
   r2(1)=(ellipcE-cacae(phi(1),tanht))*cosht/2.0_qp+xi(1)*r3(1)/cosht
  
! Partie croissante
   r1(2)=(ellipcF-cacaf(phi(2),tanht))/cosht/4.0_qp
   r2(2)=(ellipcE-cacae(phi(2),tanht))*cosht/2.0_qp+xi(2)*r3(2)/cosht

   open(67,file="suivi_sect2.dat",POSITION='APPEND')
   write(67,*) real(zc),imag(zc),& !2
               real(x(1)),imag(x(1)),real(phi(1)),imag(phi(1)),& !6
               real(x(2)),imag(x(2)),real(phi(2)),imag(phi(2)),& !10
               real(xic(1)),imag(xic(1)),real(xic(2)),imag(xic(2)),real(xic(3)),imag(xic(3)),& !16
               real(r3(1)),imag(r3(1)),real(r3(2)),imag(r3(2)),& !20
               real(r1(1)),imag(r1(1)),real(r2(1)),imag(r2(1)),& !24
               real(r1(2)),imag(r1(2)),real(r2(2)),imag(r2(2))   !28
   close(67)
  endif

 enddo

elseif(sect==3)then

 if(axer)then
!On part de l'axe réel, il faut initialiser xi, r et phi
  if(abs(imag(zdep))>1.e-15_qp) stop "Le point de départ de suivi n'est pas sur l axe réel"
  call bornesbrutal(real(zdep),3,xir)
  xi(1)=cmplx(xir(1),0.0_qp,QPC)
  call racinescubc(zdep,xia,bidon)

  cosht=real(zdep)/2.0_qp
  sinht= sqrt(cosht**2.0_qp-1.0_qp)
  tanht=sinht/cosht
  x(1)=xi(1)/sinht
  phi(1)=asin(x(1))
  r3(1)= sinht*sqrt((1.0_qp-x(1)**2.0_qp)/(1.0_qp-x(1)**2.0_qp*tanht**2.0_qp))/2.0_qp
 endif

 do it=1,nt
  t=it*dt

! suivi de xi
  zc=(1.0_qp-t)*zdep+t*zfin
  call racinescubc(zc,xic,bidon)

!  mindistrac=min(abs((/xic(1)-xic(2),xic(1)-xic(3),xic(3)-xic(2)/)))
  xia(1)=xic(minloc(abs(xic-xia(1)),DIM=1))
  xia(2)=xic(minloc(abs(xic-xia(2)),DIM=1))
  xia(3)=xic(minloc(abs(xic-xia(3)),DIM=1))

  xic=xic-x0+xq**2.0_qp/4.0_qp
  xi(1)=xic(minloc(abs(xic-xi(1)),DIM=1))
!  write(6,*)"xi(1),xic(1:2)=",real(xi(1)),real(xic)
!  write(6,*)"xi(1),xic(1:2)=",imag(xi(1)),imag(xic)

! suivi de r et phi
  cosht=zc/2.0_qp
  if(real(zc).GE.0.0_qp) sinht= sqrt(cosht**2.0_qp-1.0_qp)
  if(real(zc)<0.0_qp) sinht=-sqrt(cosht**2.0_qp-1.0_qp)
  tanht=sinht/cosht
! Partie croissante
  x(1)=xi(1)/sinht
  p=minloc(abs((/phi(1)-asin(x(1)),phi(1)-(-PI-asin(x(1))),phi(1)-(PI-asin(x(1)))/)),DIM=1)
  if(p==1)then
   phi(1)=asin(x(1))
  elseif(p==2)then
   phi(1)=-PI-asin(x(1))
  elseif(p==3)then
   phi(1)=PI-asin(x(1))
  endif
  r3c(1)=  sinht*sqrt((1.0_qp-x(1)**2.0_qp)/(1.0_qp-x(1)**2.0_qp*tanht**2.0_qp))/2.0_qp

! Compensation des sauts de signe de r3 et phi
  r3(1)=r3c(1)*(-1.0_qp)**(minloc(abs((/r3(1)+r3c(1),r3(1)-r3c(1)/)),DIM=1))
!  write(6,*)"r3(1)=",r3(1)
  if(ecrire)then

   ellipcE=cacae(pis2c,tanht)
   ellipcF=cacaf(pis2c,tanht)
   r1(1)=(ellipcF-cacaf(phi(1),tanht))/cosht/4.0_qp
   r2(1)=(ellipcE-cacae(phi(1),tanht))*cosht/2.0_qp+xi(1)*r3(1)/cosht

   open(67,file="suivi_sect3.dat",POSITION='APPEND')
   write(67,*) real(zc),imag(zc),& !2
               real(x(1)),imag(x(1)),real(phi(1)),imag(phi(1)),& !6
               real(xic(1)),imag(xic(1)),real(xic(2)),imag(xic(2)),real(xic(3)),imag(xic(3)),& !12
               real(r3(1)),imag(r3(1)),& !14
               real(r1(1)),imag(r1(1)),& !16
               real(r2(1)),imag(r2(1))   !18
   close(67)
  endif

 enddo
 xi(2)=-1e200_qp
else
 STOP 'Mauvais secteur dans suivi'
endif
END SUBROUTINE suivi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE emax(x,arg,ene,den)
!Renvoie E(x,u=1)-om et sa dérivée par rapport à x
USE nrtype

IMPLICIT NONE
REAL(QP), INTENT(IN) :: x
REAL(QP), DIMENSION(:), INTENT(IN) :: arg 
REAL(QP), INTENT(OUT) :: ene,den
REAL(QP)  xiP,xiM,eP,eM,om

om=arg(1)

xiP=x**2.0_qp+x*xq+xq**2.0_qp/4.0_qp-x0
xiM=x**2.0_qp-x*xq+xq**2.0_qp/4.0_qp-x0

eP=sqrt(xiP**2.0_qp+1.0_qp)
eM=sqrt(xiM**2.0_qp+1.0_qp)
ene=eP+eM-om
den=xiP/eP*(2.0_qp*x+xq)+xiM/eM*(2.0_qp*x-xq)
END SUBROUTINE emax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION En(x,u)
USE nrtype

IMPLICIT NONE
REAL(QP), INTENT(IN) :: x,u
REAL(QP)  En 
REAL(QP)  xiP,xiM,eP,eM

xiP=x**2.0_qp+x*xq*u+xq**2.0_qp/4.0_qp-x0
xiM=x**2.0_qp-x*xq*u+xq**2.0_qp/4.0_qp-x0

eP=sqrt(xiP**2.0_qp+1.0_qp)
eM=sqrt(xiM**2.0_qp+1.0_qp)
En=eP+eM
END FUNCTION En 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE deriv(x,arg,dEdx,ddEddx)
!Renvoie dE(x,u=1)/dx et d2E(x,u=1)/dx2
USE nrtype

IMPLICIT NONE
REAL(QP), INTENT(IN)  :: x
REAL(QP), DIMENSION(:), INTENT(IN)  :: arg
REAL(QP), INTENT(OUT) :: dEdx,ddEddx
REAL(QP) xiP,xiM,eP,eM

xiP=x**2.0_qp+x*xq+xq**2.0_qp/4.0_qp-x0
xiM=x**2.0_qp-x*xq+xq**2.0_qp/4.0_qp-x0

eP=sqrt(xiP**2.0_qp+1.0_qp)
eM=sqrt(xiM**2.0_qp+1.0_qp)

dEdx=xiP/eP*(2.0_qp*x+xq)+xiM/eM*(2.0_qp*x-xq)

ddEddx=(2.0_qp*x+xq)**2.0_qp/eP+(2.0_qp*x-xq)**2.0_qp/eM &
&       +2.0_qp*xiP/eP+2.0_qp*xiM/eM &
&       -xiP**2.0_qp/eP**3.0_qp*(2.0_qp*x+xq)**2.0_qp &
&       -xiM**2.0_qp/eM**3.0_qp*(2.0_qp*x-xq)**2.0_qp
END SUBROUTINE deriv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE oang
!Calcule la position des trois points anguleux o(1),o(2),o(3) de la densité spectrale
!ainsi que la valeur de x, alias k, où min_x E(u=1) est atteint (valeur inscrite dans o(4)) 
USE nrtype
USE recettes
IMPLICIT NONE
REAL(QP) d2en0,x2
o(3)=2.0_qp*sqrt(1.0_qp+(xq**2.0_qp/4.0_qp-x0)**2.0_qp)
if(xq**2.0_qp/4.0_qp>x0)then
 ptbranchmt=1
 o(1)=o(3)
 o(2)=o(3)
 o(4)=-1e-31_qp
else
 o(1)=2.0_qp !min_k,uE=2Delta
 d2en0=-64*x0*(1+x0**2)+48*(1+x0**2)*xq**2-12*x0*xq**4+xq**6    !derivee seconde en k=0
 if(d2en0>0.0_qp)then           !Si elle est positive, le min est en zero
  ptbranchmt=2
  o(2)=o(3)
  o(4)=-1e-31_qp
 else
  ptbranchmt=3
  x2=rtsafeq(deriv,(/bidon/),1.e-14_qp,sqrt(x0)+2.0_qp,1.e-19_qp)!Le min est entre 0 et rac(x0)
  o(2)=En(x2,1.0_qp)
  o(4)=x2
 endif
endif
END SUBROUTINE oang
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE alpha(al)
! Calcule le point de branchement dans le prolongement analytique
!par les fenetres 2 et 3. S'obtient comme la solution de Re>0 et Im<0 à discr=0 (c-a-d, il existe une racine double) où discr
!est le discriminant de l'équation cubique donnant les bornes d'intégration
USE recettes, ONLY : zroots
IMPLICIT NONE
COMPLEX(QPC), INTENT(OUT) :: al
COMPLEX(QPC) :: coeff(1:5),racines(1:4),bidon(1:3),discr
REAL(QP) :: coeffr(1:5)

coeffr(1)= 256*x0**3*xq**6 - 192*x0**2*xq**8 + 256*x0**4*xq**8 + 48*x0*xq**10 - 256*x0**3*xq**10 - 4*xq**12 & 
     +  96*x0**2*xq**12 - 16*x0*xq**14 + xq**16
coeffr(2)=-432*xq**4 - 192*x0**2*xq**4 - 480*x0*xq**6 - 256*x0**3*xq**6 + 132*xq**8 + 64*x0**2*xq**8 & 
     +  16*x0*xq**10 - 4*xq**12
coeffr(3)=48*x0*xq**2 + 132*xq**4 + 96*x0**2*xq**4 + 16*x0*xq**6 + 6*xq**8
coeffr(4)=-4 - 16*x0*xq**2 - 4*xq**4
coeffr(5)=1

coeff=cmplx(coeffr,0.0_qp,kind=qpc)

call zroots(coeff,racines,.TRUE.)
al=sqrt(racines(1))

END SUBROUTINE alpha
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE racinescubbrut(zc,r,discr)
USE recettes, ONLY : zroots
IMPLICIT NONE
COMPLEX(QPC), INTENT(IN) :: zc
COMPLEX(QPC), INTENT(OUT) :: r(1:3),discr
COMPLEX(QPC) aa,bb,cc,dd,coeff(1:4)

aa=16.0_qp*xq**2.0_qp
bb=-4.0_qp*(zc**2.0_qp+8.0_qp*x0*xq**2.0_qp-2.0_qp*xq**4.0_qp)
cc=zc**2.0_qp*(8.0_qp*x0-6.0_qp*xq**2.0_qp)+(xq**3.0_qp-4.0_qp*x0*xq)**2.0_qp
dd=zc**4.0_qp-zc**2.0_qp*(16.0_qp+(xq**2.0_qp-4.0_qp*x0)**2.0_qp)/4.0_qp

coeff=(/dd,cc,bb,aa/)
call zroots(coeff,r,.TRUE.)

END SUBROUTINE racinescubbrut
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE racinescubc(zc,r,discr)
IMPLICIT NONE
COMPLEX(QPC), INTENT(IN) :: zc
COMPLEX(QPC), INTENT(OUT) :: r(1:3),discr
COMPLEX(QPC) aa,bb,cc,dd,pp,qq,uu,vv,j

j =cmplx(-0.5_qp,sqrt(3.0_qp)/2.0_qp,QPC)

aa=16.0_qp*xq**2.0_qp
bb=-4.0_qp*(zc**2.0_qp+8.0_qp*x0*xq**2.0_qp-2.0_qp*xq**4.0_qp)
cc=zc**2.0_qp*(8.0_qp*x0-6.0_qp*xq**2.0_qp)+(xq**3.0_qp-4.0_qp*x0*xq)**2.0_qp
dd=zc**4.0_qp-zc**2.0_qp*(16.0_qp+(xq**2.0_qp-4.0_qp*x0)**2.0_qp)/4.0_qp

pp=-bb**2.0_qp/(3.0_qp*aa**2.0_qp)+cc/aa
qq=bb/(27.0_qp*aa)*(2.0_qp*bb**2.0_qp/aa**2.0_qp-9.0_qp*cc/aa)+dd/aa

discr=-4.0_qp*pp**3.0_qp-27.0_qp*qq**2.0_qp

uu =((-qq +sqrt(-discr /27.0_qp))/2.0_qp)**(1.0_qp/3.0_qp)
vv = -pp /(3.0_qp*uu )

r(1) =uu    +vv    -bb/(3.0_qp*aa)
r(2) =j*uu  +vv/j  -bb/(3.0_qp*aa)
r(3) =j*j*uu+vv/j/j-bb/(3.0_qp*aa)

END SUBROUTINE racinescubc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION ffr(y,om)
IMPLICIT NONE
REAL(QP), INTENT(IN) :: y,om
REAL(QP) ffr
ffr=om**2.0_qp-2.0_qp*((y+xq**2.0_qp/4.0_qp-x0)**2.0_qp+y*xq**2.0_qp+1.0_qp)
END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION ffc(y,om)
IMPLICIT NONE
COMPLEX(QPC), INTENT(IN) :: y,om
COMPLEX(QPC) ffc
ffc=om**2.0_qp-2.0_qp*((y+xq**2.0_qp/4.0_qp-x0)**2.0_qp+y*xq**2.0_qp+1.0_qp)
END FUNCTION
END MODULE modpro
