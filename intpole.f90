MODULE intpole
USE nrtype
USE modsim
USE vars
USE intldc
IMPLICIT NONE
REAL(QP) :: ccheck,Theta,Xx
CHARACTER(len=90) fichier
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION selfEpole(k,zk)
USE modsim
USE intldc
IMPLICIT NONE
REAL(QP), INTENT(IN) :: k,zk
REAL(QP) :: selfEpole(1:6)

REAL(QP) :: qmin,qmax,EPS,dq
INTEGER iq,nq,nqeff
INTEGER nn,taille
 
open(11,file=trim(fichier)//".info")
 read(11,*)x0,ccheck,Theta,Xx,qmin,qmax,nq,nn
 write(6,*)"x0,ccheck,Theta,Xx,qmin,qmax,nq,nn=",x0,ccheck,Theta,Xx,qmin,qmax,nq,nn
close(11)

dq=(qmax-qmin)/nq

inquire(file=trim(fichier)//".dat", size=taille)
nqeff=taille/nn

write(6,*)"dq,nqeff,nqeff*dq=",dq,nqeff,qmin+nqeff*dq

EPS =1.0e-8_qp
selfEpole=qromovq(integrandeq  ,0.0_qp,nqeff*dq,6,(/bidon/),midpntvq,EPS)
selfEpole=selfEpole*2.0_qp*PI
!Comparison with the perturbative energy correction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
 FUNCTION integrandeq(q,arg,m) 
 IMPLICIT NONE

 INTEGER,  INTENT(IN) :: m
 REAL(QP), INTENT(IN), DIMENSION(:) :: q,arg
 REAL(QP), DIMENSION(size(q),m) :: integrandeq

 REAL(QP)  ptq(1:5),ptom(1:5),ptM(1:3,1:5),ptdM(1:3,1:5)
 REAL(QP)  qs,om,Ma(1:3),dM(1:3),MatCat(1:3),errom,errM(1:3),errdM(1:3)
 REAL(QP)  ddet
 REAL(QP)  IuM(1:3),IuP(1:3)
 INTEGER is,ib,iq

 do is=1,size(q)
  qs=q(is)
  if(qs<0.1)then
   om=sqrt(2.0_qp*ccheck)*qs
   Ma(1)=om**2*Xx**2/Theta/8.0_qp
   Ma(2)=Theta/2.0_qp
   Ma(3)=-Xx*om/4.0_qp
   dM(1)=-Theta*om/4.0_qp
   dM(2)=0.0_qp
   dM(3)=-Xx/4.0_qp
  else
   open(10,file=trim(fichier)//".dat",action="read",access="direct",form="unformatted",recl=nn)
   
   iq=floor((qs-qmin)/dq)+1
   read(10,rec=iq-2)ptq(1),ptom(1),ptM(:,1),ptdM(:,1)
   read(10,rec=iq-1)ptq(2),ptom(2),ptM(:,2),ptdM(:,2)
   read(10,rec=iq  )ptq(3),ptom(3),ptM(:,3),ptdM(:,3)
   read(10,rec=iq+1)ptq(4),ptom(4),ptM(:,4),ptdM(:,4)
   read(10,rec=iq+2)ptq(5),ptom(5),ptM(:,5),ptdM(:,5)
  
!   write(6,*)qs,ptq(3:4)
   call polint(ptq,ptom,qs,om,errom)
   do ib=1,3
    call polint(ptq, ptM(ib,:),qs, Ma(ib),errM(ib))
    call polint(ptq,ptdM(ib,:),qs,dM(ib),errdM(ib))
   enddo
   close(10) 
  endif

  ddet=Ma(1)*dM(2)+Ma(2)*dM(1)-2.0_qp*Ma(3)*dM(3)

  MatCat(1)=(Ma(1)+Ma(2))/2.0_qp+Ma(3)
  MatCat(2)=(Ma(1)+Ma(2))/2.0_qp-Ma(3)
  MatCat(3)=(Ma(2)-Ma(1))/2.0_qp

  xi0=k**2+qs**2-x0
  xiP=k**2+qs**2+2.0_qp*k*qs-x0
  xiM=k**2+qs**2-2.0_qp*k*qs-x0
  epsP=sqrt(xiP**2+1.0_qp)
  epsM=sqrt(xiM**2+1.0_qp)

  xmin=epsP-xiP
  xmax=epsM-xiM

  IuP=real(Iuanaly(om+zk,k,qs))
  IuM=real(Iuanaly(om-zk,k,qs))

  integrandeq(is,1)=-MatCat(2)*IuM(1)/ddet
  integrandeq(is,2)=-MatCat(1)*IuM(2)/ddet
  integrandeq(is,3)= MatCat(3)*IuM(3)/ddet

  integrandeq(is,4)= MatCat(1)*IuP(2)/ddet
  integrandeq(is,5)= MatCat(2)*IuP(1)/ddet
  integrandeq(is,6)= MatCat(3)*IuP(3)/ddet
  write(6,*)"qs,integrandeq=",qs,integrandeq(is,1)
  integrandeq(is,:)=qs**2*integrandeq(is,:)
 enddo

 END FUNCTION integrandeq 
END FUNCTION selfEpole 
END MODULE intpole
