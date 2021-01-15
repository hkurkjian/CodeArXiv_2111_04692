MODULE intpole
USE nrtype
USE modsim
USE vars
USE intldc
IMPLICIT NONE
REAL(QP) :: ccheck,Theta,Xx
CHARACTER(len=90) fichier
LOGICAL blaPole
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION selfEpole(k,zk)
USE modsim
USE intldc
IMPLICIT NONE
REAL(QP), INTENT(IN) :: k,zk
REAL(QP) :: selfEpole(1:6)

REAL(QP) :: q,qmin,qmax,EPS,kmin,kmax,dk,dq,bq(1:7),qThr
INTEGER ix,nx,iq,nq,nqeff,nbq
INTEGER nn,taille
 
open(11,file=trim(fichier)//".info")
 read(11,*)x0,ccheck,Theta,Xx,qmin,qmax,nq,nn
 if (blaPole)then
  write(6,*)"x0,ccheck,Theta,Xx,qmin,qmax,nq,nn,k,zk=",x0,ccheck,Theta,Xx,qmin,qmax,nq,nn,k,zk
 endif
close(11)

! Select threshhold q below which omq is approximated by c q
qThr=0.0205_qp

dq=(qmax-qmin)/nq

inquire(file=trim(fichier)//".dat", size=taille)
nqeff=taille/nn

! Calculate q bounds
call boundsQ(k,zk,bq,nbq)

if (blaPole)then
  write(6,*)"dq,nqeff,nqeff*dq=",dq,nqeff,qmin+nqeff*dq
  write(6,*)"Calculating all q-bounds for k,zk=",k,zk
  write(6,*)"q bounds found:",nbq,bq(1:nbq)
endif

! open(16,file="testBounds"//trim(fichier)//".dat",position="append")
! write(16,*)k,zk,bq(1:nbq)
! close(16)

EPS =1.0e-8_qp
! selfEpole=qromovq(integrandeq  ,0.0_qp,nqeff*dq,6,(/bidon/),midpntvq,EPS)
selfEpole=0.1_qp
selfEpole=selfEpole*2.0_qp*PI
!Comparison with the perturbative energy correction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
 FUNCTION integrandeq(q,arg,m) 
 IMPLICIT NONE

 INTEGER,  INTENT(IN) :: m
 REAL(QP), INTENT(IN), DIMENSION(:) :: q,arg
 COMPLEX(QP), DIMENSION(size(q),m) :: integrandeq

 REAL(QP)  ptq(1:5),ptom(1:5),ptM(1:3,1:5),ptdM(1:3,1:5)
 REAL(QP)  qs,om,Ma(1:3),dM(1:3),MatCat(1:3),errom,errM(1:3),errdM(1:3)
 REAL(QP)  ddet
 REAL(QP)  IuM(1:3),IuP(1:3)
 INTEGER is,ib,iq

 do is=1,size(q)
  qs=q(is)
  if(qs<0.1)then
  ! For small q, use analytic formulas
   om=sqrt(2.0_qp*ccheck)*qs
   Ma(1)=om**2*Xx**2/Theta/8.0_qp
   Ma(2)=Theta/2.0_qp
   Ma(3)=-Xx*om/4.0_qp
   dM(1)=-Theta*om/4.0_qp
   dM(2)=0.0_qp
   dM(3)=-Xx/4.0_qp
  else
  ! Else, use interpolation from data file
   open(10,file=trim(fichier)//".dat",action="read",access="direct",form="unformatted",recl=nn)
   
  ! Select 5 data points around q
   iq=floor((qs-qmin)/dq)+1
   read(10,rec=iq-2)ptq(1),ptom(1),ptM(:,1),ptdM(:,1)
   read(10,rec=iq-1)ptq(2),ptom(2),ptM(:,2),ptdM(:,2)
   read(10,rec=iq  )ptq(3),ptom(3),ptM(:,3),ptdM(:,3)
   read(10,rec=iq+1)ptq(4),ptom(4),ptM(:,4),ptdM(:,4)
   read(10,rec=iq+2)ptq(5),ptom(5),ptM(:,5),ptdM(:,5)
  
!   write(6,*)qs,ptq(3:4)
  ! Interpolation scheme for omega_q
   call polint(ptq,ptom,qs,om,errom)
   ! Interpolation scheme for each matrix component and derivative
   do ib=1,3
    call polint(ptq, ptM(ib,:),qs, Ma(ib),errM(ib))
    call polint(ptq,ptdM(ib,:),qs,dM(ib),errdM(ib))
   enddo
   close(10) 
  endif

  ! Derivative of the determinant
  ddet=Ma(1)*dM(2)+Ma(2)*dM(1)-2.0_qp*Ma(3)*dM(3)

  ! Cartesian basis
  MatCat(1)=(Ma(1)+Ma(2))/2.0_qp+Ma(3)
  MatCat(2)=(Ma(1)+Ma(2))/2.0_qp-Ma(3)
  MatCat(3)=(Ma(2)-Ma(1))/2.0_qp

  ! Energy functions for analytic angle integration
  xi0=k**2+qs**2-x0
  xiP=k**2+qs**2+2.0_qp*k*qs-x0
  xiM=k**2+qs**2-2.0_qp*k*qs-x0
  epsP=sqrt(xiP**2+1.0_qp)
  epsM=sqrt(xiM**2+1.0_qp)
  xmin=epsP-xiP
  xmax=epsM-xiM

  ! Angle integration
  IuP=Iuanaly(om+zk,k,qs)
  IuM=Iuanaly(om-zk,k,qs)

  ! 1->2 process integrand
  integrandeq(is,1)=-MatCat(2)*IuM(1)/ddet
  integrandeq(is,2)=-MatCat(1)*IuM(2)/ddet
  integrandeq(is,3)=-MatCat(3)*IuM(3)/ddet

  ! 0->3 process integrand
  integrandeq(is,4)= MatCat(1)*IuP(2)/ddet
  integrandeq(is,5)= MatCat(2)*IuP(1)/ddet
  integrandeq(is,6)=-MatCat(3)*IuP(3)/ddet

  ! multiply by Jacobian q^2
  write(6,*)"qs,integrandeq=",qs,integrandeq(is,1)
  integrandeq(is,:)=qs**2*integrandeq(is,:)
 enddo

 END FUNCTION integrandeq 
 SUBROUTINE intOmQ(qVal,omq)
  IMPLICIT NONE
  REAL(QP), INTENT(IN) :: qVal
  REAL(QP), INTENT(OUT) :: omq

  INTEGER iq
  REAL(QP)  ptq(1:5),ptom(1:5),ptM(1:3,1:5),ptdM(1:3,1:5), errom

  if(qVal<qThr)then
    ! For small q, use analytic formulas
     omq=sqrt(2.0_qp*ccheck)*qVal
  else

    ! write(6,*)"Calculating omega for q=",qVal

    ! Open file for interpolation of omega_q
    open(12,file=trim(fichier)//".dat",action="read",access="direct",form="unformatted",recl=nn)
    
    ! Read 5 points around qVal
    iq=floor((qVal-qmin)/dq)+1
    read(12,rec=iq-2)ptq(1),ptom(1),ptM(:,1),ptdM(:,1)
    read(12,rec=iq-1)ptq(2),ptom(2),ptM(:,2),ptdM(:,2)
    read(12,rec=iq  )ptq(3),ptom(3),ptM(:,3),ptdM(:,3)
    read(12,rec=iq+1)ptq(4),ptom(4),ptM(:,4),ptdM(:,4)
    read(12,rec=iq+2)ptq(5),ptom(5),ptM(:,5),ptdM(:,5)

    ! Close file
    close(12)

    ! Interpolation scheme for omega_q
    call polint(ptq,ptom,qVal,omq,errom)
  endif

 END SUBROUTINE intOmQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE boundsQ(k,zk,bq,nbq)
 ! Calculate all bounds on q for integration above emission continuum
  USE recettes, ONLY : rtsafe
  IMPLICIT NONE
  REAL(QP), INTENT(IN) :: k,zk
  REAL(QP), INTENT(OUT) :: bq(1:7) ! q bounds, max lengths = 7
  INTEGER, INTENT(OUT) :: nbq ! number of q bounds (max 7)

  REAL(QP) qMax, omqMax, mMax(1:3), dmMax(1:3)
  REAL(QP) qT, omqT, mT(1:3), dmT(1:3)
  REAL(QP) x1, x2, y1P, y2P, y1M, y2M
  REAL(QP) q0, qM(1:3), qX(1:3)
  INTEGER iq, nqM, nqP, iP, iM, iT
  ! REAL(QP) omqTEST, yTEST,y1TEST,y2TEST, dyTEST,dy1TEST,dy2TEST
  
  ! Initialize bounds (1e50 if not important)
  q0=   1.0e50_qp
  qM(1)=1.0e50_qp
  qM(2)=1.0e50_qp
  qM(3)=1.0e50_qp
  qX(1)=1.0e50_qp
  qX(2)=1.0e50_qp
  qX(3)=1.0e50_qp

  ! Open file for interpolation of omega_q
  open(11,file=trim(fichier)//".dat",action="read",access="direct",form="unformatted",recl=nn)
  ! Read last value for omqMax
  read(11,rec=nq)qMax,omqMax,mMax(:),dmMax(:)

  ! Initialize nqM and nqP (number of bounds)
  nqM=0
  nqP=0

  ! First data point (q=0)
  x1=0.0_qp
  y1P=-zk+sqrt((k**2-x0)**2+1.0_qp)
  y1M=-zk+sqrt((k**2-x0)**2+1.0_qp)

  ! Loop over all data points to intervals with roots
  do iq=1,nq
   
   ! Second data point
   read(11,rec=iq)qT,omqT,mT(:),dmT(:)
   if(qT<qThr)then
    ! Below threshhold, use linear approximation for omq
    call intOmQ(qT,omqT)
   endif
   x2=qT
   y2P=-zk+sqrt((k**2+qT**2+2.0_qp*k*qT-x0)**2+1.0_qp)+omqT
   y2M=-zk+sqrt((k**2+qT**2-2.0_qp*k*qT-x0)**2+1.0_qp)+omqT

   ! Find roots in interval for z=e_{k+q}+omq
   if(y1P*y2P<0.0_qp)then
    if (blaPole)then
     write(6,*)"Pole found for eP in interval",x1,x2
    endif
    nqP=nqP+1
    qX(nqP)=rtsafe(rootFunQEps,(/1.0_qp/),x1,x2,1.e-18_qp)

    ! write(6,*)"pole is ",qX(nqP)
    ! call intOmQ(qX(nqP),omqTEST)
    ! yTEST=-zk+sqrt((k**2+qX(nqP)**2+2.0_qp*k*qX(nqP)-x0)**2+1.0_qp)+omqTEST
    ! write(6,*)"function values:",y1P,yTEST,y2P

    ! call rootFunQEps(x1,(/1.0_qp/),y1TEST,dy1TEST)
    ! call rootFunQEps(x2,(/1.0_qp/),y2TEST,dy2TEST)
    ! call rootFunQEps(qX(nqP),(/1.0_qp/),yTEST,dyTEST)
    ! write(6,*)"function values:",y1TEST,yTEST,y2TEST
    ! write(6,*)"function values derivatives:",dy1TEST,dyTEST,dy2TEST
   endif

   ! Find roots in interval for z=e_{k-q}+omq
   if(y1M*y2M<0.0_qp)then
    if (blaPole)then
     write(6,*)"Pole found for eM in interval",x1,x2
    endif
    nqM=nqM+1
    qM(nqM)=rtsafe(rootFunQEps,(/-1.0_qp/),x1,x2,1.e-18_qp)
   endif

   ! Store second data point in first for next iteration
   x1=x2
   y1P=y2P
   y1M=y2M

  enddo

  ! Close file
  close(11)

  ! Sort all q bounds
  iP=1
  iM=1
  iT=1
  ! Store sorted elements
  do while((iP<=nqP).AND.(iM<=nqM))
   if(qX(iP)<qM(iM))then
    bq(iT)=qX(iP)
    iP=iP+1
    iT=iT+1
   else
    bq(iT)=qM(iM)
    iM=iM+1
    iT=iT+1
   endif
  end do

  ! Store remaining elements
  do while(iP<=nqP)
   bq(iT)=qX(iP)
   iP=iP+1
   iT=iT+1
  end do
  do while(iM<=nqM)
   bq(iT)=qM(iM)
   iM=iM+1
   iT=iT+1
  end do

  ! Calculate q0 if it exists and store
  if((1.0_qp<zk).AND.(zk<1.0_qp+omqMax))then
   if (blaPole)then
    write(6,*)"zk crosses 1+omq, zk,omqMax=",zk,omqMax
   endif
   q0=rtsafe(rootFunQ0,(/0.0_qp/),qmin,qmax,1.e-18_qp)
   bq(iT)=q0
   nbq=nqP+nqM+1
  else
   nbq=nqP+nqM
  endif

 END SUBROUTINE boundsQ
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE rootFunQ0(qIn,arg,x,dx)
  IMPLICIT NONE
  REAL(QP), INTENT(IN) :: qIn
  REAL(QP), DIMENSION(:), INTENT(IN) :: arg
  REAL(QP), INTENT(OUT) :: x,dx

  REAL(QP) om, omM, omP, dom, qInM, qInP, h

  ! Step size for derivative
  h=qIn*0.0001_qp
  ! write(6,*)"Looking for q0 with qIn,h=",qIn,h

  ! Calculate omega values from interpolation
  call intOmQ(qIn  ,om )
  qInM=qIn-h
  qInP=qIn+h
  call intOmQ(qInM,omM)
  call intOmQ(qInP,omP)
  dom=(omP-omM)/(2.0_qp*h)
  ! write(6,*)"Calculated om and dom=",om,dom

  ! Output function and derivative
  x=1.0_qp+om-zk
  dx=dom
  ! write(6,*)"Calculated x and dx=",x,dx

 END SUBROUTINE rootFunQ0
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE rootFunQEps(qIn,arg,x,dx)
  IMPLICIT NONE
  REAL(QP), INTENT(IN) :: qIn
  REAL(QP), DIMENSION(:), INTENT(IN) :: arg
  REAL(QP), INTENT(OUT) :: x,dx

  REAL(QP) om, omM, omP, dom, h
  REAL(QP) us, xi, eps

  ! Plus or minus sign
  us=arg(1)

  ! Step size for derivative
  h=qIn*0.0001_qp

  ! Calculate omega values from interpolation
  call intOmQ(qIn  ,om )
  call intOmQ(qIn-h,omM)
  call intOmQ(qIn+h,omP)
  dom=(omP-omM)/(2.0_qp*h)

  xi=k**2.0_qp+qIn**2.0_qp+2.0_qp*us*k*qIn-x0
  eps=sqrt(xi**2.0_qp+1.0_qp)

  ! Output function and derivative
  x=eps+om-zk
  dx=2*us*k*xi/eps+dom

 END SUBROUTINE rootFunQEps
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END FUNCTION selfEpole 
END MODULE intpole
