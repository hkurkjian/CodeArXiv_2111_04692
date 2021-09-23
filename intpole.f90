MODULE intpole
USE nrtype
USE modsim
USE vars
USE angularint
IMPLICIT NONE
REAL(QP) :: c0,g0
REAL(QP) :: lMpp2,lMpp4,lMmm0,lMmm2,lMmm4,lMpm1,lMpm3
REAL(QP) :: ldMpp1,ldMpp3,ldMmm1,ldMmm3,ldMpm0,ldMpm2
REAL(QP) :: kMM,kMP
REAL(QP) :: qThr,qmin,qmax,dq
REAL(QP), ALLOCATABLE, DIMENSION(:,:) :: donpol
INTEGER nq,nn,nqeff
LOGICAL blaPole
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION selfEpole(k,zk,EPSpole)
USE modsim
REAL(QP), INTENT(IN) :: k,zk,EPSpole
COMPLEX(QPC) :: selfEpole(1:6), SEint(1:6)

REAL(QP) :: q,kmin,kmax,dk,bq(1:7)
INTEGER ix,nx,iq,nbq
INTEGER nbounds, ibound
REAL(QP) bounds(1:9)
 
! Read info-file for low-q and interpolation variables
! Contains: 
!   c0,g0
!   lMpp2,lMpp4,lMmm0,lMmm2,lMmm4,lMpm1,lMpm3
!   ldMpp1,ldMpp3,ldMmm1,ldMmm3,ldMpm0,ldMpm2
!   kMM, kMP
!   qThr,qmin,qmax,dq
!   nqeff

if (blaPole)then
  write(6,*)"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  write(6,*)
  write(6,*)"++++++++++++++++++++++++++++ selfEpole ++++++++++++++++++++++++++++"
  write(6,*)
  write(6,*)"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
endif

! Open file for interpolation of omega_q
!open(32,file=trim(fichpol)//".dat",action="read",access="direct",form="unformatted",recl=nn)

! Calculate q bounds

if (blaPole)then
  write(6,*)
  write(6,*)"dq,nqeff,qmax=",dq,nqeff,qmin+nqeff*dq
  write(6,*)"Calculating all q-bounds for k,zk=",k,zk
  write(6,*)
endif
call boundsQ(k,zk,bq,nbq)
if (blaPole)then
  write(6,*)
  write(6,*)"            Coefficients of the low q expansion"
  write(6,*)
  write(6,FMT="(A6 ,2G20.10)")"c0,g0=",c0,g0
  write(6,FMT="(A42,7G20.10)")"lMpp2,lMpp4,lMmm0,lMmm2,lMmm4,lMpm1,lMpm3=",lMpp2,lMpp4,lMmm0,lMmm2,lMmm4,lMpm1,lMpm3
  write(6,FMT="(A42,6G20.10)")"ldMpp1,ldMpp3,ldMmm1,ldMmm3,ldMpm0,ldMpm2=",ldMpp1,ldMpp3,ldMmm1,ldMmm3,ldMpm0,ldMpm2
  write(6,FMT="(A8 ,2G20.10)")"kMM,kMP=",kMM,kMP
  write(6,*)
  write(6,*)"---------------------------------------------"
endif
if (blaPole)then
  write(6,*)
  write(6,*)nbq," q bounds found:",bq(1:nbq)
endif

! Add 0.0 and qMax to bounds
bounds(1)=0.0_qp
if(nbq>0)then
 bounds(2:1+nbq)=bq(1:nbq)
 bounds(2+nbq)=qmin+nqeff*dq
 nbounds=2+nbq
else
 bounds(2)=qmin+nqeff*dq
 nbounds=2
endif

if (blaPole)then
 write(6,*)
 write(6,*)"All integration boundaries:",bounds(1:nbounds)
 write(6,*)"---------------------------------------------"
endif

! Calculate q-integral
selfEpole(:)=cmplx(0.0_qp,0.0_qp,kind=qpc)
do ibound=1,nbounds-1
 if (blaPole)then
  write(6,*)"Integration from ",bounds(ibound)," to ",bounds(ibound+1)
 endif
 SEint=qromovcq(integrandeq  ,bounds(ibound),bounds(ibound+1),6,(/bidon/),midpntvcq,EPSpole)
 selfEpole(:)=selfEpole(:)+SEint(:)
 if (blaPole)then
  write(6,*)"re selfEpole (1<->2)= ",real(SEint(1:3))
  write(6,*)"im selfEpole (1<->2)= ",imag(SEint(1:3))
  write(6,*)"re selfEpole (3<->0)= ",real(SEint(4:6))
  write(6,*)"im selfEpole (3<->0)= ",imag(SEint(4:6))
  write(6,*)
  write(6,*)"---------------------------------------------"
 endif
enddo

! Calculate phi-integral
selfEpole=selfEpole*2.0_qp*PI

if (blaPole)then
  write(6,*)"Total"
  write(6,*)"re selfEpole (1<->2)= ",real(selfEpole(1:3))
  write(6,*)"im selfEpole (1<->2)= ",imag(selfEpole(1:3))
  write(6,*)"re selfEpole (3<->0)= ",real(selfEpole(4:6))
  write(6,*)"im selfEpole (3<->0)= ",imag(selfEpole(4:6))
  write(6,*)
  write(6,*)"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
endif

!close(32)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
 FUNCTION integrandeq(q,arg,m) 

 INTEGER,  INTENT(IN) :: m
 REAL(QP), INTENT(IN), DIMENSION(:) :: q,arg
 COMPLEX(QPC), DIMENSION(size(q),m) :: integrandeq

 REAL(QP)  ptq(1:5),ptom(1:5),ptM(1:3,1:5),ptdM(1:3,1:5),sdonpol(1:8,1:5)
 REAL(QP)  qs,om,Ma(1:3),dM(1:3),MatCat(1:3),errom,errM(1:3),errdM(1:3)
 REAL(QP) xiP,xiM,epsP,epsM,xmin,xmax
 REAL(QP)  ddet
 COMPLEX(QPC)  IuM(1:3),IuP(1:3)
 INTEGER is,ib,iq

 do is=1,size(q)
  qs=q(is)
  if(qs<qThr)then
  ! For small q, use analytic formulas
   om=c0*qs*(1.0_qp+g0*(qs/c0)**2.0_qp)
   Ma(1)=lMpp2*qs**2.0_qp+lMpp4*qs**4.0_qp
   Ma(2)=lMmm0+lMmm2*qs**2.0_qp+lMmm4*qs**4.0_qp
   Ma(3)=lMpm1*qs+lMpm3*qs**3.0_qp
   dM(1)=ldMpp1*qs+ldMpp3*qs**3.0_qp
   dM(2)=ldMmm1*qs+ldMmm3*qs**3.0_qp
   dM(3)=ldMpm0+ldMpm2*qs**2.0_qp
  else
   
  ! Select 5 data points around q
   iq=floor((qs-qmin)/dq)+1

   sdonpol=donpol(:,iq-2:iq+2)
   ptq =sdonpol(1  ,:)
   ptom=sdonpol(2  ,:)
   ptM =sdonpol(3:5,:)
   ptdM=sdonpol(6:8,:)
!   read(32,rec=iq-2)ptq(1),ptom(1),ptM(:,1),ptdM(:,1)
!   read(32,rec=iq-1)ptq(2),ptom(2),ptM(:,2),ptdM(:,2)
!   read(32,rec=iq  )ptq(3),ptom(3),ptM(:,3),ptdM(:,3)
!   read(32,rec=iq+1)ptq(4),ptom(4),ptM(:,4),ptdM(:,4)
!   read(32,rec=iq+2)ptq(5),ptom(5),ptM(:,5),ptdM(:,5)
  
!   write(6,*)qs,ptq(3:4)
  ! Interpolation scheme for omega_q
   call polint(ptq,ptom,qs,om,errom)
   ! Interpolation scheme for each matrix component and derivative
   do ib=1,3
    call polint(ptq, ptM(ib,:),qs, Ma(ib),errM(ib))
    call polint(ptq,ptdM(ib,:),qs,dM(ib),errdM(ib))
   enddo
  endif

  ! Derivative of the determinant
  ddet=Ma(1)*dM(2)+Ma(2)*dM(1)-2.0_qp*Ma(3)*dM(3)

  ! Cartesian basis
  MatCat(1)=(Ma(1)+Ma(2))/2.0_qp+Ma(3)
  MatCat(2)=(Ma(1)+Ma(2))/2.0_qp-Ma(3)
  MatCat(3)=(Ma(2)-Ma(1))/2.0_qp

  ! Energy functions for analytic angle integration
  xiP=k**2+qs**2+2.0_qp*k*qs-x0
  xiM=k**2+qs**2-2.0_qp*k*qs-x0
  epsP=sqrt(xiP**2+1.0_qp)
  epsM=sqrt(xiM**2+1.0_qp)
  xmin=epsP-xiP
  xmax=epsM-xiM

  ! Angle integration
  IuP=Iuanaly(om+zk,k,qs,xmin,xmax)
  IuM=conjg(Iuanaly(om-zk,k,qs,xmin,xmax))

  ! 1->2 process integrand
  integrandeq(is,1)=-MatCat(2)*IuM(1)/ddet
  integrandeq(is,2)=-MatCat(1)*IuM(2)/ddet
  integrandeq(is,3)= MatCat(3)*IuM(3)/ddet

  ! 0->3 process integrand
  integrandeq(is,4)= MatCat(1)*IuP(2)/ddet
  integrandeq(is,5)= MatCat(2)*IuP(1)/ddet
  integrandeq(is,6)= MatCat(3)*IuP(3)/ddet

  ! multiply by Jacobian q^2
  ! if(blaPole)then
  !  write(6,*)"qs,integrandeq=",qs,integrandeq(is,1)
  ! endif
  integrandeq(is,:)=qs**2*integrandeq(is,:)
 enddo

 END FUNCTION integrandeq 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE boundsQ(k,zk,bq,nbq)
 ! Calculate all bounds on q for integration above emission continuum
  USE recettes, ONLY : rtsafe
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

  ! Read last value for omqMax
  qMax  =donpol(1,nqeff)
  omqMax=donpol(2,nqeff)
  
!  read(32,rec=nqeff)qMax,omqMax,mMax(:),dmMax(:)

  ! Initialize nqM and nqP (number of bounds)
  nqM=0
  nqP=0

  ! First data point (q=0)
  x1=0.0_qp
  y1P=-zk+sqrt((k**2-x0)**2+1.0_qp)
  y1M=-zk+sqrt((k**2-x0)**2+1.0_qp)

  ! Loop over all data points to intervals with roots
  do iq=1,nqeff
   
   ! Second data point
   qT  =donpol(1,iq)
   omqT=donpol(2,iq)
!   read(32,rec=iq)qT,omqT,mT(:),dmT(:)
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
  close(31)

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
   if (blaPole) write(6,*)"zk crosses 1+omq, zk,omqMax=",zk,omqMax
   q0=rtsafe(rootFunQ0,(/0.0_qp/),qmin,nqeff*dq,1.e-18_qp)
   bq(iT)=q0
   nbq=nqP+nqM+1
  else
   nbq=nqP+nqM
  endif

 END SUBROUTINE boundsQ
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE rootFunQ0(qIn,arg,x,dx)
  REAL(QP), INTENT(IN) :: qIn
  REAL(QP), DIMENSION(:), INTENT(IN) :: arg
  REAL(QP), INTENT(OUT) :: x,dx

  REAL(QP) om, omM, omP, dom, qInM, qInP, h

  ! Step size for derivative
  h=qIn*0.00001_qp

  ! Calculate omega values from interpolation
  call intOmQ(qIn  ,om )
  qInM=qIn-h
  qInP=qIn+h
  call intOmQ(qInM,omM)
  call intOmQ(qInP,omP)
  dom=(omP-omM)/(2.0_qp*h)

  ! Output function and derivative
  x=1.0_qp+om-zk
  dx=dom

 END SUBROUTINE rootFunQ0
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE rootFunQEps(qIn,arg,x,dx)
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
  dx=2.0_qp*(q+us*k)*xi/eps+dom

 END SUBROUTINE rootFunQEps
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END FUNCTION selfEpole 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE intOmQ(qVal,omq)
 REAL(QP), INTENT(IN) :: qVal
 REAL(QP), INTENT(OUT) :: omq

 INTEGER iq
 REAL(QP)  ptq(1:5),ptom(1:5),ptM(1:3,1:5),ptdM(1:3,1:5), sdonpol(1:8,1:5),errom

 if(qVal<qThr)then
   ! For small q, use analytic formulas
   omq=c0*qVal*(1.0_qp+g0*(qVal/c0)**2.0_qp)
 else

   ! Read 5 points around qVal
   iq=floor((qVal-qmin)/dq)+1

   if(iq.GE.nqeff)then 
     write(6,*) "Attention qVal=",qVal," dépasse la valeur maximale de q=",donpol(1,nqeff)," du fichier"
     write(6,*) "valeur de omq obtenue par extrapolation polynomiale"
     sdonpol=donpol(:,nqeff-4:nqeff)
   elseif(iq.GE.nqeff-1)then
    sdonpol=donpol(:,iq-3:iq+1)
   else
    sdonpol=donpol(:,iq-2:iq+2)
   endif
   ptq =sdonpol(1  ,:)
   ptom=sdonpol(2  ,:)
   ptM =sdonpol(3:5,:)
   ptdM=sdonpol(6:8,:)
!   read(32,rec=iq-2)ptq(1),ptom(1),ptM(:,1),ptdM(:,1)
!   read(32,rec=iq-1)ptq(2),ptom(2),ptM(:,2),ptdM(:,2)
!   read(32,rec=iq  )ptq(3),ptom(3),ptM(:,3),ptdM(:,3)
!   read(32,rec=iq+1)ptq(4),ptom(4),ptM(:,4),ptdM(:,4)
!   read(32,rec=iq+2)ptq(5),ptom(5),ptM(:,5),ptdM(:,5)

   ! Interpolation scheme for omega_q
   call polint(ptq,ptom,qVal,omq,errom)
 endif

END SUBROUTINE intOmQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE rdInfo(fichpol)
  ! Read the fichpol.info file for all needed variables for interpolation
  CHARACTER(len=90), INTENT(IN) :: fichpol
  INTEGER taille,iq

  open(31,file=trim(fichpol)//".info")
    read(31,*)x0,qmin,qmax,nq,nn
    read(31,*)c0,g0
    read(31,*)lMpp2,lMpp4,lMmm0,lMmm2,lMmm4,lMpm1,lMpm3
    read(31,*)ldMpp1,ldMpp3,ldMmm1,ldMmm3,ldMpm0,ldMpm2
    read(31,*)kMM,kMP
  close(31)

  ! q step
  dq=(qmax-qmin)/nq

  ! Effective q length
  inquire(file=trim(fichpol)//".dat", size=taille)
  nqeff=taille/nn

  ! Select threshhold q below which omq is approximated by c q
  qThr=0.0405_qp
  
  !load fichpol.dat into donpol
  allocate(donpol(1:8,1:nqeff))
  open(32,file=trim(fichpol)//".dat",action="read",access="direct",form="unformatted",recl=nn)
    do iq=1,nqeff
      read(32,rec=iq)donpol(1:8,iq)!q,om,M(1:3),dM(1:3)
!      if(blaPole) write(6,*)"iq,q=",iq,donpol(1,iq)
    enddo
  close(32)

  nqeff=nqeff-2
END SUBROUTINE rdInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION contPole(k) 
  ! Calculate the lower continuum edge for the 1->2 process
  USE recettes, ONLY : rtsafe
  REAL(QP), INTENT(IN) :: k
  REAL(QP) :: contPole(1:3)

  REAL(QP) qMax, omqMax, mMax(1:3), dmMax(1:3)
  REAL(QP) ptq(1:5),ptom(1:5),ptM(1:3,1:5),ptdM(1:3,1:5), sdonpol(1:8,1:5), errom
  REAL(QP) qT, omqT, domqT, qRt, omqRt
  REAL(QP) xikq, epkq, h, domqM, domqP
  REAL(QP) k0, uC, contTMP(1:9), qTMP(1:9), omTMP(1:9)
  REAL(QP) x1, x2, y1, y2
  INTEGER nCM, iq, mini

  ! Initialize
  contTMP(:)=1.0e50_qp
  nCM=0
  ! Minimum
  if(x0<=0.0)then
    k0=0.0_qp
  else
    k0=sqrt(x0)
  endif
  
  ! Open file for interpolation of omega_q
!  open(32,file=trim(fichpol)//".dat",action="read",access="direct",form="unformatted",recl=nn)

  ! Read last value for omqMax
  qMax  =donpol(1,nqeff)
  omqMax=donpol(2,nqeff)
!  read(32,rec=nqeff)qMax,omqMax,mMax(:),dmMax(:)
  
  if((kMM<=k).AND.(k<=kMP))then
    ! Close to minimum: contPole = eps_k
    contPole(1)=sqrt((k*k-x0)**2.0_qp+1.0_qp)
    contPole(2)=0.0_qp
    contPole(3)=0.0_qp
  else if((k0+qMax < k).AND.(qMax<20.0))then
    ! At high k in BCS regime, contPole = eps_{k-qMax)+om_{qMax}
    contPole(1)=sqrt((k*k+qMax*qMax-2.0_qp*k*qMax-x0)**2.0_qp+1.0_qp)+omqMax
    contPole(2)=qMax
    contPole(3)=omqMax
  else if((k>kMP).AND.((k-kMP)<0.1_qp))then
    ! Avoid the bug at low-q by computing analytically the edge when k is close to kPM or kMM 
          contPole(2)=k-kMP
          contPole(3)=c0*contPole(2)
          contPole(1)=sqrt(((k-contPole(2))**2-x0)**2+1.0_qp)+contPole(3)
  else if((k<kMM).AND.((kMM-k)<0.1_qp))then
          contPole(2)=kMM-k
          contPole(3)=c0*contPole(2)
          contPole(1)=sqrt(((k+contPole(2))**2-x0)**2+1.0_qp)+contPole(3)
  else
    if(k<kMM)then
      ! Descending part: minimum of e_{k+q}+omq
      uC=1.0_qp;
    else
      ! Ascending part: minimum of e_{k-q}+omq
      uC=-1.0_qp;
    endif

    ! First data point (q=0)
    x1=0.0_qp
    xikq=k**2.0_qp-x0
    epkq=sqrt(xikq*xikq+1.0_qp)
    y1=2.0_qp*uC*k*xikq/epkq+c0

    ! Loop over all data points to intervals with roots
    !   (start from iq=3 for interpolation scheme)
    do iq=30,nqeff-2
! BUG ??

      ! Calculate derivative at q point
      sdonpol=donpol(:,iq-2:iq+2)
      ptq =sdonpol(1  ,:)
      ptom=sdonpol(2  ,:)
      ptM =sdonpol(3:5,:)
      ptdM=sdonpol(6:8,:)
!      read(32,rec=iq-2)ptq(1),ptom(1),ptM(:,1),ptdM(:,1)
!      read(32,rec=iq-1)ptq(2),ptom(2),ptM(:,2),ptdM(:,2)
!      read(32,rec=iq  )ptq(3),ptom(3),ptM(:,3),ptdM(:,3)
!      read(32,rec=iq+1)ptq(4),ptom(4),ptM(:,4),ptdM(:,4)
!      read(32,rec=iq+2)ptq(5),ptom(5),ptM(:,5),ptdM(:,5)
      qT=ptq(3)         ! q value
      omqT=ptom(3)      ! omq at qT

      ! Noticed that the brute-force technique breaks down at low q, so skip these values
      if(qT>0.02)then

        ! Derivative of omq
        h=qT*0.0001_qp   ! step size for derivative
        call polint(ptq,ptom,qT-h,domqM,errom)
        call polint(ptq,ptom,qT+h,domqP,errom)
        domqT=(domqP-domqM)/(2.0_qp*h) ! Derivative at qT

        ! Second data point
        x2=qT
        xikq=k*k+qT*qT+2.0_qp*k*qT*uC-x0;
        epkq=sqrt(xikq*xikq+1.0_qp);
        y2=2.0_qp*(qT+uC*k)*xikq/epkq+domqT

        ! Look for root in interval if sign change
        if(y1*y2<0.0_qp)then

          ! q value at extremum
          qRt=rtsafe(rootFunContQ,(/uC/),x1,x2,1.e-18_qp)

          ! Calculate function value
          call intOmQ(qRt,omqRt)
          nCM=nCM+1
          contTMP(nCM)=sqrt((k*k+qRt*qRt+2.0_qp*k*qRt*uC-x0)**2.0_qp+1.0_qp)+omqRt
          qTMP(nCM)=qRt
          omTMP=omqRt

          if (blaPole)then
            write(6,*)"Root found in interval",x1,x2,", for k=",k
            write(6,*)"function value=",contTMP(nCM)
          endif
        endif

        ! Store second data point in first for next iteration
        x1=x2
        y1=y2
      endif

    enddo

    ! From all extrema, take the minimum
    mini=iminloc(contTMP)
    contPole(1)=contTMP(mini)
    contPole(2)=qTMP   (mini)
    contPole(3)=omTMP  (mini)
    
  endif
  
  ! Close file
!  close(32)
  
CONTAINS
  SUBROUTINE rootFunContQ(qIn,arg,x,dx)
    REAL(QP), INTENT(IN) :: qIn
    REAL(QP), DIMENSION(:), INTENT(IN) :: arg
    REAL(QP), INTENT(OUT) :: x,dx

    REAL(QP) om, omM, omP, dom, ddom, h
    REAL(QP) us, xi, eps

    ! Plus or minus sign
    us=arg(1)

    ! Step size for derivative
    h=qIn*0.0001_qp

    ! Calculate omega values from interpolation
    call intOmQ(qIn  ,om )
    call intOmQ(qIn-h,omM)
    call intOmQ(qIn+h,omP)
    dom =(omP-omM)/(2.0_qp*h)
    ddom=(omP-2.0_qp*om+omM)/h/h

    xi=k*k+qIn*qIn+2.0_qp*us*k*qIn-x0
    eps=sqrt(xi*xi+1.0_qp)

    ! Output function and derivative
    x=2.0_qp*(qIn+us*k)*xi/eps+dom
    dx=4.0_qp*(qIn+us*k)**2.0_qp*(eps*eps-xi*xi)/eps/eps/eps+2.0_qp*xi/eps+ddom

  END SUBROUTINE rootFunContQ
END FUNCTION contPole
END MODULE intpole
