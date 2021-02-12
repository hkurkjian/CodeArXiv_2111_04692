MODULE selfcons
  USE modsim
  USE vars
  USE selftot
  IMPLICIT NONE
  LOGICAL blaSC

CONTAINS
  SUBROUTINE SCenergy(mu,k,zk,fichiers,EPS)
    USE recettes, ONLY : rtsafe
    USE recettes, ONLY : mnewt
    IMPLICIT NONE
    REAL(QP), INTENT(IN) :: k
    REAL(QP), INTENT(INOUT) :: zk
    REAL(QP), INTENT(IN) :: EPS(1:3)
    CHARACTER(len=90), INTENT(IN) fichiers(1:6) !1:2 fichldc, 3 fichlec, 4:5 fichom2, 6 fichpol

    REAL(QP) contK, xik, e0, Uk, Vk, OS
    REAL(QP) SigZero(1:6), SigCont(1:6), eZero, eCont
    REAL(QP) tst, tst2

    ! Calculate continuum
    contK=contPole(k) 

    ! Mean-field functions
    xik=k*k-x0
    e0=sqrt(xik*xik+1.0_qp)
    Uk=sqrt((1.0_qp+xik/e0)/2.0_qp)
    Vk=sqrt((1.0_qp-xik/e0)/2.0_qp)
    
    ! Look for a solution below the continuum
    OS=1.e-6_qp
    SigZero=detG(k,      OS,fich(3),EPS)
    SigCont=detG(k,contK-OS,fich(3),EPS)
    eZero=-OS+e0-(Uk*Uk*SigZero(1)+Vk*Vk*SigZero(2)+2.0_qp*Uk*Vk*SigZero(3))
    eCont=-contK+OS+e0-(Uk*Uk*SigCont(1)+Vk*Vk*SigCont(2)+2.0_qp*Uk*Vk*SigCont(3))


    if(eZero*eCont<0.0_qp)then
      if(blaSC)then
        write(6,*)"Looking for solution below continuum"
      endif
      zk=rtsafe(rootFunSC,(/bidon/),OS,contK-OS,1.e-5_qp)
    else
      write(6,*)"Solution is inside the continuum"
    endif

    if(routres)then
     call(unload_data)
    endif

  CONTAINS 
    SUBROUTINE rootFunSC(zkIn,arg,f,df)
      IMPLICIT NONE
      REAL(QP), INTENT(IN) :: zkIn
      REAL(QP), DIMENSION(:), INTENT(IN) :: arg
      REAL(QP), INTENT(OUT) :: f,df

      REAL(QP) sig(1:6), sigP(1:6), h, dsig(1:6)
      ! REAL(QP) sigM(1:6)

      ! Self-energy
      h=zkIn*0.0001_qp
      sig= selfEpole(k,zkIn)
      sigP=selfEpole(k,zkIn+h)
      ! sigM=selfEpole(k,zkIn-h)
      ! dsig=(sigP-sigM)/(2.0_qp*h)
      dsig=(sigP-sig)/h

      f=-zkIn+e0-(Uk*Uk* sig(1)+Vk*Vk* sig(2)+2.0_qp*Uk*Vk* sig(3))
      df=-1     -(Uk*Uk*dsig(1)+Vk*Vk*dsig(2)+2.0_qp*Uk*Vk*dsig(3))

      if(blaSC)then
        write(6,*)"zk=",zkIn
        write(6,*)"f,df=",f,df
      endif

    END SUBROUTINE
  END SUBROUTINE SCenergy
END MODULE selfcons
