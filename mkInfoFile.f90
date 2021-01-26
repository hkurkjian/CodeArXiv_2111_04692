program mkInfoFile
    USE eqdetat
    USE nrtype
    USE modsim
    IMPLICIT NONE
    
    REAL(QP) x0,qmin,qmax
    INTEGER nq,nn
    CHARACTER(len=90) fichier
    REAL(QP) y,c0,g0,D0
    REAL(QP) lMpp2,lMpp4,lMmm0,lMmm2,lMmm4,lMpm1,lMpm3
    REAL(QP) ldMpp1,ldMpp3,ldMmm1,ldMmm3,ldMpm0,ldMpm2
    ! REAL(QP) xC,a,b,c,d,p,q,t0,t1,kMM,kMP

    ! fichier="specUni"
    open(10,file="mkInfoFile.inp")
        read(10,*)fichier
    close(10)
    open(11,file=trim(fichier)//".info")
        read(11,*)x0,qmin,qmax,nq,nn
        write(6,*)x0,qmin,qmax,nq,nn
        write(6,*)" "
    close(11)

    y=I6(x0)/I5(x0)
    D0=1.0_qp/(x0*I5(x0)+I6(x0))**(2.0_qp/3.0_qp)

    ! Dispersion parameters
    c0=sqrt(4.0_qp/3.0_qp*(x0+y)/(1.0_qp+y*y))
    g0=(-54.0_qp-13.0_qp*y**2.0_qp+56.0_qp*y**4.0_qp+35.0_qp*y**6.0_qp &
        +2.0_qp*x0*y*(71.0_qp+32.0_qp*y**2.0_qp+y**4.0_qp) &
        -4.0_qp*x0**4.0_qp*(8.0_qp+16.0_qp*y**2.0_qp+13.0_qp*y**4.0_qp) &
        +4.0_qp*x0**3.0_qp*y*(8.0_qp+41.0_qp*y**2.0_qp+13.0_qp*y**4.0_qp) &
        +x0**2.0_qp*(-61.0_qp-252.0_qp*y**2.0_qp-21.0_qp*y**4.0_qp+50.0_qp*y**6.0_qp)) &
        /(135.0_qp*(1.0_qp+x0**2.0_qp)*(1.0_qp+y**2.0_qp)**3.0_qp)

    ! Mpp 
    lMpp2=y**2.0_qp/(4.0_qp*(1.0_qp+y**2.0_qp))
    lMpp4=(-15.0_qp*c0**4.0_qp*(3.0_qp+x0*(4.0_qp*x0+y)) &
        -24.0_qp*(1+x0**2.0_qp)*(6.0_qp+15.0_qp*g0+2.0_qp*x0*(4.0_qp*x0+y)) &
        +10.0_qp*c0**2.0_qp*(7.0_qp*y+x0*(19.0_qp+4.0_qp*x0*(4.0_qp*x0+y)))) &
        /(1920.0_qp*(1+x0**2.0_qp)*(x0+y))
    
    ! Mmm
    lMmm0=3.0_qp/4.0_qp/(x0+y)
    lMmm2=(20+(3*(-2+c0**2.0_qp*x0))/(1+x0**2.0_qp) &
        -(3*(3*c0**2.0_qp+4*x0))/(x0+y))/96
    lMmm4=-(15.0_qp*c0**4.0_qp*( &
            9.0_qp+x0*(-2.0_qp*y+x0*(13.0_qp+2.0_qp*x0*(4.0_qp*x0+y)))) &
        +12.0_qp*( &
            39.0_qp+50.0_qp*g0*(1.0_qp+x0**2.0_qp)*(3.0_qp+2.0_qp*x0**2.0_qp-x0*y) &
            +x0*(18.0_qp*y+x0*(135.0_qp+6.0_qp*x0*y &
            +4.0_qp*x0**2.0_qp*(27.0_qp+2.0_qp*x0*(4.0_qp*x0+y))))) &
        -20.0_qp*c0**2.0_qp*(10.0_qp*y+x0*( &
            13.0_qp+x0*(26.0_qp*y+x0*(41.0_qp+4.0_qp*x0*(4.0_qp*x0+y)))))) &
        /(19200.0_qp*(1.0_qp+x0**2.0_qp)**2.0_qp*(x0+y))
    
    ! Mpm
    lMpm1=(-3.0_qp*c0*y)/(8.0_qp*(x0+y))
    lMpm3=(-12.0_qp*g0*(1.0_qp+x0**2.0_qp)*y &
        -c0**4.0_qp*(x0+y) &
        +c0**2.0_qp*(6.0_qp+4.0_qp*x0**2.0_qp-2.0_qp*x0*y)) &
        /(64.0_qp*c0*(1.0_qp+x0**2.0_qp)*(x0+y))

    ! dMpp
    ldMpp1=-3.0_qp*c0/(8.0_qp*(x0+y))
    ldMpp3=(c0*(7.0_qp-3.0_qp*c0**2.0_qp*x0+4.0_qp*x0**2.0_qp)) &
            /(96.0_qp*(1.0_qp+x0**2.0_qp)) &
        +(-3.0_qp*c0**4.0_qp-6.0_qp*g0+4.0_qp*c0**2.0_qp*x0)/(32.0_qp*c0*(x0+y))

    ! dMmm
    ldMmm1=(c0*(x0/(1.0_qp+x0**2.0_qp)-3.0_qp/(x0+y)))/16.0_qp
    ldMmm3=(-30.0_qp*g0*(1.0_qp+x0**2.0_qp)*(3.0_qp+2.0_qp*x0**2.0_qp-x0*y) &
        -3.0_qp*c0**4.0_qp*( &
            9.0_qp+x0*(-2.0_qp*y+x0*(13.0_qp+2.0_qp*x0*(4.0_qp*x0+y)))) &
        +c0**2.0_qp*(20.0_qp*y+2.0_qp*x0*( &
            13.0_qp+x0*(26.0_qp*y+x0*(41.0_qp+4.0_qp*x0*(4.0_qp*x0+y)))))) &
        /(960.0_qp*c0*(1.0_qp+x0**2.0_qp)**2.0_qp*(x0+y))

    ! dMpm
    ldMpm0=(-3.0_qp*y)/(8.0_qp*(x0+y))
    ldMpm2=(-((3.0_qp*c0**2.0_qp+2.0_qp*x0)/(1.0_qp+x0**2.0_qp)) &
        +6.0_qp/(x0+y))/64.0_qp

    ! Correct units
    lMpp2=lMpp2/D0**(3.0_qp/2.0_qp)*8.0_qp*PI
    lMpp4=lMpp4/D0**(3.0_qp/2.0_qp)*8.0_qp*PI
    lMmm0=lMmm0/D0**(3.0_qp/2.0_qp)*8.0_qp*PI
    lMmm2=lMmm2/D0**(3.0_qp/2.0_qp)*8.0_qp*PI
    lMmm4=lMmm4/D0**(3.0_qp/2.0_qp)*8.0_qp*PI
    lMpm1=lMpm1/D0**(3.0_qp/2.0_qp)*8.0_qp*PI
    lMpm3=lMpm3/D0**(3.0_qp/2.0_qp)*8.0_qp*PI
    ldMpp1=ldMpp1/D0**(3.0_qp/2.0_qp)*8.0_qp*PI
    ldMpp3=ldMpp3/D0**(3.0_qp/2.0_qp)*8.0_qp*PI
    ldMmm1=ldMmm1/D0**(3.0_qp/2.0_qp)*8.0_qp*PI
    ldMmm3=ldMmm3/D0**(3.0_qp/2.0_qp)*8.0_qp*PI
    ldMpm0=ldMpm0/D0**(3.0_qp/2.0_qp)*8.0_qp*PI
    ldMpm2=ldMpm2/D0**(3.0_qp/2.0_qp)*8.0_qp*PI

    ! ! Interval around minimum
    ! xC=2.449775705506458_qp

    ! ! Calculate polynomial coefficients a k^6 + b k^4 + c k^2 + d
    ! a=4.0_qp
    ! b=-(c0**2.0_qp+8.0_qp*x0)
    ! c=2.0_qp*x0*(c0**2.0_qp+2.0_qp*x0)
    ! d=-c0**2.0_qp*(x0**2.0_qp+1.0_qp)

    ! ! Depressed cubic t^3 + p t + q
    ! p=(3.0_qp*a*c-b**2.0_qp)/(3.0_qp*a**2.0_qp)
    ! q=(2.0_qp*b**3.0_qp-9.0_qp*a*b*c+27.0_qp*a**2.0_qp*d)/(27.0_qp*a**3.0_qp)

    ! if(x0<xC)then ! kMM=0
    !     kMM=0.0_qp
    !     t0=(-q/2.0_qp+sqrt(q**2.0_qp/4.0_qp+p**3.0_qp/27.0_qp))**(1.0_qp/3.0_qp) &
    !      +(-q/2.0_qp-sqrt(q**2.0_qp/4.0_qp+p**3.0_qp/27.0_qp))**(1.0_qp/3.0_qp)
    !     kMP=sqrt(t0-b/(3.0_qp*a))
    ! else
    !     t0=2.0_qp*sqrt(-p/3.0_qp)*cos(acos(3.0_qp/2.0_qp*q/p*sqrt(-3.0_qp/p))/3.0_qp)
    !     t1=2.0_qp*sqrt(-p/3.0_qp)*cos(acos(3.0_qp/2.0_qp*q/p*sqrt(-3.0_qp/p))/3.0_qp &
    !                                     -2.0_qp*PI/3.0_qp)
    !     kMM=sqrt(t1-b/(3.0_qp*a));
    !     kMP=sqrt(t0-b/(3.0_qp*a));
    ! endif

    write(6,*)y,c0,g0,D0
    write(6,*)" "
    write(6,*)lMpp2,lMpp4,lMmm0,lMmm2,lMmm4,lMpm1,lMpm3
    write(6,*)" "
    write(6,*)ldMpp1,ldMpp3,ldMmm1,ldMmm3,ldMpm0,ldMpm2
    ! write(6,*)" "
    ! write(6,*)kMM,kMP

    open(11,file=trim(fichier)//".info")
        write(11,*)x0,qmin,qmax,nq,nn
        write(11,*)c0,g0
        write(11,*)lMpp2,lMpp4,lMmm0,lMmm2,lMmm4,lMpm1,lMpm3
        write(11,*)ldMpp1,ldMpp3,ldMmm1,ldMmm3,ldMpm0,ldMpm2
        ! write(11,*)kMM,kMP
    close(11)
    
end program