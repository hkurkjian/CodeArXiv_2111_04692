program test
USE modsim
USE recettes
USE dspec
USE vars
USE Zerom
USE intpole
USE eqdetat
IMPLICIT NONE

REAL(QP) om,dom,M(1:2,1:2),dM(1:2,1:2),A(1:6),dMm(1:3),q
COMPLEX(QPC) Mm2(1:6)
REAL(QP) vec(1:2000,1:10000)
REAL(QP) vec2(1:2000,1:10000)
REAL(QP) :: xik,epsk,Uk2,Vk2,k,zk,zkmax,le(1:8),bq(1:10),arr(1:8),don(1:7,1:1000,0:100),don2(1:7,1:1000),est(1:6),dest(1:6)
CHARACTER(len=90) fichdep,fich,fich2
CHARACTER(len=2)  reg,regvieux
INTEGER izk,taille,config(1:7),pos(1:8),nn,nn2,ixq,ixqbis,compteur,nxq,iom
COMPLEX(QPC) Gamm(1:2,1:2),Matt(1:2,1:2),MatCat(1:2,1:2),det
COMPLEX(QPC) SigPole(1:6)

REAL(QP) ptq,ptom,ptM(1:3),ptdM(1:3)
INTEGER iq, ik, iz, ix
REAL(QP) nnn,mmm,intell

LOGICAL errtype1,errtype2,interpol
REAL(QP) ThetaT, XxT, ccheckT

! x0=0.860436686125678599999999999999999961_qp
fichier="./datSpectre/specUni"
! blaPole=.TRUE.

open(1, file = 'unisd.dat', status = 'new')
do ik=1,50
    k=ik/25.0_qp
    do iz=1,50
        zk=1.0_qp+iz/25.0_qp
        SigPole=selfEpole(k,zk)
        write(1,*)k,zk,real(SigPole(1:3)),imag(SigPole(1:3)),real(SigPole(4:6))
    enddo
    write(6,*)"Calculation succesful for k=",k
enddo
close(11)
! k=1.1_qp
! zk=sqrt((k**2.0_qp-x0)**2.0_qp+1.0_qp)


! call system("rm tstUniDAT.dat")
! open(1, file = 'tstUniDAT.dat', status = 'new')
! do ik=1,250
!     k=ik/50.0_qp
!     zk=sqrt((k**2.0_qp-x0)**2.0_qp+1.0_qp)
!     SigPole=selfEpole(k,zk)
!     write(1,*)k,zk,real(SigPole(1:3)),imag(SigPole(1:3)),real(SigPole(4:6))
! enddo
! close(1)


! k=1.14_qp
! zk=sqrt((k**2.0_qp-x0)**2.0_qp+1.0_qp)
! SigPole=selfEpole(k,zk)
! write(6,*)"Sigma=",SigPole(:)

! open(10,file="testCIJ.dat")
! write(10,*)"x0, Theta, Xx, ccheck"
! do ix=1,101
!     x0=ix/10.0_qp-5.1_qp
!     ThetaT =4.0_qp*PI*I5(x0)
!     XxT    =4.0_qp*PI*I6(x0)
!     ccheckT=mc2sDelta(x0)
!     write(10,*)x0,ThetaT,XxT,ccheckT
! enddo
! close(10)

! nn=128
! inquire(file="./datSpectre/"//trim(fichier)//".dat", size=taille)
! nxq=taille/nn
! write(6,*)nxq,taille,nn

! open(10,file="./datSpectre/"//trim(fichier)//".dat",action="read",access="direct",form="unformatted",recl=nn)
! open(11,file="dat21_read.dat")
! ! write(11,*)"q    omq    Mpp    Mmm    Mpm    dMpp    dMmm    dMpm"
! write(11,*)"q    omq"
! do iq=1,nxq
!     read(10,rec=iq)ptq,ptom,ptM(:),ptdM(:)
!     write(11,*)ptq,ptom
! enddo
! close(10)
! close(11)


end program