PROGRAM example
USE nrtype
USE recettes
USE modpro

 IMPLICIT NONE
 COMPLEX(QPC) I(1:3)
 REAL(QP) Ir(1:3),rhor(1:3)

! function intsuromega has no argument but shares the values 
! of xq,x0 and z with the main program (for a list of
! shared variables see the preamble of modpro in intsuromega.f90).
! It returns a length 3 vector with I(1)=M++, I(2)=M-- and I(3)=M+-
! All matrix elements are in units of Delta or k_Delta and multiplied
! by (2pi)^3.
!
! function dspec takes the value the frequency om as sole argument
! but shares the values of xq,x0 with the main program.
! outputs the spectral density  the spectral density Im(M(om+i0^+))
!
! if axereel=.FALSE. real(intsuromega()) and imag(intsuromega()) correspond to Re(M(z)) and Im(M(z))
! if axereel=.TRUE. real(intsuromega) is Re(M(om+i0^+)) and dspec(real(z)) give the spectral density Im(M(om+i0^+))
!
! bla1 and bla2 control the level of blabla on stdout

 xq=1.0_qp
 x0=1.0_qp

!Example on the real axis
 axereel=.TRUE.
 z=cmplx(3.1_qp,0.0_qp)

 Ir=real(intsuromega())
 rhor=dspec(real(z))
 I=cmplx(Ir,rhor,kind=qpc)

 write(6,*)"axereel=",axereel
 write(6,*)"z=",z
 write(6,*)"I=",I(:)

!Example outside the real axis
 axereel=.FALSE.
 z=cmplx(3.1_qp,1.0_qp)
 I=intsuromega()

 write(6,*)
 write(6,*)"axereel=",axereel
 write(6,*)"z=",z
 write(6,*)"I=",I(:)

END PROGRAM example
