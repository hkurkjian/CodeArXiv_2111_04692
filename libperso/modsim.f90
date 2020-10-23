MODULE modsim
 USE nrtype; USE nrutil, ONLY :arth
       INTERFACE iminloc
        MODULE PROCEDURE iminloc_s,iminloc_d,iminloc_q
       END INTERFACE iminloc
       INTERFACE polint 
        MODULE PROCEDURE polints,polintd,polintq,polintc,polintcq
       END INTERFACE polint
       INTERFACE trapz 
        MODULE PROCEDURE trapzs,trapzd,trapzq,trapzc,trapzcq
       END INTERFACE trapz
       INTERFACE midpnt
        MODULE PROCEDURE midpnts,midpntd,midpntq,midpntc,midpntcq
       END INTERFACE midpnt
       INTERFACE midinf
        MODULE PROCEDURE midinfs,midinfd,midinfq,midinfc,midinfcq
       END INTERFACE midinf
       INTERFACE midsqu
        MODULE PROCEDURE midsqus,midsqud,midsquq,midsquc,midsqucq
       END INTERFACE midsqu
       INTERFACE midsql
        MODULE PROCEDURE midsqls,midsqld,midsqlq,midsqlc,midsqlcq
       END INTERFACE midsql
       INTERFACE racinf
        MODULE PROCEDURE racinfq,racinfcq
       END INTERFACE racinf
       INTERFACE qromb
        MODULE PROCEDURE qrombs,qrombd,qrombq,qrombc,qrombcq
       END INTERFACE qromb
       INTERFACE qromo
        MODULE PROCEDURE qromos,qromod,qromoq,qromoc,qromocq
       END INTERFACE qromo
!       INTERFACE func
!        FUNCTION funcs(x,arg)
!        USE nrtype
!        REAL(SP), DIMENSION(:), INTENT(IN) :: x
!        REAL(SP), DIMENSION(:), INTENT(IN) :: arg
!        REAL(SP), DIMENSION(size(x)) :: funcs
!        END FUNCTION funcs
!        FUNCTION funcd(x,arg)
!        USE nrtype
!        REAL(DP), DIMENSION(:), INTENT(IN) :: x
!        REAL(DP), DIMENSION(:), INTENT(IN) :: arg
!        REAL(DP), DIMENSION(size(x)) :: funcd
!        END FUNCTION funcd
!       END INTERFACE
 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE trapzs(func,a,b,ar,s,n)
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), INTENT(IN), DIMENSION(:) :: ar
       REAL(SP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(SP), DIMENSION(:), INTENT(IN) :: x
         REAL(SP), DIMENSION(:), INTENT(IN) :: arg
         REAL(SP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended trapezoidal rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of b a f(x)dx. Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2n-2
!       additional interior points. s should not be modified between sequential calls.
       REAL(SP) :: del,fsum
       INTEGER(I4B) :: it
       if (n == 1) then
       s=0.5_sp*(b-a)*sum(func( (/ a,b /) ,ar))
       else
       it=2**(n-2)
       del=(b-a)/it !This is the spacing of the points to be added.
       fsum=sum(func(arth(a+0.5_sp*del,del,it),ar))
       s=0.5_sp*(s+del*fsum) !This replaces s by its refined value.
       end if
       END SUBROUTINE trapzs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE trapzd(func,a,b,ar,s,n)
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), INTENT(IN), DIMENSION(:) :: ar
       REAL(DP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(DP), DIMENSION(:), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         REAL(DP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended trapezoidal rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of b a f(x)dx. Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2n-2
!       additional interior points. s should not be modified between sequential calls.
       REAL(DP) :: del,fsum
       INTEGER(I4B) :: it
       if (n == 1) then
       s=0.5_dp*(b-a)*sum(func( (/ a,b /) ,ar))
       else
       it=2**(n-2)
       del=(b-a)/it !This is the spacing of the points to be added.
       fsum=sum(func(arth(a+0.5_dp*del,del,it),ar))
       s=0.5_dp*(s+del*fsum) !This replaces s by its refined value.
       end if
       END SUBROUTINE trapzd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE trapzq(func,a,b,ar,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: ar
       REAL(QP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended trapezoidal rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of b a f(x)dx. Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2n-2
!       additional interior points. s should not be modified between sequential calls.
       REAL(QP) :: del,fsum
       INTEGER(I4B) :: it
       if (n == 1) then
       s=0.5_qp*(b-a)*sum(func( (/ a,b /) ,ar))
       else
       it=2**(n-2)
       del=(b-a)/it !This is the spacing of the points to be added.
       fsum=sum(func(arth(a+0.5_qp*del,del,it),ar))
       s=0.5_qp*(s+del*fsum) !This replaces s by its refined value.
       end if
       END SUBROUTINE trapzq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE trapzc(func,a,b,ar,s,n)
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), INTENT(IN), DIMENSION(:) :: ar
       COMPLEX(DPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(DP), DIMENSION(:), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(DPC), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended trapezoidal rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of b a f(x)dx. Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2n-2
!       additional interior points. s should not be modified between sequential calls.
       REAL(DP) :: del
       COMPLEX(DPC) :: fsum
       INTEGER(I4B) :: it
       if (n == 1) then
       s=0.5_dp*(b-a)*sum(func( (/ a,b /) ,ar))
       else
       it=2**(n-2)
       del=(b-a)/it !This is the spacing of the points to be added.
       fsum=sum(func(arth(a+0.5_dp*del,del,it),ar))
       s=0.5_dp*(s+del*fsum) !This replaces s by its refined value.
       end if
       END SUBROUTINE trapzc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE trapzcq(func,a,b,ar,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: ar
       COMPLEX(QPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended trapezoidal rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of b a f(x)dx. Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2n-2
!       additional interior points. s should not be modified between sequential calls.
       REAL(QP) :: del
       COMPLEX(QPC) :: fsum
       INTEGER(I4B) :: it
       if (n == 1) then
       s=0.5_qp*(b-a)*sum(func( (/ a,b /) ,ar))
       else
       it=2**(n-2)
       del=(b-a)/it !This is the spacing of the points to be added.
       fsum=sum(func(arth(a+0.5_qp*del,del,it),ar))
       s=0.5_qp*(s+del*fsum) !This replaces s by its refined value.
       end if
       END SUBROUTINE trapzcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midpnts(func,a,b,arg,s,n)
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), INTENT(IN), DIMENSION(:) :: arg
       REAL(SP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(SP), DIMENSION(:), INTENT(IN) :: x
         REAL(SP), DIMENSION(:), INTENT(IN) :: arg
         REAL(SP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended midpoint rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding
!       (2/3)×3n-1 additional interior points. s should not be modified between sequential calls.
       REAL(SP) :: del
       INTEGER(I4B) :: it
       REAL(SP), DIMENSION(2*3**(n-2)) :: x
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_sp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_sp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it) 
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
        s=s/3.0_sp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       END SUBROUTINE midpnts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midpntd(func,a,b,arg,s,n)
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), INTENT(IN), DIMENSION(:) :: arg
       REAL(DP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(DP), DIMENSION(:), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         REAL(DP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended midpoint rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding
!       (2/3)×3n-1 additional interior points. s should not be modified between sequential calls.
       REAL(DP) :: del
       INTEGER(I4B) :: it
       REAL(DP), DIMENSION(2*3**(n-2)) :: x
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_dp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_dp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_dp*del,3.0_dp*del,it) 
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_dp*del
        s=s/3.0_dp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       END SUBROUTINE midpntd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midpntq(func,a,b,arg,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       REAL(QP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended midpoint rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding
!       (2/3)×3n-1 additional interior points. s should not be modified between sequential calls.
       REAL(QP) :: del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it) 
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s=s/3.0_qp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       END SUBROUTINE midpntq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midpntvq(func,a,b,m,arg,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       REAL(QP), INTENT(INOUT), DIMENSION(m) :: s
       INTEGER(I4B), INTENT(IN) :: m,n
       INTERFACE
        FUNCTION func(x,arg,mm)
         USE nrtype
         INTEGER, INTENT(IN) ::  mm
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), DIMENSION(size(x),mm) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended midpoint rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding
!       (2/3)×3n-1 additional interior points. s should not be modified between sequential calls.
       REAL(QP) :: del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       if (n == 1) then
        s(:)=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg,m ),dim=1)
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it) 
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s(:)=s(:)/3.0_qp+del*sum(func(x,arg,m),dim=1) !The new sum is combined with the old integral
       end if !to give a refined integral.
       END SUBROUTINE midpntvq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midpntc(func,a,b,arg,s,n)
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(DPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(DP), DIMENSION(:), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(DPC), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended midpoint rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding
!       (2/3)×3n-1 additional interior points. s should not be modified between sequential calls.
       REAL(DP) :: del
       INTEGER(I4B) :: it
       REAL(DP), DIMENSION(2*3**(n-2)) :: x
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_dp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_dp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_dp*del,3.0_dp*del,it) 
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_dp*del
        s=s/3.0_dp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       END SUBROUTINE midpntc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midpntcq(func,a,b,arg,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(QPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended midpoint rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding
!       (2/3)×3n-1 additional interior points. s should not be modified between sequential calls.
       REAL(QP) :: del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it) 
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s=s/3.0_qp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       END SUBROUTINE midpntcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midinfs(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: aa,bb
       REAL(SP), INTENT(IN), DIMENSION(:) :: arg
       REAL(SP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(SP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(SP), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(SP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(SP), DIMENSION(2*3**(n-2)) :: x
       if(aa*bb <= 0.0) STOP 'bornes dans midinf'
       b=1.0_sp/aa !These two statements change the limits of integration accordingly
       a= 1.0_sp/bb
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_sp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_sp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
        s=s/3.0_sp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(SP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(SP), DIMENSION(size(x)) :: func
        func=funk(1.0_sp/x,arg)/x**2
        END FUNCTION func
       END SUBROUTINE midinfs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midinfd(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: aa,bb
       REAL(DP), INTENT(IN), DIMENSION(:) :: arg
       REAL(DP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(DP), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(DP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(DP), DIMENSION(2*3**(n-2)) :: x
       if(aa*bb <= 0.0) STOP 'bornes dans midinf'
       b=1.0_dp/aa !These two statements change the limits of integration accordingly
       a= 1.0_dp/bb
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_dp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_dp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_dp*del,3.0_dp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_dp*del
        s=s/3.0_dp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(DP), DIMENSION(size(x)) :: func
        func=funk(1.0_dp/x,arg)/x**2
        END FUNCTION func
       END SUBROUTINE midinfd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midinfq(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: aa,bb
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       REAL(QP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for miqpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(QP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       if(aa*bb <= 0.0) STOP 'bornes dans midinf'
       b=1.0_qp/aa !These two statements change the limits of integration accordingly
       a= 1.0_qp/bb
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s=s/3.0_qp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x)) :: func
        func=funk(1.0_qp/x,arg)/x**2
        END FUNCTION func
       END SUBROUTINE midinfq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midinfvq(funk,aa,bb,m,arg,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: aa,bb
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       REAL(QP), INTENT(INOUT), DIMENSION(m) :: s
       INTEGER(I4B), INTENT(IN) :: m,n
       INTERFACE
        FUNCTION funk(x,arg,mm)
         USE nrtype
         INTEGER, INTENT(IN) ::  mm
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), DIMENSION(size(x),mm) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended midpoint rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding
!       (2/3)×3n-1 additional interior points. s should not be modified between sequential calls.
       REAL(QP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       if(aa*bb <= 0.0) STOP 'bornes dans midinf'
       b=1.0_qp/aa !These two statements change the limits of integration accordingly
       a= 1.0_qp/bb
       if (n == 1) then
        s(:)=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg,m ),dim=1)
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it) 
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s(:)=s(:)/3.0_qp+del*sum(func(x,arg,m),dim=1) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg,mm) !This internal function effects the change of variable.
        INTEGER, INTENT(IN) ::  mm
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x),mm) :: func
        INTEGER is
        func=funk(1.0_qp/x,arg,mm)
        do is=1,size(x)
         func(is,:)=func(is,:)/x(is)**2
        enddo
        END FUNCTION func
       END SUBROUTINE midinfvq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midinfc(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: aa,bb
       REAL(DP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(DPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(DPC), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(DP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(DP), DIMENSION(2*3**(n-2)) :: x
       if(aa*bb <= 0.0) STOP 'bornes dans midinfc'
       b=1.0_dp/aa !These two statements change the limits of integration accordingly
       a= 1.0_dp/bb
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_dp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_dp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_dp*del,3.0_dp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_dp*del
        s=s/3.0_dp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(DPC), DIMENSION(size(x)) :: func
        func=funk(1.0_dp/x,arg)/x**2
        END FUNCTION func
       END SUBROUTINE midinfc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midinfcq(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: aa,bb
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(QPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(QPC), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(QP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       if(aa*bb <= 0.0) STOP 'bornes dans midinfcq'
       b=1.0_qp/aa !These two statements change the limits of integration accordingly
       a= 1.0_qp/bb
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s=s/3.0_qp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(QPC), DIMENSION(size(x)) :: func
        func=funk(1.0_qp/x,arg)/x**2
        END FUNCTION func
       END SUBROUTINE midinfcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE racinfq(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: aa,bb
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       REAL(QP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for miqpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(QP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       if(aa*bb <= 0.0) STOP 'bornes dans racinf'
       b=1.0_qp/sqrt(aa) !These two statements change the limits of integration accordingly
       a= 1.0_qp/sqrt(bb)
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s=s/3.0_qp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x)) :: func
        func=2.0_qp*funk(1.0_qp/x**2,arg)/x**3
        END FUNCTION func
       END SUBROUTINE racinfq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE racinfvq(funk,aa,bb,m,arg,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: aa,bb
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       REAL(QP), INTENT(INOUT), DIMENSION(m) :: s
       INTEGER(I4B), INTENT(IN) :: m,n
       INTERFACE
        FUNCTION funk(x,arg,mm)
         USE nrtype
         INTEGER, INTENT(IN) ::  mm
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), DIMENSION(size(x),mm) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended midpoint rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding
!       (2/3)×3n-1 additional interior points. s should not be modified between sequential calls.
       REAL(QP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       if(aa*bb <= 0.0) STOP 'bornes dans racinf'
       b=1.0_qp/sqrt(aa) !These two statements change the limits of integration accordingly
       a= 1.0_qp/sqrt(bb)
       if (n == 1) then
        s(:)=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg,m ),dim=1)
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it) 
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s(:)=s(:)/3.0_qp+del*sum(func(x,arg,m),dim=1) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg,mm) !This internal function effects the change of variable.
        INTEGER, INTENT(IN) ::  mm
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x),mm) :: func
        INTEGER is
        func=funk(1.0_qp/x**2,arg,mm)
        do is=1,size(x)
         func(is,:)=2.0_qp*func(is,:)/x(is)**3
        enddo
        END FUNCTION func
       END SUBROUTINE racinfvq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE racinfcq(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: aa,bb
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(QPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(QPC), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for miqpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(QP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       if(aa*bb <= 0.0) STOP 'bornes dans racinf'
       b=1.0_qp/sqrt(aa) !These two statements change the limits of integration accordingly
       a= 1.0_qp/sqrt(bb)
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s=s/3.0_qp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(QPC), DIMENSION(size(x)) :: func
        func=2.0_qp*funk(1.0_qp/x**2,arg)/x**3
        END FUNCTION func
       END SUBROUTINE racinfcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqus(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: aa,bb
       REAL(SP), INTENT(IN), DIMENSION(:) :: arg
       REAL(SP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(SP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(SP), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(SP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(SP), DIMENSION(2*3**(n-2)) :: x
       b=sqrt(bb-aa) !These two statements change the limits of integration accordingly
       a= 0.0
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_sp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_dp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_dp*del,3.0_dp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_dp*del
        s=s/3.0_dp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(SP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(SP), DIMENSION(size(x)) :: func
        func=2.0_dp*x*funk(bb-x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqud(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: aa,bb
       REAL(DP), INTENT(IN), DIMENSION(:) :: arg
       REAL(DP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(DP), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(DP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(DP), DIMENSION(2*3**(n-2)) :: x
       b=sqrt(bb-aa) !These two statements change the limits of integration accordingly
       a= 0.0
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_dp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_dp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_dp*del,3.0_dp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_dp*del
        s=s/3.0_dp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(DP), DIMENSION(size(x)) :: func
        func=2.0_dp*x*funk(bb-x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqud
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsquq(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: aa,bb
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       REAL(QP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for miqpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(QP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       b=sqrt(bb-aa) !These two statements change the limits of integration accordingly
       a= 0.0
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s=s/3.0_qp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x)) :: func
        func=2.0_qp*x*funk(bb-x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsquq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsquc(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: aa,bb
       REAL(DP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(DPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(DPC), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(DP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(DP), DIMENSION(2*3**(n-2)) :: x
       b=sqrt(bb-aa) !These two statements change the limits of integration accordingly
       a= 0.0
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_dp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_dp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_dp*del,3.0_dp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_dp*del
        s=s/3.0_dp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(DPC), DIMENSION(size(x)) :: func
        func=2.0_dp*x*funk(bb-x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsquc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqucq(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: aa,bb
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(QPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(QPC), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(QP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       b=sqrt(bb-aa) !These two statements change the limits of integration accordingly
       a= 0.0
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s=s/3.0_qp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(QPC), DIMENSION(size(x)) :: func
        func=2.0_qp*x*funk(bb-x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqucq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqls(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: aa,bb
       REAL(SP), INTENT(IN), DIMENSION(:) :: arg
       REAL(SP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(SP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(SP), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(SP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(SP), DIMENSION(2*3**(n-2)) :: x
       b=sqrt(bb-aa) !These two statements change the limits of integration accordingly
       a= 0.0
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_sp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_sp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
        s=s/3.0_sp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(SP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(SP), DIMENSION(size(x)) :: func
        func=2.0_sp*x*funk(aa+x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqld(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: aa,bb
       REAL(DP), INTENT(IN), DIMENSION(:) :: arg
       REAL(DP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(DP), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(DP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(DP), DIMENSION(2*3**(n-2)) :: x
       b=sqrt(bb-aa) !These two statements change the limits of integration accordingly
       a= 0.0
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_dp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_dp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_dp*del,3.0_dp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_dp*del
        s=s/3.0_dp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(DP), DIMENSION(size(x)) :: func
        func=2.0_dp*x*funk(aa+x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqld
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqlq(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: aa,bb
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       REAL(QP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(QP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       b=sqrt(bb-aa) !These two statements change the limits of integration accordingly
       a= 0.0
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s=s/3.0_qp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x)) :: func
        func=2.0_qp*x*funk(aa+x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqlq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqlc(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: aa,bb
       REAL(DP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(DPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(DPC), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(DP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(DP), DIMENSION(2*3**(n-2)) :: x
       b=sqrt(bb-aa) !These two statements change the limits of integration accordingly
       a= 0.0
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_dp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_dp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_dp*del,3.0_dp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_dp*del
        s=s/3.0_dp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(DPC), DIMENSION(size(x)) :: func
        func=2.0_dp*x*funk(aa+x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqlc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqlcq(funk,aa,bb,arg,s,n)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: aa,bb
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(QPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION funk(x,arg)
        USE nrtype
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(QPC), DIMENSION(size(x)) :: funk
        END FUNCTION funk
       END INTERFACE
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       REAL(QP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       b=sqrt(bb-aa) !These two statements change the limits of integration accordingly
       a= 0.0
       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s=s/3.0_qp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(QPC), DIMENSION(size(x)) :: func
        func=2.0_qp*x*funk(aa+x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqlcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qrombs(func,a,b,ar)
       IMPLICIT NONE
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(SP), DIMENSION(:), INTENT(IN) :: x
         REAL(SP), DIMENSION(:), INTENT(IN) :: arg
         REAL(SP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), INTENT(IN), DIMENSION(:) :: ar
       REAL(SP) :: qrombs
       INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(SP), PARAMETER :: EPS=1.0e-7_sp
!       Returns the integral of the function func from a to b. Integration is performed by Romberg’s
!       method of order 2K, where, e.g., K=2 is Simpson’s rule.
!       Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation error
!       estimate; JMAX limits the total number of steps; K is the number of points used in the
!       extrapolation.
       REAL(SP), DIMENSION(JMAXP) :: h,s !These store the successive trapezoidal approximations and their relative stepsizes
       REAL(SP) :: dqromb 
       INTEGER(I4B) :: j
!       write(6,*)'qrombs'
       h(1)=1.0
       do j=1,JMAX
        call trapz(func,a,b,ar,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qrombs,dqromb)
         if (abs(dqromb) <= EPS*abs(qrombs)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=0.25_sp*h(j) !This is a key step: The factor is 0.25 even
!        though the stepsize is decreased by only
!        0.5. This makes the extrapolation a polynomial
!        in h2 as allowed by equation (4.2.1),
!        not just a polynomial in h.
       end do
       STOP 'qromb:nbr d itérations dépassé'
       END FUNCTION qrombs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qrombd(func,a,b,ar)
       IMPLICIT NONE
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(DP), DIMENSION(:), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         REAL(DP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), INTENT(IN), DIMENSION(:) :: ar
       REAL(DP) :: qrombd
       INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(DP), PARAMETER :: EPS=1.0e-6_dp
!       Returns the integral of the function func from a to b. Integration is performed by Romberg’s
!       method of order 2K, where, e.g., K=2 is Simpson’s rule.
!       Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation error
!       estimate; JMAX limits the total number of steps; K is the number of points used in the
!       extrapolation.
       REAL(DP), DIMENSION(JMAXP) :: h,s !These store the successive trapezoidal approximations and their relative stepsizes
       REAL(DP) :: dqromb 
       INTEGER(I4B) :: j
!       write(6,*)'qrombd'
       h(1)=1.0
       do j=1,JMAX
        call trapz(func,a,b,ar,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qrombd,dqromb)
         if (abs(dqromb) <= EPS*abs(qrombd)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=0.25_dp*h(j) !This is a key step: The factor is 0.25 even
!        though the stepsize is decreased by only
!        0.5. This makes the extrapolation a polynomial
!        in h2 as allowed by equation (4.2.1),
!        not just a polynomial in h.
       end do
       STOP 'qromb:nbr d itérations dépassé'
       END FUNCTION qrombd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qrombq(func,a,b,ar)
       IMPLICIT NONE
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: ar
       REAL(QP) :: qrombq
       INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(QP), PARAMETER :: EPS=1.0e-6_qp
!       Returns the integral of the function func from a to b. Integration is performed by Romberg’s
!       method of order 2K, where, e.g., K=2 is Simpson’s rule.
!       Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation error
!       estimate; JMAX limits the total number of steps; K is the number of points used in the
!       extrapolation.
       REAL(QP), DIMENSION(JMAXP) :: h,s !These store the successive trapezoidal approximations and their relative stepsizes
       REAL(QP) :: dqromb 
       INTEGER(I4B) :: j
!       write(6,*)'qrombq'
       h(1)=1.0
       do j=1,JMAX
        call trapz(func,a,b,ar,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_qp,qrombq,dqromb)
         if (abs(dqromb) <= EPS*abs(qrombq)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=0.25_qp*h(j) !This is a key step: The factor is 0.25 even
!        though the stepsize is decreased by only
!        0.5. This makes the extrapolation a polynomial
!        in h2 as allowed by equation (4.2.1),
!        not just a polynomial in h.
       end do
       STOP 'qromb:nbr d itérations dépassé'
       END FUNCTION qrombq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qrombc(func,a,b,ar)
       IMPLICIT NONE
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(DP), DIMENSION(:), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(DPC), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), INTENT(IN), DIMENSION(:) :: ar
       COMPLEX(DPC) :: qrombc
       INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(DP), PARAMETER :: EPS=1.0e-8_dp
!       Returns the integral of the function func from a to b. Integration is performed by Romberg’s
!       method of order 2K, where, e.g., K=2 is Simpson’s rule.
!       Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation error
!       estimate; JMAX limits the total number of steps; K is the number of points used in the
!       extrapolation.
       REAL(DP), DIMENSION(JMAXP) :: h !These store the successive trapezoidal approximations and their relative stepsizes
       COMPLEX(DPC), DIMENSION(JMAXP) :: s !These store the successive trapezoidal approximations and their relative stepsizes
       COMPLEX(DPC) :: dqromb 
       INTEGER(I4B) :: j
!       write(6,*)'qrombc'
       h(1)=1.0
       do j=1,JMAX
        call trapz(func,a,b,ar,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qrombc,dqromb)
         if (abs(dqromb) <= EPS*abs(qrombc)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=0.25_dp*h(j) !This is a key step: The factor is 0.25 even
!        though the stepsize is decreased by only
!        0.5. This makes the extrapolation a polynomial
!        in h2 as allowed by equation (4.2.1),
!        not just a polynomial in h.
       end do
       STOP 'qromb:nbr d itérations dépassé'
       END FUNCTION qrombc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qrombcq(func,a,b,ar)
       IMPLICIT NONE
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: ar
       COMPLEX(QPC) :: qrombcq
       INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(QP), PARAMETER :: EPS=1.0e-8_qp
!       Returns the integral of the function func from a to b. Integration is performed by Romberg’s
!       method of order 2K, where, e.g., K=2 is Simpson’s rule.
!       Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation error
!       estimate; JMAX limits the total number of steps; K is the number of points used in the
!       extrapolation.
       REAL(QP), DIMENSION(JMAXP) :: h !These store the successive trapezoidal approximations and their relative stepsizes
       COMPLEX(QPC), DIMENSION(JMAXP) :: s !These store the successive trapezoidal approximations and their relative stepsizes
       COMPLEX(QPC) :: dqromb 
       INTEGER(I4B) :: j
!       write(6,*)'qrombcq'
       h(1)=1.0
       do j=1,JMAX
        call trapz(func,a,b,ar,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_qp,qrombcq,dqromb)
         if (abs(dqromb) <= EPS*abs(qrombcq)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=0.25_qp*h(j) !This is a key step: The factor is 0.25 even
!        though the stepsize is decreased by only
!        0.5. This makes the extrapolation a polynomial
!        in h2 as allowed by equation (4.2.1),
!        not just a polynomial in h.
       end do
       STOP 'qromb:nbr d itérations dépassé'
       END FUNCTION qrombcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromos(func,a,b,arg,choose,EPS)
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), DIMENSION(:), INTENT(IN) :: arg
       REAL(SP) :: qromos
       INTERFACE
         FUNCTION func(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(SP), DIMENSION(:), INTENT(IN) :: x
         REAL(SP), DIMENSION(:), INTENT(IN) :: arg
         REAL(SP), DIMENSION(size(x)) :: func
         END FUNCTION func
         SUBROUTINE choose(funk,aa,bb,arg,s,n)
         USE nrtype
         IMPLICIT NONE
         REAL(SP), INTENT(IN) :: aa,bb
         REAL(SP), DIMENSION(:), INTENT(IN) :: arg
         REAL(SP), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         INTERFACE
          FUNCTION funk(x,arg)
          USE nrtype
          IMPLICIT NONE
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(:), INTENT(IN) :: arg
          REAL(SP), DIMENSION(size(x)) :: funk
          END FUNCTION funk
         END INTERFACE
         END SUBROUTINE choose
       END INTERFACE
       INTEGER(I4B), PARAMETER :: JMAX=14,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(SP) :: EPS
!         Romberg integration on an open interval. Returns the integral of the function func from a
!         to b, using any specified integrating subroutine choose and Romberg’s method. Normally
!         choose will be an open formula, not evaluating the function at the endpoints. It is assumed
!         that choose triples the number of steps on each call, and that its error series contains only
!         even powers of the number of steps. The routines midpnt, midinf, midsql, midsqu,
!         and midexp are possible choices for choose. The parameters have the same meaning as in qromb.
       REAL(SP), DIMENSION(JMAXP) :: h,s
       REAL(SP) :: dqromo
       INTEGER(I4B) :: j
       h(1)=1.0
       do j=1,JMAX
        call choose(func,a,b,arg,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qromos,dqromo)
         if (abs(dqromo) <= EPS*abs(qromos)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=h(j)/9.0_sp !This is where the assumption of step tripling and an even error series is used.
       end do
       STOP 'Nombre d itération dépassé dans qromos'
       END FUNCTION qromos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromod(func,a,b,arg,choose,EPS)
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), DIMENSION(:), INTENT(IN) :: arg
       REAL(DP) :: qromod
       INTERFACE
         FUNCTION func(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(DP), DIMENSION(:), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         REAL(DP), DIMENSION(size(x)) :: func
         END FUNCTION func
         SUBROUTINE choose(funk,aa,bb,arg,s,n)
         USE nrtype
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: aa,bb
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         REAL(DP), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         INTERFACE
          FUNCTION funk(x,arg)
          USE nrtype
          IMPLICIT NONE
          REAL(DP), DIMENSION(:), INTENT(IN) :: x
          REAL(DP), DIMENSION(:), INTENT(IN) :: arg
          REAL(DP), DIMENSION(size(x)) :: funk
          END FUNCTION funk
         END INTERFACE
         END SUBROUTINE choose
       END INTERFACE
       INTEGER(I4B), PARAMETER :: JMAX=15,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(DP) :: EPS
!         Romberg integration on an open interval. Returns the integral of the function func from a
!         to b, using any specified integrating subroutine choose and Romberg’s method. Normally
!         choose will be an open formula, not evaluating the function at the endpoints. It is assumed
!         that choose triples the number of steps on each call, and that its error series contains only
!         even powers of the number of steps. The routines midpnt, midinf, midsql, midsqu,
!         and midexp are possible choices for choose. The parameters have the same meaning as in qromb.
       REAL(DP), DIMENSION(JMAXP) :: h,s
       REAL(DP) :: dqromo
       INTEGER(I4B) :: j
       h(1)=1.0
       do j=1,JMAX
        call choose(func,a,b,arg,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qromod,dqromo)
         if (abs(dqromo) <= EPS*abs(qromod)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=h(j)/9.0_dp !This is where the assumption of step tripling and an even error series is used.
       end do
!       STOP 'Nombre d itération dépassé dans qromo'
       write(6,*) 'Nombre d itération dépassé dans qromod'
       END FUNCTION qromod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromoq(func,a,b,arg,choose,EPS)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), DIMENSION(:), INTENT(IN) :: arg
       REAL(QP) :: qromoq
       INTERFACE
         FUNCTION func(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), DIMENSION(size(x)) :: func
         END FUNCTION func
         SUBROUTINE choose(funk,aa,bb,arg,s,n)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), INTENT(IN) :: aa,bb
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         INTERFACE
          FUNCTION funk(x,arg)
          USE nrtype
          IMPLICIT NONE
          REAL(QP), DIMENSION(:), INTENT(IN) :: x
          REAL(QP), DIMENSION(:), INTENT(IN) :: arg
          REAL(QP), DIMENSION(size(x)) :: funk
          END FUNCTION funk
         END INTERFACE
         END SUBROUTINE choose
       END INTERFACE
       INTEGER(I4B), PARAMETER :: JMAX=16,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(QP) :: EPS
!         Romberg integration on an open interval. Returns the integral of the function func from a
!         to b, using any specified integrating subroutine choose and Romberg’s method. Normally
!         choose will be an open formula, not evaluating the function at the endpoints. It is assumed
!         that choose triples the number of steps on each call, and that its error series contains only
!         even powers of the number of steps. The routines midpnt, midinf, midsql, midsqu,
!         and midexp are possible choices for choose. The parameters have the same meaning as in qromb.
       REAL(QP), DIMENSION(JMAXP) :: h,s
       REAL(QP) :: dqromo
       INTEGER(I4B) :: j
       h(1)=1.0
       do j=1,JMAX
        call choose(func,a,b,arg,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_qp,qromoq,dqromo)
         if (abs(dqromo) <= EPS*abs(qromoq)) RETURN
         if (abs(qromoq) <= EPS*10e-10_qp) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=h(j)/9.0_qp !This is where the assumption of step tripling and an even error series is used.
       end do
!       STOP 'Nombre d itération dépassé dans qromo'
       write(6,*) 'Nombre d itération dépassé dans qromoq'
       END FUNCTION qromoq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromovq(func,a,b,m,arg,choose,EPS)
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: m !Dimension du vecteur à intégrer
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), DIMENSION(:), INTENT(IN) :: arg
       REAL(QP), DIMENSION(m) :: qromovq
       INTERFACE
         FUNCTION func(x,arg,mm)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) ::  mm
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), DIMENSION(size(x),mm) :: func
         END FUNCTION func
         SUBROUTINE choose(funk,aa,bb,mm,arg,s,n)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: mm
         REAL(QP), INTENT(IN) :: aa,bb
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), DIMENSION(mm), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         INTERFACE
          FUNCTION funk(x,arg,mmm)
          USE nrtype
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: mmm
          REAL(QP), DIMENSION(:), INTENT(IN) :: x
          REAL(QP), DIMENSION(:), INTENT(IN) :: arg
          REAL(QP), DIMENSION(size(x),mmm) :: funk
          END FUNCTION funk
         END INTERFACE
         END SUBROUTINE choose
       END INTERFACE
       INTEGER(I4B), PARAMETER :: JMAX=16,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(QP) :: EPS
!         Romberg integration on an open interval. Returns the integral of the function func from a
!         to b, using any specified integrating subroutine choose and Romberg’s method. Normally
!         choose will be an open formula, not evaluating the function at the endpoints. It is assumed
!         that choose triples the number of steps on each call, and that its error series contains only
!         even powers of the number of steps. The routines midpnt, midinf, midsql, midsqu,
!         and midexp are possible choices for choose. The parameters have the same meaning as in qromb.
       REAL(QP), DIMENSION(JMAXP,m) :: h,s
       REAL(QP) :: dqromo(m)
       INTEGER(I4B) :: j,im
       LOGICAL conv(m)
       h(1,:)=1.0
       do j=1,JMAX
        call choose(func,a,b,m,arg,s(j,:),j)
        if (j >= K) then
         do im=1,m
          call polint(h(j-KM:j,im),s(j-KM:j,im),0.0_qp,qromovq(im),dqromo(im))
          conv(im)=abs(dqromo(im)) <= EPS*abs(qromovq(im))
         enddo
         if (all(conv)) RETURN
        end if
        s(j+1,:)=s(j,:)
        h(j+1,:)=h(j,:)/9.0_qp !This is where the assumption of step tripling and an even error series is used.
       end do
!       STOP 'Nombre d itération dépassé dans qromovq'
       write(6,*) 'Nombre d itération dépassé dans qromovq'
       END FUNCTION qromovq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromoc(func,a,b,arg,choose,EPS)
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), DIMENSION(:), INTENT(IN) :: arg
       COMPLEX(DPC) :: qromoc
       INTERFACE
         FUNCTION func(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(DP), DIMENSION(:), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(DPC), DIMENSION(size(x)) :: func
         END FUNCTION func
         SUBROUTINE choose(funk,aa,bb,arg,s,n)
         USE nrtype
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: aa,bb
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(DPC), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         INTERFACE
          FUNCTION funk(x,arg)
          USE nrtype
          IMPLICIT NONE
          REAL(DP), DIMENSION(:), INTENT(IN) :: x
          REAL(DP), DIMENSION(:), INTENT(IN) :: arg
          COMPLEX(DPC), DIMENSION(size(x)) :: funk
          END FUNCTION funk
         END INTERFACE
         END SUBROUTINE choose
       END INTERFACE
       INTEGER(I4B), PARAMETER :: JMAX=14,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(DP) :: EPS
!         Romberg integration on an open interval. Returns the integral of the function func from a
!         to b, using any specified integrating subroutine choose and Romberg’s method. Normally
!         choose will be an open formula, not evaluating the function at the endpoints. It is assumed
!         that choose triples the number of steps on each call, and that its error series contains only
!         even powers of the number of steps. The routines midpnt, midinf, midsql, midsqu,
!         and midexp are possible choices for choose. The parameters have the same meaning as in qromb.
       REAL(DP), DIMENSION(JMAXP) :: h
       COMPLEX(DPC), DIMENSION(JMAXP) :: s
       COMPLEX(DPC) :: dqromo
       INTEGER(I4B) :: j
       h(1)=1.0
       do j=1,JMAX
        call choose(func,a,b,arg,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qromoc,dqromo)
         if (abs(dqromo) <= EPS*abs(qromoc)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=h(j)/9.0_dp !This is where the assumption of step tripling and an even error series is used.
       end do
       STOP 'Nombre d itération dépassé dans qromoc'
       END FUNCTION qromoc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromocq(func,a,b,arg,choose,EPS)
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), DIMENSION(:), INTENT(IN) :: arg
       COMPLEX(QPC) :: qromocq
       INTERFACE
         FUNCTION func(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), DIMENSION(size(x)) :: func
         END FUNCTION func
         SUBROUTINE choose(funk,aa,bb,arg,s,n)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), INTENT(IN) :: aa,bb
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         INTERFACE
          FUNCTION funk(x,arg)
          USE nrtype
          IMPLICIT NONE
          REAL(QP), DIMENSION(:), INTENT(IN) :: x
          REAL(QP), DIMENSION(:), INTENT(IN) :: arg
          COMPLEX(QPC), DIMENSION(size(x)) :: funk
          END FUNCTION funk
         END INTERFACE
         END SUBROUTINE choose
       END INTERFACE
       INTEGER(I4B), PARAMETER :: JMAX=17,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(QP) :: EPS
!         Romberg integration on an open interval. Returns the integral of the function func from a
!         to b, using any specified integrating subroutine choose and Romberg’s method. Normally
!         choose will be an open formula, not evaluating the function at the endpoints. It is assumed
!         that choose triples the number of steps on each call, and that its error series contains only
!         even powers of the number of steps. The routines midpnt, midinf, midsql, midsqu,
!         and midexp are possible choices for choose. The parameters have the same meaning as in qromb.
       REAL(QP), DIMENSION(JMAXP) :: h
       COMPLEX(QPC), DIMENSION(JMAXP) :: s
       COMPLEX(QPC) :: dqromo
       INTEGER(I4B) :: j
       h(1)=1.0
       do j=1,JMAX
        call choose(func,a,b,arg,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_qp,qromocq,dqromo)
         if (abs(dqromo) <= EPS*abs(qromocq)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=h(j)/9.0_qp !This is where the assumption of step tripling and an even error series is used.
       end do
       STOP 'Nombre d itération dépassé dans qromocq'
       END FUNCTION qromocq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE polints(xa,ya,x,y,dy)
       IMPLICIT NONE
       REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
       REAL(SP), INTENT(IN) :: x
       REAL(SP), INTENT(OUT) :: y,dy
!       Given arrays xa and ya of length N, and given a value x, this routine returns a value y,
!       and an error estimate dy. If P(x) is the polynomial of degree N − 1 such that P(xai) =
!       yai, i = 1, . . . ,N, then the returned value y = P(x).
       INTEGER(I4B) :: m,n,ns
       REAL(SP), DIMENSION(size(xa)) :: c,d,den,ho
       n=size(xa)
       if(size(xa).NE.size(ya)) STOP 'Taille des entrees dans polint'
       c=ya !Initialize the tableau of c’s and d’s.
       d=ya
       ho=xa-x
       ns=iminloc(abs(x-xa)) !Find index ns of closest table entry.
       y=ya(ns) !This is the initial approximation to y.
       ns=ns-1
       do m=1,n-1 !For each column of the tableau,
        den(1:n-m)=ho(1:n-m)-ho(1+m:n) !we loop over the current c’s and d’s and upidate them.
        if(any(den(1:n-m) == 0.0)) STOP 'polint: entrées dégénérées'
!       This error can occur only if two input xa’s are (to within roundoff) identical.
        den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
        d(1:n-m)=ho(1+m:n)*den(1:n-m) !Here the c’s and d’s are updated.
        c(1:n-m)=ho(1:n-m)*den(1:n-m)
        if (2*ns < n-m) then 
!       After each column in the tableau is completed, we decide
!       which correction, c or d, we want to add to our accumulating
!       value of y, i.e., which path to take through
!       the tableau—forking up or down. We do this in such a
!       way as to take the most “straight line” route through the
!       tableau to its apex, updating ns accordingly to keep track
!       of where we are. This route keeps the partial approximations
!       centered (insofar as possible) on the target x. The
!       last dy added is thus the error indication.
         dy=c(ns+1)
        else
         dy=d(ns)
         ns=ns-1
        end if
        y=y+dy
       end do
       END SUBROUTINE polints
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE polintd(xa,ya,x,y,dy)
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(OUT) :: y,dy
!       Given arrays xa and ya of length N, and given a value x, this routine returns a value y,
!       and an error estimate dy. If P(x) is the polynomial of degree N − 1 such that P(xai) =
!       yai, i = 1, . . . ,N, then the returned value y = P(x).
       INTEGER(I4B) :: m,n,ns
       REAL(DP), DIMENSION(size(xa)) :: c,d,den,ho
       n=size(xa)
       if(size(xa).NE.size(ya)) STOP 'Taille des entrees dans polint'
       c=ya !Initialize the tableau of c’s and d’s.
       d=ya
       ho=xa-x
       ns=iminloc(abs(x-xa)) !Find index ns of closest table entry.
       y=ya(ns) !This is the initial approximation to y.
       ns=ns-1
       do m=1,n-1 !For each column of the tableau,
        den(1:n-m)=ho(1:n-m)-ho(1+m:n) !we loop over the current c’s and d’s and upidate them.
        if(any(den(1:n-m) == 0.0)) STOP 'polint: entrées dégénérées'
!       This error can occur only if two input xa’s are (to within roundoff) identical.
        den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
        d(1:n-m)=ho(1+m:n)*den(1:n-m) !Here the c’s and d’s are updated.
        c(1:n-m)=ho(1:n-m)*den(1:n-m)
        if (2*ns < n-m) then 
!       After each column in the tableau is completed, we decide
!       which correction, c or d, we want to add to our accumulating
!       value of y, i.e., which path to take through
!       the tableau—forking up or down. We do this in such a
!       way as to take the most “straight line” route through the
!       tableau to its apex, updating ns accordingly to keep track
!       of where we are. This route keeps the partial approximations
!       centered (insofar as possible) on the target x. The
!       last dy added is thus the error indication.
         dy=c(ns+1)
        else
         dy=d(ns)
         ns=ns-1
        end if
        y=y+dy
       end do
       END SUBROUTINE polintd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE polintq(xa,ya,x,y,dy)
       IMPLICIT NONE
       REAL(QP), DIMENSION(:), INTENT(IN) :: xa,ya
       REAL(QP), INTENT(IN) :: x
       REAL(QP), INTENT(OUT) :: y,dy
!       Given arrays xa and ya of length N, and given a value x, this routine returns a value y,
!       and an error estimate dy. If P(x) is the polynomial of degree N − 1 such that P(xai) =
!       yai, i = 1, . . . ,N, then the returned value y = P(x).
       INTEGER(I4B) :: m,n,ns
       REAL(QP), DIMENSION(size(xa)) :: c,d,den,ho
       n=size(xa)
       if(size(xa).NE.size(ya)) STOP 'Taille des entrees dans polint'
       c=ya !Initialize the tableau of c’s and d’s.
       d=ya
       ho=xa-x
       ns=iminloc(abs(x-xa)) !Find index ns of closest table entry.
       y=ya(ns) !This is the initial approximation to y.
       ns=ns-1
       do m=1,n-1 !For each column of the tableau,
        den(1:n-m)=ho(1:n-m)-ho(1+m:n) !we loop over the current c’s and d’s and upidate them.
        if(any(den(1:n-m) == 0.0)) STOP 'polint: entrées dégénérées'
!       This error can occur only if two input xa’s are (to within roundoff) identical.
        den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
        d(1:n-m)=ho(1+m:n)*den(1:n-m) !Here the c’s and d’s are updated.
        c(1:n-m)=ho(1:n-m)*den(1:n-m)
        if (2*ns < n-m) then 
!       After each column in the tableau is completed, we decide
!       which correction, c or d, we want to add to our accumulating
!       value of y, i.e., which path to take through
!       the tableau—forking up or down. We do this in such a
!       way as to take the most “straight line” route through the
!       tableau to its apex, updating ns accordingly to keep track
!       of where we are. This route keeps the partial approximations
!       centered (insofar as possible) on the target x. The
!       last dy added is thus the error indication.
         dy=c(ns+1)
        else
         dy=d(ns)
         ns=ns-1
        end if
        y=y+dy
       end do
       END SUBROUTINE polintq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE polintc(xa,ya,x,y,dy)
       IMPLICIT NONE
       COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: ya
       REAL(DP), DIMENSION(:), INTENT(IN) :: xa
       REAL(DP), INTENT(IN) :: x
       COMPLEX(DPC), INTENT(OUT) :: y,dy
!       Given arrays xa and ya of length N, and given a value x, this routine returns a value y,
!       and an error estimate dy. If P(x) is the polynomial of degree N − 1 such that P(xai) =
!       yai, i = 1, . . . ,N, then the returned value y = P(x).
       INTEGER(I4B) :: m,n,ns
       COMPLEX(DPC), DIMENSION(size(xa)) :: c,d,den
       REAL(DP), DIMENSION(size(xa)) :: ho
       n=size(xa)
       if(size(xa).NE.size(ya)) STOP 'Taille des entrees dans polint'
       c=ya !Initialize the tableau of c’s and d’s.
       d=ya
       ho=xa-x
       ns=iminloc(abs(x-xa)) !Find index ns of closest table entry.
       y=ya(ns) !This is the initial approximation to y.
       ns=ns-1
       do m=1,n-1 !For each column of the tableau,
        den(1:n-m)=ho(1:n-m)-ho(1+m:n) !we loop over the current c’s and d’s and upidate them.
        if(any(abs(den(1:n-m)) == 0.0)) STOP 'polint: entrées dégénérées'
!       This error can occur only if two input xa’s are (to within roundoff) identical.
        den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
        d(1:n-m)=ho(1+m:n)*den(1:n-m) !Here the c’s and d’s are updated.
        c(1:n-m)=ho(1:n-m)*den(1:n-m)
        if (2*ns < n-m) then 
!       After each column in the tableau is completed, we decide
!       which correction, c or d, we want to add to our accumulating
!       value of y, i.e., which path to take through
!       the tableau—forking up or down. We do this in such a
!       way as to take the most “straight line” route through the
!       tableau to its apex, updating ns accordingly to keep track
!       of where we are. This route keeps the partial approximations
!       centered (insofar as possible) on the target x. The
!       last dy added is thus the error indication.
         dy=c(ns+1)
        else
         dy=d(ns)
         ns=ns-1
        end if
        y=y+dy
       end do
       END SUBROUTINE polintc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE polintcq(xa,ya,x,y,dy)
       IMPLICIT NONE
       COMPLEX(QPC), DIMENSION(:), INTENT(IN) :: ya
       REAL(QP), DIMENSION(:), INTENT(IN) :: xa
       REAL(QP), INTENT(IN) :: x
       COMPLEX(QPC), INTENT(OUT) :: y,dy
!       Given arrays xa and ya of length N, and given a value x, this routine returns a value y,
!       and an error estimate dy. If P(x) is the polynomial of degree N − 1 such that P(xai) =
!       yai, i = 1, . . . ,N, then the returned value y = P(x).
       INTEGER(I4B) :: m,n,ns
       COMPLEX(QPC), DIMENSION(size(xa)) :: c,d,den
       REAL(QP), DIMENSION(size(xa)) :: ho
       n=size(xa)
       if(size(xa).NE.size(ya)) STOP 'Taille des entrees dans polint'
       c=ya !Initialize the tableau of c’s and d’s.
       d=ya
       ho=xa-x
       ns=iminloc(abs(x-xa)) !Find index ns of closest table entry.
       y=ya(ns) !This is the initial approximation to y.
       ns=ns-1
       do m=1,n-1 !For each column of the tableau,
        den(1:n-m)=ho(1:n-m)-ho(1+m:n) !we loop over the current c’s and d’s and upidate them.
        if(any(abs(den(1:n-m)) == 0.0)) STOP 'polint: entrées dégénérées'
!       This error can occur only if two input xa’s are (to within roundoff) identical.
        den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
        d(1:n-m)=ho(1+m:n)*den(1:n-m) !Here the c’s and d’s are updated.
        c(1:n-m)=ho(1:n-m)*den(1:n-m)
        if (2*ns < n-m) then 
!       After each column in the tableau is completed, we decide
!       which correction, c or d, we want to add to our accumulating
!       value of y, i.e., which path to take through
!       the tableau—forking up or down. We do this in such a
!       way as to take the most “straight line” route through the
!       tableau to its apex, updating ns accordingly to keep track
!       of where we are. This route keeps the partial approximations
!       centered (insofar as possible) on the target x. The
!       last dy added is thus the error indication.
         dy=c(ns+1)
        else
         dy=d(ns)
         ns=ns-1
        end if
        y=y+dy
       end do
       END SUBROUTINE polintcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION iminloc_s(arr)
       REAL(SP), DIMENSION(:), INTENT(IN) :: arr
       INTEGER(I4B), DIMENSION(1) :: imin
       INTEGER(I4B) :: iminloc_s
       imin=minloc(arr(:))
       iminloc_s=imin(1)
       END FUNCTION iminloc_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION iminloc_d(arr)
       REAL(DP), DIMENSION(:), INTENT(IN) :: arr
       INTEGER(I4B), DIMENSION(1) :: imin
       INTEGER(I4B) :: iminloc
       imin=minloc(arr(:))
       iminloc_d=imin(1)
       END FUNCTION iminloc_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION iminloc_q(arr)
       REAL(QP), DIMENSION(:), INTENT(IN) :: arr
       INTEGER(I4B), DIMENSION(1) :: imin
       INTEGER(I4B) :: iminloc
       imin=minloc(arr(:))
       iminloc_q=imin(1)
       END FUNCTION iminloc_q
END MODULE modsim
