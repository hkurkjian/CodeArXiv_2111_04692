MODULE recettes 
        USE nrtype
        USE nrutil, ONLY : assert
        IMPLICIT NONE
        INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
        INTERFACE tri
                MODULE PROCEDURE tri_s,tri_d,tri_q
        END INTERFACE
        INTERFACE indexx
                MODULE PROCEDURE indexx_sp,indexx_dp,indexx_qp
        END INTERFACE
        INTERFACE laguer
                MODULE PROCEDURE laguer_d,laguer_q
        END INTERFACE
        INTERFACE zroots
                MODULE PROCEDURE zroots_d,zroots_q
        END INTERFACE
        INTERFACE rf
                 MODULE PROCEDURE rf_s,rf_d,rf_q,rf_c,rf_cq
        END INTERFACE
        INTERFACE rd
                 MODULE PROCEDURE rd_s,rd_d,rd_q,rd_c,rd_cq
        END INTERFACE
        INTERFACE rc
                 MODULE PROCEDURE rc_q
        END INTERFACE
        INTERFACE rj
                 MODULE PROCEDURE rj_q
        END INTERFACE
        INTERFACE elle
                 MODULE PROCEDURE elle_s,elle_d,elle_q,elle_c,elle_rc,elle_cq,elle_rcq
        END INTERFACE
        INTERFACE ellf
                 MODULE PROCEDURE ellf_s,ellf_d,ellf_q,ellf_c,ellf_rc,ellf_cq,ellf_rcq
        END INTERFACE
        INTERFACE ellpi
                 MODULE PROCEDURE ellpi_q
        END INTERFACE
        INTERFACE argum
                 MODULE PROCEDURE argum_d,argum_q
        END INTERFACE
!        INTERFACE arth
!                MODULE PROCEDURE arth_r, arth_d, arth_q, arth_i
!        END INTERFACE
        INTERFACE rtsafe
                MODULE PROCEDURE rtsafed, rtsafeq
        END INTERFACE
        INTERFACE mnewt
                MODULE PROCEDURE mnewt_s, mnewt_d, mnewt_q
        END INTERFACE
        INTERFACE ludcmp
                MODULE PROCEDURE ludcmp_s, ludcmp_d, ludcmp_q, ludcmp_cq
        END INTERFACE
        INTERFACE lubksb
                MODULE PROCEDURE lubksb_s, lubksb_d, lubksb_q
        END INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
        SUBROUTINE mnewt_s(ntrial,x,tolx,tolf,usrfun)
        INTEGER(I4B), INTENT(IN) :: ntrial
        REAL(SP), INTENT(IN) :: tolx,tolf
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
        INTERFACE
        SUBROUTINE usrfun(x,fvec,fjac)
        USE nrtype
        IMPLICIT NONE
        REAL(SP), DIMENSION(:), INTENT(IN) :: x
        REAL(SP), DIMENSION(:), INTENT(OUT) :: fvec
        REAL(SP), DIMENSION(:,:), INTENT(OUT) :: fjac
        END SUBROUTINE usrfun
        END INTERFACE
!        Given an initial guess x for a root in N dimensions, take ntrial Newton-Raphson steps to
!        improve the root. Stop if the root converges in either summed absolute variable increments
!        tolx or summed absolute function values tolf.
        INTEGER(I4B) :: i
        INTEGER(I4B), DIMENSION(size(x)) :: indx
        REAL(SP) :: d
        REAL(SP), DIMENSION(size(x)) :: fvec,p
        REAL(SP), DIMENSION(size(x),size(x)) :: fjac
        do i=1,ntrial
         call usrfun(x,fvec,fjac)
!        User subroutine supplies function values at x in fvec and Jacobian matrix in fjac.
         if (sum(abs(fvec)) <= tolf) RETURN !Check function convergence.
         p=-fvec !Right-hand side of linear equations.
         call ludcmp(fjac,indx,d) !Solve linear equations using LU decomposition.
         call lubksb(fjac,indx,p) 
         x=x+p !Update solution.
         if (sum(abs(p)) <= tolx) RETURN !Check root convergence.
        end do
        END SUBROUTINE mnewt_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE mnewt_d(ntrial,x,tolx,tolf,usrfun)
        INTEGER(I4B), INTENT(IN) :: ntrial
        REAL(DP), INTENT(IN) :: tolx,tolf
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
        INTERFACE
        SUBROUTINE usrfun(x,fvec,fjac)
        USE nrtype
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: x
        REAL(DP), DIMENSION(:), INTENT(OUT) :: fvec
        REAL(DP), DIMENSION(:,:), INTENT(OUT) :: fjac
        END SUBROUTINE usrfun
        END INTERFACE
!        Given an initial guess x for a root in N dimensions, take ntrial Newton-Raphson steps to
!        improve the root. Stop if the root converges in either summed absolute variable increments
!        tolx or summed absolute function values tolf.
        INTEGER(I4B) :: i
        INTEGER(I4B), DIMENSION(size(x)) :: indx
        REAL(DP) :: d
        REAL(DP), DIMENSION(size(x)) :: fvec,p
        REAL(DP), DIMENSION(size(x),size(x)) :: fjac
        do i=1,ntrial
         call usrfun(x,fvec,fjac)
!        User subroutine supplies function values at x in fvec and Jacobian matrix in fjac.
         if (sum(abs(fvec)) <= tolf) RETURN !Check function convergence.
         p=-fvec !Right-hand side of linear equations.
         call ludcmp(fjac,indx,d) !Solve linear equations using LU decomposition.
         call lubksb(fjac,indx,p) 
         x=x+p !Update solution.
         if (sum(abs(p)) <= tolx) RETURN !Check root convergence.
        end do
        END SUBROUTINE mnewt_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE mnewt_q(ntrial,x,tolx,tolf,usrfun)
        INTEGER(I4B), INTENT(IN) :: ntrial
        REAL(QP), INTENT(IN) :: tolx,tolf
        REAL(QP), DIMENSION(:), INTENT(INOUT) :: x
        INTERFACE
        SUBROUTINE usrfun(x,fvec,fjac)
        USE nrtype
        IMPLICIT NONE
        REAL(QP), DIMENSION(:), INTENT(IN) :: x
        REAL(QP), DIMENSION(:), INTENT(OUT) :: fvec
        REAL(QP), DIMENSION(:,:), INTENT(OUT) :: fjac
        END SUBROUTINE usrfun
        END INTERFACE
!        Given an initial guess x for a root in N dimensions, take ntrial Newton-Raphson steps to
!        improve the root. Stop if the root converges in either summed absolute variable increments
!        tolx or summed absolute function values tolf.
        INTEGER(I4B) :: i
        INTEGER(I4B), DIMENSION(size(x)) :: indx
        REAL(QP) :: d
        REAL(QP), DIMENSION(size(x)) :: fvec,p
        REAL(QP), DIMENSION(size(x),size(x)) :: fjac
        do i=1,ntrial
         call usrfun(x,fvec,fjac)
!        User subroutine supplies function values at x in fvec and Jacobian matrix in fjac.
         if (sum(abs(fvec)) <= tolf) RETURN !Check function convergence.
         p=-fvec !Right-hand side of linear equations.
         call ludcmp(fjac,indx,d) !Solve linear equations using LU decomposition.
         call lubksb(fjac,indx,p) 
         x=x+p !Update solution.
         if (sum(abs(p)) <= tolx) RETURN !Check root convergence.
        end do
        END SUBROUTINE mnewt_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE indexx_sp(arr,index)
        USE nrutil, ONLY :arth,assert_eq,nrerror,swap
        IMPLICIT NONE
        REAL(SP), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
        INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
!        Indexes an array arr, i.e., outputs the array index of length N such that arr(index(j))
!        is in ascending order for j = 1, 2, . . . ,N. The input quantity arr is not changed.
        REAL(SP) :: a
        INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
        INTEGER(I4B), DIMENSION(NSTACK) :: istack
        n=assert_eq(size(index),size(arr),'indexx_sp')
        index=arth(1,1,n)
        jstack=0
        l=1
        r=n
        do
         if (r-l < NN) then
          do j=l+1,r
           indext=index(j)
           a=arr(indext)
           do i=j-1,l,-1
            if (arr(index(i)) <= a) exit
            index(i+1)=index(i)
           end do
           index(i+1)=indext
          end do
          if (jstack == 0) RETURN
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
         else
          k=(l+r)/2
          call swap(index(k),index(l+1))
          call icomp_xchg(index(l),index(r))
          call icomp_xchg(index(l+1),index(r))
          call icomp_xchg(index(l),index(l+1))
          i=l+1
          j=r
          indext=index(l+1)
          a=arr(indext)
          do
           do
            i=i+1
            if (arr(index(i)) >= a) exit
           end do
           do
            j=j-1
            if (arr(index(j)) <= a) exit
           end do
           if (j < i) exit
           call swap(index(i),index(j))
          end do
          index(l+1)=index(j)
          index(j)=indext
          jstack=jstack+2
          if (jstack > NSTACK) then
           write(6,*)'indexx_sp:NSTACK too small'
           STOP
          endif
          if (r-i+1 >= j-l) then
           istack(jstack)=r
           istack(jstack-1)=i
           r=j-1
          else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
          end if
         end if
        end do
        CONTAINS
        SUBROUTINE icomp_xchg(i,j)
        INTEGER(I4B), INTENT(INOUT) :: i,j
        INTEGER(I4B) :: swp
        if (arr(j) < arr(i)) then
        swp=i
        i=j
        j=swp
        end if
        END SUBROUTINE icomp_xchg
        END SUBROUTINE indexx_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE indexx_dp(arr,index)
        USE nrutil, ONLY :arth,assert_eq,nrerror,swap
        REAL(DP), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
        INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
!        Indexes an array arr, i.e., outputs the array index of length N such that arr(index(j))
!        is in ascending order for j = 1, 2, . . . ,N. The input quantity arr is not changed.
        REAL(DP) :: a
        INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
        INTEGER(I4B), DIMENSION(NSTACK) :: istack
        n=assert_eq(size(index),size(arr),'indexx_dp')
        index=arth(1,1,n)
        jstack=0
        l=1
        r=n
        do
         if (r-l < NN) then
          do j=l+1,r
           indext=index(j)
           a=arr(indext)
           do i=j-1,l,-1
            if (arr(index(i)) <= a) exit
            index(i+1)=index(i)
           end do
           index(i+1)=indext
          end do
          if (jstack == 0) RETURN
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
         else
          k=(l+r)/2
          call swap(index(k),index(l+1))
          call icomp_xchg(index(l),index(r))
          call icomp_xchg(index(l+1),index(r))
          call icomp_xchg(index(l),index(l+1))
          i=l+1
          j=r
          indext=index(l+1)
          a=arr(indext)
          do
           do
            i=i+1
            if (arr(index(i)) >= a) exit
           end do
           do
            j=j-1
            if (arr(index(j)) <= a) exit
           end do
           if (j < i) exit
           call swap(index(i),index(j))
          end do
          index(l+1)=index(j)
          index(j)=indext
          jstack=jstack+2
          if (jstack > NSTACK) then
           write(6,*)'indexx_dp:NSTACK too small'
           STOP
          endif
          if (r-i+1 >= j-l) then
           istack(jstack)=r
           istack(jstack-1)=i
           r=j-1
          else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
          end if
         end if
        end do
        CONTAINS
        SUBROUTINE icomp_xchg(i,j)
        INTEGER(I4B), INTENT(INOUT) :: i,j
        INTEGER(I4B) :: swp
        if (arr(j) < arr(i)) then
        swp=i
        i=j
        j=swp
        end if
        END SUBROUTINE icomp_xchg
        END SUBROUTINE indexx_dp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE indexx_qp(arr,index)
        USE nrutil, ONLY :arth,assert_eq,nrerror,swap
        IMPLICIT NONE
        REAL(QP), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
        INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
!        Indexes an array arr, i.e., outputs the array index of length N such that arr(index(j))
!        is in ascending order for j = 1, 2, . . . ,N. The input quantity arr is not changed.
        REAL(QP) :: a
        INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
        INTEGER(I4B), DIMENSION(NSTACK) :: istack
        n=assert_eq(size(index),size(arr),'indexx_qp')
        index=arth(1,1,n)
        jstack=0
        l=1
        r=n
        do
         if (r-l < NN) then
          do j=l+1,r
           indext=index(j)
           a=arr(indext)
           do i=j-1,l,-1
            if (arr(index(i)) <= a) exit
            index(i+1)=index(i)
           end do
           index(i+1)=indext
          end do
          if (jstack == 0) RETURN
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
         else
          k=(l+r)/2
          call swap(index(k),index(l+1))
          call icomp_xchg(index(l),index(r))
          call icomp_xchg(index(l+1),index(r))
          call icomp_xchg(index(l),index(l+1))
          i=l+1
          j=r
          indext=index(l+1)
          a=arr(indext)
          do
           do
            i=i+1
            if (arr(index(i)) >= a) exit
           end do
           do
            j=j-1
            if (arr(index(j)) <= a) exit
           end do
           if (j < i) exit
           call swap(index(i),index(j))
          end do
          index(l+1)=index(j)
          index(j)=indext
          jstack=jstack+2
          if (jstack > NSTACK) then
           write(6,*)'indexx_qp:NSTACK too small'
           STOP
          endif
          if (r-i+1 >= j-l) then
           istack(jstack)=r
           istack(jstack-1)=i
           r=j-1
          else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
          end if
         end if
        end do
        CONTAINS
        SUBROUTINE icomp_xchg(i,j)
        INTEGER(I4B), INTENT(INOUT) :: i,j
        INTEGER(I4B) :: swp
        if (arr(j) < arr(i)) then
        swp=i
        i=j
        j=swp
        end if
        END SUBROUTINE icomp_xchg
        END SUBROUTINE indexx_qp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE laguer_d(a,x,its)
        USE nrutil, ONLY :poly,poly_term
        INTEGER(I4B), INTENT(OUT) :: its
        COMPLEX(DPC), INTENT(INOUT) :: x
        COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: a
        REAL(DP), PARAMETER :: EPS=epsilon(1.0_dp)
        INTEGER(I4B), PARAMETER :: MR=8,MT=10,MAXIT=MT*MR
!        Given an array of M + 1 complex coefficients a of the polynomial
!        and
!        given a complex value x, this routine improves x by Laguerre's method until it converges,
!        within the achievable roundoff limit, to a root of the given polynomial. The number of
!        iterations taken is returned as its.
!        Parameters: EPS is the estimated fractional roundoff error. We try to break (rare) limit
!        cycles with MR different fractional values, once every MT steps, for MAXIT total allowed
!        iterations.
        INTEGER(I4B) :: iter,m
        REAL(DP) :: abx,abp,abm,err
        COMPLEX(DPC) :: dx,x1,f,g,h,sq,gp,gm,g2
        COMPLEX(DPC), DIMENSION(size(a)) :: b,d
        REAL(DP), DIMENSION(MR) :: frac = &
        (/ 0.5_dp,0.25_dp,0.75_dp,0.13_dp,0.38_dp,0.62_dp,0.88_dp,1.0_dp /)
!        Fractions used to break a limit cycle.
        m=size(a)-1
        do iter=1,MAXIT !Loop over iterations up to allowed maximum.
         its=iter
         abx=abs(x)
         b(m+1:1:-1)=poly_term(a(m+1:1:-1),x) !Efficient computation of the polynomial
         d(m:1:-1)=poly_term(b(m+1:2:-1),x) ! and its first two derivatives.
         f=poly(x,d(2:m))
         err=EPS*poly(abx,abs(b(1:m+1))) !Esimate of roundoff in evaluating polynomial.
         if (abs(b(1)) <= err) RETURN !We are on the root.
         g=d(1)/b(1) !The generic case: Use Laguerre's formula.
         g2=g*g
         h=g2-2.0_dp*f/b(1)
         sq=sqrt((m-1)*(m*h-g2))
         gp=g+sq
         gm=g-sq
         abp=abs(gp)
         abm=abs(gm)
         if (abp < abm) gp=gm
         if (max(abp,abm) > 0.0) then
          dx=m/gp
         else
          dx=exp(cmplx(log(1.0_dp+abx),iter,kind=dpc))
         end if
         x1=x-dx
         if (x == x1) RETURN !Converged.
         if (mod(iter,MT) /= 0) then
          x=x1
         else !Every so often we take a fractional step, to break any limit cycle (itself a rare occurrence).
          x=x-dx*frac(iter/MT)
         end if
        end do
        write(6,*) 'Nbr d itérations dépassé dans laguer_d'
        STOP
!        Very unusual — can occur only for complex roots. Try a different starting guess for the root.
        END SUBROUTINE laguer_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE laguer_q(a,x,its)
        USE nrutil, ONLY :poly,poly_term
        INTEGER(I4B), INTENT(OUT) :: its
        COMPLEX(QPC), INTENT(INOUT) :: x
        COMPLEX(QPC), DIMENSION(:), INTENT(IN) :: a
        REAL(QP), PARAMETER :: EPS=epsilon(1.0_qp)
        INTEGER(I4B), PARAMETER :: MR=8,MT=10,MAXIT=MT*MR
!        Given an array of M + 1 complex coefficients a of the polynomial
!        and
!        given a complex value x, this routine improves x by Laguerre's method until it converges,
!        within the achievable roundoff limit, to a root of the given polynomial. The number of
!        iterations taken is returned as its.
!        Parameters: EPS is the estimated fractional roundoff error. We try to break (rare) limit
!        cycles with MR different fractional values, once every MT steps, for MAXIT total allowed
!        iterations.
        INTEGER(I4B) :: iter,m
        REAL(QP) :: abx,abp,abm,err
        COMPLEX(QPC) :: dx,x1,f,g,h,sq,gp,gm,g2
        COMPLEX(QPC), DIMENSION(size(a)) :: b,d
        REAL(QP), DIMENSION(MR) :: frac = &
        (/ 0.5_qp,0.25_qp,0.75_qp,0.13_qp,0.38_qp,0.62_qp,0.88_qp,1.0_qp /)
!        Fractions used to break a limit cycle.
        m=size(a)-1
        do iter=1,MAXIT !Loop over iterations up to allowed maximum.
         its=iter
         abx=abs(x)
         b(m+1:1:-1)=poly_term(a(m+1:1:-1),x) !Efficient computation of the polynomial
         d(m:1:-1)=poly_term(b(m+1:2:-1),x) ! and its first two derivatives.
         f=poly(x,d(2:m))
         err=EPS*poly(abx,abs(b(1:m+1))) !Esimate of roundoff in evaluating polynomial.
         if (abs(b(1)) <= err) RETURN !We are on the root.
         g=d(1)/b(1) !The generic case: Use Laguerre's formula.
         g2=g*g
         h=g2-2.0_qp*f/b(1)
         sq=sqrt((m-1)*(m*h-g2))
         gp=g+sq
         gm=g-sq
         abp=abs(gp)
         abm=abs(gm)
         if (abp < abm) gp=gm
         if (max(abp,abm) > 0.0) then
          dx=m/gp
         else
          dx=exp(cmplx(log(1.0_qp+abx),iter,kind=qpc))
         end if
         x1=x-dx
         if (x == x1) RETURN !Converged.
         if (mod(iter,MT) /= 0) then
          x=x1
         else !Every so often we take a fractional step, to break any limit cycle (itself a rare occurrence).
          x=x-dx*frac(iter/MT)
         end if
        end do
        write(6,*) 'Nbr d itérations dépassé dans laguer_d'
        STOP
!        Very unusual — can occur only for complex roots. Try a different starting guess for the root.
        END SUBROUTINE laguer_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE zroots_d(a,roots,polish)
        USE nrutil, ONLY :assert_eq, poly_term
        COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: a
        COMPLEX(DPC), DIMENSION(:), INTENT(OUT) :: roots
        LOGICAL(LGT), INTENT(IN) :: polish
        REAL(DP), PARAMETER :: EPS=1.0e-6_dp
!        Given the array of M + 1 complex coefficients a of the polynomial
!        this
!        routine successively calls laguer and finds all M complex roots. The logical variable
!        polish should be input as .true. if polishing (also by Laguerre's method) is desired,
!       .false. if the roots will be subsequently polished by other means.
!        Parameter: EPS is a small number.
        INTEGER(I4B) :: j,its,m
        INTEGER(I4B), DIMENSION(size(roots)) :: indx
        COMPLEX(DPC) :: x
        COMPLEX(DPC), DIMENSION(size(a)) :: ad
        m=assert_eq(size(roots),size(a)-1,'zroots')
        ad(:)=a(:) !Copy of coefficients for successive deflation.
        do j=m,1,-1 !Loop over each root to be found.
         x=cmplx(0.0_dp,kind=dpc)
!        Start at zero to favor convergence to smallest remaining root.
         call laguer(ad(1:j+1),x,its) !Find the root.
         if (abs(aimag(x)) <= 2.0_dp*EPS**2*abs(real(x))) &
         x=cmplx(real(x),kind=dpc)
         roots(j)=x
         ad(j:1:-1)=poly_term(ad(j+1:2:-1),x) !Forward deflation.
        end do
        if (polish) then
         do j=1,m !Polish the roots using the undeflated coefficients.
          call laguer(a(:),roots(j),its)
         end do
        end if
        call indexx(real(roots),indx) !Sort roots by their real parts.
        roots=roots(indx)
        END SUBROUTINE zroots_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE zroots_q(a,roots,polish)
        USE nrutil, ONLY :assert_eq, poly_term
        COMPLEX(QPC), DIMENSION(:), INTENT(IN) :: a
        COMPLEX(QPC), DIMENSION(:), INTENT(OUT) :: roots
        LOGICAL(LGT), INTENT(IN) :: polish
        REAL(QP), PARAMETER :: EPS=1.0e-14_qp
!        Given the array of M + 1 complex coefficients a of the polynomial
!        this
!        routine successively calls laguer and finds all M complex roots. The logical variable
!        polish should be input as .true. if polishing (also by Laguerre's method) is desired,
!       .false. if the roots will be subsequently polished by other means.
!        Parameter: EPS is a small number.
        INTEGER(I4B) :: j,its,m
        INTEGER(I4B), DIMENSION(size(roots)) :: indx
        COMPLEX(QPC) :: x
        COMPLEX(QPC), DIMENSION(size(a)) :: ad
        m=assert_eq(size(roots),size(a)-1,'zroots')
        ad(:)=a(:) !Copy of coefficients for successive deflation.
        do j=m,1,-1 !Loop over each root to be found.
         x=cmplx(0.0_qp,kind=qpc)
!        Start at zero to favor convergence to smallest remaining root.
         call laguer(ad(1:j+1),x,its) !Find the root.
         if (abs(aimag(x)) <= 2.0_qp*EPS**2*abs(real(x))) &
         x=cmplx(real(x),kind=qpc)
         roots(j)=x
         ad(j:1:-1)=poly_term(ad(j+1:2:-1),x) !Forward deflation.
        end do
        if (polish) then
         do j=1,m !Polish the roots using the undeflated coefficallcients.
          call laguer(a(:),roots(j),its)
         end do
        end if
        call indexx(aimag(roots),indx) !Sort roots by their real parts.
        roots=roots(indx)
        END SUBROUTINE zroots_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE ludcmp_s(a,indx,d)
        USE nrutil, ONLY : assert_eq,imaxloc,outerprod,swap
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
        REAL(SP), INTENT(OUT) :: d
!        Given an N × N input matrix a, this routine replaces it by the LU decomposition of a
!        rowwise permutation of itself. On output, a is arranged as in equation (2.3.14); indx is an
!        output vector of length N that records the row permutation effected by the partial pivoting;
!        d is output as ±1 depending on whether the number of row interchanges was even or odd,
!        respectively. This routine is used in combination with lubksb to solve linear equations or
!        invert a matrix.
        REAL(SP), DIMENSION(size(a,1)) :: vv !vv stores the implicit scaling of each row.
        REAL(SP), PARAMETER :: TINY=1.0e-20_sp !A small number.
        INTEGER(I4B) :: j,n,imax
        n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
        d=1.0 !No row interchanges yet.
        vv=maxval(abs(a),dim=2) !Loop over rows to get the implicit scaling information.
        if (any(vv == 0.0))then 
         write(*,*) 'singular matrix in ludcmp'
         STOP
        endif
!        There is a row of zeros.
        vv=1.0_sp/vv !Save the scaling.
        do j=1,n
        imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j))) !Find the pivot row.
        if (j /= imax) then !Do we need to interchange rows?
        call swap(a(imax,:),a(j,:)) !Yes, do so...
        d=-d !...and change the parity of d.
        vv(imax)=vv(j) !Also interchange the scale factor.
        end if
        indx(j)=imax
        if (a(j,j) == 0.0) a(j,j)=TINY
!        If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
!        For some applications on singular matrices, it is desirable to substitute TINY
!        for zero.
        a(j+1:n,j)=a(j+1:n,j)/a(j,j) !Divide by the pivot element.
        a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
!        Reduce remaining submatrix.
        end do
        END SUBROUTINE ludcmp_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE lubksb_s(a,indx,b)
        USE nrutil, ONLY : assert_eq
        REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
!        Solves the set of N linear equations A · X = B. Here the N × N matrix a is input, not
!        as the original matrix A, but rather as its LU decomposition, determined by the routine
!        ludcmp. indx is input as the permutation vector of length N returned by ludcmp. b is
!        input as the right-hand-side vector B, also of length N, and returns with the solution vector
!        X. a and indx are not modified by this routine and can be left in place for successive calls
!        with different right-hand sides b. This routine takes into account the possibility that b will
!        begin with many zero elements, so it is efficient for use in matrix inversion.
        INTEGER(I4B) :: i,n,iii,ll
        REAL(SP) :: summ
        n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
        iii=0 !When iii is set to a positive value, it will become the index
!        of the first nonvanishing element of b. We now do
!        the forward substitution, equation (2.3.6). The only new
!        wrinkle is to unscramble the permutation as we go.
        do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
        if (iii /= 0) then
        summ=summ-dot_product(a(i,iii:i-1),b(iii:i-1))
        else if (summ /= 0.0) then
        iii=i !A nonzero element was encountered, so from now on we will
        end if !have to do the dot product above.
        b(i)=summ
        end do
        do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
        b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
        end do
        END SUBROUTINE lubksb_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE ludcmp_d(a,indx,d)
        USE nrutil, ONLY : assert_eq,imaxloc,outerprod,swap
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
        REAL(DP), INTENT(OUT) :: d
!        Given an N × N input matrix a, this routine replaces it by the LU decomposition of a
!        rowwise permutation of itself. On output, a is arranged as in equation (2.3.14); indx is an
!        output vector of length N that records the row permutation effected by the partial pivoting;
!        d is output as ±1 depending on whether the number of row interchanges was even or odd,
!        respectively. This routine is used in combination with lubksb to solve linear equations or
!        invert a matrix.
        REAL(DP), DIMENSION(size(a,1)) :: vv !vv stores the implicit scaling of each row.
        REAL(DP), PARAMETER :: TINY=1.0e-20_dp !A small number.
        INTEGER(I4B) :: j,n,imax
        n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
        d=1.0 !No row interchanges yet.
        vv=maxval(abs(a),dim=2) !Loop over rows to get the implicit scaling information.
        if (any(vv == 0.0))then 
         write(*,*) 'singular matrix in ludcmp'
         STOP
        endif
!        There is a row of zeros.
        vv=1.0_dp/vv !Save the scaling.
        do j=1,n
        imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j))) !Find the pivot row.
        if (j /= imax) then !Do we need to interchange rows?
        call swap(a(imax,:),a(j,:)) !Yes, do so...
        d=-d !...and change the parity of d.
        vv(imax)=vv(j) !Also interchange the scale factor.
        end if
        indx(j)=imax
        if (a(j,j) == 0.0) a(j,j)=TINY
!        If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
!        For some applications on singular matrices, it is desirable to substitute TINY
!        for zero.
        a(j+1:n,j)=a(j+1:n,j)/a(j,j) !Divide by the pivot element.
        a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
!        Reduce remaining submatrix.
        end do
        END SUBROUTINE ludcmp_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE lubksb_d(a,indx,b)
        USE nrutil, ONLY : assert_eq
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
!        Solves the set of N linear equations A · X = B. Here the N × N matrix a is input, not
!        as the original matrix A, but rather as its LU decomposition, determined by the routine
!        ludcmp. indx is input as the permutation vector of length N returned by ludcmp. b is
!        input as the right-hand-side vector B, also of length N, and returns with the solution vector
!        X. a and indx are not modified by this routine and can be left in place for successive calls
!        with different right-hand sides b. This routine takes into account the possibility that b will
!        begin with many zero elements, so it is efficient for use in matrix inversion.
        INTEGER(I4B) :: i,n,iii,ll
        REAL(DP) :: summ
        n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
        iii=0 !When iii is set to a positive value, it will become the index
!        of the first nonvanishing element of b. We now do
!        the forward substitution, equation (2.3.6). The only new
!        wrinkle is to unscramble the permutation as we go.
        do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
        if (iii /= 0) then
        summ=summ-dot_product(a(i,iii:i-1),b(iii:i-1))
        else if (summ /= 0.0) then
        iii=i !A nonzero element was encountered, so from now on we will
        end if !have to do the dot product above.
        b(i)=summ
        end do
        do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
        b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
        end do
        END SUBROUTINE lubksb_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE ludcmp_q(a,indx,d)
        USE nrutil, ONLY : assert_eq,imaxloc,outerprod,swap
        REAL(QP), DIMENSION(:,:), INTENT(INOUT) :: a
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
        REAL(QP), INTENT(OUT) :: d
!        Given an N × N input matrix a, this routine replaces it by the LU decomposition of a
!        rowwise permutation of itself. On output, a is arranged as in equation (2.3.14); indx is an
!        output vector of length N that records the row permutation effected by the partial pivoting;
!        d is output as ±1 depending on whether the number of row interchanges was even or odd,
!        respectively. This routine is used in combination with lubksb to solve linear equations or
!        invert a matrix.
        REAL(QP), DIMENSION(size(a,1)) :: vv !vv stores the implicit scaling of each row.
        REAL(QP), PARAMETER :: TINY=1.0e-20_qp !A small number.
        INTEGER(I4B) :: j,n,imax
        n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
        d=1.0 !No row interchanges yet.
        vv=maxval(abs(a),dim=2) !Loop over rows to get the implicit scaling information.
        if (any(vv == 0.0))then
         write(*,*) 'singular matrix in ludcmp'
         STOP
        endif
!        There is a row of zeros.
        vv=1.0_qp/vv !Save the scaling.
        do j=1,n
        imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j))) !Find the pivot row.
        if (j /= imax) then !Do we need to interchange rows?
        call swap(a(imax,:),a(j,:)) !Yes, do so...
        d=-d !...and change the parity of d.
        vv(imax)=vv(j) !Also interchange the scale factor.
        end if
        indx(j)=imax
        if (a(j,j) == 0.0) a(j,j)=TINY
!        If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
!        For some applications on singular matrices, it is desirable to substitute TINY
!        for zero.
        a(j+1:n,j)=a(j+1:n,j)/a(j,j) !Divide by the pivot element.
        a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
!        Reduce remaining submatrix.
        end do
        END SUBROUTINE ludcmp_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE lubksb_q(a,indx,b)
        USE nrutil, ONLY : assert_eq
        REAL(QP), DIMENSION(:,:), INTENT(IN) :: a
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
        REAL(QP), DIMENSION(:), INTENT(INOUT) :: b
!        Solves the set of N linear equations A · X = B. Here the N × N matrix a is input, not
!        as the original matrix A, but rather as its LU decomposition, determined by the routine
!        ludcmp. indx is input as the permutation vector of length N returned by ludcmp. b is
!        input as the right-hand-side vector B, also of length N, and returns with the solution vector
!        X. a and indx are not modified by this routine and can be left in place for successive calls
!        with different right-hand sides b. This routine takes into account the possibility that b will
!        begin with many zero elements, so it is efficient for use in matrix inversion.
        INTEGER(I4B) :: i,n,iii,ll
        REAL(QP) :: summ
        n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
        iii=0 !When iii is set to a positive value, it will become the index
!        of the first nonvanishing element of b. We now do
!        the forward substitution, equation (2.3.6). The only new
!        wrinkle is to unscramble the permutation as we go.
        do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
        if (iii /= 0) then
        summ=summ-dot_product(a(i,iii:i-1),b(iii:i-1))
        else if (summ /= 0.0) then
        iii=i !A nonzero element was encountered, so from now on we will
        end if !have to do the dot product above.
        b(i)=summ
        end do
        do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
        b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
        end do
        END SUBROUTINE lubksb_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE ludcmp_cq(a,indx,d)
        USE nrutil, ONLY : assert_eq,imaxloc,outerprod,swap
        COMPLEX(QPC), DIMENSION(:,:), INTENT(INOUT) :: a
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
        REAL(QP), INTENT(OUT) :: d
!        Given an N × N input matrix a, this routine replaces it by the LU decomposition of a
!        rowwise permutation of itself. On output, a is arranged as in equation (2.3.14); indx is an
!        output vector of length N that records the row permutation effected by the partial pivoting;
!        d is output as ±1 depending on whether the number of row interchanges was even or odd,
!        respectively. This routine is used in combination with lubksb to solve linear equations or
!        invert a matrix.
        REAL(QP), DIMENSION(size(a,1)) :: vv !vv stores the implicit scaling of each row.
        REAL(QP), PARAMETER :: TINY=1.0e-80_qp !A small number.
        INTEGER(I4B) :: j,n,imax
        n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
        d=1.0 !No row interchanges yet.
        vv=maxval(abs(a),dim=2) !Loop over rows to get the implicit scaling information.
        if (any(vv == 0.0))then
         write(*,*) 'singular matrix in ludcmp'
         STOP
        endif
!        There is a row of zeros.
        vv=1.0_qp/vv !Save the scaling.
        do j=1,n
         imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j))) !Find the pivot row.
         if (j /= imax) then !Do we need to interchange rows?
          call swap(a(imax,:),a(j,:)) !Yes, do so...
          d=-d !...and change the parity of d.
          vv(imax)=vv(j) !Also interchange the scale factor.
         end if
         indx(j)=imax
         if (a(j,j) == 0.0) a(j,j)=TINY
!        If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
!        For some applications on singular matrices, it is desirable to substitute TINY
!        for zero.
         a(j+1:n,j)=a(j+1:n,j)/a(j,j) !Divide by the pivot element.
         a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
!        Reduce remaining submatrix.
        end do
        END SUBROUTINE ludcmp_cq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rd_s(x,y,z)
        REAL(SP), INTENT(IN) :: x,y,z
        REAL(SP) :: rd_s
        REAL(SP), PARAMETER :: ERRTOL=0.0015_sp,TINY=1.0e-38_sp,BIG=4.5e37_sp,&
        C1=3.0_sp/14.0_sp,C2=1.0_sp/6.0_sp,C3=9.0_sp/22.0_sp,&
        C4=3.0_sp/26.0_sp,C5=0.25_sp*C3,C6=1.5_sp*C4
        !Computes Carlson's elliptic integral of the second kind, RD(x, y, z). x and y must be
        !nonnegative, and at most one can be zero. z must be positive. TINY must be at least twice
        !the negative 2/3 power of the machine overflow limit. BIG must be at most 0.1×ERRTOL
        !times the negative 2/3 power of the machine underflow limit.
        REAL(SP) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
        ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
        call assert(min(x,y) >= 0.0, min(x+y,z) >= TINY, max(x,y,z) <= BIG, &
        'rd_s args')
        xt=x
        yt=y
        zt=z
        sum=0.0
        fac=1.0
        do
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=0.25_sp*fac
        xt=0.25_sp*(xt+alamb)
        yt=0.25_sp*(yt+alamb)
        zt=0.25_sp*(zt+alamb)
        ave=0.2_sp*(xt+yt+3.0_sp*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        ea=delx*dely
        eb=delz*delz
        ec=ea-eb
        ed=ea-6.0_sp*eb
        ee=ed+ec+ec
        rd_s=3.0_sp*sum+fac*(1.0_sp+ed*(-C1+C5*ed-C6*delz*ee)&
        +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
        END FUNCTION rd_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rd_d(x,y,z)
        REAL(DP), INTENT(IN) :: x,y,z
        REAL(DP) :: rd_d
        REAL(DP), PARAMETER :: ERRTOL=0.0015_dp,TINY=1.0e-200_dp,BIG=4.5e200_dp,&
        C1=3.0_dp/14.0_dp,C2=1.0_dp/6.0_dp,C3=9.0_dp/22.0_dp,&
        C4=3.0_dp/26.0_dp,C5=0.25_dp*C3,C6=1.5_dp*C4
        !Computes Carlson's elliptic integral of the second kind, RD(x, y, z). x and y must be
        !nonnegative, and at most one can be zero. z must be positive. TINY must be at least twice
        !the negative 2/3 power of the machine overflow limit. BIG must be at most 0.1×ERRTOL
        !times the negative 2/3 power of the machine underflow limit.
        REAL(DP) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
        ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
        call assert(min(x,y) >= 0.0, min(x+y,z) >= TINY, max(x,y,z) <= BIG, &
        'rd_d args')
        xt=x
        yt=y
        zt=z
        sum=0.0d0
        fac=1.0d0
        do
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=0.25_dp*fac
        xt=0.25_dp*(xt+alamb)
        yt=0.25_dp*(yt+alamb)
        zt=0.25_dp*(zt+alamb)
        ave=0.2_dp*(xt+yt+3.0_dp*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        ea=delx*dely
        eb=delz*delz
        ec=ea-eb
        ed=ea-6.0_dp*eb
        ee=ed+ec+ec
        rd_d=3.0_dp*sum+fac*(1.0_dp+ed*(-C1+C5*ed-C6*delz*ee)&
        +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
        END FUNCTION rd_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rd_q(x,y,z)
        REAL(QP), INTENT(IN) :: x,y,z
        REAL(QP) :: rd_q
        REAL(QP), PARAMETER :: ERRTOL=0.000025_qp,TINY=1.0e-200_qp,BIG=4.5e200_qp,&
        C1=3.0_qp/14.0_qp,C2=1.0_qp/6.0_qp,C3=9.0_qp/22.0_qp,&
        C4=3.0_qp/26.0_qp,C5=0.25_qp*C3,C6=1.5_qp*C4
        !Computes Carlson's elliptic integral of the second kind, RD(x, y, z). x and y must be
        !nonnegative, and at most one can be zero. z must be positive. TINY must be at least twice
        !the negative 2/3 power of the machine overflow limit. BIG must be at most 0.1×ERRTOL
        !times the negative 2/3 power of the machine underflow limit.
        REAL(QP) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
        ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
        call assert(min(x,y) >= 0.0, min(x+y,z) >= TINY, max(x,y,z) <= BIG, &
        'rd_q args')
        xt=x
        yt=y
        zt=z
        sum=0.0d0
        fac=1.0d0
        do
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=0.25_qp*fac
        xt=0.25_qp*(xt+alamb)
        yt=0.25_qp*(yt+alamb)
        zt=0.25_qp*(zt+alamb)
        ave=0.2_qp*(xt+yt+3.0_qp*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        ea=delx*dely
        eb=delz*delz
        ec=ea-eb
        ed=ea-6.0_qp*eb
        ee=ed+ec+ec
        rd_q=3.0_qp*sum+fac*(1.0_qp+ed*(-C1+C5*ed-C6*delz*ee)&
        +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
        END FUNCTION rd_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rd_c(x,y,z)
        COMPLEX*16, INTENT(IN) :: x,y,z
        COMPLEX*16 :: rd_c
        REAL(DP), PARAMETER :: ERRTOL=0.0015_dp,TINY=1.0e-200_dp,BIG=4.5e200_dp,&
        C1=3.0_dp/14.0_dp,C2=1.0_dp/6.0_dp,C3=9.0_dp/22.0_dp,&
        C4=3.0_dp/26.0_dp,C5=0.25_dp*C3,C6=1.5_dp*C4
        !Computes Carlson's elliptic integral of the second kind, RD(x, y, z). x and y must be
        !nonnegative (ligne de coupure ]-inf,0]), and at most one can be zero. z must be positive. TINY must be at least twice
        !the negative 2/3 power of the machine overflow limit. BIG must be at most 0.1×ERRTOL
        !times the negative 2/3 power of the machine underflow limit.
        COMPLEX*16 :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
        ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
        call assert((/.NOT.((real(x) <= 0.0).AND.(abs(imag(x))<TINY)), &
                      .NOT.((real(y) <= 0.0).AND.(abs(imag(y))<TINY)), &
                      .NOT.((real(z) <= 0.0).AND.(abs(imag(z))<TINY)), &
                      min(abs(x)+abs(y),abs(z)) >= TINY, &
                      max(abs(x),abs(y),abs(z)) <= BIG/),'rd_c args')
        xt=x
        yt=y
        zt=z
        sum=dcmplx(0.0d0,0.d0)
        fac=dcmplx(1.0d0,0.d0)
        do
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=0.25_dp*fac
        xt=0.25_dp*(xt+alamb)
        yt=0.25_dp*(yt+alamb)
        zt=0.25_dp*(zt+alamb)
        ave=0.2_dp*(xt+yt+3.0_dp*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        ea=delx*dely
        eb=delz*delz
        ec=ea-eb
        ed=ea-6.0_dp*eb
        ee=ed+ec+ec
        rd_c=3.0_dp*sum+fac*(1.0_dp+ed*(-C1+C5*ed-C6*delz*ee)&
        +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
        END FUNCTION rd_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rd_cq(x,y,z)
        COMPLEX(QPC), INTENT(IN) :: x,y,z
        COMPLEX(QPC) :: rd_cq
        REAL(QP), PARAMETER :: ERRTOL=0.00005_qp,TINY=1.0e-200_qp,BIG=4.5e200_qp,&
        C1=3.0_qp/14.0_qp,C2=1.0_qp/6.0_qp,C3=9.0_qp/22.0_qp,&
        C4=3.0_qp/26.0_qp,C5=0.25_qp*C3,C6=1.5_qp*C4
        !Computes Carlson's elliptic integral of the second kind, RD(x, y, z). x and y must be
        !nonnegative (ligne de coupure ]-inf,0]), and at most one can be zero. z must be positive. TINY must be at least twice
        !the negative 2/3 power of the machine overflow limit. BIG must be at most 0.1×ERRTOL
        !times the negative 2/3 power of the machine underflow limit.
        COMPLEX(QPC) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
        ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
        call assert((/.NOT.((real(x) <= 0.0).AND.(abs(imag(x))<TINY)), &
                      .NOT.((real(y) <= 0.0).AND.(abs(imag(y))<TINY)), &
                      .NOT.((real(z) <= 0.0).AND.(abs(imag(z))<TINY)), &
                      min(abs(x)+abs(y),abs(z)) >= TINY, &
                      max(abs(x),abs(y),abs(z)) <= BIG/),'rd_cq args')
        xt=x
        yt=y
        zt=z
        sum=cmplx(0.0,0.0,kind=QPC)
        fac=cmplx(1.0,0.0,kind=QPC)
        do
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=0.25_qp*fac
        xt=0.25_qp*(xt+alamb)
        yt=0.25_qp*(yt+alamb)
        zt=0.25_qp*(zt+alamb)
        ave=0.2_qp*(xt+yt+3.0_qp*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        ea=delx*dely
        eb=delz*delz
        ec=ea-eb
        ed=ea-6.0_qp*eb
        ee=ed+ec+ec
        rd_cq=3.0_qp*sum+fac*(1.0_qp+ed*(-C1+C5*ed-C6*delz*ee)&
        +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
        END FUNCTION rd_cq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rf_s(x,y,z)
        REAL(SP), INTENT(IN) :: x,y,z
        REAL(SP) :: rf_s
        REAL(SP), PARAMETER :: ERRTOL=0.08_sp,TINY=1.5e-38_sp,BIG=3.0e37_sp,&
        THIRD=1.0_sp/3.0_sp,&
        C1=1.0_sp/24.0_sp,C2=0.1_sp,C3=3.0_sp/44.0_sp,C4=1.0_sp/14.0_sp
        !Computes Carlson's elliptic integral of the first kind, RF (x, y, z). x, y, and z must be
        !nonnegative, and at most one can be zero. TINY must be at least 5 times the machine
        !underflow limit, BIG at most one-fifth the machine overflow limit.
        REAL(SP) :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
        call assert(min(x,y,z) >= 0.0, min(x+y,x+z,y+z) >= TINY, &
        max(x,y,z) <= BIG, 'rf_s args')
        xt=x
        yt=y
        zt=z
        do
         sqrtx=sqrt(xt)
         sqrty=sqrt(yt)
         sqrtz=sqrt(zt)
         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
         xt=0.25_sp*(xt+alamb)
         yt=0.25_sp*(yt+alamb)
         zt=0.25_sp*(zt+alamb)
         ave=THIRD*(xt+yt+zt)
         delx=(ave-xt)/ave
         dely=(ave-yt)/ave
         delz=(ave-zt)/ave
         if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        e2=delx*dely-delz**2
        e3=delx*dely*delz
        rf_s=(1.0_sp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
        END FUNCTION rf_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rf_d(x,y,z)
        REAL(DP), INTENT(IN) :: x,y,z
        REAL(DP) :: rf_d
        REAL(DP), PARAMETER :: ERRTOL=0.0025_dp,TINY=1.5e-100_dp,BIG=3.0e100_dp,&
        THIRD=1.0_dp/3.0_dp,&
        C1=1.0_dp/24.0_dp,C2=0.1_dp,C3=3.0_dp/44.0_dp,C4=1.0_dp/14.0_dp
        !Computes Carlson's elliptic integral of the first kind, RF (x, y, z). x, y, and z must be
        !nonnegative, and at most one can be zero. TINY must be at least 5 times the machine
        !underflow limit, BIG at most one-fifth the machine overflow limit.
        REAL(DP) :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
        if(.NOT.((min(x,y,z) >= 0.d0).AND.(min(x+y,x+z,y+z)>=TINY).AND.(max(x,y,z)<=BIG)))then
          write(6,*) 'erreur ds les arguments de rf_d'
          write(6,*) 'x,y,z=',x,y,z
          STOP
        endif
        xt=x
        yt=y
        zt=z
        do
         sqrtx=sqrt(xt)
         sqrty=sqrt(yt)
         sqrtz=sqrt(zt)
         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
         xt=0.25_dp*(xt+alamb)
         yt=0.25_dp*(yt+alamb)
         zt=0.25_dp*(zt+alamb)
         ave=THIRD*(xt+yt+zt)
         delx=(ave-xt)/ave
         dely=(ave-yt)/ave
         delz=(ave-zt)/ave
         if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        e2=delx*dely-delz**2
        e3=delx*dely*delz
        rf_d=(1.0_dp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
        END FUNCTION rf_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rf_q(x,y,z)
        REAL(QP), INTENT(IN) :: x,y,z
        REAL(QP) :: rf_q
        REAL(QP), PARAMETER :: ERRTOL=0.000025_qp,TINY=1.5e-100_qp,BIG=3.0e100_qp,&
        THIRD=1.0_qp/3.0_qp,&
        C1=1.0_qp/24.0_qp,C2=0.1_qp,C3=3.0_qp/44.0_qp,C4=1.0_qp/14.0_qp
        !Computes Carlson's elliptic integral of the first kind, RF (x, y, z). x, y, and z must be
        !nonnegative, and at most one can be zero. TINY must be at least 5 times the machine
        !underflow limit, BIG at most one-fifth the machine overflow limit.
        REAL(QP) :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
        if(.NOT.((min(x,y,z) >= 0.d0).AND.(min(x+y,x+z,y+z)>=TINY).AND.(max(x,y,z)<=BIG)))then
          write(6,*) 'erreur ds les arguments de rf_q'
          write(6,*) 'x,y,z=',x,y,z
          STOP
        endif
        xt=x
        yt=y
        zt=z
        do
         sqrtx=sqrt(xt)
         sqrty=sqrt(yt)
         sqrtz=sqrt(zt)
         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
         xt=0.25_qp*(xt+alamb)
         yt=0.25_qp*(yt+alamb)
         zt=0.25_qp*(zt+alamb)
         ave=THIRD*(xt+yt+zt)
         delx=(ave-xt)/ave
         dely=(ave-yt)/ave
         delz=(ave-zt)/ave
         if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        e2=delx*dely-delz**2
        e3=delx*dely*delz
        rf_q=(1.0_qp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
        END FUNCTION rf_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rf_c(x,y,z)
        COMPLEX*16, INTENT(IN) :: x,y,z
        COMPLEX*16 :: rf_c
        REAL(DP), PARAMETER :: ERRTOL=0.0025_dp,TINY=1.5e-100_dp,BIG=3.0e100_dp,&
        THIRD=1.0_dp/3.0_dp,&
        C1=1.0_dp/24.0_dp,C2=0.1_dp,C3=3.0_dp/44.0_dp,C4=1.0_dp/14.0_dp
        !Computes Carlson's elliptic integral of the first kind, RF (x, y, z). x, y, and z must be
        !nonnegative, and at most one can be zero. TINY must be at least 5 times the machine
        !underflow limit, BIG at most one-fifth the machine overflow limit.
        COMPLEX*16 :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
        call assert((/.NOT.((real(x) <= 0.0).AND.(abs(imag(x))<TINY)), &
                      .NOT.((real(y) <= 0.0).AND.(abs(imag(y))<TINY)), &
                      .NOT.((real(z) <= 0.0).AND.(abs(imag(z))<TINY)), &
                      min(abs(x)+abs(y),abs(x)+abs(z),abs(y)+abs(z)) >= TINY, &
                      max(abs(x),abs(y),abs(z)) <= BIG/),'rf_c args')
        xt=x
        yt=y
        zt=z
        do
         sqrtx=sqrt(xt)
         sqrty=sqrt(yt)
         sqrtz=sqrt(zt)
         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
         xt=0.25_dp*(xt+alamb)
         yt=0.25_dp*(yt+alamb)
         zt=0.25_dp*(zt+alamb)
         ave=THIRD*(xt+yt+zt)
         delx=(ave-xt)/ave
         dely=(ave-yt)/ave
         delz=(ave-zt)/ave
         if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        e2=delx*dely-delz**2
        e3=delx*dely*delz
        rf_c=(1.0_dp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
        END FUNCTION rf_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rf_cq(x,y,z)
        COMPLEX(QPC), INTENT(IN) :: x,y,z
        COMPLEX(QPC) :: rf_cq
        REAL(QP), PARAMETER :: ERRTOL=0.000025_qp,TINY=1.5e-100_qp,BIG=3.0e100_qp,&
        THIRD=1.0_qp/3.0_qp,&
        C1=1.0_qp/24.0_qp,C2=0.1_qp,C3=3.0_qp/44.0_qp,C4=1.0_qp/14.0_qp
        !Computes Carlson's elliptic integral of the first kind, RF (x, y, z). x, y, and z must be
        !nonnegative, and at most one can be zero. TINY must be at least 5 times the machine
        !underflow limit, BIG at most one-fifth the machine overflow limit.
        COMPLEX(QPC) :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
        call assert((/.NOT.((real(x) <= 0.0).AND.(abs(imag(x))<TINY)), &
                      .NOT.((real(y) <= 0.0).AND.(abs(imag(y))<TINY)), &
                      .NOT.((real(z) <= 0.0).AND.(abs(imag(z))<TINY)), &
                      min(abs(x)+abs(y),abs(x)+abs(z),abs(y)+abs(z)) >= TINY, &
                      max(abs(x),abs(y),abs(z)) <= BIG/),'rf_cq args')
        xt=x
        yt=y
        zt=z
        do
         sqrtx=sqrt(xt)
         sqrty=sqrt(yt)
         sqrtz=sqrt(zt)
         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
         xt=0.25_qp*(xt+alamb)
         yt=0.25_qp*(yt+alamb)
         zt=0.25_qp*(zt+alamb)
         ave=THIRD*(xt+yt+zt)
         delx=(ave-xt)/ave
         dely=(ave-yt)/ave
         delz=(ave-zt)/ave
         if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        e2=delx*dely-delz**2
        e3=delx*dely*delz
        rf_cq=(1.0_qp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
        END FUNCTION rf_cq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rj_q(x,y,z,p)
        USE nrutil, ONLY : assert
        REAL(QP), INTENT(IN) :: x,y,z,p
        REAL(QP) :: rj_q
        REAL(QP), PARAMETER :: ERRTOL=0.0005_qp,TINY=2.5e-33_qp,BIG=9.0e41_qp,&
        C1=3.0_qp/14.0_qp,C2=1.0_qp/3.0_qp,C3=3.0_qp/22.0_qp,&
        C4=3.0_qp/26.0_qp,C5=0.75_qp*C3,C6=1.5_qp*C4,C7=0.5_qp*C2,&
        C8=C3+C3
        !Computes Carlson’s elliptic integral of the third kind, RJ(x, y, z, p). x, y, and z must be
        !nonnegative, and at most one can be zero. p must be nonzero. If p < 0, the Cauchy
        !principal value is returned. TINY must be at least twice the cube root of the machine
        !underflow limit, BIG at most one-fifth the cube root of the machine overflow limit.
        REAL(QP) :: a,alamb,alpha,ave,b,bet,delp,delx,&
        dely,delz,ea,eb,ec,ed,ee,fac,pt,rho,sqrtx,sqrty,sqrtz,&
        sm,tau,xt,yt,zt
        call assert(min(x,y,z) >= 0.0, min(x+y,x+z,y+z,abs(p)) >= TINY, &
         max(x,y,z,abs(p)) <= BIG, "rj_q args")
        sm=0.0
        fac=1.0
        if (p > 0.0) then
         xt=x
         yt=y
         zt=z
         pt=p
        else
         xt=min(x,y,z)
         zt=max(x,y,z)
         yt=x+y+z-xt-zt
         a=1.0_qp/(yt-p)
         b=a*(zt-yt)*(yt-xt)
         pt=yt+b
         rho=xt*zt/yt
         tau=p*pt/yt
        end if
        do
         sqrtx=sqrt(xt)
         sqrty=sqrt(yt)
         sqrtz=sqrt(zt)
         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
         alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
         bet=pt*(pt+alamb)**2
         sm=sm+fac*rc(alpha,bet)
         fac=0.25_qp*fac
         xt=0.25_qp*(xt+alamb)
         yt=0.25_qp*(yt+alamb)
         zt=0.25_qp*(zt+alamb)
         pt=0.25_qp*(pt+alamb)
         ave=0.2_qp*(xt+yt+zt+pt+pt)
         delx=(ave-xt)/ave
         dely=(ave-yt)/ave
         delz=(ave-zt)/ave
         delp=(ave-pt)/ave
         if (max(abs(delx),abs(dely),abs(delz),abs(delp)) <= ERRTOL) exit
        end do
        ea=delx*(dely+delz)+dely*delz
        eb=delx*dely*delz
        ec=delp**2
        ed=ea-3.0_qp*ec
        ee=eb+2.0_qp*delp*(ea-ec)
        rj_q=3.0_qp*sm+fac*(1.0_qp+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8&
        +delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
        if (p <= 0.0) rj_q=a*(b*rj_q+3.0_qp*(rc(rho,tau)-rf(xt,yt,zt)))
        END FUNCTION rj_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rc_q(x,y)
        USE nrutil, ONLY : assert
        REAL(QP), INTENT(IN) :: x,y
        REAL(QP) :: rc_q
        REAL(QP), PARAMETER :: ERRTOL=0.0004_qp,TINY=1.69e-68_qp,&
        SQRTNY=1.3e-34_qp,BIG=3.0e37_qp,TNBG=TINY*BIG,&
        COMP1=2.236_qp/SQRTNY,COMP2=TNBG*TNBG/25.0_qp,&
        THIRD=1.0_qp/3.0_qp,&
        C1=0.3_qp,C2=1.0_qp/7.0_qp,C3=0.375_qp,C4=9.0_qp/22.0_qp
!        Computes Carlson’s degenerate elliptic integral, RC(x, y). x must be nonnegative and y
!        must be nonzero. If y < 0, the Cauchy principal value is returned. TINY must be at least
!        5 times the machine underflow limit, BIG at most one-fifth the machine maximum overflow
!        limit.
        REAL(QP) :: alamb,ave,s,w,xt,yt
        call assert( (/x >= 0.0,y /= 0.0,x+abs(y) >= TINY,x+abs(y) <= BIG, &
         y >= -COMP1 .or. x <= 0.0 .or. x >= COMP2/),"rc_q")
        if (y > 0.0) then
         xt=x
         yt=y
         w=1.0
        else
         xt=x-y
         yt=-y
         w=sqrt(x)/sqrt(xt)
        end if
        do
         alamb=2.0_qp*sqrt(xt)*sqrt(yt)+yt
         xt=0.25_qp*(xt+alamb)
         yt=0.25_qp*(yt+alamb)
         ave=THIRD*(xt+yt+yt)
         s=(yt-ave)/ave
         if (abs(s) <= ERRTOL) exit
        end do
        rc_q=w*(1.0_qp+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
        END FUNCTION rc_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION elle_s(phi,ak)
        REAL(SP), INTENT(IN) :: phi,ak
        REAL(SP) :: elle_s
        !Legendre elliptic integral of the 2nd kind E(φ, k), evaluated using Carlson's functions RD
        !and RF . The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
        REAL(SP) :: cc,q,s
        s=sin(phi)
        cc=cos(phi)**2
        q=(1.0_sp-s*ak)*(1.0_sp+s*ak)
        if(isnan(cc).OR.isnan(q))then
         write(6,*)'Arguments dans elle, phi,ak=',phi,ak
         write(6,*)'cc,q=',cc,q
         STOP
        endif
        elle_s=s*(rf(cc,q,1.0_sp)-((s*ak)**2)*rd(cc,q,1.0_sp)/3.0_sp)
        END FUNCTION elle_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION elle_d(phi,ak)
        REAL(DP), INTENT(IN) :: phi,ak
        REAL(DP) :: elle_d
        !Legendre elliptic integral of the 2nd kind E(φ, k), evaluated using Carlson's functions RD
        !and RF . The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
        REAL(DP) :: cc,q,s
        s=sin(phi)
        cc=cos(phi)**2
        q=(1.0_dp-s*ak)*(1.0_dp+s*ak)
        if(isnan(cc).OR.isnan(q))then
         write(6,*)'Arguments dans elle, phi,ak=',phi,ak
         write(6,*)'cc,q=',cc,q
         STOP
        endif
        elle_d=s*(rf(cc,q,1.0_dp)-((s*ak)**2)*rd(cc,q,1.0_dp)/3.0_dp)
        END FUNCTION elle_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION elle_q(phi,ak)
        REAL(QP), INTENT(IN) :: phi,ak
        REAL(QP) :: elle_q
        !Legendre elliptic integral of the 2nd kind E(φ, k), evaluated using Carlson's functions RD
        !and RF . The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
        REAL(QP) :: cc,q,s
        s=sin(phi)
        cc=cos(phi)**2
        q=(1.0_qp-s*ak)*(1.0_qp+s*ak)
        if(isnan(cc).OR.isnan(q))then
         write(6,*)'Arguments dans elle, phi,ak=',phi,ak
         write(6,*)'cc,q=',cc,q
         STOP
        endif
        elle_q=s*(rf(cc,q,1.0_qp)-((s*ak)**2)*rd(cc,q,1.0_qp)/3.0_qp)
        END FUNCTION elle_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION elle_c(phi,ak)
        COMPLEX(DPC), INTENT(IN) :: phi,ak
        COMPLEX(DPC) :: elle_c
        !Legendre elliptic integral of the 2nd kind E(φ, k), evaluated using Carlson's functions RD
        !and RF . The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
        COMPLEX(DPC) :: cc,q,s,un
        un=dcmplx(1.d0,0.d0)
        s=sin(phi)
        cc=cos(phi)**2
        q=(1.0_dp-s*ak)*(1.0_dp+s*ak)
        elle_c=s*(rf(cc,q,un)-((s*ak)**2)*rd(cc,q,un)/3.0_dp)
        END FUNCTION elle_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION elle_rc(phi,ak)
        COMPLEX(DPC), INTENT(IN) :: ak
        REAL(DP), INTENT(IN) :: phi
        COMPLEX(DPC) :: elle_rc
        !Legendre elliptic integral of the 2nd kind E(φ, k), evaluated using Carlson's functions RD
        !and RF . The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
        COMPLEX(DPC) :: cc,q,s,un
        un=dcmplx(1.d0,0.d0)
        s=sin(phi)
        cc=cos(phi)**2
        q=(1.0_dp-s*ak)*(1.0_dp+s*ak)
        elle_rc=s*(rf(cc,q,un)-((s*ak)**2)*rd(cc,q,un)/3.0_dp)
        END FUNCTION elle_rc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION elle_cq(phi,ak)
        COMPLEX(QPC), INTENT(IN) :: phi,ak
        COMPLEX(QPC) :: elle_cq
        !Legendre elliptic integral of the 2nd kind E(φ, k), evaluated using Carlson's functions RD
        !and RF . The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
        COMPLEX(QPC) :: cc,q,s,un
        un=cmplx(1.0_qp,0.0_qp)
        s=sin(phi)
        cc=cos(phi)**2
        q=(1.0_qp-s*ak)*(1.0_qp+s*ak)
        elle_cq=s*(rf(cc,q,un)-((s*ak)**2)*rd(cc,q,un)/3.0_qp)
        END FUNCTION elle_cq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION elle_rcq(phi,ak)
        COMPLEX(QPC), INTENT(IN) :: ak
        REAL(QP), INTENT(IN) :: phi
        COMPLEX(QPC) :: elle_rcq
        !Legendre elliptic integral of the 2nd kind E(φ, k), evaluated using Carlson's functions RD
        !and RF . The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
        COMPLEX(QPC) :: cc,q,s,un
        un=cmplx(1.0_qp,0.0_qp)
        s=sin(phi)
        cc=cos(phi)**2
        q=(1.0_qp-s*ak)*(1.0_qp+s*ak)
        elle_rcq=s*(rf(cc,q,un)-((s*ak)**2)*rd(cc,q,un)/3.0_qp)
        END FUNCTION elle_rcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellf_s(phi,ak)
        REAL(SP), INTENT(IN) :: phi,ak
        REAL(SP) :: ellf_s
!        Legendre elliptic integral of the 1st kind F(φ, k), evaluated using Carlson's function RF .
!        The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
        REAL(SP) :: s
        s=sin(phi)
        if(isnan(s).OR.isnan(ak))then
         write(6,*)'Arguments dans ellf, phi,ak=',phi,ak
         write(6,*)'s=',s
         STOP
        endif
        ellf_s=s*rf(cos(phi)**2,(1.0_sp-s*ak)*(1.0_sp+s*ak),1.0_sp)
        END FUNCTION ellf_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellf_d(phi,ak)
        REAL(DP), INTENT(IN) :: phi,ak
        REAL(DP) :: ellf_d
!        Legendre elliptic integral of the 1st kind F(φ, k), evaluated using Carlson's function RF .
!        The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
        REAL(DP) :: s
        s=sin(phi)
        if(isnan(s).OR.isnan(ak))then
         write(6,*)'Arguments dans ellf, phi,ak=',phi,ak
         write(6,*)'s=',s
         STOP
        endif
        ellf_d=s*rf(cos(phi)**2,(1.0_dp-s*ak)*(1.0_dp+s*ak),1.0_dp)
        END FUNCTION ellf_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellf_q(phi,ak)
        REAL(QP), INTENT(IN) :: phi,ak
        REAL(QP) :: ellf_q
!        Legendre elliptic integral of the 1st kind F(φ, k), evaluated using Carlson's function RF .
!        The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
        REAL(QP) :: s
        s=sin(phi)
        if(isnan(s).OR.isnan(ak))then
         write(6,*)'Arguments dans ellf, phi,ak=',phi,ak
         write(6,*)'s=',s
         STOP
        endif
        ellf_q=s*rf(cos(phi)**2,(1.0_qp-s*ak)*(1.0_qp+s*ak),1.0_qp)
        END FUNCTION ellf_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellf_c(phi,ak)
        COMPLEX(DPC), INTENT(IN) :: phi,ak
        COMPLEX(DPC) :: ellf_c
!        Legendre elliptic integral of the 1st kind F(φ, k), evaluated using Carlson's function RF .
!        The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
        COMPLEX(DPC) :: s
        s=sin(phi)
        ellf_c=s*rf(cos(phi)**2,(1.0_dp-s*ak)*(1.0_dp+s*ak),dcmplx(1.0_dp,0.0_dp))
        END FUNCTION ellf_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellf_rc(phi,ak)
        REAL(DP), INTENT(IN) :: phi
        COMPLEX(DPC), INTENT(IN) :: ak
        COMPLEX(DPC) :: ellf_rc
!        Legendre elliptic integral of the 1st kind F(φ, k), evaluated using Carlson's function RF .
!        The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
        COMPLEX(DPC) :: s
        s=dcmplx(sin(phi),0.d0)
        ellf_rc=s*rf(dcmplx(cos(phi)**2,0.d0),(1.0_dp-s*ak)*(1.0_dp+s*ak),dcmplx(1.0_dp,0.0_dp))
        END FUNCTION ellf_rc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellf_cq(phi,ak)
        COMPLEX(QPC), INTENT(IN) :: phi,ak
        COMPLEX(QPC) :: ellf_cq
!        Legendre elliptic integral of the 1st kind F(φ, k), evaluated using Carlson's function RF .
!        The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
        COMPLEX(QPC) :: s
        s=sin(phi)
        ellf_cq=s*rf(cos(phi)**2,(1.0_qp-s*ak)*(1.0_qp+s*ak),cmplx(1.0_qp,0.0_qp,QPC))
        END FUNCTION ellf_cq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellf_rcq(phi,ak)
        REAL(QP), INTENT(IN) :: phi
        COMPLEX(QPC), INTENT(IN) :: ak
        COMPLEX(QPC) :: ellf_rcq
!        Legendre elliptic integral of the 1st kind F(φ, k), evaluated using Carlson's function RF .
!        The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
        COMPLEX(QPC) :: s
        s=cmplx(sin(phi),0.0_qp)
        ellf_rcq=s*rf(cmplx(cos(phi)**2,0.0_qp,QPC),(1.0_qp-s*ak)*(1.0_qp+s*ak),cmplx(1.0_qp,0.0_qp,QPC))
        END FUNCTION ellf_rcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellpi_q(phi,en,ak)
        REAL(QP), INTENT(IN) :: phi,en,ak
        REAL(QP) :: ellpi_q
!        Legendre elliptic integral of the 3rd kind Π(φ, n, k), evaluated using Carlson’s functions RJ
!        and RF . (Note that the sign convention on n is opposite that of Abramowitz and Stegun.)
!        The ranges of φ and k are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
!        definition mathematica int_0^phi dt/(1-n*sin^2 t)/sqrt(1-k*sin^2 t)
        REAL(QP) :: cc,enss,q,s
        s=sin(phi)
        enss=-en*s*s
        cc=cos(phi)**2
        q=(1.0_qp-s**2*ak)
        ellpi_q=s*(rf(cc,q,1.0_qp)-enss*rj(cc,q,1.0_qp,1.0_qp+enss)/3.0_qp)
        END FUNCTION ellpi_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rtnewt(funcd,arg,x1,x2,xacc)
        REAL(DP), INTENT(IN) :: x1,x2,xacc
        REAL(DP), DIMENSION(:), INTENT(IN) :: arg
        REAL(DP) :: rtnewt
        INTERFACE
         SUBROUTINE funcd(x,arg,fval,fderiv)
         USE nrtype
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x !x: variable sur laquelle effectuer la recherche de racine. arg: tous les autres variables muettes
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg !x: variable sur laquelle effectuer la recherche de racine. arg: tous les autres variables muettes
         REAL(DP), INTENT(OUT) :: fval,fderiv
         END SUBROUTINE funcd
        END INTERFACE
        INTEGER(I4B), PARAMETER :: MAXIT=25
!        Using the Newton-Raphson method, find the root of a function known to lie in the interval
!        [x1, x2]. The root rtnewt will be refined until its accuracy is known within ±xacc. funcd
!        is a user-supplied subroutine that returns both the function value and the first derivative of
!        the function.
!        Parameter: MAXIT is the maximum number of iterations.
        INTEGER(I4B) :: j
        REAL(DP) :: df,dx,f
        rtnewt=0.5_dp*(x1+x2) !Premier essai
        do j=1,MAXIT
        call funcd(rtnewt,arg,f,df)
        dx=f/df
        rtnewt=rtnewt-dx
        if ((x1-rtnewt)*(rtnewt-x2) < 0.0)then
         write (6,*) 'nrerror: rtnewt valeur hors de lintervalle'
         STOP
        endif
        if (abs(dx) < xacc) RETURN !Convergence.
        end do
        write (6,*) 'nrerreur: nombre ditérations dépassé dans rtnewt '
        STOP 
        END FUNCTION rtnewt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rtsafed(funcd,arg,x1,x2,xacc)
        REAL(DP), INTENT(IN) :: x1,x2,xacc
        REAL(DP), DIMENSION(:), INTENT(IN) :: arg
        REAL(DP) :: rtsafed
        INTERFACE
         SUBROUTINE funcd(x,arg,fval,fderiv)
         USE nrtype
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x !x: variable sur laquelle effectuer la recherche de racine. arg: tous les autres variables muettes
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg !x: variable sur laquelle effectuer la recherche de racine. arg: tous les autres variables muettes
         REAL(DP), INTENT(OUT) :: fval,fderiv
         END SUBROUTINE funcd
        END INTERFACE
        INTEGER(I4B), PARAMETER :: MAXIT=1800
!       Using a combination of Newton-Raphson and bisection, find the root of a function bracketed
!       between x1 and x2. The root, returned as the function value rtsafed, will be refined until
!       its accuracy is known within ±xacc. funcd is a user-supplied subroutine that returns both
!       the function value and the first derivative of the function.
!       Parameter: MAXIT is the maximum allowed number of iterations.
        INTEGER(I4B) :: j
        REAL(DP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
        call funcd(x1,arg,fl,df)
        call funcd(x2,arg,fh,df)
        if ((fl > 0.0 .and. fh > 0.0) .or. &
            (fl < 0.0 .and. fh < 0.0)) then
           write(6,*)'nerreur: racine pas encadrée dans rtsafed'
           STOP
        endif
        if (fl == 0.0) then
           rtsafed=x1
           RETURN
        else if (fh == 0.0) then
           rtsafed=x2
           RETURN
        else if (fl < 0.0) then !Orient the search so that f(xl) < 0.
           xl=x1
           xh=x2
        else
           xh=x1
           xl=x2
        end if
        rtsafed=0.5_dp*(x1+x2) !Initialize the guess for root,
        dxold=abs(x2-x1)      !the “stepsize before last,”
        dx=dxold              !and the last step.
        call funcd(rtsafed,arg,f,df)
        do j=1,MAXIT !Loop over allowed iterations.
           if (((rtsafed-xh)*df-f)*((rtsafed-xl)*df-f) > 0.0 .or. &
              abs(2.0_dp*f) > abs(dxold*df) ) then !Bisect if Newton out of range, or not decreasing fast enough.
              dxold=dx
              dx=0.5_dp*(xh-xl)
              rtsafed=xl+dx
              if (xl == rtsafed) RETURN !Change in root is negligible.
           else !Newton step acceptable. Take it.
           dxold=dx
           dx=f/df
           temp=rtsafed
           rtsafed=rtsafed-dx
           if (temp == rtsafed) RETURN
           end if
           if (abs(dx) < xacc) RETURN !Convergence criterion.
           call funcd(rtsafed,arg,f,df) !One new function evaluation per iteration.
           !write(6,*)'x,f=',rtsafed,f
           if (f < 0.0) then !Maintain the bracket on the root.
              xl=rtsafed
           else
              xh=rtsafed
           end if
        end do
        write (6,*) 'nrerreur: nombre ditérations dépassé dans rtsafed '
        STOP
        END FUNCTION rtsafed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rtsafeq(funcd,arg,x1,x2,xacc)
        REAL(QP), INTENT(IN) :: x1,x2,xacc
        REAL(QP), DIMENSION(:), INTENT(IN) :: arg
        REAL(QP) :: rtsafeq
        INTERFACE
         SUBROUTINE funcd(x,arg,fval,fderiv)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), INTENT(IN) :: x !x: variable sur laquelle effectuer la recherche de racine. arg: tous les autres variables muettes
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg !x: variable sur laquelle effectuer la recherche de racine. arg: tous les autres variables muettes
         REAL(QP), INTENT(OUT) :: fval,fderiv
         END SUBROUTINE funcd
        END INTERFACE
        INTEGER(I4B), PARAMETER :: MAXIT=1800
!       Using a combination of Newton-Raphson and bisection, find the root of a function bracketed
!       between x1 and x2. The root, returned as the function value rtsafe, will be refined until
!       its accuracy is known within ±xacc. funcd is a user-supplied subroutine that returns both
!       the function value and the first derivative of the function.
!       Parameter: MAXIT is the maximum allowed number of iterations.
        INTEGER(I4B) :: j
        REAL(QP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
        call funcd(x1,arg,fl,df)
        call funcd(x2,arg,fh,df)
        if ((fl > 0.0 .and. fh > 0.0) .or. &
            (fl < 0.0 .and. fh < 0.0)) then
           write(6,*)'nerreur: racine pas encadrée dans rtsafeq'
           STOP
        endif
        if (fl == 0.0) then
           rtsafeq=x1
           RETURN
        else if (fh == 0.0) then
           rtsafeq=x2
           RETURN
        else if (fl < 0.0) then !Orient the search so that f(xl) < 0.
           xl=x1
           xh=x2
        else
           xh=x1
           xl=x2
        end if
        rtsafeq=0.5_qp*(x1+x2) !Initialize the guess for root,
        dxold=abs(x2-x1)      !the “stepsize before last,”
        dx=dxold              !and the last step.
        call funcd(rtsafeq,arg,f,df)
        do j=1,MAXIT !Loop over allowed iterations.
           if (((rtsafeq-xh)*df-f)*((rtsafeq-xl)*df-f) > 0.0 .or. &
              abs(2.0_qp*f) > abs(dxold*df) ) then !Bisect if Newton out of range, or not decreasing fast enough.
              dxold=dx
              dx=0.5_qp*(xh-xl)
              rtsafeq=xl+dx
              if (xl == rtsafeq) RETURN !Change in root is negligible.
           else !Newton step acceptable. Take it.
           dxold=dx
           dx=f/df
           temp=rtsafeq
           rtsafeq=rtsafeq-dx
           if (temp == rtsafeq) RETURN
           end if
           if (abs(dx) < xacc) RETURN !Convergence criterion.
           call funcd(rtsafeq,arg,f,df) !One new function evaluation per iteration.
           !write(6,*)'x,f=',rtsafeq,f
           if (f < 0.0) then !Maintain the bracket on the root.
              xl=rtsafeq
           else
              xh=rtsafeq
           end if
        end do
        write (6,*) 'nrerreur: nombre ditérations dépassé dans rtsafeq '
        STOP
        END FUNCTION rtsafeq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE tri_s(arr)
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
!        Sorts an array arr into ascending numerical order, by straight insertion. arr is replaced
!        on output by its sorted rearrangement.
        INTEGER(I4B) :: i,j,n
        REAL(SP) :: a
        n=size(arr)
        do j=2,n !Pick out each element in turn.
          a=arr(j)
          do i=j-1,1,-1 !Look for the place to insert it.
            if (arr(i) <= a) exit
            arr(i+1)=arr(i)
          end do
          arr(i+1)=a !Insert it.
        end do
        END SUBROUTINE tri_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE tri_d(arr)
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
!        Sorts an array arr into ascending numerical order, by straight insertion. arr is replaced
!        on output by its sorted rearrangement.
        INTEGER(I4B) :: i,j,n
        REAL(DP) :: a
        n=size(arr)
        do j=2,n !Pick out each element in turn.
          a=arr(j)
          do i=j-1,1,-1 !Look for the place to insert it.
            if (arr(i) <= a) exit
            arr(i+1)=arr(i)
          end do
          arr(i+1)=a !Insert it.
        end do
        END SUBROUTINE tri_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE tri_q(arr)
        REAL(QP), DIMENSION(:), INTENT(INOUT) :: arr
!        Sorts an array arr into ascending numerical order, by straight insertion. arr is replaced
!        on output by its sorted rearrangement.
        INTEGER(I4B) :: i,j,n
        REAL(QP) :: a
        n=size(arr)
        do j=2,n !Pick out each element in turn.
          a=arr(j)
          do i=j-1,1,-1 !Look for the place to insert it.
            if (arr(i) <= a) exit
            arr(i+1)=arr(i)
          end do
          arr(i+1)=a !Insert it.
        end do
        END SUBROUTINE tri_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION argum_d(z)
        COMPLEX(DPC), INTENT(IN) :: z
        REAL(DP) argum_d
        argum_d=atan2(imag(z),real(z))
        END FUNCTION argum_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION argum_q(z)
        COMPLEX(QPC), INTENT(IN) :: z
        REAL(QP) argum_q
        argum_q=atan2(imag(z),real(z))
        END FUNCTION argum_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE frenel(x,s,c)
        USE nrutil, ONLY : nrerror
        REAL(QP), INTENT(IN) :: x
        REAL(QP), INTENT(OUT) :: s,c
        INTEGER(I4B), PARAMETER :: MAXIT=1000
        REAL(QP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x),BIG=huge(x)*EPS,XMIN=1.5_qp
        !Computes the Fresnel integrals S(x) and C(x) for all real x.
        !Parameters: MAXIT is the maximum number of iterations allowed; EPS is the relative error;
        !FPMIN is a number near the smallest representable floating-point number; BIG is a number
        !near the machine overflow limit; XMIN is the dividing line between using the series and
        !continued fraction.
        INTEGER(I4B) :: k,n
        REAL(QP) :: a,ax,fact,pix2,sign,sum,sumc,sums,term,test
        COMPLEX(QPC) :: b,cc,d,h,del,cs
        LOGICAL(LGT) :: odd
        ax=abs(x)
        if (ax < sqrt(FPMIN)) then !Special case: avoid failure of convergence test because of underflow.
         s= 0.0 
         c=ax
        else if (ax <= XMIN) then !Evaluate both series simultaneously.
         sum=0.0
         sums=0.0
         sumc=ax
         sign=1.0_qp
         fact=PIs2*ax*ax
         odd=.true.
         term=ax
         n=3
         do k=1,MAXIT
            term=term*fact/k
            sum=sum+sign*term/n
            test=abs(sum)*EPS
            if (odd) then
               sign=-sign
               sums=sum
               sum=sumc
            else
               sumc=sum
               sum=sums
            end if
            if (term < test) exit
            odd=.not. odd
            n=n+2
         end do
         if (k > MAXIT) call nrerror('frenel: series failed')
         s=sums
         c=sumc
        else !Evaluate continued fraction by modified Lentz’s method (§5.2).
         pix2=PI*ax*ax 
         b=cmplx(1.0_qp,-pix2,kind=qpc)
         cc=BIG
         d=1.0_qp/b
         h=d
         n=-1
         do k=2,MAXIT
          n=n+2
          a=-n*(n+1)
          b=b+4.0_qp
          d=1.0_qp/(a*d+b) !Denominators cannot be zero.
          cc=b+a/cc
          del=cc*d
          h=h*del
          if (absc(del-1.0_qp) <= EPS) exit
         end do
         if (k > MAXIT) call nrerror('cf failed in frenel')
         h=h*cmplx(ax,-ax,kind=qpc)
         cs=cmplx(0.5_qp,0.5_qp,kind=qpc)*(1.0_qp-cmplx(cos(0.5_qp*pix2),sin(0.5_qp*pix2),kind=qpc)*h)
         c=real(cs)
         s=aimag(cs)
        end if
        if (x < 0.0) then !Use antisymmetry.
         c=-c
         s=-s
        end if
        CONTAINS
         FUNCTION absc(z)
         IMPLICIT NONE
         COMPLEX(QPC), INTENT(IN) :: z
         REAL(QP) :: absc
         absc=abs(real(z))+abs(aimag(z))
         END FUNCTION absc
        END SUBROUTINE frenel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION cacaf(phi,ak)
! Calcule l'intégrale elliptique de première espèce prolongée au plan complexe
! Déplace les lignes de coupure par relèvement adiabatique de la racine de l'intégrande
! Des points de brancherment demeurent là où l'argument de la racine carrée s'annulle
        COMPLEX(QPC), INTENT(IN) :: phi,ak
        COMPLEX(QPC) :: cacaf
        COMPLEX(QPC) :: int1,int2
        REAL(QP) t,dt,pref,arg(1:4)
        INTEGER it,nt
        
        cacaf=cmplx(0.0_qp,0.0_qp,kind=qpc)
        nt=54000
        dt=1.0_qp/real(nt)
        arg=(/real(phi),imag(phi),real(ak),imag(ak)/)
        int1=1.0_qp
! Règle de Simpson
        do it=0,nt
         if((it==0).OR.(it==nt))then
          pref=1.0_qp
         elseif(modulo(it,2)==1)then
          pref=4.0_qp
         elseif(modulo(it,2)==0)then
          pref=2.0_qp
         endif
         t=it*dt
         int2=int1
         int1=intef(t,arg)
         int1=int1*(-1.0_qp)**(minloc(abs((/int1+int2,int1-int2/)),DIM=1))
!Suivi adiabatique de l'intégrande. Échoue aux points de branchement où intef,intee=0
!         write(6,*)'t,int1=',t,int1
         cacaf=cacaf+pref*int1
        enddo
        cacaf=phi*cacaf*dt/3.0_qp
        
        CONTAINS
         FUNCTION intef(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC) :: intef
         COMPLEX(QPC) phi,ak
        
         phi=arg(1)+iiq*arg(2)
         ak =arg(3)+iiq*arg(4)
        
         intef=1.0_qp/sqrt(1.0_qp-ak**2.0_qp*sin(phi*x)**2.0_qp)
!         intef=1.0_qp
        
         END FUNCTION
        END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION cacae(phi,ak)
! Calcule l'intégrale elliptique de seconde espèce prolongée au plan complexe
! Déplace les lignes de coupure par relèvement adiabatique de la racine de l'intégrande
! Des points de brancherment demeurent là où l'argument de la racine carrée s'annulle
        COMPLEX(QPC), INTENT(IN) :: phi,ak
        COMPLEX(QPC) :: cacae
        COMPLEX(QPC) :: int1,int2
        REAL(QP) t,dt,pref,arg(1:4)
        INTEGER it,nt
        
        cacae=cmplx(0.0_qp,0.0_qp,kind=qpc)
        nt=54000
        dt=1.0_qp/real(nt)
        arg=(/real(phi),imag(phi),real(ak),imag(ak)/)
        int1=1.0_qp
        do it=0,nt
! Règle de Simpson
         if((it==0).OR.(it==nt))then
          pref=1.0_qp
         elseif(modulo(it,2)==1)then
          pref=4.0_qp
         elseif(modulo(it,2)==0)then
          pref=2.0_qp
         endif
         t=it*dt
         int2=int1
         int1=intee(t,arg)
         int1=int1*(-1.0_qp)**(minloc(abs((/int1+int2,int1-int2/)),DIM=1))
!Suivi adiabatique de l'intégrande. Échoue aux points de branchement où intef,intee=0
!         write(6,*)'t,int1=',t,int1
         cacae=cacae+pref*int1
        enddo
        cacae=phi*cacae*dt/3.0_qp
        
        CONTAINS
         FUNCTION intee(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC) :: intee
         COMPLEX(QPC) phi,ak
        
         phi=arg(1)+iiq*arg(2)
         ak =arg(3)+iiq*arg(4)
        
         intee=sqrt(1.0_qp-ak**2.0_qp*sin(phi*x)**2.0_qp)
!         intee=1.0_qp
        
         END FUNCTION
        END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION cacaf2(phi,ak)
! Calcule l'intégrale elliptique de première espèce prolongée au plan complexe
! Déplace les lignes de coupure par relèvement adiabatique de la racine de l'intégrande
! Des points de branchement demeurent là où l'argument de la racine carrée s'annule
        USE modsim
        IMPLICIT NONE
        COMPLEX(QPC), INTENT(IN) :: phi,ak
        COMPLEX(QPC) :: cacaf2
        COMPLEX(QPC) :: int1,int2
        REAL(QP) t,dt,pref,arg(1:4)
        INTEGER it,nt
        
        cacaf2=cmplx(0.0_qp,0.0_qp,kind=qpc)
        arg=(/real(phi),imag(phi),real(ak),imag(ak)/)
        cacaf2=qromocq(intef,0.0_qp,1.0_qp,arg,midsqlcq,1.0e-14_qp)
        cacaf2=phi*cacaf2
        
        CONTAINS
         FUNCTION intef(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), DIMENSION(size(x)) :: intef
         COMPLEX(QPC) phi,ak
        
         phi=arg(1)+iiq*arg(2)
         ak =arg(3)+iiq*arg(4)
        
         intef=1.0_qp/sqrt(1.0_qp-ak**2.0_qp*sin(phi*x)**2.0_qp)
        
         END FUNCTION
        END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION cacae2(phi,ak)
! Calcule l'intégrale elliptique de première espèce prolongée au plan complexe
! Déplace les lignes de coupure par relèvement adiabatique de la racine de l'intégrande
! Des points de branchement demeurent là où l'argument de la racine carrée s'annule
        USE modsim
        IMPLICIT NONE
        COMPLEX(QPC), INTENT(IN) :: phi,ak
        COMPLEX(QPC) :: cacae2
        COMPLEX(QPC) :: int1,int2
        REAL(QP) t,dt,pref,arg(1:4)
        INTEGER it,nt
        
        cacae2=cmplx(0.0_qp,0.0_qp,kind=qpc)
        arg=(/real(phi),imag(phi),real(ak),imag(ak)/)
        cacae2=qromocq(intee,0.0_qp,1.0_qp,arg,midpntcq,1.0e-14_qp)
        cacae2=phi*cacae2
        
        CONTAINS
         FUNCTION intee(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), DIMENSION(size(x)) :: intee
         COMPLEX(QPC) phi,ak
        
         phi=arg(1)+iiq*arg(2)
         ak =arg(3)+iiq*arg(4)
        
         intee=sqrt(1.0_qp-ak**2.0_qp*sin(phi*x)**2.0_qp)
        
         END FUNCTION
        END FUNCTION
END MODULE recettes
