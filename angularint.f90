MODULE angularint
CONTAINS
FUNCTION Iuanaly(en,k,q,xmin,xmax) !Computes analytically the integrals int_-1^1 num/(en+eps(\vec{k}-\vec{q})+i*0^+)
 USE recettes            !with num=U^2 in Iuanaly(1), V^2 in Iuanaly(2) and UV in Iuanaly(3)
 IMPLICIT NONE           !Use conjg the get the integrals with a -i*0^+
 REAL(QP), INTENT(IN) :: en,k,q,xmin,xmax
 COMPLEX(QPC) Iuanaly(1:3)
 COMPLEX(QPC) x1,x2

 REAL(QP) I1,I2,I3,imI1,imI2,rac

 Iuanaly(:)=0.0_qp
 imI1=0.0_qp
 imI2=0.0_qp

!xmin/xmax are en(ergy) independent. Precompute them as:
!   xiP=k**2+q**2+2.0_qp*k*q-x0
!   xiM=k**2+q**2-2.0_qp*k*q-x0
!   epsP=sqrt(xiP**2+1.0_qp)
!   epsM=sqrt(xiM**2+1.0_qp)
!   xmin=epsP-xiP
!   xmax=epsM-xiM

 if(abs(en)>1.0_qp)then
  rac=sqrt(en**2-1.0_qp)
  x1=cmplx(-en+rac,0.0_qp,kind=qpc)
  x2=cmplx(-en-rac,0.0_qp,kind=qpc)
 else
  rac=sqrt(1.0_qp-en**2)
  x1=cmplx(-en,+rac,kind=qpc)
  x2=cmplx(-en,-rac,kind=qpc)
 endif
 I1=log(abs((xmax-x1)*(xmax-x2)/(xmin-x1)/(xmin-x2)))
 I3=log(xmax/xmin)
 if(abs(en)<1.0_qp)then
  I2=2.0_qp*(argum(xmax-x1)-argum(xmin-x1))/rac
 else
  I2=log(abs((xmax-x1)*(xmin-x2)/(xmax-x2)/(xmin-x1)))/rac
  if((xmax>real(x1)).AND.(xmin<real(x1)))then
   imI1=imI1+PI*sign(1.0_qp,en)
   imI2=imI2+PI*sign(1.0_qp,en)
  endif
  if((xmax>real(x2)).AND.(xmin<real(x2)))then
   imI1=imI1-PI*sign(1.0_qp,en)
   imI2=imI2+PI*sign(1.0_qp,en)
  endif
  imI2=imI2/rac
 endif
 Iuanaly(1)=( (I1+iiq*imI1)-en*(I2+iiq*imI2))          /(4*k*q) !I_V^2
 Iuanaly(2)=(-(I1+iiq*imI1)-en*(I2+iiq*imI2)+2.0_qp*I3)/(4*k*q) !I_U^2
 Iuanaly(3)=  (I2+iiq*imI2)                            /(4*k*q) !I_UV
END FUNCTION Iuanaly
END MODULE angularint
