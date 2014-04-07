!> @brief Integrands to be used by numerical integration procedures
!> @author Daniel Menendez Crespo

module functions 

    implicit none

!==========================================================
contains
!==========================================================

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = x^2 e^{-x} \f]
!----------------------------------------------------------
pure function f1(x) result(f)

    real(kind=8), intent(in) :: x
    real(kind=8) :: f

    f = x*x*exp(-x)

end function f1

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = e^{3x} \sin(3x) \f]
!----------------------------------------------------------
pure function f2(x) result(f)

    real(kind=8), intent(in) :: x
    real(kind=8) :: f

    f = exp(3*x)*sin(3*x)

end function f2

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = e^{-x} \f]
!----------------------------------------------------------
pure function f3(x) result(f)

    real(kind=8), intent(in) :: x
    real(kind=8) :: f

    f = exp(-x)

end function f3

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = \frac{x^2}{1 + x^3} \f]
!----------------------------------------------------------
pure function f4(x) result(f)

    real(kind=8), intent(in) :: x
    real(kind=8) :: f

    f = x*x/(1+x*x*x)

end function f4

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = x^2 \sin(nx) \f]
!----------------------------------------------------------
pure function f5(x) result(f)

    real(kind=8), intent(in) :: x
    real(kind=8) :: f
    real(kind=8) :: n

    n = 4.d0

    f = x*x*sin(n*x)

end function f5

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = \sqrt{4 - x^2 - y^2 - z^2} \f]
!----------------------------------------------------------
pure function f6(x,y,z) result(f)

    real(kind=8), intent(in) :: x,y,z
    real(kind=8) :: f

    f = sqrt(4.d0-x*x-y*y-z*z)

end function f6

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ r = x^2 \f]
!----------------------------------------------------------
pure function f(x) result(r)   ! f(x) = x**2

    real(kind=8), intent(in) :: x
    real(kind=8)             :: r

    r = x*x

end function f

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ r = 3y \f]
!----------------------------------------------------------
pure function g(y)  result(r)  ! g(y) = 1/(sqrt(xn**2- y**2))

    real(kind=8), intent(in) :: y
    real(kind=8)             :: r

    r = 3* y

end function g

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ du = u\cos(t) \f]
!----------------------------------------------------------
pure subroutine f_1(t, u, du)

    real(kind=8), intent(in)  :: t, u
    real(kind=8), intent(out) :: du
    intrinsic :: cos

    du = u*cos(t)

end subroutine f_1

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ r = \frac{tx - x^2}{t^2} \f]
!----------------------------------------------------------
pure function f_2(t,x) result(r)

    real(kind=8), intent(in)  :: t, x
    real(kind=8)              :: r

    r = (t*x -x**2.0)/t**2.0

    return

end function f_2 

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ r = \frac{2}{\sqrt{\pi\exp(-x^2)}} \f]
!----------------------------------------------------------
pure function f_3(x) result(r)

    use constants, only : c => dml_const

    real(kind=8), intent(in) :: x
    real(kind=8)             :: r
    intrinsic                :: sqrt, exp

    !f_2 = 1.0d0/x
    r = 2.d0/sqrt(c%pi)*exp(-x*x)

end function f_3

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ r = \frac{t}{0.5 + \ln(t)}\f]
!----------------------------------------------------------
pure function u(t) result(r)

    real(kind=8), intent(in) :: t
    real(kind=8)             :: r
    intrinsic                :: log

    r = t/(0.5d0 + log(t))

    return

end function u 

!----------------------------------------------------------
!> @author Numerical Recipes
!----------------------------------------------------------
SUBROUTINE GAULEG(X1,X2,x,w,n)

    IMPLICIT NONE

    INTEGER(KIND=4) :: N,M,I,J
    REAL(KIND=8), dimension(n) :: W,X
    REAL(KIND=8) :: XM,XL,YN,X1,X2,EPS,P1,P2,P3,YJ,PI,Z1,Z,PP,YI
    YN=DFLOAT(N)
    EPS=1.D-14
    M=(N+1)/2
    XM=0.5D0*(X1+X2)
    XL=0.5D0*(X2-X1)
    PI=DACOS(-1.D0)
    DO 12 I=1,M
    YI=DFLOAT(I)
    Z=DCOS(PI*(YI-0.25D0)/(YN+0.5D0))
1   CONTINUE
    P1=1.D0
    P2=0.D0
    DO 11 J=1,N
    P3=P2
    P2=P1
    YJ=DFLOAT(J)
    P1=((2.D0*YJ-1.D0)*Z*P2-(YJ-1.D0)*P3)/YJ
    11    CONTINUE
    PP=YN*(Z*P1-P2)/(Z*Z-1.D0)
    Z1=Z
    Z=Z1-P1/PP
    IF(DABS(Z-Z1).GT.EPS)GO TO 1
    X(I)=XM-XL*Z
    X(N+1-I)=XM+XL*Z
    W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
    W(N+1-I)=W(I)
12  CONTINUE

END SUBROUTINE GAULEG


end module functions
