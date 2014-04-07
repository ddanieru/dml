!> @brief Functions to be used by numerical optimization procedures
!> @author Daniel Menendez Crespo

module functions

use interfaces

    implicit none

!==========================================================
contains
!==========================================================

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = \exp(x) - x^2 + 3x - 2 \f]
!----------------------------------------------------------
pure function f1(x) result(f)  

    real(kind=8)              :: f
    real(kind=8), intent (in) :: x

    f = exp(x) - x*x + 3.0d0*x - 2.0d0

end function f1 

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = x^5 - \frac{4x}{5} + 0.5 \f]
!----------------------------------------------------------
pure function f2(x) result(f)

    real(kind=8)             :: f
    real(kind=8), intent(in) :: x

    f = x**(5.d0) - (4.d0*x)/5.d0 + 0.5d0

end function f2

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = \exp(-x) - 2 \f]
!----------------------------------------------------------
pure function f3(x) result(f)

    real(kind=8)             :: f
    real(kind=8), intent(in) :: x

    f = exp(-x) - 2.d0 

end function f3
  
!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = 3\exp(x) - 4\cos(x) \f]
!----------------------------------------------------------
pure function f4(x) result(f)

    real(kind=8)             :: f
    real(kind=8), intent(in) :: x

    f = 3*exp(x) - 4.d0*cos(x)

end function f4

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = x^2 - 3xy + (y+2)^4 \f]
!----------------------------------------------------------
pure function f6(x) result(f)

    real(kind=8)             :: f
    real(kind=8), dimension(2), intent(in) :: x

    f = x(1)*x(1) - 3.d0*x(1)*x(2) + (x(2)+2.d0)**(4.d0)

end function f6

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = \exp(x^2 + (y+2)^2) \f]
!----------------------------------------------------------
pure function f7(x) result(f) 

    real(kind=8)             :: f 
    real(kind=8), dimension(2), intent(in) :: x

    f = exp(x(1)*x(1) + (x(2)+2.d0)**2.d0)

end function f7

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = 3x - \cos(yz) - 0.5 \f]
!----------------------------------------------------------
pure function f8(x) result(f)

    real(kind=8)             :: f
    real(kind=8), dimension(3), intent(in) :: x

    f = 3.d0*x(1) - cos(x(2)*x(3)) - 0.5d0

end function f8

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = x^2 - 81(y+0.1)^2 + \sin(z) + 1.06 \f]
!----------------------------------------------------------
pure function f9(x) result(f) 

    real(kind=8)             :: f
    real(kind=8), dimension(3), intent(in) :: x 

    f = x(1)**2.d0 - 81.d0*(x(2)+0.1d0)**2 + sin(x(3)) + 1.06

end function f9

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f = \exp(-xy) + 20z + \frac{10\pi-3}{3} \f]
!----------------------------------------------------------
pure function f10(x) result(f)

    real(kind=8)             :: f   
    real(kind=8), dimension(3), intent(in) :: x 
    real(kind=8), parameter  :: pi=acos(-1.d0)

    f = exp(-x(1)*x(2)) + 20.d0*x(3) + (10.d0*pi-3.d0)/3.d0

end function f10

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  
!> Evaluates the functions at a point and returns its values
!> in a vector
!----------------------------------------------------------
pure subroutine funs(xyz, fval)

    real(kind=8), dimension(3), intent(in)  :: xyz
    real(kind=8), dimension(3), intent(out) :: fval

    fval(1) = f8(xyz) 
    fval(2) = f9(xyz) 
    fval(3) = f10(xyz) 

end subroutine funs

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ f'(x) \approx \frac{f(x+h)-f(x-h)}{2h} \f]
!> @details Numeric finite central first derivative 
!----------------------------------------------------------
pure subroutine derivative_1d(f, x, h, fp)

    procedure(func) :: f
    real(kind=8), intent(in)  :: x, h
    real(kind=8), intent(out) :: fp

    fp = 0.d0
    fp = f(x+h) - f(x-h)
    fp = fp/(2.d0*h)

end subroutine derivative_1d

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ \nabla f(\mathbf{x}) \f]
!> @details Gradient vector for a ND function
!> @param[out] grad \f$ \nabla f(\mathbf{x}) \f$
!> @param[in]  x    Coordinates vector
!> @param[in]  f    Function to evaluate
!----------------------------------------------------------
subroutine gradient(x, f, grad)

    real(kind=8)            :: f
    !procedure(f3D)          :: f
    real(kind=8), dimension(:), intent(in)  :: x
    real(kind=8), dimension(:), intent(out) :: grad
    real(kind=8), dimension(:), allocatable :: x1, x2
    real(kind=8), parameter :: h=1.d-4
    integer                 :: i

    i = size(x)
    allocate(x1(i))
    allocate(x2(i))

    do i=1, size(x)
        x1 = x
        x2 = x
        x1(i) = x(i) + h
        x2(i) = x(i) - h
        grad(i) = f(x1) - f(x2)
    end do

    grad(:) = grad(:)/(2.d0*h)               

end subroutine gradient

!----------------------------------------------------------
!> @author  Daniel Menendez Crespo
!> @brief   \f[ J(x) \f]
!> @details Jacobian matrix from a set of functions
!> @param[out] jac \f$ J(x) \f$
!> @param[in]  x   Coordinates vector
!----------------------------------------------------------
subroutine jacobian(x, jac)

    real(kind=8), dimension(3),     intent(in) :: x
    real(kind=8), dimension(3,3), intent(out)  :: jac
    real(kind=8), parameter  :: h=1.d-3
    integer                  :: i, j


    call gradient(x, f8 , jac(:,1))
    call gradient(x, f9 , jac(:,2))
    call gradient(x, f10, jac(:,3))

end subroutine jacobian

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  \f[ H(x) \f]
!> @details Central difference approximation to the
!>          Hessian matrix from a set of functions
!> @param[out] hess \f$ H(x) \f$
!> @param[in]  y    Coordinates vector
!> @param[in]  f    Function to evaluate
!> @todo generalize to any m x m matrix
!----------------------------------------------------------
subroutine hessian(y, f, hess)

    real(kind=8)            :: f
    !procedure(f2D)          :: f
    real(kind=8), dimension(2),   intent(in)  :: y
    real(kind=8), dimension(2,2), intent(out) :: hess
    real(kind=8), dimension(2)                :: x
    real(kind=8), parameter :: h=1.d-4
    integer                 :: i


    hess = 0.d0
    x = y

    !> Diagonal elements
    do i=1, 2
        x(i)      = y(i) + h
        hess(i,i) = f(x)
        x(i)      = y(i) - h
        hess(i,i) = hess(i,i) + f(x)
        hess(i,i) = hess(i,i) - 2.d0*f(y)
    end do

    !> Non diagonal elements
    x(1) = y(1) + h/2.d0
    x(2) = y(2) + h/2.d0
    hess(1,2) = f(x) 
    x(1) = y(1) - h/2.d0
    x(2) = y(2) + h/2.d0
    hess(1,2) = hess(1,2) - f(x)
    x(1) = y(1) + h/2.d0
    x(2) = y(2) - h/2.d0
    hess(1,2) = hess(1,2) - f(x)
    x(1) = y(1) - h/2.d0
    x(2) = y(2) - h/2.d0
    hess(1,2) = hess(1,2) + f(x)

    hess      = hess/(h*h)
    hess(2,1) = hess(1,2)

end subroutine hessian

!----------------------------------------------------------
!> @author  Daniel Menendez Crespo
!> @brief   Inversion of a matrix by Gauss-Jordan method
!> @details Subroutine that calculates the inverse of a matrix 
!>          within the Gauss-Jordan method
!----------------------------------------------------------
subroutine inverse(mat)

    real(kind=8), dimension(:,:), intent(inout) :: mat
    real(kind=8), dimension(:,:), allocatable :: a
    logical, dimension(:), allocatable :: done

    real(kind=8), parameter :: thres = 1.0d-13
    real(kind=8) :: max_dia, m
    integer :: i, j, k, l, k_max, n

    n = size(mat, dim=1)
    allocate(a(n,2*n))
    allocate(done(n))

    a(:,:n) = mat

    done = .false.

! mat = ones()
    mat = 0.d0
    do i=1, n
        mat(i,i) = 1.d0
    end do

! Extend the matrix
    a(:,(n+1):) = mat

    do k=1, n

!   Choose maximum diagonal
        max_dia = 0.0d0
        do l=1, n
            if (done(l)) cycle
            if (abs(a(l,l)).gt.max_dia) then
                k_max = l
                max_dia = abs(a(l,l))
            end if
        end do

        if (max_dia.le.thres) then
            write(*,*) 'Error: singular matrix given'
            return 
        end if

!   Scale k_max row
        m = a(k_max,k_max)
        do j=1, 2*n
            a(k_max,j) = a(k_max,j) / m
        end do

!   Treat rows before
        do i=1, k_max-1
            m = a(i,k_max)
            do j=1, 2*n
                a(i,j) = a(i,j) - m*a(k_max,j)
            end do
        end do

!   Treat rows after
        do i=k_max+1, n
            m = a(i,k_max)
            do j=1, 2*n
                a(i,j) = a(i,j) - m*a(k_max,j)
            end do
        end do

!   Mark kmax as done
        done(k_max) = .true.
    end do

    mat = a(:n,n+1:)


end subroutine inverse

end module functions
