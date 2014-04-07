!> @brief Numerical integration procedures
!> @author Daniel Menendez Crespo

module dml_int

use interfaces
use functions

    implicit none

!==========================================================
contains
!==========================================================

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  Gaussian quadrature with Gauss-Legendre polinomials
!> @warning Specific to only one function since error is computed
!>          against the real value of one function
!----------------------------------------------------------
subroutine quadgleg(f, a, b, n, ounit)

    procedure(func)          :: f
    real(kind=8), intent(in) :: a, b
    integer     , intent(in) :: n, ounit

    integer      :: i, k
    real(kind=8) :: suma, c, m

    real(kind=8), dimension(n) :: res
    real(kind=8), dimension(n)   :: w, t
    
    write(ounit,*) 'n    integral      Error'
    write(ounit,*) '----------------------------'

    res = 0.d0

    c = (a+b)*0.5d0
    m = (b-a)*0.5d0

    do i=2, n
        call gauleg(-1.d0, 1.d0, t, w, i)

        suma = 0.d0
        do k=1, i
            suma = suma + w(k)*f(c+m*t(k))*m
        end do
        res(i) = suma
        write(ounit,8) i, res(i), abs(res(i)-0.160603d0)
        !write(ounit,8) i, res(i), abs(res(i)-res(i-1))
    end do

8   format(i2, 20f13.8)

end subroutine quadgleg

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  Monte-Carlo integration with the hit and miss method 
!----------------------------------------------------------
subroutine hitmiss(f, a, b, h, n, res)

    procedure(func)           :: f
    real(kind=8), intent(in)  :: a, b, h
    real(kind=8), intent(out) :: res
    integer, intent(in)       :: n

    integer      :: i, ntrue
    real(kind=8) :: x, y, ratio, ran


    do i=1, n
        call random_number(ran)
        x = b*ran
        call random_number(ran)
        y = h*ran
        if (y.le.f(x)) then
            ntrue = ntrue + 1
        end if
        ratio = real(ntrue, kind=8)/real(n, kind=8)
        res = ratio*b*h
    end do

end subroutine hitmiss

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  Monte-Carlo integration with the sample mean method
!----------------------------------------------------------
subroutine samplemean(f, a, b, h, n, res)
    procedure(func) :: f
    real(kind=8), intent(in) :: a, b, h
    real(kind=8), intent(out) :: res
    integer, intent(in) :: n

    integer :: i
    real(kind=8) :: x, area, suma, ran


    area = (b-a)*h
    suma = 0.d0
    do i=1, n
        call random_number(ran)
        x = b*ran
        suma = suma + f(x)
    end do
    res = area*suma/n

end subroutine samplemean

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  Monte-Carlo integration with the sample mean method
!>         in 3D
!----------------------------------------------------------
subroutine samplemean3d(f, x0, x1, y0, y1, z0, z1, n, res)

    procedure(f3D)            :: f
    real(kind=8), intent(in)  :: x0, x1, y0, y1, z0, z1
    real(kind=8), intent(out) :: res
    integer, intent(in)       :: n

    integer :: i
    real(kind=8) :: x, y, z,  volume, suma, ran


    volume = (x1-x0)*(y1-y0)*(z1-z0)
    suma = 0.d0
    do i=1, n
        call random_number(ran)
        x = x1*ran + x0
        call random_number(ran)
        y = y1*ran + y0
        call random_number(ran)
        z = z1*ran + z0
        suma = suma + f(x,y,z)
    end do
    res = volume*suma/n

end subroutine samplemean3d

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  Trapezoid composite rule integration
!----------------------------------------------------------
subroutine trapecio(f, x0, xn, n, integral)
!                   in  in in  in  out
!                   fun re re  int re

    procedure(func)           :: f
    real(kind=8), intent(in)  :: x0, xn
    integer,      intent(in)  :: n
    real(kind=8), intent(out) :: integral

    integer                   :: k
    real(kind=8)              :: s, h
    intrinsic                 :: real

    s = 0.d0
    h = (xn-x0)/real(n, kind=8)

    do k=1, n-1
      s = s + f(x0 + h*k)
    end do

    integral = h * (f(x0)/2.d0 + f(xn)/2.d0 + s)

    return

end subroutine trapecio

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  Simpson composite rule integration
!----------------------------------------------------------
subroutine simpson(f, x0, xn, n, integral)

    procedure(func)           :: f
    real(kind=8), intent(in)  :: x0, xn   ! valores de contorno
    integer,      intent(in)  :: n        ! numero de puntos
    real(kind=8), intent(out) :: integral ! el resultado

    integer                   :: k        ! contador
    real(kind=8)              :: s, h

    s = 0.d0
    h = (xn-x0)/real(n, kind=8)  ! incremento

    do k=1, n-1, 2 ! impares
      s = s + 4.d0 * f(x0 + h*k)
    end do

    do k=2, n-2, 2  ! pares
       s = s + 2.d0 * f(x0 + h*k)
    end do

    integral = h * (f(x0) + f(xn) + s) / 3.d0

    return

end subroutine simpson

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  Runge-Kutta order 4 integration
!> @todo use abstract interface
!----------------------------------------------------------
pure subroutine rk4(t0, u0, dt, f_1, u1)

    real(kind=8), intent(in)  :: t0, u0, dt
    real(kind=8), intent(out) :: u1
    real(kind=8)              :: f1,f2,f3,f4
    interface f
        pure subroutine f_1(t,u,du)
            real(kind=8), intent(in)  :: t, u
            real(kind=8), intent(out) :: du
        end subroutine f_1
    end interface f

!
!  Cuatro puntos
!
    call f_1(t0, u0, f1)
    call f_1(t0 + dt / 2.0D+00, u0 + dt * f1 / 2.0D+00, f2)
    call f_1(t0 + dt / 2.0D+00, u0 + dt * f2 / 2.0D+00, f3)
    call f_1(t0 + dt, u0 + dt * f3, f4)
!
!  Los combina  para estimar para t1 = t0 + dt.
!
    u1 = u0 + dt*(f1 + 2.0D+00 * f2 + 2.0D+00 * f3 + f4) / 6.0D+00

    return

end subroutine rk4 

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  Romberg method integration
!> @details
!> - convergence check
!> - optional iteration values
!----------------------------------------------------------
pure subroutine romberg(f, a, b, maxi, eps, val, res)

    procedure(func)                           :: f
    integer,      intent(in)                  :: maxi
    real(kind=8), intent(in)                  :: a, b, eps
    real(kind=8), intent(out)                 :: val
    real(kind=8), dimension(:,:), allocatable, &
                  intent(out),    optional    :: res

    real(kind=8), dimension(0:maxi,0:maxi)    :: r
    real(kind=8)                              :: h, suma
    integer                                   :: i,n,m

    h = b-a  
    r(0,0) = 0.5d0*h*(f(a) + f(b))    
    !
    do n=1, maxi - 1 
        h = 0.5d0*h   
        
        suma = 0.0d0   
        do i=1, 2**(n-1)      
            suma = suma + f(a + real(2*i-1, kind=8)*h)   
        end do
        
        r(n,0) = 0.5d0*r(n-1,0) + h*suma   
        
        do m=1, n
            !r(n,m) = r(n,m-1) + (r(n,m-1) - r(n-1,m-1)) / ((4.d0**m) - 1.0d0)  
            r(n,m) = r(n,m-1) + (r(n,m-1) - r(n-1,m-1)) / real((4**m) - 1, kind=8)  
        end do
        ! convergence check
        if (dabs(r(n,n)-r(n,n-1)).le.eps) then
            exit
        end if
    end do
 
    val = r(n-1, n-1)
    if (present(res)) then
        allocate(res(n,n))
        res = r(0:n-1, 0:n-1)
    end if

    return

end subroutine romberg

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  Takes Trapezoid or Simpson and increases points
!>          until convergence
!----------------------------------------------------------
subroutine multpts(method, f, a, b, thres, ounit)

    procedure(func)                :: f
    real(kind=8), intent(in)       :: a,b,thres
    integer, intent(in)            :: ounit
    integer, parameter             :: k_max=8
    real(kind=8), dimension(k_max) :: res
    integer                        :: i, k, n

    n = 1
    k = 1
    call method(f, a, b, n, res(k))
    write(ounit,*) '     Integral                    error                 No. points'
    write(ounit,*) '------------------------------------------------------------------------'
    do k=1, k_max-1
        n = 10*n
        call method(f, a, b, n, res(k+1))
        write(ounit,*) res(k+1), abs(res(k+1) - res(k)), n
        if (abs(res(k+1)-res(k)).lt.thres) then
            write(ounit,*)'Convergence reached'
            write(ounit,*) 'Result: ', res(k)
            exit 
        end if
    end do 

end subroutine multpts

end module dml_int
