!> @brief Numerical Optimization procedures
!> @author Daniel Menendez Crespo

module dml_opt

use functions

    implicit none

!==========================================================
contains
!==========================================================

!----------------------------------------------------------
!> @author  Daniel Menendez Crespo
!> @brief   Bisection method
!> @details Finds the roots of a function. Boltzano based.
!----------------------------------------------------------
subroutine bisect(f, a, b, m, thres, outunit)                                        

    procedure(func)          :: f
    integer, intent (in)     :: m, outunit
    real(kind=8), intent(in) :: thres
    real(kind=8) :: a, b, c, fa, fb, fc, error                              
    integer :: n

    fa = f(a)                                                       
    fb = f(b)
    if ((abs(fa).le.thres).or.(abs(fb).le.thres)) then
        write(outunit,*) 'Either a or b are roots  ','  f(a)=',fa,'  f(b)=',fb
        return
    end if
    if (sign(1.0d0,fa) == sign(1.0d0,fb)) then                       
        write(outunit,*),"function has same sign at",a,"and",b                       
        return
    end if                                                            
    error = b - a
    do n = 0, m                                                    
        error = error*0.5d0
        c  = a + error
        fc = f(c)
        if (abs(fc).lt.thres) then
            write(outunit,*)'Root: c = ', c
            write(outunit,*)'No. of iterations:', n
            return
        else
            if (sign(1.0d0,fa) /= sign(1.0d0,fc)) then           
                b = c                  
                fb = fc                                                
            else                                                
                a = c    
                fa = fc  
            end if                                               
        end if
    end do                                                 

    close(outunit)


end subroutine bisect 

!----------------------------------------------------------
!> @author  Daniel Menendez Crespo
!> @brief   Regula falsi method
!> @details Finds the roots of a function. Boltzano based.
!----------------------------------------------------------
subroutine regfalsi(f, a, b, m, thres, outunit)                        

    procedure(func)          :: f
    integer, intent (in)     :: m, outunit
    real(kind=8), intent(in) :: thres
    real(kind=8) :: a, b, c, fa, fb, fc, error               
    integer :: n

    fa = f(a)                                                       
    fb = f(b)
    if ((abs(fa).le.thres).or.(abs(fb).le.thres)) then
        write(7,*) 'Either a or b are roots  ','  fa=',fa,'  fb=',fb
        return
    end if
    if (sign(1.0d0,fa) == sign(1.0d0,fb)) then                       
        write(*,*) "function has same sign at",a,"and",b                       
        return
    end if                                                            
    do n = 0, m                                                    
        error = b - a
        c = a - (fa*error)/(fb - fa)
        fc = f(c)
        if (abs(fc).lt.thres) then
            write(outunit,*)'Root: c = ', c
            write(outunit,*)'No. of iterations:', n
            return
        else
            if (sign(1.0d0,fa) /= sign(1.0d0,fc)) then           
                b = c                  
                fb = fc                                                
            else                                                
                a = c    
                fa = fc  
            end if                                               
        end if
    end do                                                 

end subroutine regfalsi 

!----------------------------------------------------------
!> @author  Daniel Menendez Crespo
!> @brief   Quadratic interpolation method with equidistant points
!> @details Finds the minimum of a function. 
!----------------------------------------------------------
subroutine qinterp(f, x0, h0, k_max, thres)

    procedure(func)          :: f
    integer, intent(in)      :: k_max
    real(kind=8), intent(in) :: thres, h0
    real(kind=8), intent(inout) :: x0

    real(kind=8), dimension(k_max) :: x, h, xp, xpp, m
    integer                  :: k

    x(1) = x0 
    h(1) = h0
    k = 1

    do while (abs(x(k+1)-x(k)).gt.thres)
        m(k+1) = x(k+1)
        xp(k+1) = x(k+1) - h(k+1)
        xpp(k+1) = x(k+1) + h(k+1)
           
        if ((f(m(k+1)).ge.f(xp(k+1))).or.(f(m(k+1)).ge.f(xpp(k+1)))) then
            k = 0
            call random_number(x(1))
            write(*,*) 'New starting point: ', x(1)
            cycle
        end if

        k = k + 1
        if ((k+1).eq.k_max) then
            write(*,*) 'Max No. iterations'
            return
        end if
        x(k+1) = x(k) - h(k)*(0.5d0) *              &
               (f(xpp(k)) - f(xp(k)))/(f(xpp(k))  &
               -2.d0*f(m(k)) + f(xp(k)))
        h(k+1) = h(k)*(0.5d0)**k
    end do

    x0 = x(k+1)


end subroutine qinterp

!----------------------------------------------------------
!> @author  Daniel Menendez Crespo
!> @brief   Newton-Raphson method 1D
!> @details Finds the roots of a 1D function. 
!----------------------------------------------------------
subroutine nroots1D(start, fun, derivative_1d)

    procedure(func)          :: fun
    procedure(fp1D)          :: derivative_1d
    real(kind=8), intent(in) :: start

    integer     ,    parameter :: nmax  = 1e5
    real(kind=8),    parameter :: thres = 1.d-4 !> for convergence
    real(kind=8),    parameter :: h     = 1.d-4 !> for derivative
    real(kind=8), dimension(nmax) :: x, f, fp
    real(kind=8) :: t
    integer      :: j

    x = 0.d0
    f  = 0.d0
    fp = 0.d0

    open(7,file='nroots.out',status='unknown')
    write(7,*)

    x(1) = start

    f(1) = fun(x(1))
    call derivative_1d(fun, x(1), h, fp(1))

    do j=1, nmax
        x(j+1) = x(j) - f(j)/fp(j)
        t = x(j+1)
        f(j+1) = fun(t)
        call derivative_1d(fun, t, h, fp(j+1))
        if (abs(x(j+1)-x(j)).le.thres) then
            write(7,*)'Root: x = ', x(j+1)
            write(7,*)'No. of iterations:', j
            return
        end if
    end do

end subroutine nroots1D

!----------------------------------------------------------
!> @author  Daniel Menendez Crespo
!> @brief   Newton-Raphson method 2D
!> @details This program minimizes a 2d-function using the 
!>          Newton-Raphson method
!> @warning hess is inverted while the name is still the same
!----------------------------------------------------------
subroutine nrmin(x,fun, itmax)

    real(kind=8)                                :: fun
    integer, intent(in)                         :: itmax
    real(kind=8), dimension(:,:), intent(inout) :: x
    integer, parameter                          :: dim=2
    real(kind=8), parameter                     :: thres=1.d-4
    real(kind=8), dimension(dim,dim)            :: hess
    real(kind=8), dimension(dim)                :: grad,add
    integer :: it,i,j

    open(7,file='nrmin.out',status='unknown')

    it = 0.d0

    call gradient(x(:,it),fun,grad)

    do while (maxval(abs(grad)).gt.thres)
        it = it + 1

        if (it.eq.itmax) then
            write(*,*) 'Max No. iterations'
            return
        end if

        call hessian(x(:,it-1),fun,hess)

        call inverse(hess)

        add = matmul(hess, grad)

        x(:,it) = x(:,it-1) - add(:)

        call gradient(x(:,it),fun,grad)

    end do

    write(7,*)
    write(7,*) ' Minimum:'
    write(7,9) '  x               y              f(x,y)'
    write(7,*) '-------------------------------------------------------'
    write(7,99)  x(1,it-1), x(2,it-1), fun(x(1,it-1),x(2,it-1))

9 format (10x, a, 50x, a, 50x, a)
99 format (3f18.13)
    
end subroutine nrmin

!----------------------------------------------------------
!> @author  Daniel Menendez Crespo
!> @brief   Newton-Raphson method 3D roots
!> @details Solves a set of non-linear equations of 
!>          3D-functions using the Newton-Raphson method
!> @warning jac is inverted while the name is still the same
!----------------------------------------------------------
subroutine nreqs(xyz_0, thres, itmax)

    use functions

    integer, parameter                       :: dim=3
    integer, intent(in)                      :: itmax
    real(kind=8), intent(in)                 :: thres
    real(kind=8), dimension(dim), intent(in) :: xyz_0
    real(kind=8), dimension(dim, 0:itmax)    :: xyz
    real(kind=8), dimension(dim,dim)         :: jac
    real(kind=8), dimension(dim)             :: add, fval
    integer                                  :: it, i, j

    open(7,file='nreqs.out',status='unknown')

    it = 0
    xyz(:,it) = xyz_0
      
    call funs(xyz(:,it), fval)
    
    do while (maxval(abs(fval)).gt.thres)

        call jacobian(xyz(:,it), jac)

        call inverse(jac)

        call funs(xyz(:,it), fval)

        add = matmul(jac, fval)

        xyz(:,it+1) = xyz(:,it) - add(:)
    
        call funs(xyz(:,it+1), fval)

        it = it + 1
        if (it.gt.itmax) then
            write(*,*) 'Max No. iterations'
            return
        end if
    end do

    write(7,99) "Solution: ", xyz(:,it)
    write(7,999) 'No. of iterations: ', it
     
99 format (a13,3f15.8)
999 format (a25,i3)

end subroutine nreqs

end module dml_opt
