!> @author Daniel Menendez Crespo
!> @brief 
!>     Evaluate an integral using quadleg
!>     with 2-7 number of points.


program gaussquad
   
use functions, only: f1
use dml_int  , only: quadgleg
   
    implicit none

    integer, parameter      :: n = 7, ounit=7
    real(kind=8), parameter :: a = 0.d0
    real(kind=8), parameter :: b = 1.d0

    open(ounit,file='gleg.out',status='unknown')

    write(ounit,*) '==================================================='
    write(ounit,*) 'Gaussian quadrature with Gauss-Legendre polinomials'
    write(ounit,*) '==================================================='
    write(ounit,*) 
    write(ounit,*) '\int_0^1 x^2 e^{-x} \mathrm{d}x \approx 0.16'
    write(ounit,*) 
    call quadgleg(f1, a, b, n, ounit)


end program gaussquad
