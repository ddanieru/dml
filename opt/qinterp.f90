!> @file:   qinterp.f90
!> @brief   Quadratic interpolation method test
!> @author  Daniel Menendez Crespo
!> @details Minimizes a function using a parabolic or quadratic 
!>          interpolation method with equidistant points     

program quadinterp

    use functions
    use dml_opt

    implicit none

    integer, parameter             :: k_max  = 20
    real(kind=8), parameter        :: thres = 1.0d-4
    real(kind=8), parameter        :: h0 = 0.1d0
    real(kind=8)                   :: x0
    real(kind=8), dimension(k_max) :: x

    open(7,file='qinterp.out',status='unknown')

    x0 = 0.4d0

    write(7,*) 'Find x that satisfies:'
    write(7,*) '\frac{\mathrm{d}}{\mathrm{d}x} x^5 - \frac{4x}{5} + 0.5 = 0'
    write(7,*) "starting at x_0 = ", x0
    call qinterp(f2, x0, h0, k_max, thres)

    write(7,99) 'Found x = ', x0
    write(*,*) x0-0.63245552

99 format (a10,es13.7)

end program quadinterp
