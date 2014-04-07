!> @brief  Romberg integration method
!> @author Daniel Menendez Crespo

program romberg_f

use constants , only: c => dml_const
use functions , only: f4 , f5
use dml_int   , only: romberg

    implicit none

    integer      :: i,j, maxi
    real(kind=8) :: a, b, thres
    real(kind=8) :: integral

    real(kind=8), dimension(:,:), allocatable :: r


!----------------------------------------------------------
! Romberg
!----------------------------------------------------------
    open(7,file='romberg.out',status='unknown')

    write(7,*) '============='
    write(7,*) 'Romberg'
    write(7,*) '============='

    a = 0.0d0
    b = 1.0d0
    maxi = 20
    thres = 1.0e-10
!    write (*,*) 'Integracion en el intervalo', a, b
!    print*, 'de la funcion error'
!    print*, 'Con no maximo de puntos:', maxi
!    print '(a20, es20.10)', 'Con error tolerado:', thres
!    print '(a38, es20.10)', 'Valor dado por la funcion intrinseca:', erf(b)
!!
!    print *
!    print *,' Extrapolacion '
!    print *
!    print 99,'n','R(n,0)','R(n,1)','R(n,2)','R(n,3)','R(n,4)'
!    print *
!!
!    call romberg(f_3, a, b, maxi, thres, integral, r)
!!
!    do i=1, size(r, dim=1) 
!        print 999, i, (r(i,j), j=1, i)
!    end do
!    print '(A6, 2X, es20.10, 2X, A2, 2X, i5, 2X, A10)', &
!          'Result:', integral, 'in', size(r, dim=1) - 1, 'iterations'
!    print '(A2, 1X, i5, 1X, 1A, 1X, i5, A4, 2X, es20.10)', &
!          'R(', size(r, dim=1), ',', size(r, dim=1), ')  =', integral
!
!
!    call romberg(f_3, a, b, maxi, thres, integral)
!    print '(A6, 2X, es20.10)', 'Result', integral
!
!    call romberg(func4, a, b, maxi, thres, integral)
!    print '(A6, 2X, es20.10)', 'Result', integral

    write(7,*)
    write(7,*) ' Extrapolation '
    write(7,*)
    write(7,*) '\int_0^1 \frac{x^2}{1 + x^3} \mathrm{d}x'
    write(7,*) 
    write(7,99) 'n','R(n,0)','R(n,1)','R(n,2)','R(n,3)','R(n,4)'
    write(7,*)

    call romberg(f4, a, b, maxi, thres, integral, r)
    write(7,*) 

    do i=1, size(r, dim=1) 
        write(7,999) i, (r(i,j), j=1, i)
    end do
    write(7,'(A6, 2X, es20.10, 2X, A2, 2X, i5, 2X, A10)') &
          'Result:', integral, 'in', size(r, dim=1) - 1, 'iterations'
    write(7,'(A2, 1X, i5, 1X, 1A, 1X, i5, A4, 2X, es20.10)') &
          'R(', size(r, dim=1), ',', size(r, dim=1), ')  =', integral

    write(7,*)
    write(7,*) '\int_0^{\pi/4} x^2 \sin(4x) \mathrm{d}x'
    write(7,*) 
    write(7,99) 'n','R(n,0)','R(n,1)','R(n,2)','R(n,3)','R(n,4)'
    write(7,*)

    write(7,*) '\int_0^{\pi/4} x^2 \sin(4x) \mathrm{d}x'
    write(7,*) 
    call romberg(f5, a, c%pi*0.25d0, maxi, thres, integral, r)

    do i=1, size(r, dim=1) 
        write(7,999) i, (r(i,j), j=1, i)
    end do
    write(7,'(A6, 2X, es20.10, 2X, A2, 2X, i5, 2X, A10)') &
          'Result:', integral, 'in', size(r, dim=1) - 1, 'iterations'
    write(7,'(A2, 1X, i5, 1X, 1A, 1X, i5, A4, 2X, es20.10)') &
          'R(', size(r, dim=1), ',', size(r, dim=1), ')  =', integral


    if (allocated(r)) then
        deallocate(r)
    end if

!     los formatos
99 format(a6,a13,4(a22))
!999  format(1x,i5,2x,8(e13.6,2x))
999 format(1x,i5,2x,8(es20.10,2x))
!999  format(1x,i5,2x,8(e16.9,2x))

    stop

end program romberg_f

