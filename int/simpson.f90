!> @brief  Trapezoid and Simpson composite rules tested
!> @author Daniel Menendez Crespo

program simpson_f

use constants, only: c => dml_const
use functions, only: f1, f2, f3
use dml_int  , only: simpson, multpts

    implicit none

    integer, parameter :: k_max=8, ounit=7
    real(kind=8), parameter :: thres=1.d-8
    !real(kind=8) :: integral


!----------------------------------------------------------
! Simpson
!----------------------------------------------------------
    open(ounit,file='simpson.out',status='unknown')

    write(ounit,*) '============='
    write(ounit,*) 'Simpson'
    write(ounit,*) '============='
    write(ounit,*) 

    write(ounit,*) '\int_0^1 x^2 e^{-x} \mathrm{d}x \approx 0.16'
    write(ounit,*) 
    call multpts(simpson, f1, 0.0d0, 1.0d0, thres, ounit)
    write(ounit,*) 
    write(ounit,*) '\int_0^{\pi/4} e^{-x} \mathrm{d}x \approx 2.65'
    write(ounit,*) 
    call multpts(simpson, f2, 0.0d0, c%pi*0.25d0, thres, ounit)
    write(ounit,*) 
    write(ounit,*) '\int_0^10 e^{-x} \mathrm{d}x \approx 1'
    write(ounit,*) 
    call multpts(simpson, f3, 0.0d0, 10.d0, thres, ounit)

!    call simpson(f1, 0.0d0, 1.0d0, 1000, integral)
!    write (8,*) 'Integral con formula simpson : ', integral
!
!    call simpson(f2, 0.0d0, c%pi*0.25d0, 1000, integral)
!    write (8,*) 'Integral con formula simpson : ', integral
!
!    call simpson(f3, 0.0d0, 10.d0, 1000, integral)
!    write (8,*) 'Integral con formula simpson : ', integral
!
!    ! integra f con simpson
!    call simpson(f, 0.0d0, 1.0d0, 1000, integral)
!    write (*,*) 'Integral con formula simpson : ', integral
!
!    ! integra g con simpson
!    call simpson(g, 0.0d0, 1.0d0, 1000, integral)
!    write (*,*) 'Integral con formula simpson : ', integral

    close(ounit)

end program simpson_f
