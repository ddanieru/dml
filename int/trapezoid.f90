!> @brief  Trapezoid and Simpson composite rules tested
!> @author Daniel Menendez Crespo

program trap_simpson

use constants, only: c => dml_const
use functions, only: f1, f2, f3
use dml_int  , only: trapecio, multpts

    implicit none

    integer, parameter :: k_max=8, ounit=7
    real(kind=8), parameter :: thres=1.d-8
    !real(kind=8) :: integral


!----------------------------------------------------------
! Trapecio
!----------------------------------------------------------
    open(ounit,file='trap.out',status='unknown')

    write(ounit,*) '============='
    write(ounit,*) 'Trapezoid'
    write(ounit,*) '============='
    write(ounit,*) 

    write(ounit,*) '\int_0^1 x^2 e^{-x} \mathrm{d}x \approx 0.16'
    write(ounit,*) 
    call multpts(trapecio, f1, 0.0d0, 1.0d0, thres, ounit)
    write(ounit,*) 
    write(ounit,*) '\int_0^{\pi/4} e^{-x} \mathrm{d}x \approx 2.65'
    write(ounit,*) 
    call multpts(trapecio, f2, 0.0d0, c%pi*0.25d0, thres, ounit)
    write(ounit,*) 
    write(ounit,*) '\int_0^10 e^{-x} \mathrm{d}x \approx 1'
    write(ounit,*) 
    call multpts(trapecio, f3, 0.0d0, 10.d0, thres, ounit)

!    call trapecio(f1, 0.0d0, 1.0d0, 1000, integral)
!    write (7,*) 'Integral con formula trapecio : ', integral
!
!    call trapecio(f2, 0.0d0, c%pi*0.25d0, 1000, integral)
!    write (7,*) 'Integral con formula trapecio : ', integral
!
!    call trapecio(f3, 0.0d0, 10.d0, 1000, integral)
!    write (7,*) 'Integral con formula trapecio : ', integral
!
!    ! integra f con trapecio
!    call trapecio(f, 0.0d0, 1.0d0, 1000, integral)
!    write (*,*) 'Integral con formula trapecio : ', integral
!
!    ! integra g con trapecio
!    call trapecio(g, 0.0d0, 1.0d0, 1000, integral)
!    write (*,*) 'Integral con formula trapecio : ', integral


    close(ounit)

end program trap_simpson
