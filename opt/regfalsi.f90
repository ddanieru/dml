!> @file: regfalsi.f90
!> @brief  Regula falsi method test
!> @author Daniel Menendez Crespo

program regulafalsi

    use functions
    use dml_opt

    implicit none

    integer, parameter :: m=60, outunit=7
    real(kind=8) :: fa = 0.0d0, fb = 1.0d0 
    real(kind=8) :: f3a = -2.d0, f3b = 2.0d0 
    real(kind=8), parameter :: thres = 1.0d-5

    open(outunit,file='regfalsi.out',status='unknown')

    write(outunit,*) 'Find x that satisfies \exp(x) - x^2 + 3x - 2 = 0'
    write(outunit,*) "for x \in [", fa, ",", fb, "]"
    call regfalsi(f1, fa, fb, m, thres, outunit)
    write(outunit,*) 
    write(outunit,*) 'Find x that satisfies \exp(-x) - 2 = 0'
    write(outunit,*) "for x \in [", f3a, ",", f3b, "]"
    call regfalsi(f3, f3a, f3b, m, thres, outunit)


end program regulafalsi
