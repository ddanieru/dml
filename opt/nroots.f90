!> @brief   Newton-Raphson roots of onedimensional functions
!> @author  Daniel Menendez Crespo
!> @details This program finds the roots of  1D-functions using the 
!>          Newton-Raphson method

program nroots

    use functions
    use dml_opt

    implicit none
    real(kind=8), parameter :: start = 0.4d0

    call nroots1D(start, f4, derivative_1d)

end program nroots
