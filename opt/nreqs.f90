!> @brief   Newton-Raphson roots of sistem of non linear equations
!> @author  Daniel Menendez Crespo
!> @details This program solves a set of non-linear equations of 
!>          nd-functions using the Newton-Raphson method

program nr_eqs

    use dml_opt

    implicit none

    real(kind=8), parameter :: thres = 1.d-8
    integer,      parameter :: itmax=10000
    real(kind=8), parameter :: x0 = 0.1d0, y0 = 0.1d0, z0 = -0.1d0
    real(kind=8), dimension(3) :: xyz_0

    xyz_0(:) = [x0, y0, z0]

    call nreqs(xyz_0, thres, itmax)

end program nr_eqs


