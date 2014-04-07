!> @brief   Newton-Raphson minimization of multidimensional functions
!> @author  Daniel Menendez Crespo
!> @details This program minimizes a function using a 
!>          Newton-Raphson method            

program nr_min

    use functions
    use dml_opt

    implicit none

    integer, parameter :: dim=2,maxit=30
    real(kind=8), dimension(dim, 0:maxit) :: x

    x(:,0) = [0.d0, 0.d0]

    call nrmin(x,f6,maxit)  

    call nrmin(x,f7,maxit)

end program nr_min
