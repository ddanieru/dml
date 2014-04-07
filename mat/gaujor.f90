!----------------------------------------------------------
!> @author  Daniel Menendez Crespo
!> @brief   Inversion of a matrix by Gauss-Jordan method
!> @details Subroutine that calculates the inverse of a matrix 
!>          within the Gauss-Jordan method
!----------------------------------------------------------
program gaujor

use dml_mat

    implicit none

    integer, parameter :: maxdim = 77
    real(kind=8), dimension(:,:), allocatable :: a
    integer :: n, i, j

    n = 3
    allocate(a(n,n))
    a = reshape([2.d0, -1.d0, 0.d0, &
                 -1.d0, 2.d0,-1.d0, &
                 0.d0, -1.d0, 2.d0], [n,n])

    call inverse(a)

    write(6,*) 'Inverse matrix'
    do i=1, n
        write(6,'(13f13.7)') (a(i,j),j=1,n)
    end do

    deallocate(a)


end program gaujor
