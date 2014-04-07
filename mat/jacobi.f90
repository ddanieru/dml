!====================================================================
!> @author
!> Daniel Menendez Crespo
!> @brief
!> Jacobi method 
!> @details 
!> eigenvalues and eigenvectors of a real symmetric matrix
!====================================================================
program main

use dml_mat

    implicit none

    integer, parameter :: n=2
    real(kind=8) :: a(n,n), x(n,n)
    real(kind=8), parameter:: thres=1.0d-13
    integer i, j

    !a(:,1) = [1,2,3]
    !a(:,2) = [2,2,-2]
    !a(:,3) = [3,-2,4]
    a(:,1) = [2,4]
    a(:,2) = [4,5]

    ! print a header and the original matrix
    write (*,*) ' Cyclic Jacobi method for matrix:'
    do i=1,n
        write (*,99) (a(i,j),j=1,n)
    end do

    call Jacobi(a,n,x,thres)

    ! print solutions
    write (*,*) 
    write (*,*) ' Eigenvalues'
    write (*,99) (a(i,i),i=1,n)
    write (*,*) 
    write (*,*) ' Eigenvectors'
    do i=1, n
        write (*,99)  (x(i,j),j=1,n)
    end do

99 format (x, 6f14.8)

end program main

