!---------------------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief Gaussian elimination with partial pivoting
!---------------------------------------------------------------------
program gauss

use dml_mat

    implicit none

    real(kind=8), dimension(:,:), allocatable :: a
    integer      :: i, j, n

    n = 3
    allocate(a(n,n+1))
    a = reshape([2.d0, -3.d0, -2.d0, &
                 1.d0, -1.d0, 1.d0,  &
                 -1.d0, 2.d0, 2.d0,  &
                 8.d0, -11.d0, -3.d0], shape(a))

    write(6,*) 'Initial expanded matrix'
    do i=1, n
       write(6,'(13f13.7)') (a(i,j),j=1,n+1)
    end do
    call gaussel(a)

    write(6,*) 'U matrix'
    do i=1, n
        write(6,'(13f13.7)') (a(i,j),j=1,n+1)
    end do

    deallocate(a)


end program gauss
