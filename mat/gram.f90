!---------------------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief Gram-Schmidt orthogonalization
!---------------------------------------------------------------------
program gram

    use dml_mat

    implicit none

    integer, parameter           :: n = 4
    real(kind=8), parameter      :: thres = 1.d-13
    real(kind=8), dimension(n,n) :: a, v
    integer                      :: i,j,k
    real(kind=8)                 :: l, xdia, xoff

    v = 0.d0
    l = 0.d0
    a = reshape([1.d0,2.d0,3.d0,4.d0,&
                 4.d0,3.d0,2.d0,1.d0,&
                 4.d0,2.d0,3.d0,1.d0,&
                 1.d0,2.d0,-6.d0,1.d0],shape(a))
    !a = reshape([1.d0,2.d0,3.d0,4.d0,&
    !             5.d0,6.d0,7.d0,8.d0,&
    !             9.d0,10.d0,11.d0,12.d0],shape(a))

    do i=1,n
        print*, a(:,i)
    end do

    call gram_schmidt(a,v)

    print*, 'Solution'
    do i=1,n
        print*, v(:,i)
    end do

    print*, 'Check solution'

    do j=1, n
        do k=1,j-1
            l = dot_product(v(:,j),v(:,k))
            if (abs(l) .gt. xoff) xoff = abs(l)
        end do
        l = 0.d0
        !l = norm2(v(:,j))**2
        do i=1, n
            l = l + v(i,j)**2
        end do
        l = l - 1.d0
        if (abs(l) .gt. xdia) xdia = abs(l)
    end do

    print*, 'Deviation in diagonal elements: ', xdia
    print*, 'Deviation in non-diagonal elements: ', xoff


    stop

end program gram
