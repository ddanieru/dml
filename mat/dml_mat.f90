
!> @brief Numerical Linear Algebra procedures
!> @author Daniel Menendez Crespo

module dml_mat

    implicit none

!==========================================================
contains
!==========================================================

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  Gram-Schmidt orthogonalization
!----------------------------------------------------------
subroutine gram_schmidt(a, v)

    real(kind=8), dimension(:,:), intent(in)  :: a
    real(kind=8), dimension(:,:), intent(out) :: v
    real(kind=8), parameter                   :: thres = 1.d-13
    integer                                   :: i,j,k,n
    real(kind=8)                              :: l

    n = size(a, dim=1)

    v(:,1) = a(:,1)/sqrt(dot_product(a(:,1),a(:,1)))
    print*, dot_product(a(:,1),a(:,1))

    do j= 2, n
        v(:,j) = a(:,j)
        do k=1, j-1
            l = dot_product(a(:,j),v(:,k))
            v(:,j) = v(:,j) - l*v(:,k)
        end do
        v(:,j) = v(:,j)/sqrt(dot_product(v(:,j),v(:,j)))
        if (sqrt(dot_product(v(:,j),v(:,j))) .le. thres) then
            write(*,*) 'Error: Linear dependency found.'
            return
        end if
    end do

end subroutine gram_schmidt

!----------------------------------------------------------
!> @author Daniel Menendez Crespo
!> @brief  Gaussian elimination with partial pivoting
!----------------------------------------------------------
subroutine gaussel(a)

    real(kind=8), dimension(:,:), allocatable, intent(inout) :: a
    real(kind=8), dimension(:),   allocatable :: x
    real(kind=8), parameter            :: thres = 1.0d-13

    real(kind=8) :: m, piv, copy
    integer      :: i, j, k, kpiv, n

    n = size(a, dim=1)
    allocate(x(n))

    do k=1, n-1

!   Choose pivot
        piv = a(k,k)
        kpiv  = k
        do j = k+1,n
            if (abs(piv).ge.abs(a(j,k))) cycle
            piv = a(j,k)
            kpiv  = j
        end do
      
        if (abs(piv).le.thres) then
            write(*,*) 'Error: Singular matrix given'
            return
        end if
      
!   Swap rows
        if (kpiv.ne.k) then
            do j=k, n+1
                copy = a(kpiv,j)
                a(kpiv,j) = a(k,j)
                a(k,j) = copy 
            end do
        end if
        kpiv = k
      
!   Eliminate following (n-k) equations
        do i=k+1, n
            m = a(i,k) / piv ! calculate factor
            do j=k, n+1
                a(i,j) = a(i,j) - m*a(kpiv,j) ! Update matrix
            end do
        end do
    end do

    deallocate(x)

end subroutine gaussel

!----------------------------------------------------------
!> @author  Daniel Menendez Crespo
!> @brief   Inversion of a matrix by Gauss-Jordan method
!> @details Subroutine that calculates the inverse of a matrix 
!>          within the Gauss-Jordan method
!----------------------------------------------------------
subroutine inverse(mat)

    real(kind=8), dimension(:,:), intent(inout) :: mat
    real(kind=8), dimension(:,:), allocatable :: a
    logical, dimension(:), allocatable :: done

    real(kind=8), parameter :: thres = 1.0d-13
    real(kind=8) :: max_dia, m
    integer :: i, j, k, l, k_max, n

    n = size(mat, dim=1)
    allocate(a(n,2*n))
    allocate(done(n))

    a(:,:n) = mat

    done = .false.

! mat = ones()
    mat = 0.d0
    do i=1, n
        mat(i,i) = 1.d0
    end do

! Extend the matrix
    a(:,(n+1):) = mat

    do k=1, n

!   Choose maximum diagonal
        max_dia = 0.0d0
        do l=1, n
            if (done(l)) cycle
            if (abs(a(l,l)).gt.max_dia) then
                k_max = l
                max_dia = abs(a(l,l))
            end if
        end do

        if (max_dia.le.thres) then
            write(*,*) 'Error: singular matrix given'
            return 
        end if

!   Scale k_max row
        m = a(k_max,k_max)
        do j=1, 2*n
            a(k_max,j) = a(k_max,j) / m
        end do

!   Treat rows before
        do i=1, k_max-1
            m = a(i,k_max)
            do j=1, 2*n
                a(i,j) = a(i,j) - m*a(k_max,j)
            end do
        end do

!   Treat rows after
        do i=k_max+1, n
            m = a(i,k_max)
            do j=1, 2*n
                a(i,j) = a(i,j) - m*a(k_max,j)
            end do
        end do

!   Mark kmax as done
        done(k_max) = .true.
    end do

    mat = a(:n,n+1:)


end subroutine inverse

!----------------------------------------------------------
!> @author  Daniel Menendez Crespo
!> @brief   Eigenvalues and eigenvectors by Cyclic Jacobian
!> @details 'Matrix Computations' Golub & van Loan 1996 section 8.4.4
!> @param[in,out] a  real symmetric matrix that is returned diagonalized
!> @param[in] n     size of the matrix
!> @param[in] thres square of non diagonal elements 
!> @param[out] v    eigenvectors
!----------------------------------------------------------
subroutine jacobi(a,n,v,thres)

    implicit none

    integer, parameter :: maxit = 200
    integer :: p, q, r, n, it
    real(kind=8), dimension(n,n) :: a, v
    real(kind=8) :: thres, off, eps
    real(kind=8) :: tau,t, c, s, cc, ss, rho

! Equivalent to v = ones()
    v = 0.0d0
    do p=1, n
        v(p,p) = 1.0d0
    end do

! Square of sum of  all off-diagonal elements (squared)
! The idea behind Jacobi method is to systematically reduce this
    off = 0.0d0
    do p=1, n
        do q=1, n
            if (p.ne.q) then
                off = off + a(p,q)**2
            end if
        end do
    end do

! If it is already diagonal
    if (off.le.thres) return

! half of average for off-diagonal elements 
    eps = 0.1d0*off/real(n*n)

    it = 1
    do while (off.gt.thres)
        do p=1, n-1
            do q=p+1, n
                if (a(q,p)**2.le.eps) cycle  ! skip small ones
                off = off - 2.0d0*a(q,p)**2 ! (8.4.2) 'Matrix Computations'
                eps = 0.1d0*off/real(n*n)

            ! Calculate coefficient c and s for Jacobi/Givens rotations
            ! 2-by-2 Symmetric Schur Decomposition

                tau = (a(p,p)-a(q,q)) / (2.d0*a(q,p))
                t = sign(1.d0,tau) / (abs(tau) + sqrt(1.d0+tau**2))
                c = 1.d0/sqrt(1.d0 + t**2)
                s = t*c

            ! Refresh p and q rows

                do r=1, n
                    cc = c*a(p,r) + s*a(q,r)
                    ss = c*a(q,r) - s*a(p,r)
                    a(p,r) = cc
                    a(q,r) = ss
                end do

            ! Refresh p and q columns

                do r=1, n
                    cc = c*a(r,p) + s*a(r,q)
                    ss = c*a(r,q) - s*a(r,p) 
                    a(r,p) = cc
                    a(r,q) = ss
                end do
!                rho = s / (1.d0 + c)
!                do r=1, p-1
!                   a(r,p) = a(r,p) - s * (a(r,q) + rho * a(r,p))
!                   a(p,r) = a(r,p)
!                end do
!                do r=p+1, n
!                   a(r,p) = a(r,p) - s * (a(r,q) + rho * a(r,p))
!                   a(p,r) = a(r,p)
!                end do
!
!                do r=1, q-1
!                   a(r,q) = a(r,q) + s * (a(r,p) - rho * a(r,q))
!                   a(q,r) = a(r,q)
!                end do
!                do r=q+1, n
!                   a(r,q) = a(r,q) + s * (a(r,p) - rho * a(r,q))
!                   a(q,r) = a(r,q)
!                end do
!
!                a(p,p) = a(p,p) - t*a(p,q)
!                a(q,q) = a(q,q) + t*a(p,q)
!                a(p,q) = 0.0d0
!                a(q,p) = 0.0d0

            ! Update eigenvectors matrix

                do r=1, n
                    cc =  c*v(r,p) + s*v(r,q)
                    ss =  c*v(r,q) - s*v(r,p)
                    v(r,p) = cc
                    v(r,q) = ss
                end do
            end do
        end do
        if (it.eq.maxit) then
            write(*,*) 'Maximum number of iterations reached'
            return
        end if
        it = it + 1
    end do

end subroutine jacobi

end module dml_mat

