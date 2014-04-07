!> @author Daniel Menendez Crespo
!> @brief
!>     This program uses 'hit and miss', 'sample mean' in 1D and 3D
!>      Monte Carlo integration to do integration. 

program mc

    use functions, only: f3, f6
    use dml_int, only: hitmiss, samplemean, samplemean3d
    
    implicit none

    real(kind=8), parameter :: ej3=1.70579270d0 !> real value for ej3
    real(kind=8), parameter :: ej5=0.999955d0 !> real value for ej5 
    real(kind=8), parameter :: a= 0.d0, b= 10.d0 !> limits for ej5
    real(kind=8), parameter :: x0= 0.d0, x1= 9.d0/10.d0 !> limit for ej3
    real(kind=8), parameter :: y0= 0.d0, y1= 1.d0 !> limit for ej3
    real(kind=8), parameter :: z0= 0.d0, z1= 11.d0/10.d0 !> limit for ej3
    integer      :: n, k
    real(kind=8) :: h !> height of the box
    real(kind=8), dimension(6) :: hm, sm, sm3d !> values obtained in each method

    open(unit=7, file='mc.out', status='unknown')
    write(7,*)'HM=hit and miss'
    write(7,*)'SM=sample mean'
    write(7,*) 
    write(7,*) 'Mono-dimensional'
    write(7,*) 
    write(7,*) '\int_0^10 e^{-x} \mathrm{d}x \approx 1'
    write(7,*) 
    write(7,*)'        n         HM              SM             err(HM)         err(SM)'

    n = 10
    hm = 0.d0
    sm = 0.d0

    do k=1, 6
        h = 1.d0
        call hitmiss(f3, a, b, h, n, hm(k))
        call samplemean(f3, a, b, h, n, sm(k))

        write(7,99) n, hm(k), sm(k), abs(hm(k)-ej5), abs(sm(k)-ej5)

        n = n*10
    end do

    write(7,*) 
    write(7,*) 'Multi-dimensional'
    write(7,*) 
    write(7,*) '\int_0^{9/10}\int_0^1\int_0^{11/10} \sqrt{4 - x^2 - y^2 - z^2} &
                \mathrm{d}x \, \mathrm{d}y \, \mathrm{d}z \approx 1.70'
    write(7,*) 
    write(7,*)'        n         SM              err(SM)'
    n = 10

    do k=1, 6
        call samplemean3d(f6, x0, x1, y0, y1, z0, z1, n, sm3d(k))
        write(7,99) n, sm3d(k), abs(sm3d(k)-ej3)
        n = n*10
    end do

99   format(i10, 4f16.8)

end program mc
