!> @brief  Interfaces module
!> @author Daniel Menendez Crespo
!> @details  It defines abstract interfaces to
!>           do callback with any one-dimensional function, ...

module interfaces

    implicit none

    abstract interface                       
        !> One dimensional functions
        pure real(kind=8) function func(x)       
            real(kind=8), intent(in) :: x  
        end function func                                                    
        !> Three dimensional functions
        pure real(kind=8) function f3D(x,y,z)   
            real(kind=8), intent(in) :: x,y,z 
        end function f3D                                                    
        !> One dimensional derivatives
        pure subroutine fp1D(f, x, h, fp)     
            procedure(func) :: f
            real(kind=8), intent(in) :: x, h   
            real(kind=8), intent(out) :: fp   
        end subroutine fp1D                                                    
    end interface
end module interfaces
