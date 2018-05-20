! Модуль, печатающий любую двумерную матрицу
module exercise2
    implicit none
    contains
    subroutine printMatrix(a) 
        ! Входной параметр
        integer :: nRows, nColumns, i, j
        real, dimension(:, :), intent(in) :: a 
        
        nRows = size(a, 1)
        nColumns = size(a, 2)
        
        do i = 1, nRows
            print *, (a(i, j), j = 1, nColumns)
            print *, '' 
        end do  
    end subroutine printMatrix
end module exercise2