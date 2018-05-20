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

program loop
    use exercise2
    implicit none
    
    ! Раздел описания переменных
    real, dimension(1 : 4, 1 : 4) :: a
    
    ! Раздел операторов
    data a/1, 0, 25, 0,    57, 1, 133, 0,    0, 0, 1, 0,    0, 0, 0, 1/
    call printMatrix(a)
    pause
end program loop