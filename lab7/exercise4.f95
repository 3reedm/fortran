! Вычисление определённого интеграла методом трапеций
program exercise4
    implicit none
    
    !-------------------
    real, dimension(1 : 2, 1 : 2) :: A 
    real, dimension(1 : 2) :: F, X
    interface
        recursive subroutine jacobi(A, F, X, n)
            real, dimension(:,:), intent(in) :: A
            real, dimension(:), intent(in) :: F
            integer, intent(in) :: n
            real, dimension(:), intent(inout) :: X
        end subroutine jacobi
        subroutine printMatrix(a) 
            real, dimension(:, :), intent(inout) :: a  
        end subroutine printMatrix
    end interface
    
    data A/2, 5,  1, 7/, F/11, 13/, X/1, 1/
    
    !-------------------
    write(*, *) '>>>>>>Решение СЛАУ методом Якоби'
    write(*, *) ''
    
    write(*, *) '>>>>>>Матрица коэффициентов: '
    call printMatrix(A)
    write(*, *) ''
    write(*, *) '>>>>>>Столбец свободных членов: '
    write(*, "(f5.2)") F
    write(*, *) ''
    
    call jacobi(A,F,X,2)
    write(*, *) '>>>>>>Вектор решений: '
    write(*, "(f5.2)") X
end program exercise4

! Метод Якоби
subroutine jacobi(A, F, X, n)
    implicit none
    
    !---------------------------------------------------------------------------
    real, dimension(:,:), intent(in) :: A
    real, dimension(:), intent(in) :: F
    integer, intent(in) :: n
    real, dimension(:), intent(inout) :: X
    real, dimension(:), allocatable :: T
    integer :: i,j
    real :: norm
    real, parameter :: eps = 0.01

    allocate(T(1:n))
    
    norm = 0.1
	do while (norm .gt. eps)
		do i=1,n
			T(i) = F(i)
			do j=1,n 
				if (i .ne. j) then
					T(i) = T(i) - A(i,j)*X(j)
                end if
			end do
			T(i) = T(i) / A(i,i)
		end do
        
        norm = abs(X(1)-T(1))
		do i=1,n 
			if (abs(X(i)-T(i)) .gt. norm) then
				norm = abs(X(i)-T(i))
            end if
			X(i) = T(i)
		end do
	end do

    deallocate(T)
    return
end subroutine jacobi

! Подпрограмма, печатающая любую двумерную матрицу
subroutine printMatrix(a) 
    ! Входной параметр
    integer :: nRows, nColumns, i, j
    real, dimension(:, :), intent(inout) :: a 
    character(len = 20) :: fmt ! переменная для задания форматирования
       
    nRows = size(a, 1)
    nColumns = size(a, 2)
    
    do i = 1, nRows
        do j = 1, nColumns
            if(a(i, j) .eq. 0) then
                a(i, j) = abs(a(i, j))
            end if
        end do
    end do 
    
    do i = 1, nRows
        write(fmt, *) nColumns
        write(*, "("//adjustl(fmt)//"(f8.3), t1)") (a(i, j), j = 1, nColumns)
    end do  
end subroutine printMatrix
