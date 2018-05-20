! Умножение матриц с помощью стандартного метода, метода Винограда, функции matmul. Сравнение их производительности
program exercise3
    use exercise2  ! Из предыдущей задачи (для печати матрицы)
    implicit none
    
    ! Раздел описания переменных
    real, dimension(1 : 3, 1 : 2) :: A 
    real, dimension(1 : 2, 1 : 3) :: B
    real, dimension(1 : 3, 1 : 3) :: C
    ! integer :: count, count1
    
    ! Раздел операторов
    ! call random_number(A)
    ! call random_number(B)
    data A/1, 4, 5,   6, 9, 3/, B/2, 7,   3, 5,   6, 1/
    
    ! print *, "Проверка производительности алгоритмов для больших матриц (1000 x 1000)"
    ! print *, ""
    print *, "Стандартный алгоритм умножения матриц"
    print *, ""
    ! call system_clock(count)
    call stdMM(A, B, C, 3, 3, 2)
    ! call system_clock(count1)
    call printMatrix(C)
    print *, ""
    ! print *, "time: ", count1 - count
    ! print *, ""
    
    print *, "Алгоритм Винограда умножения матриц"
    print *, ""
    ! call system_clock(count)
    call winogradMM(A, B, C, 3, 3, 2)
    ! call system_clock(count1)
    call printMatrix(C)
    print *, ""
    ! print *, "time: ", count1 - count
    ! print *, ""
    
    print *, "Использование функции matmul для умножения матриц"
    print *, ""
    ! call system_clock(count)
    C = matmul(A, B)
    ! call system_clock(count1)
    call printMatrix(C)
    ! print *, ""
    ! print *, "time: ", count1 - count
    
    pause
end program exercise3

subroutine stdMM(A, B, C, nRows, nColumns, nAB)
    ! Описание переменных
    integer, intent(in) :: nRows, nColumns, nAB
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: A
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: B
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: C
    integer :: l, i, j
    
    ! Действия 
    do i = 1, nRows
        do j = 1, nColumns
            C(i, j) = A(i, 1) * B(1, j)
            do l = 2, nAB
                C(i, j) = C(i, j) + A(i, l) * B(l, j)
            end do
        end do
    end do
end subroutine stdMM

subroutine winogradMM(A, B, C, nRows, nColumns, nAB)
    ! Описание переменных
    integer, intent(in) :: nRows, nColumns, nAB
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: A
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: B
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: C
    real, dimension(1 : nRows) :: rowFactor
    real, dimension(1 : nColumns) :: columnFactor
    integer :: l, i, j, d 
    
    ! Действия   
    d = nAB / 2
    
    ! Вычисление rowFactor для матрицы A
    do i = 1, nRows
        rowFactor(i) = A(i, 1) * A(i, 2)
        do j = 2, d
            rowFactor(i) = rowFactor(i) + A(i, 2 * j - 1) * A(i, 2 * j)
        end do
    end do
    
    ! Вычисление columnFactor для матрицы B
    do i = 1, nColumns
        columnFactor(i) = B(1, i) * B(2, i)
        do j = 2, d
            columnFactor(i) = columnFactor(i) + B(2 * j - 1, i) * B(2 * j, i)
        end do
    end do
    
    ! Вычисление матрицы C 
    do i = 1, nRows
        do j = 1, nColumns
            C(i, j) = - rowFactor(i) - columnFactor(j)
            do l = 1, d
                C(i, j) = C(i, j) + (A(i, 2 * l - 1) + B(2 * l, j)) * (A(i, 2 * l) + B(2 * l - 1, j))
            end do
        end do
    end do
    
    ! Прибавление членов в случае нечётной разметки
    if (mod(nAB, 2) .eq. 1) then
        do i = 1, nRows
            do j = 1, nColumns
                C(i, j) = C(i, j) + A(i, nAB) * B(nAB, j)
            end do
        end do
    end if
end subroutine winogradMM

