! Параллельная версия умножения матриц с помощью стандартного метода, метода Винограда. Сравнение их производительности с производительностью функции matmul
program exercise2
    implicit none
    include "omp_lib.h" 
    real, allocatable :: A(:, :), B(:, :), C(:, :)
    integer :: m, n, k, th_num, t1, t2
    real :: t3
    interface
        subroutine printMatrix(B)
            real, dimension(:, :), intent(inout) :: B 
        end subroutine printMatrix
    end interface
    
    write(*,  *) "                          Умножение матриц"
    write(*,  *) ""
    write(*,  *) "        Введите размерности матриц (A = m * n, B = n * k, C = m * k)"
    write(*,  *) ""
    write(*,  *) "Введите m ="   
    read(*, *) m
    write(*,  *) "Введите n ="
    read(*, *) n
    write(*,  *) "Введите k ="
    read(*, *) k
    write(*,  *) ""
    write(*, *) "Введите th_num (количество нитей) =" 
    read(*, *) th_num
    
    ! Выделение памяти под массивы
    allocate(A(1 : m, 1 : n), B(1 : n, 1 : k), C(1 : m, 1 : k))
    
    ! Заполнение массивов A и B (числами от 0 до 1)
    call random_number(A)
    call random_number(B)
    
    ! Печать матриц (A, B)
    !print *, ""
    !print *, "Матрица A"
    !call printMatrix(A)
    !print *, ""
    !print *, "Матрица B"
    !call printMatrix(B)
    print *, ""
    print *, ""
    
    print *, "                  Проверка производительности алгоритмов "
    print *, ""
    
    print *, "Стандартный алгоритм умножения матриц"
    call system_clock(t1)
    call stdMM(A, B, C, m, k, n)
    call system_clock(t2)
    !call printMatrix(C)
    t3 = real(t2 - t1) / 1000
    write(*, '(A, f8.3)') " последовательное время (секунды): ", t3
    call system_clock(t1)
    call stdPMM(A, B, C, m, k, n, th_num)
    call system_clock(t2)
    !call printMatrix(C)
    t3 = real(t2 - t1) / 1000 - t3
    write(*, '(A, f8.3)') " параллельное время (секунды): ", real(t2 - t1) / 1000
    write(*, '(A, f8.3)') " выигрыш параллельного варианта (секунды): ", -t3
    print *, ""
    
    print *, "Алгоритм Винограда умножения матриц"
    call system_clock(t1)
    call winogradMM(A, B, C, m, k, n)
    call system_clock(t2)
    t3 = real(t2 - t1) / 1000
    !call printMatrix(C)
    write(*, '(A, f8.3)') " последовательное время (секунды): ", real(t2 - t1) / 1000
    call system_clock(t1)
    call winogradPMM(A, B, C, m, k, n, th_num)
    call system_clock(t2)
    t3 = real(t2 - t1) / 1000 - t3
    !call printMatrix(C)
    write(*, '(A, f8.3)') " параллельное время (секунды): ", real(t2 - t1) / 1000
    write(*, '(A, f8.3)') " выигрыш параллельного варианта (секунды): ", -t3
    print *, ""
    
    print *, "Использование функции matmul для умножения матриц"
    call system_clock(t1)
    C = matmul(A, B)
    call system_clock(t2)
    !call printMatrix(C)
    write(*, '(A, f8.3)') " время (секунды): ", real(t2 - t1) / 1000
    print *, ""
end program exercise2

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

subroutine stdPMM(A, B, C, nRows, nColumns, nAB, th_num)
    ! Описание переменных
    integer, intent(in) :: nRows, nColumns, nAB, th_num
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: A
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: B
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: C
    integer :: i, j
    
    ! Действия 
    !$omp parallel do num_threads(th_num)
        do i = 1, nRows
            do j = 1, nColumns
                C(i, j) = A(i, 1) * B(1, j)
                do l = 2, nAB
                    C(i, j) = C(i, j) + A(i, l) * B(l, j)
                end do
            end do
        end do
    !$omp end parallel do
end subroutine stdPMM

subroutine winogradPMM(A, B, C, nRows, nColumns, nAB, th_num)
    ! Описание переменных
    integer, intent(in) :: nRows, nColumns, nAB, th_num
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: A
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: B
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: C
    real, dimension(1 : nRows) :: rowFactor
    real, dimension(1 : nColumns) :: columnFactor
    integer :: l, i, j, d 
    
    ! Действия   
    d = nAB / 2
    
    !$omp parallel num_threads(th_num)
        ! Вычисление rowFactor для матрицы A
        !$omp do
            do i = 1, nRows
                rowFactor(i) = A(i, 1) * A(i, 2)
                do j = 2, d
                    rowFactor(i) = rowFactor(i) + A(i, 2 * j - 1) * A(i, 2 * j)
                end do
            end do
        !$omp end do
    
        ! Вычисление columnFactor для матрицы B
        !$omp do
            do i = 1, nColumns
                columnFactor(i) = B(1, i) * B(2, i)
                do j = 2, d
                    columnFactor(i) = columnFactor(i) + B(2 * j - 1, i) * B(2 * j, i)
                end do
            end do
        !$omp end do
        
        ! Вычисление матрицы C
        !$omp do
            do i = 1, nRows
                do j = 1, nColumns
                    C(i, j) = - rowFactor(i) - columnFactor(j)
                    do l = 1, d
                        C(i, j) = C(i, j) + (A(i, 2 * l - 1) + B(2 * l, j)) * (A(i, 2 * l) + B(2 * l - 1, j))
                    end do
                end do
            end do
        !$omp end do
        
        ! Прибавление членов в случае нечётной разметки
        if (mod(nAB, 2) .eq. 1) then
            !$omp do
                do i = 1, nRows
                    do j = 1, nColumns
                        C(i, j) = C(i, j) + A(i, nAB) * B(nAB, j)
                    end do
                end do
            !$omp end do
        end if
    !$omp end parallel
end subroutine winogradPMM

! Подпрограмма для печати любой двумерной матрицы
subroutine printMatrix(A) 
    integer :: nRows, nColumns, i, j
    real, dimension(:, :), intent(inout) :: A 
    character(len = 20) :: fmt
       
    nRows = size(A, 1)
    nColumns = size(A, 2)
    
    do i = 1, nRows
        do j = 1, nColumns
            if (A(i, j) .eq. 0) then
                A(i, j) = abs(A(i, j))
            end if
        end do
    end do 
    
    do i = 1, nRows
        write(fmt, *) nColumns
        write(*, "("//adjustl(fmt)//"(f8.3), t1)") (A(i, j), j = 1, nColumns)
    end do  
end subroutine printMatrix