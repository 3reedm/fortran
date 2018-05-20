! Проверка кода
program exercise1
    write(*, *) '>>>>>>Exercise 1'
    call lapack_ex1
    write(*, *) ''
    
    write(*, *) '>>>>>>Exercise 2'
    call lapack_ex2
    write(*, *) ''
    
    write(*, *) '>>>>>>Exercise 3'
    call lapack_ex3
    write(*, *) ''
    
    write(*, *) '>>>>>>Exercise 4'
    call lapack_ex4
    write(*, *) ''
end program exercise1

! Решение системы уравнений
subroutine lapack_ex1
    implicit none
    
    !---------------------------------------------------------------------------
    real, allocatable :: A(:, :), b(:)
    integer, allocatable :: pivot(:)
    integer :: n, i, j, ok
    character(len = 20) :: fmt ! переменная для задания форматирования
    
    !---------------------------------------------------------------------------
    !write(*, *) "Введите размерность СЛАУ"
    !read(*, *) n
    write(*, *) "Размерность СЛАУ = 3"
    n = 3
    allocate(A(1 : n, 1 : n))
    allocate(b(1 : n))
    allocate(pivot(1 : n))
    
    ! Заполняем числами от 0 до 1
    !call random_number(A)
    !call random_number(b)
    A(1, 1) = 1
    A(1, 2) = 0
    A(1, 3) = 0
    A(2, 1) = 0
    A(2, 2) = 0
    A(2, 3) = 3
    A(3, 1) = 0
    A(3, 2) = 4
    A(3, 3) = 0
    
    b(1) = 3
    b(2) = 2
    b(3) = 2
    
    !!
    write(fmt, *) n
    write(*, *) "----- Матрица A -----"
    write(*, "("//adjustl(fmt)//"(f8.3), t1)") ((A(i, j), j = 1, n), i = 1, n)
    write(*, *) ''
    write(*, *) "----- Вектор b -----"
    write(*, "("//adjustl(fmt)//"(f8.3), t1)") (b(i), i = 1, n)
    write(*, *) ''
    
    !!
    call SGESV(n, 1, A, n, pivot, b, n, ok)
    
    !!
    write(*, *) "--- Вектор перестановок ---"
    write(*,"("//adjustl(fmt)//"(i3), t1)") (pivot(i), i = 1, n)
    write(*, *) ''
    
    write(*, *) "--- Решение системы ---"
    write(*, "("//adjustl(fmt)//"(f8.3), t1)") (b(i), i = 1, n)
    write(*, *) ''
    
    if (ok .eq. 0) then
        write(*, *) "Успешное выполнение"
    else if(ok .gt. 0) then
        write(*, *) "Ошибка!"
    else if(ok .lt. 0) then
        write(*, *) "Неправильно задан параметр №", abs(ok)
    end if
end subroutine lapack_ex1

! LUP разложение
subroutine lapack_ex2
    implicit none
    
    !---------------------------------------------------------------------------
    real, allocatable :: A(:, :)
    integer, allocatable :: pivot(:)
    integer :: m, n, i, j, lda, info
    character(len = 20) :: fmt ! переменная для задания форматирования
    interface
        subroutine printMatrix(B)
            integer :: nRows, nColumns, k, l
            real, dimension(:, :), intent(in) :: B        
        end subroutine printMatrix
    end interface
    
    !---------------------------------------------------------------------------
    !write(*, *) "Введите размерность матрицы m и n"
    !read(*, *) m, n
    write(*, *) "Размерность матрицы 3x3"
    m = 3
    n = 3
    lda = n
    allocate(A(1 : m, 1 : n))
    allocate(pivot(1 : m))
    
    ! Заполняем числами от 0 до 1
    !call random_number(A)
    A(1, 1) = 1
    A(1, 2) = 0
    A(1, 3) = 0
    A(2, 1) = 0
    A(2, 2) = 0
    A(2, 3) = 3
    A(3, 1) = 0
    A(3, 2) = 4
    A(3, 3) = 0
    
    !!
    write(*, *) "----- Матрица A -----"
    call printMatrix(A)
    write(*, *) ''
    
    !!
    call sgetrf(m, n, A, m, pivot, info)
    
    !!
    write(*, *) "----- Матрица LU -----"
    call printMatrix(A)
    write(*, *) ''
    
    write(fmt, *) n
    write(*, *) "--- Вектор перестановок ---"
    write(*,"("//adjustl(fmt)//"(i3), t1)") (pivot(i), i = 1, n)
    write(*, *) ''
    
    if (info .eq. 0) then
        write(*, *) "Успешное выполнение"
    else if(info .gt. 0) then
        write(*, *) "Ошибка!"
    else if(info .lt. 0) then
        write(*, *) "Неправильно задан параметр №", abs(info)
    end if
end subroutine lapack_ex2

! Нахождение собственных значений
subroutine lapack_ex3
    implicit none
    
    !---------------------------------------------------------------------------
    real, allocatable :: A(:, :), Vl(:, :), Vr(:, :)
    real, allocatable :: wr(:), wi(:), work(:)
    character(len = 1) :: jobvl, jobvr
    integer :: n, i, j, lda, ldvl, ldvr, lwork, info
    character(len = 20) :: fmt ! переменная для задания форматирования
    interface
        subroutine printMatrix(B)
            integer :: nRows, nColumns, k, l
            real, dimension(:, :), intent(in) :: B        
        end subroutine printMatrix
    end interface
    
    !---------------------------------------------------------------------------
    !write(*, *) "Введите размерность квадратной матрицы n"
    !read(*, *) n
    write(*, *) "Размерность матрицы 3x3"
    n = 3
    lda = n
    ldvl = n
    ldvr = n
    jobvl = 'V'
    jobvr = 'V'
    lwork = 5 * n
    allocate(A(1 : n, 1 : n))
    allocate(Vl(1 : n, 1 : n))
    allocate(Vr(1 : n, 1 : n))
    allocate(work(1 : lwork))
    allocate(wr(1 : n))
    allocate(wi(1 : n))
    
    ! Заполняем числами от 0 до 1
    !call random_number(A)
    A(1, 1) = 1
    A(1, 2) = 0
    A(1, 3) = 0
    A(2, 1) = 0
    A(2, 2) = 0
    A(2, 3) = 3
    A(3, 1) = 0
    A(3, 2) = 4
    A(3, 3) = 0
    
    !!
    write(*, *) "----- Матрица A -----"
    call printMatrix(A)
    write(*, *) ''
    
    !!
    call sgeev(jobvl, jobvr, n, A, lda, wr, wi, Vl, ldvl, Vr, ldvr, work, lwork, info)
    
    !!    
    write(*, *) "Собственные значения"
    write(fmt, *) n
    do i = 1, n
        write(*, "(f8.3, A3, f8.3, A1)") wr(i), ' + ', wi(i), 'i'        
    end do
    write(*, *) ''
    
    write(*, *) "Левые собственные векторы"
    call printMatrix(Vl)
    write(*, *) ''
    
    write(*, *) "Правые собственные векторы"
    call printMatrix(Vr)
    write(*, *) ''
    
    if (info .eq. 0) then
        write(*, *) "Успешное выполнение"
    else if(info .gt. 0) then
        write(*, *) "Ошибка!"
    else if(info .lt. 0) then
        write(*, *) "Неправильно задан параметр №", abs(info)
    end if
end subroutine lapack_ex3

!! Использование процедуры для SVD разложения
!! действительной матрицы одинарной точности
subroutine lapack_ex4
    implicit none
    
    !---------------------------------------------------------------------------
    character(len = 1) :: jobu, jobvt
    real, allocatable :: A(:, :), U(:, :), Vt(:, :), s(:), work(:)
    integer :: m, n, lda, ldu, ldvt, lwork, info, i, j
    character(len = 20) :: fmt ! переменная для задания форматирования
    interface
        subroutine printMatrix(B)
            integer :: nRows, nColumns, k, l
            real, dimension(:, :), intent(in) :: B        
        end subroutine printMatrix
    end interface
    
    !---------------------------------------------------------------------------
    !write(*, *) "Введите размерности матрицы m и n"
    !read(*, *) m, n
    write(*, *) "Размерность матрицы 3x3"
    m = 3
    n = 3
    lda = m
    ldu = m
    ldvt = n
    jobu = 'A'
    jobvt = 'A'
    lwork = 20 * max(m, n) ! Служебная переменная, чем больше, тем лучше
    allocate(A(1 : m, 1 : n))
    allocate(U(1 : m, 1 : m))
    allocate(Vt(1 : n, 1 : n))
    
    ! Диагональные элементы сингулярной матрицы
    allocate(s(1 : min(n, m)))
    
    ! Служебный массив
    allocate(work(1 : lwork))
    
    ! Заполняем числами от 0 до 1
    !call random_number(A)
    A(1, 1) = 1
    A(1, 2) = 0
    A(1, 3) = 0
    A(2, 1) = 0
    A(2, 2) = 0
    A(2, 3) = 3
    A(3, 1) = 0
    A(3, 2) = 4
    A(3, 3) = 0
    
    !!
    write(*, *) "----- Матрица A -----"
    call printMatrix(A)
    write(*, *) ''
    
    !!
    call sgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, Vt, ldvt, work, lwork, info)
    
    !!
    write(fmt, *) min(n, m)
    write(*, *) "--- Диагональные элементы сингулярной матрицы ---"
    write(*, "("//adjustl(fmt)//"(f8.3), t1)") (s(i), i = 1, min(n, m))
    write(*, *) ''
    
    write(*, *) "----- Матрица U -----"
    call printMatrix(U)
    write(*, *) ''
    
    write(*, *) "----- Матрица Vt -----"
    call printMatrix(Vt)
    write(*, *) ''
    
    if (info .eq. 0) then
        write(*, *) "Успешное выполнение"
    else if(info .gt. 0) then
        write(*, *) "Ошибка!"
    else if(info .lt. 0) then
        write(*, *) "Неправильно задан параметр №", abs(info)
    end if
end subroutine lapack_ex4

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
