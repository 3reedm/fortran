! Нахождение детерминанта матрицы
program exercise2
    implicit none
    
    integer :: timeStart, timeStop, n, i, j
    real, allocatable :: A(:, :)
    character(len = 10) :: fmt
    
    write(*, *) '      Вычисление определителя матрицы (не учитывая погрешность)'
    write(*, *) ''
    
    write(*, *) "Введите размерность квадратной матрицы [n]"
    read(*, *) n
    write(*, *) ''
    
    ! Заполняем числами от 0 до 10
    allocate(A(1 : n, 1 : n))
    call random_number(A)
    do i = 1, n
        do j = 1, n
            A(i, j) = anint(A(i, j) * 10)
        end do
    end do
    
    write(*, *) '>>>>>>C помощью собственных векторов'
    call system_clock(timeStart)
    call methodEV(A, n)
    call system_clock(timeStop)
    write(fmt, '(f8.3)') real(timeStop - timeStart) / 1000
    write(*, *) ''
    write(*, *) '     Время выполнения' // fmt // ' секунд' 
    write(*, *) ''
    write(*, *) ''
    
    write(*, *) '>>>>>>C помощью LUP-разложения'
    call system_clock(timeStart)
    call methodLUP(A, n)
    call system_clock(timeStop)
    write(fmt, '(f8.3)') real(timeStop - timeStart) / 1000
    write(*, *) ''
    write(*, *) '     Время выполнения' // fmt // ' секунд' 
    write(*, *) ''
    write(*, *) ''
    
    write(*, *) '>>>>>>С помощью SVD-разложения'
    call system_clock(timeStart)
    call methodSVD(A, n)
    call system_clock(timeStop)
    write(fmt, '(f8.3)') real(timeStop - timeStart) / 1000
    write(*, *) ''
    write(*, *) '     Время выполнения' // fmt // ' секунд' 
end program exercise2

! Нахождение детерминанта с помощью собственных векторов
subroutine methodEV(B, n)
    implicit none
    
    !---------------------------------------------------------------------------
    integer, intent(in) :: n
    real, dimension(1 : n, 1 : n) :: B, A, Vl, Vr
    real, dimension(1 : n) :: wr, wi
    real, dimension(1 : 5 * n) :: work
    character(len = 1) :: jobvl, jobvr
    integer :: i, j, lda, ldvl, ldvr, lwork, info
    complex :: det
    character(len = 20) :: fmt ! переменная для задания форматирования
    interface
        subroutine printMatrix(B)
            integer :: nRows, nColumns, k, l
            real, dimension(:, :), intent(in) :: B        
        end subroutine printMatrix
    end interface
    
    !---------------------------------------------------------------------------
    lda = n
    ldvl = n
    ldvr = n
    jobvl = 'V'
    jobvr = 'V'
    lwork = 5 * n
    
    ! Копируем матрицу B в A, чтобы предотвратить её изменение
    A = B
    
    !!
    write(*, *) "----- Матрица A -----"
    call printMatrix(A)
    write(*, *) ''
    
    !!
    call sgeev(jobvl, jobvr, n, A, lda, wr, wi, Vl, ldvl, Vr, ldvr, work, lwork, info)
    
    det = cmplx(wr(1), wi(1))
    do i = 2, n
        det = det * cmplx(wr(i), wi(i))
    end do
    
    !!
    write(fmt, '(f15.3)') real(det)
    write(*, *) "     Детерминант матрицы равен " // fmt
end subroutine methodEV
                              
! Нахождение детерминанта с помощью LUP-разложения
subroutine methodLUP(B, n)
    implicit none
    
    !---------------------------------------------------------------------------
    integer, intent(in) :: n
    real, dimension(1 : n, 1 : n) :: B, A
    integer, dimension(1 : n) :: pivot
    integer :: i, j, lda, info, s
    real :: det
    character(len = 20) :: fmt ! переменная для задания форматирования
    interface
        subroutine printMatrix(B)
            integer :: nRows, nColumns, k, l
            real, dimension(:, :), intent(in) :: B        
        end subroutine printMatrix
    end interface
    
    !---------------------------------------------------------------------------
    lda = n
    
    ! Копируем матрицу B в A, чтобы предотвратить её изменение
    A = B
    
    !!
    write(*, *) "----- Матрица A -----"
    call printMatrix(A)
    write(*, *) ''
    
    !!
    call sgetrf(n, n, A, n, pivot, info)
    
    det = 1
    do i = 1, n
        det = det * A(i, i)
    end do
    
    s = 0
    do i = 1, n
        if (i .ne. pivot(i)) then
            s = s + 1
        end if
    end do
    if (mod(s, 2) .ne. 0) then
        det = -det
    end if
    
    !!
    write(fmt, '(f15.3)') det 
    write(*, *) "     Детерминант матрицы равен " // fmt
end subroutine methodLUP

! Нахождение детерминанта с помощью SVD-разложения
subroutine methodSVD(B, n)
    implicit none
    
    !---------------------------------------------------------------------------
    integer, intent(in) :: n
    real, dimension(1 : n, 1 : n) :: B, A, U, Vt
    real, dimension(1 : n) :: s
    real, dimension(1 : 20 * n) :: work
    integer :: lda, ldu, ldvt, lwork, info, i, j
    character(len = 1) :: jobu, jobvt
    real :: det
    character(len = 20) :: fmt ! переменная для задания форматирования
    interface
        subroutine printMatrix(B)
            integer :: nRows, nColumns, k, l
            real, dimension(:, :), intent(in) :: B        
        end subroutine printMatrix
    end interface
    
    !---------------------------------------------------------------------------
    lda = n
    ldu = n
    ldvt = n
    jobu = 'A'
    jobvt = 'A'
    lwork = 20 * n ! Служебная переменная, чем больше, тем лучше
    
    ! Копируем матрицу B в A, чтобы предотвратить её изменение
    A = B
    
    !!
    write(*, *) "----- Матрица A -----"
    call printMatrix(A)
    write(*, *) ''
    
    !!
    call sgesvd(jobu, jobvt, n, n, A, lda, s, U, ldu, Vt, ldvt, work, lwork, info)
    
    det = 1
    do i = 1, n
        det = det * s(i)
    end do
    
    !!
    write(fmt, '(f15.3)') abs(det)
    write(*, *) "     Детерминант матрицы равен (только абсолютная величина) " // fmt
end subroutine methodSVD

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