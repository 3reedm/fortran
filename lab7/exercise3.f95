! Умножение матриц с помощью стандартного метода, рекурсивного метода, метода Штрассена, функции matmul. Сравнение их производительности
program exercise3
    implicit none
    
    ! Раздел описания переменных
    real, dimension(1 : 2, 1 : 2) :: A 
    real, dimension(1 : 2, 1 : 2) :: B
    real, dimension(1 : 2, 1 : 2) :: C
    !integer :: count, count1
    interface
        recursive subroutine recurseMM(A, B, C, v)
            integer, intent(in) :: v
            real, dimension(1 : v, 1 : v), intent(in) :: A
            real, dimension(1 : v, 1 : v), intent(in) :: B
            real, dimension(1 : v, 1 : v), intent(out) :: C
        end subroutine recurseMM
        recursive subroutine shtrassenMM(A, B, C, v)
            integer, intent(in) :: v
            real, dimension(1 : v, 1 : v), intent(in) :: A
            real, dimension(1 : v, 1 : v), intent(in) :: B
            real, dimension(1 : v, 1 : v), intent(out) :: C
        end subroutine shtrassenMM
        subroutine printMatrix(a) 
            real, dimension(:, :), intent(in) :: a  
        end subroutine printMatrix
    end interface 
    
    ! Раздел операторов
    !call random_number(A)
    !call random_number(B)
    data A/1, 4,  6, 9/, B/2, 7,   3, 5/
    
    ! print *, "Проверка производительности алгоритмов для больших матриц (1000 x 1000)"
    ! print *, ""
    print *, "Стандартный алгоритм умножения матриц"
    print *, ""
    !call system_clock(count)
    call stdMM(A, B, C, 2, 2, 2)
    !call system_clock(count1)
    call printMatrix(C)
    print *, ""
    !print *, "time: ", count1 - count
    !print *, ""
    
    print *, "Рекурсивный алгоритм умножения матриц"
    print *, ""
    !call system_clock(count)
    call recurseMM(A, B, C, 2)
    !call system_clock(count1)
    call printMatrix(C)
    print *, ""
    !print *, "time: ", count1 - count
    !print *, ""
    
    print *, "Алгоритм Штрассена умножения матриц"
    print *, ""
    !call system_clock(count)
    call shtrassenMM(A, B, C, 2)
    !call system_clock(count1)
    call printMatrix(C)
    print *, ""
    !print *, "time: ", count1 - count
    !print *, ""
    
    print *, "Использование функции matmul для умножения матриц"
    print *, ""
    !call system_clock(count)
    C = matmul(A, B)
    !call system_clock(count1)
    call printMatrix(C)
    !print *, ""
    !print *, "time: ", count1 - count
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

recursive subroutine recurseMM(A, B, C, v)
    integer, intent(in) :: v
    real, dimension(1 : v, 1 : v), intent(in) :: A
    real, dimension(1 : v, 1 : v), intent(in) :: B
    real, dimension(1 : v, 1 : v), intent(out) :: C
    integer :: h
    real, dimension(:,:), allocatable :: A11, A12, A21, A22, B11, B12, B21, B22
    real, dimension(:,:), allocatable :: P1, P2, P3, P4, P5, P6, P7, P8
    
    ! Действия      
    h = v/2

    allocate(P1(h,h), P2(h,h), P3(h,h), P4(h,h), P5(h,h), P6(h,h), P7(h,h), P8(h,h))
    allocate(A11(h,h), A12(h,h), A21(h,h), A22(h,h), B11(h,h), B12(h,h), B21(h,h), B22(h,h))

    A11(1:h,1:h) = A(1:h,1:h)
    A12(1:h,1:h) = A(1:h,h+1:v)
    A21(1:h,1:h) = A(h+1:v,1:h)
    A22(1:h,1:h) = A(h+1:v,h+1:v)

    B11(1:h,1:h) = B(1:h,1:h)
    B12(1:h,1:h) = B(1:h,h+1:v)
    B21(1:h,1:h) = B(h+1:v,1:h)
    B22(1:h,1:h) = B(h+1:v,h+1:v)

    call shtrassenMM(A11,B11,P1,h)
    call shtrassenMM(A12,B21,P2,h) 
    call shtrassenMM(A11,B12,P3,h)
    call shtrassenMM(A12,B22,P4,h)
    call shtrassenMM(A21,B11,P5,h)
    call shtrassenMM(A22,B21,P6,h)
    call shtrassenMM(A21,B12,P7,h)
    call shtrassenMM(A22,B22,P8,h)

    deallocate(B11,B12,B21,B22)
    deallocate(A11,A12,A21,A22)

    C(1:h,1:h) = P1 + P2
    C(1:h,h+1:v) = P3 + P4
    C(h+1:v,1:h) = P5 + P6
    C(h+1:v,h+1:v) = P7 + P8

    deallocate(P1,P2,P3,P4,P5,P6,P7,P8)
end subroutine recurseMM

recursive subroutine shtrassenMM(A, B, C, v)
    ! Описание переменных
    integer, intent(in) :: v
    real, dimension(1 : v, 1 : v), intent(in) :: A
    real, dimension(1 : v, 1 : v), intent(in) :: B
    real, dimension(1 : v, 1 : v), intent(out) :: C
    integer :: h
    real, dimension(:,:), allocatable :: P1, P2, P3, P4, P5, P6, P7
    real, dimension(:,:), allocatable :: A11, A12, A21, A22, B11, B12, B21, B22
    
    ! Действия   
    if (v .le. 64) then 
        call stdMM(A,B,C,v,v,v)
        return
    end if
    
    h = v/2

    allocate(P1(h,h), P2(h,h), P3(h,h), P4(h,h), P5(h,h), P6(h,h), P7(h,h))
    allocate(A11(h,h), A12(h,h), A21(h,h), A22(h,h), B11(h,h), B12(h,h), B21(h,h), B22(h,h))

    A11(1:h,1:h) = A(1:h,1:h)
    A12(1:h,1:h) = A(1:h,h+1:v)
    A21(1:h,1:h) = A(h+1:v,1:h)
    A22(1:h,1:h) = A(h+1:v,h+1:v)

    B11(1:h,1:h) = B(1:h,1:h)
    B12(1:h,1:h) = B(1:h,h+1:v)
    B21(1:h,1:h) = B(h+1:v,1:h)
    B22(1:h,1:h) = B(h+1:v,h+1:v)

    call shtrassenMM(A11+A22,B11+B22,P1,h)
    call shtrassenMM(A21+A22,B11,P2,h) 
    call shtrassenMM(A11,B12-B22,P3,h)
    call shtrassenMM(A22,B21-B11,P4,h)
    call shtrassenMM(A11+A12,B22,P5,h)
    call shtrassenMM(A21-A11,B11+B12,P6,h)
    call shtrassenMM(A12-A22,B21+B22,P7,h)

    deallocate(B11,B12,B21,B22)
    deallocate(A11,A12,A21,A22)

    C(1:h,1:h) = P1 + P4 + P7 - P5
    C(1:h,h+1:v) = P3 + P5
    C(h+1:v,1:h) = P2 + P4
    C(h+1:v,h+1:v) = P1 - P2 + P3 + P6

    deallocate(P1,P2,P3,P4,P5,P6,P7)
end subroutine shtrassenMM

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