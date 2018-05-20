! ��ࠫ���쭠� ����� 㬭������ ����� � ������� �⠭���⭮�� ��⮤�, ��⮤� �����ࠤ�. �ࠢ����� �� �ந�����⥫쭮�� � �ந�����⥫쭮���� �㭪樨 matmul
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
    
    write(*,  *) "                          ��������� �����"
    write(*,  *) ""
    write(*,  *) "        ������ ࠧ��୮�� ����� (A = m * n, B = n * k, C = m * k)"
    write(*,  *) ""
    write(*,  *) "������ m ="   
    read(*, *) m
    write(*,  *) "������ n ="
    read(*, *) n
    write(*,  *) "������ k ="
    read(*, *) k
    write(*,  *) ""
    write(*, *) "������ th_num (������⢮ ��⥩) =" 
    read(*, *) th_num
    
    ! �뤥����� ����� ��� ���ᨢ�
    allocate(A(1 : m, 1 : n), B(1 : n, 1 : k), C(1 : m, 1 : k))
    
    ! ���������� ���ᨢ�� A � B (�᫠�� �� 0 �� 1)
    call random_number(A)
    call random_number(B)
    
    ! ����� ����� (A, B)
    !print *, ""
    !print *, "����� A"
    !call printMatrix(A)
    !print *, ""
    !print *, "����� B"
    !call printMatrix(B)
    print *, ""
    print *, ""
    
    print *, "                  �஢�ઠ �ந�����⥫쭮�� �����⬮� "
    print *, ""
    
    print *, "�⠭����� ������ 㬭������ �����"
    call system_clock(t1)
    call stdMM(A, B, C, m, k, n)
    call system_clock(t2)
    !call printMatrix(C)
    t3 = real(t2 - t1) / 1000
    write(*, '(A, f8.3)') " ��᫥����⥫쭮� �६� (ᥪ㭤�): ", t3
    call system_clock(t1)
    call stdPMM(A, B, C, m, k, n, th_num)
    call system_clock(t2)
    !call printMatrix(C)
    t3 = real(t2 - t1) / 1000 - t3
    write(*, '(A, f8.3)') " ��ࠫ���쭮� �६� (ᥪ㭤�): ", real(t2 - t1) / 1000
    write(*, '(A, f8.3)') " �먣��� ��ࠫ���쭮�� ��ਠ�� (ᥪ㭤�): ", -t3
    print *, ""
    
    print *, "������ �����ࠤ� 㬭������ �����"
    call system_clock(t1)
    call winogradMM(A, B, C, m, k, n)
    call system_clock(t2)
    t3 = real(t2 - t1) / 1000
    !call printMatrix(C)
    write(*, '(A, f8.3)') " ��᫥����⥫쭮� �६� (ᥪ㭤�): ", real(t2 - t1) / 1000
    call system_clock(t1)
    call winogradPMM(A, B, C, m, k, n, th_num)
    call system_clock(t2)
    t3 = real(t2 - t1) / 1000 - t3
    !call printMatrix(C)
    write(*, '(A, f8.3)') " ��ࠫ���쭮� �६� (ᥪ㭤�): ", real(t2 - t1) / 1000
    write(*, '(A, f8.3)') " �먣��� ��ࠫ���쭮�� ��ਠ�� (ᥪ㭤�): ", -t3
    print *, ""
    
    print *, "�ᯮ�짮����� �㭪樨 matmul ��� 㬭������ �����"
    call system_clock(t1)
    C = matmul(A, B)
    call system_clock(t2)
    !call printMatrix(C)
    write(*, '(A, f8.3)') " �६� (ᥪ㭤�): ", real(t2 - t1) / 1000
    print *, ""
end program exercise2

subroutine stdMM(A, B, C, nRows, nColumns, nAB)
    ! ���ᠭ�� ��६�����
    integer, intent(in) :: nRows, nColumns, nAB
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: A
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: B
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: C
    integer :: l, i, j
    
    ! ����⢨� 
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
    ! ���ᠭ�� ��६�����
    integer, intent(in) :: nRows, nColumns, nAB
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: A
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: B
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: C
    real, dimension(1 : nRows) :: rowFactor
    real, dimension(1 : nColumns) :: columnFactor
    integer :: l, i, j, d 
    
    ! ����⢨�   
    d = nAB / 2
    
    ! ���᫥��� rowFactor ��� ������ A
    do i = 1, nRows
        rowFactor(i) = A(i, 1) * A(i, 2)
        do j = 2, d
            rowFactor(i) = rowFactor(i) + A(i, 2 * j - 1) * A(i, 2 * j)
        end do
    end do
    
    ! ���᫥��� columnFactor ��� ������ B
    do i = 1, nColumns
        columnFactor(i) = B(1, i) * B(2, i)
        do j = 2, d
            columnFactor(i) = columnFactor(i) + B(2 * j - 1, i) * B(2 * j, i)
        end do
    end do
    
    ! ���᫥��� ������ C
    do i = 1, nRows
        do j = 1, nColumns
            C(i, j) = - rowFactor(i) - columnFactor(j)
            do l = 1, d
                C(i, j) = C(i, j) + (A(i, 2 * l - 1) + B(2 * l, j)) * (A(i, 2 * l) + B(2 * l - 1, j))
            end do
        end do
    end do
          
    ! �ਡ������� 童��� � ��砥 ����⭮� ࠧ��⪨
    if (mod(nAB, 2) .eq. 1) then  
        do i = 1, nRows
            do j = 1, nColumns
                C(i, j) = C(i, j) + A(i, nAB) * B(nAB, j)
            end do
        end do
    end if
end subroutine winogradMM

subroutine stdPMM(A, B, C, nRows, nColumns, nAB, th_num)
    ! ���ᠭ�� ��६�����
    integer, intent(in) :: nRows, nColumns, nAB, th_num
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: A
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: B
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: C
    integer :: i, j
    
    ! ����⢨� 
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
    ! ���ᠭ�� ��६�����
    integer, intent(in) :: nRows, nColumns, nAB, th_num
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: A
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: B
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: C
    real, dimension(1 : nRows) :: rowFactor
    real, dimension(1 : nColumns) :: columnFactor
    integer :: l, i, j, d 
    
    ! ����⢨�   
    d = nAB / 2
    
    !$omp parallel num_threads(th_num)
        ! ���᫥��� rowFactor ��� ������ A
        !$omp do
            do i = 1, nRows
                rowFactor(i) = A(i, 1) * A(i, 2)
                do j = 2, d
                    rowFactor(i) = rowFactor(i) + A(i, 2 * j - 1) * A(i, 2 * j)
                end do
            end do
        !$omp end do
    
        ! ���᫥��� columnFactor ��� ������ B
        !$omp do
            do i = 1, nColumns
                columnFactor(i) = B(1, i) * B(2, i)
                do j = 2, d
                    columnFactor(i) = columnFactor(i) + B(2 * j - 1, i) * B(2 * j, i)
                end do
            end do
        !$omp end do
        
        ! ���᫥��� ������ C
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
        
        ! �ਡ������� 童��� � ��砥 ����⭮� ࠧ��⪨
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

! ����ணࠬ�� ��� ���� �� ��㬥୮� ������
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