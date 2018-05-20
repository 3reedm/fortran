! ���᫥��� ���⭮� ������ � ������� SVD-ࠧ�������
program exercise3
    write(*, *) '>>>>>>���᫥��� ���⭮� ������ � ������� SVD-ࠧ�������'
    call methodSVD
end program exercise3

subroutine methodSVD
    implicit none
    
    !---------------------------------------------------------------------------
    character(len = 1) :: jobu, jobvt
    real, allocatable :: A(:, :), E(:, :), U(:, :), Vt(:, :), s(:), work(:)
    integer :: n, lda, ldu, ldvt, lwork, info, i, j
    logical :: flag
    real :: det
    character(len = 20) :: fmt ! ��६����� ��� ������� �ଠ�஢����
    interface
        subroutine printMatrix(B)
            integer :: nRows, nColumns, k, l
            real, dimension(:, :), intent(in) :: B        
        end subroutine printMatrix
    end interface
    
    !---------------------------------------------------------------------------
    write(*, *) "������ ࠧ��୮��� �����⭮� ������ n"
    read(*, *) n
    lda = n
    ldu = n
    ldvt = n
    jobu = 'A'
    jobvt = 'A'
    lwork = 20 * n ! ��㦥���� ��६�����, 祬 �����, ⥬ ����
    allocate(A(1 : n, 1 : n))
    allocate(E(1 : n, 1 : n))
    allocate(U(1 : n, 1 : n))
    allocate(Vt(1 : n, 1 : n))
    
    ! ���������� ������ ᨭ���୮� ������
    allocate(s(1 : n))
    
    ! ��㦥��� ���ᨢ
    allocate(work(1 : lwork))
    
    ! ������塞 �᫠�� �� 0 �� 10
    call random_number(A)
    do i = 1, n
        do j = 1, n
            A(i, j) = anint(A(i, j) * 10) 
        end do
    end do
    
    !!
    write(*, *) "----- ����� A -----"
    call printMatrix(A)
    write(*, *) ''
    
    !!
    call sgesvd(jobu, jobvt, n, n, A, lda, s, U, ldu, Vt, ldvt, work, lwork, info)
    
    do i = 1, n
        if (abs(s(i)) .lt. 1e-10) then
            flag = .true.
        else 
            flag = .false.
        end if
    end do
    
    if (flag .eqv. .false.) then
        do i = 1, n
            do j = 1, n
                if (i .eq. j) then
                    E(i, j) = 1 / s(i) 
                else
                    E(i, j) = 0
                end if
            end do
        end do
        
        A = matmul(transpose(Vt), matmul(E, transpose(U)))
        
        write(*, *) "----- �� ���⭠� ����� -----"
        call printMatrix(A)
    else
        write(*, *) "----- ���⭮� ������ �� ������� (����� ��஦������) -----"
    end if
end subroutine methodSVD

! ����ணࠬ��, ������� ���� ��㬥��� ������
subroutine printMatrix(a) 
    ! �室��� ��ࠬ���
    integer :: nRows, nColumns, i, j
    real, dimension(:, :), intent(inout) :: a 
    character(len = 20) :: fmt ! ��६����� ��� ������� �ଠ�஢����
       
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
