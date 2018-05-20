! ��宦����� ���ନ���� ������
program exercise2
    implicit none
    
    integer :: timeStart, timeStop, n, i, j
    real, allocatable :: A(:, :)
    character(len = 10) :: fmt
    
    write(*, *) '      ���᫥��� ��।���⥫� ������ (�� ���뢠� ����譮���)'
    write(*, *) ''
    
    write(*, *) "������ ࠧ��୮��� �����⭮� ������ [n]"
    read(*, *) n
    write(*, *) ''
    
    ! ������塞 �᫠�� �� 0 �� 10
    allocate(A(1 : n, 1 : n))
    call random_number(A)
    do i = 1, n
        do j = 1, n
            A(i, j) = anint(A(i, j) * 10)
        end do
    end do
    
    write(*, *) '>>>>>>C ������� ᮡ�⢥���� ����஢'
    call system_clock(timeStart)
    call methodEV(A, n)
    call system_clock(timeStop)
    write(fmt, '(f8.3)') real(timeStop - timeStart) / 1000
    write(*, *) ''
    write(*, *) '     �६� �믮������' // fmt // ' ᥪ㭤' 
    write(*, *) ''
    write(*, *) ''
    
    write(*, *) '>>>>>>C ������� LUP-ࠧ�������'
    call system_clock(timeStart)
    call methodLUP(A, n)
    call system_clock(timeStop)
    write(fmt, '(f8.3)') real(timeStop - timeStart) / 1000
    write(*, *) ''
    write(*, *) '     �६� �믮������' // fmt // ' ᥪ㭤' 
    write(*, *) ''
    write(*, *) ''
    
    write(*, *) '>>>>>>� ������� SVD-ࠧ�������'
    call system_clock(timeStart)
    call methodSVD(A, n)
    call system_clock(timeStop)
    write(fmt, '(f8.3)') real(timeStop - timeStart) / 1000
    write(*, *) ''
    write(*, *) '     �६� �믮������' // fmt // ' ᥪ㭤' 
end program exercise2

! ��宦����� ���ନ���� � ������� ᮡ�⢥���� ����஢
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
    character(len = 20) :: fmt ! ��६����� ��� ������� �ଠ�஢����
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
    
    ! �����㥬 ������ B � A, �⮡� �।������ �� ���������
    A = B
    
    !!
    write(*, *) "----- ����� A -----"
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
    write(*, *) "     ���ନ���� ������ ࠢ�� " // fmt
end subroutine methodEV
                              
! ��宦����� ���ନ���� � ������� LUP-ࠧ�������
subroutine methodLUP(B, n)
    implicit none
    
    !---------------------------------------------------------------------------
    integer, intent(in) :: n
    real, dimension(1 : n, 1 : n) :: B, A
    integer, dimension(1 : n) :: pivot
    integer :: i, j, lda, info, s
    real :: det
    character(len = 20) :: fmt ! ��६����� ��� ������� �ଠ�஢����
    interface
        subroutine printMatrix(B)
            integer :: nRows, nColumns, k, l
            real, dimension(:, :), intent(in) :: B        
        end subroutine printMatrix
    end interface
    
    !---------------------------------------------------------------------------
    lda = n
    
    ! �����㥬 ������ B � A, �⮡� �।������ �� ���������
    A = B
    
    !!
    write(*, *) "----- ����� A -----"
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
    write(*, *) "     ���ନ���� ������ ࠢ�� " // fmt
end subroutine methodLUP

! ��宦����� ���ନ���� � ������� SVD-ࠧ�������
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
    character(len = 20) :: fmt ! ��६����� ��� ������� �ଠ�஢����
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
    lwork = 20 * n ! ��㦥���� ��६�����, 祬 �����, ⥬ ����
    
    ! �����㥬 ������ B � A, �⮡� �।������ �� ���������
    A = B
    
    !!
    write(*, *) "----- ����� A -----"
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
    write(*, *) "     ���ନ���� ������ ࠢ�� (⮫쪮 ��᮫�⭠� ����稭�) " // fmt
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