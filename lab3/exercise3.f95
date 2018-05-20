! ��ᯠࠫ��������� ���஢�� ����쪮�
program exercise3
    implicit none
    include "omp_lib.h"
    integer :: i, j, n, th_num, first, t1, t2
    logical :: sorted
    real, allocatable :: a(:), b(:)
    
    ! ��砫�� �����
    write(*,  *) "                 ����஢�� ����쪮�"
    write(*, *) ""
    write(*,  *) "������ �᫮ ����⮢ ="
    read(*, *) n
    write(*, *) "������ th_num (������⢮ ��⥩) =" 
    read(*, *) th_num
    write(*, *) "" 
    
    allocate(a(1 : n))
    allocate(b(1 : n))
    call random_number(a)
    b = a
        
    ! ����� ��室���� ���ᨢ� 
    !write(*, *) "��室�� ���ᨢ"
    !write(*, '(f8.3)') a
    !write(*, *) ""
    
    ! �ਣ������ ������ ���஢�� ����쪮�
    call system_clock(t1) 
    do i = 1, n
        do j = 1, n - i
            if (a(j) .gt. a(j + 1)) then
                call swap(a(j), a(j + 1))
            end if
        end do
    end do
    call system_clock(t2)
    
    write(*, '(A, f8.3)') "����⨫� �६� � ᥪ㭤�� (��᫥����⥫쭮, 1-� ��ਠ��)", real(t2 - t1) / 1000 
    !write(*, '(f8.3)') a
    
    ! ���ࠡ�⠭�� ������ ����쪠 Odd-Even transposition sort
    a = b
    call system_clock(t1)
    sorted = .false.
    do while (.not. sorted) ! ������ ���� sorted ����
        sorted = .true. 
        
        do i = 2, n - 1, 2
            if (a(i) .gt. a(i + 1)) then
                call swap(a(i), a(i + 1)) 
                sorted = .false.
            end if
        end do
        
        do i = 1, n - 1, 2
            if (a(i) .gt. a(i + 1)) then
                call swap(a(i), a(i + 1)) 
                sorted = .false.
            end if
        end do
    end do
    call system_clock(t2)
    
    write(*, '(A, f8.3)') "����⨫� �६� � ᥪ㭤�� (��᫥����⥫쭮, 2-�� ��ਠ��)", real(t2 - t1) / 1000 
    !write(*, '(f8.3)') a
    
    ! ��ࠫ����� ���ࠡ�⠭�� ������ ����쪠
    a = b
    call system_clock(t1)
    !$omp parallel num_threads(th_num)
        do i = 1, n
            first = mod(i, 2)
            !$omp do
                do j = first, n - 1, 2
                    if (a(j) .gt. a(j + 1)) then
                        call swap(a(j), a(j + 1))
                    end if
                end do
            !$omp end do
        end do
    !$omp end parallel    
    call system_clock(t2)
    
    write(*, '(A, f8.3)') "����⨫� �६� � ᥪ㭤�� (��ࠫ���쭮)", real(t2 - t1) / 1000 
    !write(*, '(f8.3)') a
end program exercise3

! ����ணࠬ�� ��� ����⠭���� ����⮢ ���ᨢ� 
subroutine swap(a, b)
    real :: a, b, tmp
   
    tmp = a
    a = b
    b = tmp
end subroutine swap