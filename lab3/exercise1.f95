! ��ࠫ���쭠� ����� ��� ����䥭�
program exercise1
    implicit none 
    include "omp_lib.h" 
    integer :: i, j, n, p_num, th_num, t1, t2
    logical, allocatable :: A(:) 
    
    ! ���뫨 䠩�
    open(10, file = "primes.txt") 
    
    write(*,  *) "                  ���᫥��� ������ �ᥫ" 
    write(*,  *) "" 
    write(*,  *) "������ n =" 
    read(*, *) n
    write(*, *) "������ th_num =" 
    read(*, *) th_num
    write(*,  *) ""
    
    ! �뤥�塞 ������ ��� ���ᨢ
    allocate(A(2 : n)) 
    
    !-------------------------------------
    A = .true.
    
    call system_clock(t1) 
    !$omp parallel num_threads(th_num)
        !$omp do schedule(dynamic, 10)
            do i = 4, n, 2
                A(i) = .false.
            end do
        !$omp end do
        
        ! �����।�⢥��� ������
        !$omp do schedule(dynamic, 10)
            do i = 2, nint(sqrt(real(n)))
                if (A(i) .eqv. .true.) then
                    do j = i**2, n, i
                        A(j) = .false.
                    end do
                end if
            end do
        !$omp end do
    !$omp end parallel
    call system_clock(t2) 
    
    write(*, '(A, f8.3, A)') " �६� ࠡ��� �����⬠: ", real(t2 - t1) / 1000, " ᥪ㭤"
    !-------------------------------------

    ! �����뢠�� �������� ����� �᫠ � 䠩�
    p_num = 0 
    !$omp parallel do schedule(dynamic, 10) ordered num_threads(th_num) 
            do i = 2, n
                !$omp ordered
                if (A(i)) then
                    p_num = p_num + 1 
                    write(10, *) i
                end if
                !$omp end ordered
            end do
    !$omp end parallel do
    
    write(*, *), "������� ", p_num, " ������ �ᥫ"
    
    ! ����뢠�� 䠩� ��᫥ ⮣� ��� �����訫� � ��� ࠡ���
    close(10, status = "keep")
end program exercise1