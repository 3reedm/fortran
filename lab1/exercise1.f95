! ���������� ����� ��������� � ������� ����������� �������
program exercise1
    implicit none
    
    ! ������ �������� ����������
    integer :: n
    real :: m
    interface
        recursive function fibonachchi(n) result(res)
            integer :: n, res
        end function
    end interface
    
    ! ������ ����������
    print *, '������� ����� �����, ������� ������ ��������, � ������������������ ���������'
    read *, m 
    
    if ((m - 1) .ge. 1e-10) then
        n = fibonachchi(nint(m))    
        print '(A44, I6)', '����� ������ ������������� ��������� �����:', n
    else 
        print *, '������ ����� � ������������������ ��������� ���� �� �����'
    end if
    pause
end program exercise1

recursive function fibonachchi(n) result(res)
    integer :: n, res
    
    if (n .ge. 1) then
        if (n .gt. 1) then
            res = fibonachchi(n - 1) + fibonachchi(n - 2)
        else
            res = 1 
        end if
    end if
end function fibonachchi