! �������� ����⮢ ���ᨢ� �����, ४��ᨢ��, �������樮��� ��⮤���.
program exercise1
    implicit none
    
    !----------------------------------------
    integer :: n,i
    real, allocatable :: b(:) 
    real :: summ
    interface
        real function sum1(b)
            real, dimension(:) :: b
        end function sum1
        real recursive function sum2(b)
            real, dimension(:) :: b
        end function sum2
        real function sum3(b)
            real, dimension(:) :: b
        end function sum3
    end interface
    
    n = 10
    allocate(b(1:n))
    do i=1,n
        b(i) = i
    end do
    write(*, "(A)") '>>>>>>���ᨢ: '
    write(*, "(f5.2)") b
    write(*, *) ''
    
    !----------------------------------------
    summ = sum1(b)
    write(*, "(A, f20.2, f5.2)") '>>>>>>���筮� �㬬�஢����: ', summ
    write(*, *) ''
    
    summ = sum2(b)
    write(*, "(A, f20.2)") '>>>>>>�����ᨢ��� �㬬�஢����: ', summ
    write(*, *) ''
    
    summ = sum3(b)
    write(*, "(A, f20.2)") '>>>>>>�������樮���� �㬬�஢����: ', summ
    write(*, *) ''
    
    deallocate(b)
end program exercise1

! ���筮� �㬬�஢����
real function sum1(b) result(res)
    implicit none
    
    !---------------------------------------------------------------------------
    real, dimension(:) :: b
    integer :: i,n
    
    n = size(b)
    do i=1,n
        res = res + b(i)
    end do
end function sum1

! �����ᨢ��� �㬬�஢����
real recursive function sum2(b) result(res)
    implicit none
    
    !---------------------------------------------------------------------------
    real, dimension(:) :: b
    real, allocatable :: c(:), d(:)
    integer :: i,n,k
    
    n = size(b)
    if (n .ge. 2) then
        k = n / 2
        
        allocate(c(1:k), d(1:n-k))
        
        do i=1,k
            c(i) = b(i)
        end do
        do i=1,n-k
            d(i) = b(i+k)
        end do
        
        res = sum2(c) + sum2(d)
        
        deallocate(c,d)
    else
        res = b(1) 
    end if
end function sum2

! �������樮���� �㬬�஢����
real function sum3(b) result(res)
    implicit none
    
    !---------------------------------------------------------------------------
    real, dimension(:) :: b
    integer :: i,n
    real :: c,y,t
    
    n = size(b)
    c = 0.0
    do i=1,n
        y = b(i) - c        !���� �� ���: c - ����.
        t = res + y         !���, res ������, y ����, ⠪ �� ����訥 ࠧ��� y ������.
        c = (t-res) - y     !(t - res) ����⠭�������� ���訥 ࠧ��� y; ���⠭�� y ����⠭�������� -(����訥 ࠧ��� y)
        res = t             !�����ࠨ�᪨, c �ᥣ�� ������ ࠢ������ ���. ��ॣ���� ᫨誮� ��⨬�������� ��������஢!
        !� ᫥���騩 ࠧ ����ﭭ� ����訥 ࠧ��� ���� ������ �ਡ������ � y.
    end do
end function sum3