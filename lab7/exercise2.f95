! ���᫥��� ����⨢��� �㬬� � �ந�������� ����⮢ ���ᨢ� (ࠧ���� � ������)
program exercise1
    implicit none
    
    !----------------------------------------
    integer :: n,i
    real, allocatable :: b(:) 
    real :: res
    interface
        real recursive function csum(b)
            real, dimension(:) :: b
        end function csum
        real recursive function prod(b)
            real, dimension(:) :: b
        end function prod
    end interface
    
    n = 5
    allocate(b(1:n))
    do i=1,n
        b(i) = i
    end do
    write(*, "(A)") '>>>>>>���ᨢ: '
    write(*, "(f5.2)") b
    write(*, *) ''
    
    !----------------------------------------
    res = csum(b)
    write(*, "(A, f20.2, f5.2)") '>>>>>>�㬬��⨢��� �㬬�: ', res
    write(*, *) ''
    
    res = prod(b)
    write(*, "(A, f20.2)") '>>>>>>�ந��������: ', res
    write(*, *) ''
    
    deallocate(b)
end program exercise1

! �����ᨢ��� ����⨢��� �㬬�
real recursive function csum(b) result(res)
    implicit none
    
    !---------------------------------------------------------------------------
    real, dimension(:) :: b
    real, allocatable :: c(:), d(:)
    integer :: i,n,k
    
    n = size(b)
    if (n .ge. 2) then
        k = n - 1
        
        allocate(c(1:k))
        
        do i=1,k
            c(i) = b(i)
        end do
        
        res = csum(c) + sum(b)
        
        deallocate(c)
    else
        res = b(1) 
    end if
end function csum

! �����ᨢ��� �ந��������
real recursive function prod(b) result(res)
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
        
        res = prod(c) * prod(d)
        
        deallocate(c,d)
    else
        res = b(1) 
    end if
end function prod