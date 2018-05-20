! Сложение элементов массива обычным, рекурсивным, компенсационным методами.
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
    write(*, "(A)") '>>>>>>Массив: '
    write(*, "(f5.2)") b
    write(*, *) ''
    
    !----------------------------------------
    summ = sum1(b)
    write(*, "(A, f20.2, f5.2)") '>>>>>>Обычное суммирование: ', summ
    write(*, *) ''
    
    summ = sum2(b)
    write(*, "(A, f20.2)") '>>>>>>Рекурсивное суммирование: ', summ
    write(*, *) ''
    
    summ = sum3(b)
    write(*, "(A, f20.2)") '>>>>>>Компенсационное суммирование: ', summ
    write(*, *) ''
    
    deallocate(b)
end program exercise1

! Обычное суммирование
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

! Рекурсивное суммирование
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

! Компенсационное суммирование
real function sum3(b) result(res)
    implicit none
    
    !---------------------------------------------------------------------------
    real, dimension(:) :: b
    integer :: i,n
    real :: c,y,t
    
    n = size(b)
    c = 0.0
    do i=1,n
        y = b(i) - c        !Пока все хорошо: c - ноль.
        t = res + y         !Увы, res велика, y мало, так что младшие разряды y потеряны.
        c = (t-res) - y     !(t - res) восстанавливает старшие разряды y; вычитание y восстанавливает -(младшие разряды y)
        res = t             !Алгебраически, c всегда должна равняться нулю. Берегитесь слишком оптимизирующих компиляторов!
        !В следующий раз потерянные младшие разряды будут заново прибавлены к y.
    end do
end function sum3