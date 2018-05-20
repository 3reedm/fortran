! Вычисление чисел Фибоначчи с помощью трёх различных алгоритмов: рекурсивный, полиномиальный, матричный
program exercise1
    implicit none
    
    ! Раздел описания переменных
    integer :: n
    real :: m
    integer, external :: pfibonachchi, mfibonachchi 
    interface
        recursive function rfibonachchi(n) result(res)
            integer :: n, res
        end function rfibonachchi
    end interface
    
    ! Раздел операторов
    print *, 'Введите номер числа, которое хотите получить, в последовательности Фибоначчи'
    read *, m 
    
    if (m .ge. 0) then
        n = rfibonachchi(nint(m))   
        print *,''
        print *, '<--Рекурсивный алгоритм-->'   
        print '(A44, I6)', 'Этому номеру соответствует следующее число:', n
        print *,''
        
        n = pfibonachchi(nint(m))   
        print *, '<--Полиномиальный алгоритм-->' 
        print '(A44, I6)', 'Этому номеру соответствует следующее число:', n
        print *,''
        
        n = mfibonachchi(nint(m))   
        print *, '<--Матричный алгоритм-->'   
        print '(A44, I6)', 'Этому номеру соответствует следующее число:', n
    else 
        print *, 'Такого числа в последовательности Фибоначчи быть не может'
    end if
end program exercise1

recursive function rfibonachchi(n) result(res)
    integer :: n, res
        
    if (n .eq. 0) then
        res = 0
    end if
    
    if (n .ge. 1) then
        if (n .gt. 1) then
            res = rfibonachchi(n - 1) + rfibonachchi(n - 2)
        else
            res = 1 
        end if
    end if
end function rfibonachchi

integer function pfibonachchi(n) result(res)
    integer :: n, i
    integer, dimension(0:n) :: f
    
    if (n .ge. 1) then
        f(0) = 0
        f(1) = 1
        do i=2,n 
            f(i) = f(i-1) + f(i-2)
        end do
        res = f(n)
    end if
    
    if (n .eq. 0) then
        res = 0
    end if
end function pfibonachchi

integer function mfibonachchi(n) result(res)
    integer :: n, i
    integer :: f1(0:1), f2(n:n+1), P(0:1,0:1), lp(0:1,0:1)
    
    data f1/0,1/, P/0,1,1,1/
    lp = P
    
    do i=2,n
        lp = matmul(lp,P)
    end do
    
    f2 = matmul(f1,lp) 
    
    if (n .eq. 0) then
        res = 0
    else
        res = f2(n)
    end if
end function mfibonachchi