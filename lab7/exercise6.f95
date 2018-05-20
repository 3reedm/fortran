! Вычисление кратного интеграла методом Монте-Карло
program exercise6
    implicit none

    integer, parameter :: h = 100
    real, dimension(1:2) :: a, b
    real, external :: f, moka
    real :: res

    !----------------
    data a/1,-3/, b/3,5/
    
    write(*, *), '>>>>>>>>>Вычисление кратного интеграла методом Монте-Карло'
    write(*, *), ''
    write(*, '(A)'), '>>>>>>>>>Функция: sin(x)'
    res = moka(f, a, b, 2, h) 
    write(*, *), ''
    write(*, '(A,f5.2)'), '>>>>>>>>>Результат: ', res
end program exercise6

real function moka(f, a, b, m, h) result(res)
    implicit none
    
    real, external :: f
    integer :: m, h
    real, dimension(1:m) :: a, b

    real :: x, dx, summ, mu
    integer :: j,i

    mu = 1
    do i=1,m
        mu = mu*(b(i)-a(i))
    end do
    mu = mu/2
    
    dx = (b(1)-a(1)) / h
    
    summ = 0
    do i=1,m
        x = a(1)
        summ = summ + f(a(1)) + f(b(1))
        do j=1,h-1
            x = x + dx
            summ = summ + f(x)
        end do
    end do
    
    res = mu/h*summ
end function moka

real function f(x) 
    real :: x
    f = sin(x)
end function f