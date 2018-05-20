! Вычисление определённого интеграла методом трапеций
program exercise5
    integer, parameter :: h = 100
    real, parameter :: x0 = 0, x1 = 3.14159
    real, external  :: f
    real :: res

    write(*, *), '>>>>>>>>>Вычисление определённого интеграла методом трапеций'
    write(*, *), ''
    write(*, '(A,f5.2,A,f5.2)'), '>>>>>>>>>Функция: sin(x), x0 =', x0, ', x1 =', x1
    res = trap(f, x0, x1, h) 
    write(*, *), ''
    write(*, '(A,f5.2)'), '>>>>>>>>>Результат:', res
end program exercise5

real function trap(f, x0, x1, h) result(res)
    real, external ::  f
    real :: x0, x1
    integer :: h

    real :: x, dx, summ
    integer :: j

    dx = (x1-x0) / h
    summ = f(x0) + f(x1)
    x = x0
    
    do j=1,h-1 
        x = x + dx
        summ = summ + 2.0*f(x)
    end do
    
    res = dx * summ / 2.0
end

real function f(x) 
    real :: x
    f = sin(x)
end