! Вычисление определённого интеграла методом Монте-Карло (последовательный и параллельный вариант). Сравнение их производительности.
program exercise4
    implicit none
    
    ! Раздел описания переменных
    integer, parameter :: h = 1000000
    real, dimension(1:2) :: a, b
    real, external  :: f
    real :: res
    integer :: count, count1, time, time1, i
    interface 
        real function moka(f, a, b, m, h) result(res)
            real, external :: f
            integer :: m, h
            real, dimension(1 : m) :: a, b
        end function moka
        real function mokaP(f, a, b, m, h, th_num) result(res)
            real, external :: f
            integer :: m, h, th_num
            real, dimension(1 : m) :: a, b
        end function mokaP
    end interface
    
    
    ! Раздел операторов
    data a/1,-3/, b/3,5/
    
    print *, "Проверка производительности алгоритмов"
    print *, ""

    write(*, *), '>>>>>>>>>Вычисление кратного интеграла методом Монте-Карло'
    write(*, *), ''
    write(*, '(A)'), '>>>>>>>>>Функция: sin(x)'
    call system_clock(count)
    res = moka(f, a, b, 2, h)
    call system_clock(count1)    
    write(*, *), ''
    write(*, '(A,f5.2)'), '>>>>>>>>>Результат:', res
    print *, ""
    time = count1 - count
    print *, "time: ", real(time)/1000.
    print *, ""
    print *, ""
    
    write(*, *), '>>>>>>>>>Вычисление кратного интеграла '
    write(*, *), '         методом Монте-Карло (параллельный вариант)'
    do i=2,16 
        call system_clock(count)
        res = mokaP(f, a, b, 2, h, i)
        call system_clock(count1)
        write(*, *), ''
        write(*, '(A,f5.2)'), '>>>>>>>>>Результат:', res
        time1 = count1 - count
        print *, "num_threads: ", i
        print *, "time: ", real(time1)/1000.
        print *, "acceleration: ", real(time)/real(time1)
        print *, "efficiency: ", real(time)/real(i*time1)
        print *, ""
    end do
end program exercise4

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

real function mokaP(f, a, b, m, h, th_num) result(res)
    implicit none
    
    real, external :: f
    integer :: m, h, th_num
    real, dimension(1:m) :: a, b
    real, dimension(1 : h-1) :: summ

    real :: x, dx, summm, mu, tmp
    integer :: j,i

    mu = 1
    do i=1,m
        mu = mu*(b(i)-a(i))
    end do
    mu = mu/2
    
    dx = (b(1)-a(1)) / h
    
    !$omp parallel num_threads(th_num)
            do i=1,m
                !$omp single
                    summm = summm + f(a(1)) + f(b(1))
                !$omp end single
                !$omp do
                    do j=1,h-1 
                        tmp = a(1) + j*dx
                        summ(j) = f(tmp)
                    end do
                !$omp end do
                !$omp single
                    summm = summm + sum(summ)
                !$omp end single
            end do
    !$omp end parallel
    
    res = mu/h*summm
end function mokaP

real function f(x) 
    real :: x
    f = sin(x)
end