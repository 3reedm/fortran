! ���᫥��� ��।��񭭮�� ��⥣ࠫ� ��⮤�� �࠯�権 (��᫥����⥫�� � ��ࠫ����� ��ਠ��). �ࠢ����� �� �ந�����⥫쭮��.
program exercise3
    implicit none
    
    ! ������ ���ᠭ�� ��६�����
    integer, parameter :: h = 1000000
    real, parameter :: x0 = 0, x1 = 3.14159
    real, external  :: f
    real :: res
    integer :: count, count1, time, time1, i
    interface 
        real function trap(f, x0, x1, h) result(res)
            real, external ::  f
            real :: x0, x1
            integer :: h
        end function trap
        real function trapP(f, x0, x1, h, th_num) result(res)
            real, external ::  f
            real :: x0, x1
            integer :: h, th_num
        end function trapP
    end interface
    
    
    ! ������ �����஢
    print *, "�஢�ઠ �ந�����⥫쭮�� �����⬮�"
    print *, ""

    write(*, *), '>>>>>>>>>���᫥��� ��।��񭭮�� ��⥣ࠫ� ��⮤�� �࠯�権'
    write(*, *), ''
    write(*, '(A,f5.2,A,f5.2)'), '>>>>>>>>>�㭪��: sin(x), x0 =', x0, ', x1 =', x1
    call system_clock(count)
    res = trap(f, x0, x1, h)
    call system_clock(count1)    
    write(*, *), ''
    write(*, '(A,f5.2)'), '>>>>>>>>>�������:', res
    print *, ""
    time = count1 - count
    print *, "time: ", real(time)/1000.
    print *, ""
    print *, ""
    
    write(*, *), '>>>>>>>>>���᫥��� ��।��񭭮�� ��⥣ࠫ� '
    write(*, *), '         ��⮤�� �࠯�権 (��ࠫ����� ��ਠ��)'
    do i=2,16 
        call system_clock(count)
        res = trapP(f, x0, x1, h, i)
        call system_clock(count1)
        write(*, *), ''
        write(*, '(A,f5.2)'), '>>>>>>>>>�������:', res
        time1 = count1 - count
        print *, "num_threads: ", i
        print *, "time: ", real(time1)/1000.
        print *, "acceleration: ", real(time)/real(time1)
        print *, "efficiency: ", real(time)/real(i*time1)
        print *, ""
    end do
end program exercise3

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

real function trapP(f, x0, x1, h, th_num) result(res)
    real, external ::  f
    real :: x0, x1
    integer :: h, th_num
    real, dimension(1 : h-1) :: summ

    real :: x, dx, summm, tmp
    integer :: j

    dx = (x1-x0) / h
    summ = f(x0) + f(x1)
    x = x0
    
    !$omp parallel num_threads(th_num)
        !$omp do
            do j=1,h-1 
                tmp = x + j*dx
                summ(j) = 2.0*f(tmp)
            end do
        !$omp end do
    !$omp end parallel
    summm = sum(summ)
    
    res = dx * summm / 2.0
end

real function f(x) 
    real :: x
    f = sin(x)
end