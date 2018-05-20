! Параллельная версия решета Эратосфена
program exercise1
    implicit none 
    include "omp_lib.h" 
    integer :: i, j, n, p_num, th_num, t1, t2
    logical, allocatable :: A(:) 
    
    ! Открыли файл
    open(10, file = "primes.txt") 
    
    write(*,  *) "                  Вычисление простых чисел" 
    write(*,  *) "" 
    write(*,  *) "Введите n =" 
    read(*, *) n
    write(*, *) "Введите th_num =" 
    read(*, *) th_num
    write(*,  *) ""
    
    ! Выделяем память под массив
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
        
        ! Непосредственно алгоритм
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
    
    write(*, '(A, f8.3, A)') " Время работы алгоритма: ", real(t2 - t1) / 1000, " секунд"
    !-------------------------------------

    ! Записываем найденные простые числа в файл
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
    
    write(*, *), "Найдено ", p_num, " простых чисел"
    
    ! Закрываем файл после того как завершили с ним работу
    close(10, status = "keep")
end program exercise1