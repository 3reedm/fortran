! Параллельная MPI-версия вычисления определённого интеграла
program exercise3
    implicit none
    include 'mpif.h'
    
    real(kind = 8) :: resi, res, h, sum, x, f, a, b
    integer :: n, myid, numprocs, i, ierr
    
    ! Задание подынтегральной функции
    f(x) = sin(x)

    call mpi_init(ierr)
    
    ! Вывод id процесса и общее число процессов (для каждого)
    call mpi_comm_rank(mpi_comm_world, myid, ierr)
    call mpi_comm_size(mpi_comm_world, numprocs, ierr)
    print *, "Процесс ", myid, " из ", numprocs
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    ! Процесс с id = 0 запрашивает исходные данные
    if (myid .eq. 0) then
        print *, ''
        print *, 'Введите точность (число разбиейний) :'
        read *, n
        print *, ''
        
        print *, 'Введите границы интеграла: [a, b]'
        read *, a, b
        print *, ''
        
        ! Вычисление шага интервала
        h = (b - a) / n
    endif
    
    ! Рассылка процессом id = 0 данных, полученных на предыдущем шаге, всем остальным процессам группы
    call mpi_bcast(n, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(a, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
    call mpi_bcast(h, 1, mpi_double_precision, 0, mpi_comm_world, ierr)

    ! Вычисление определённого интеграла всеми процессами (для каждого свои узлы) 
    sum  = 0.0
    do i = myid + 1, n, numprocs
        x = a + h * (dble(i) - 0.5)
        sum = sum + f(x)
    end do
    resi = h * sum 

    ! Суммирование результатов всех процессов
    call mpi_reduce(resi, res, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)

    ! Процесс с id = 0 выводит результат
    if (myid .eq. 0) then
        write(*, '(A, f12.5)'), ' Интеграл = ', res
    endif

    call mpi_finalize(ierr)
end program exercise3