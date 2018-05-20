! Обмен сообщениями между процессами по кругу
program exercise2
    implicit none
    include 'mpif.h'
 
    ! Данные
    integer :: rounds, myid, numprocs, msgTag, i, j, ierr
    integer, parameter :: n = 4096
    real(kind = 8) :: t1, t2, ts
    character(len = 10) :: tmp
    character(len = n) :: msg
    
    ! Действия
    call mpi_init(ierr)
    
    ! Вывод id процесса и общее число процессов (для каждого)
    call mpi_comm_rank(mpi_comm_world, myid, ierr)
    call mpi_comm_size(mpi_comm_world, numprocs, ierr)
    
    !print *, "Процесс ", myid, " из ", numprocs
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    ! Ввод начальных данных
    if (myid .eq. 0) then
        msg = 'Good Times'
        msgTag = 10
        call get_environment_variable("ROUNDS", tmp)
        read(tmp, *) rounds
        print *, ''
        write(*, '(A, I5)') ' Длина сообщения (байт): ', n
        print *, ''
        print *, '<<<<Работа подпрограммы с блокировкой<<<<'
    end if
    
    ! Рассылка 0-ым процессом стартовых данных всем остальным процессам
    call mpi_bcast(msgTag, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(rounds, 1, mpi_integer, 0, mpi_comm_world, ierr)
    
    ! Выполнение подпрограммы с использованием блокирующих передач
    t1 = mpi_wtime(ierr)
    call subSendBlock(msg, n, msgTag, myid, numprocs, rounds) 
    t2 = mpi_wtime(ierr)
    ts = t2 - t1
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    if (myid .eq. 0) then
        write(*, '(A, f8.3)') ' >>>>Время работы (сек)>>>>: ', ts
        print *, ''
        !print *, '>>>>Сообщение, полученное каждым процессом (для проверки)>>>>'
    end if
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    !print *, msg
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    ! Выполнение подпрограммы с использованием неблокирующих передач
    if (myid .eq. 0) then
        msg = 'Bad Times'
        msgTag = 20
        print *, ''
        print *, '<<<<Работа подпрограммы без блокировки<<<<'
    end if
    
    ! Рассылка 0-ым процессом стартовых данных всем остальным процессам
    call mpi_bcast(msgTag, 1, mpi_integer, 0, mpi_comm_world, ierr)
    
    t1 = mpi_wtime(ierr)
    call subSendUnblock(msg, n, msgTag, myid, numprocs, rounds)
    t2 = mpi_wtime(ierr)
    ts = t2 - t1
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    if (myid .eq. 0) then
        write(*, '(A, f8.3)') ' >>>>Время работы (сек)>>>>: ', ts
        !print *, ''
        !print *, '>>>>Сообщение, полученное каждым процессом (для проверки)>>>>'
    end if
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    !print *, msg
    
    call mpi_finalize(ierr)
end program exercise2

subroutine subSendBlock(msg, n, msgTag, myid, numprocs, rounds)  
    include 'mpif.h'
    
    ! Данные
    integer :: i, j, ierr
    integer, intent(in) :: rounds, myid, numprocs, n, msgTag
    integer, dimension(mpi_status_size) :: status
    character(len = n), intent(inout) :: msg

    ! Выполнение
    do i = 0, rounds - 1
        do j = 0, numprocs - 1
            if ((myid .eq. j) .and. (myid .ne. (numprocs - 1))) then
                call mpi_send(msg, n, mpi_character, j + 1, msgTag, mpi_comm_world, ierr)
                !write(*, *), "Процесс №", j, "послал сообщение: ", msg
            else if ((myid .eq. j) .and. (myid .eq. (numprocs - 1))) then
                call mpi_send(msg, n, mpi_character, 0, msgTag, mpi_comm_world, ierr)
                !write(*, *), "Процесс №", j, "послал сообщение: ", msg
            end if
            call mpi_barrier(mpi_comm_world, ierr)
            
            if ((myid .eq. (j + 1)) .and. ((j + 1) .ne. numprocs)) then
                call mpi_recv(msg, n, mpi_character, j, msgTag, mpi_comm_world, status, ierr)
                !write(*, *), "Процесс №", j + 1, "принял сообщение: ", msg
            else if ((myid .eq. 0) .and. ((j + 1) .eq. numprocs)) then
                call mpi_recv(msg, n, mpi_character, j, msgTag, mpi_comm_world, status, ierr)
                !write(*, *), "Процесс №", 0, "принял сообщение: ", msg
            end if
        end do
    end do
end subroutine subSendBlock

subroutine subSendUnblock(msg, n, msgTag, myid, numprocs, rounds)
    include 'mpif.h'
    
    ! Данные
    integer :: i, j, ierr, request
    integer, dimension(mpi_status_size) :: status
    integer, intent(in) :: rounds, myid, numprocs, n, msgTag
    character(len = n), intent(inout) :: msg
    
    ! Выполнение
    request = myid
    do i = 0, rounds - 1
        do j = 0, numprocs - 1
            if ((myid .eq. j) .and. (myid .ne. (numprocs - 1))) then
                call mpi_isend(msg, n, mpi_character, j + 1, msgTag, mpi_comm_world, request, ierr)
                call mpi_wait(request, status, ierr)
                !write(*, *), "Процесс №", j, "послал сообщение: ", msg
            else if ((myid .eq. j) .and. (myid .eq. (numprocs - 1))) then
                call mpi_isend(msg, n, mpi_character, 0, msgTag, mpi_comm_world, request, ierr)
                call mpi_wait(request, status, ierr)
                !write(*, *), "Процесс №", j, "послал сообщение: ", msg
            end if
            
            if ((myid .eq. (j + 1)) .and. ((j + 1) .ne. numprocs)) then
                call mpi_irecv(msg, n, mpi_character, j, msgTag, mpi_comm_world, request, ierr)
                call mpi_wait(request, status, ierr)
                !write(*, *), "Процесс №", j + 1, "принял сообщение: ", msg
            else if ((myid .eq. 0) .and. ((j + 1) .eq. numprocs)) then
                call mpi_irecv(msg, n, mpi_character, j, msgTag, mpi_comm_world, request, ierr)
                call mpi_wait(request, status, ierr)
                !write(*, *), "Процесс №", 0, "принял сообщение: ", msg
            end if
        end do
    end do
end subroutine subSendUnblock
