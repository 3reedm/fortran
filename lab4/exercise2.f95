! ����� ᮮ�饭�ﬨ ����� ����ᠬ� �� ����
program exercise2
    implicit none
    include 'mpif.h'
 
    ! �����
    integer :: rounds, myid, numprocs, msgTag, i, j, ierr
    integer, parameter :: n = 4096
    real(kind = 8) :: t1, t2, ts
    character(len = 10) :: tmp
    character(len = n) :: msg
    
    ! ����⢨�
    call mpi_init(ierr)
    
    ! �뢮� id ����� � ��饥 �᫮ ����ᮢ (��� �������)
    call mpi_comm_rank(mpi_comm_world, myid, ierr)
    call mpi_comm_size(mpi_comm_world, numprocs, ierr)
    
    !print *, "����� ", myid, " �� ", numprocs
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    ! ���� ��砫��� ������
    if (myid .eq. 0) then
        msg = 'Good Times'
        msgTag = 10
        call get_environment_variable("ROUNDS", tmp)
        read(tmp, *) rounds
        print *, ''
        write(*, '(A, I5)') ' ����� ᮮ�饭�� (����): ', n
        print *, ''
        print *, '<<<<����� ����ணࠬ�� � �����஢���<<<<'
    end if
    
    ! ����뫪� 0-� ����ᮬ ���⮢�� ������ �ᥬ ��⠫�� ����ᠬ
    call mpi_bcast(msgTag, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(rounds, 1, mpi_integer, 0, mpi_comm_world, ierr)
    
    ! �믮������ ����ணࠬ�� � �ᯮ�짮������ ���������� ��।��
    t1 = mpi_wtime(ierr)
    call subSendBlock(msg, n, msgTag, myid, numprocs, rounds) 
    t2 = mpi_wtime(ierr)
    ts = t2 - t1
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    if (myid .eq. 0) then
        write(*, '(A, f8.3)') ' >>>>�६� ࠡ��� (ᥪ)>>>>: ', ts
        print *, ''
        !print *, '>>>>����饭��, ����祭��� ����� ����ᮬ (��� �஢�ન)>>>>'
    end if
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    !print *, msg
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    ! �믮������ ����ணࠬ�� � �ᯮ�짮������ ������������ ��।��
    if (myid .eq. 0) then
        msg = 'Bad Times'
        msgTag = 20
        print *, ''
        print *, '<<<<����� ����ணࠬ�� ��� �����஢��<<<<'
    end if
    
    ! ����뫪� 0-� ����ᮬ ���⮢�� ������ �ᥬ ��⠫�� ����ᠬ
    call mpi_bcast(msgTag, 1, mpi_integer, 0, mpi_comm_world, ierr)
    
    t1 = mpi_wtime(ierr)
    call subSendUnblock(msg, n, msgTag, myid, numprocs, rounds)
    t2 = mpi_wtime(ierr)
    ts = t2 - t1
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    if (myid .eq. 0) then
        write(*, '(A, f8.3)') ' >>>>�६� ࠡ��� (ᥪ)>>>>: ', ts
        !print *, ''
        !print *, '>>>>����饭��, ����祭��� ����� ����ᮬ (��� �஢�ન)>>>>'
    end if
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    !print *, msg
    
    call mpi_finalize(ierr)
end program exercise2

subroutine subSendBlock(msg, n, msgTag, myid, numprocs, rounds)  
    include 'mpif.h'
    
    ! �����
    integer :: i, j, ierr
    integer, intent(in) :: rounds, myid, numprocs, n, msgTag
    integer, dimension(mpi_status_size) :: status
    character(len = n), intent(inout) :: msg

    ! �믮������
    do i = 0, rounds - 1
        do j = 0, numprocs - 1
            if ((myid .eq. j) .and. (myid .ne. (numprocs - 1))) then
                call mpi_send(msg, n, mpi_character, j + 1, msgTag, mpi_comm_world, ierr)
                !write(*, *), "����� �", j, "��᫠� ᮮ�饭��: ", msg
            else if ((myid .eq. j) .and. (myid .eq. (numprocs - 1))) then
                call mpi_send(msg, n, mpi_character, 0, msgTag, mpi_comm_world, ierr)
                !write(*, *), "����� �", j, "��᫠� ᮮ�饭��: ", msg
            end if
            call mpi_barrier(mpi_comm_world, ierr)
            
            if ((myid .eq. (j + 1)) .and. ((j + 1) .ne. numprocs)) then
                call mpi_recv(msg, n, mpi_character, j, msgTag, mpi_comm_world, status, ierr)
                !write(*, *), "����� �", j + 1, "�ਭ� ᮮ�饭��: ", msg
            else if ((myid .eq. 0) .and. ((j + 1) .eq. numprocs)) then
                call mpi_recv(msg, n, mpi_character, j, msgTag, mpi_comm_world, status, ierr)
                !write(*, *), "����� �", 0, "�ਭ� ᮮ�饭��: ", msg
            end if
        end do
    end do
end subroutine subSendBlock

subroutine subSendUnblock(msg, n, msgTag, myid, numprocs, rounds)
    include 'mpif.h'
    
    ! �����
    integer :: i, j, ierr, request
    integer, dimension(mpi_status_size) :: status
    integer, intent(in) :: rounds, myid, numprocs, n, msgTag
    character(len = n), intent(inout) :: msg
    
    ! �믮������
    request = myid
    do i = 0, rounds - 1
        do j = 0, numprocs - 1
            if ((myid .eq. j) .and. (myid .ne. (numprocs - 1))) then
                call mpi_isend(msg, n, mpi_character, j + 1, msgTag, mpi_comm_world, request, ierr)
                call mpi_wait(request, status, ierr)
                !write(*, *), "����� �", j, "��᫠� ᮮ�饭��: ", msg
            else if ((myid .eq. j) .and. (myid .eq. (numprocs - 1))) then
                call mpi_isend(msg, n, mpi_character, 0, msgTag, mpi_comm_world, request, ierr)
                call mpi_wait(request, status, ierr)
                !write(*, *), "����� �", j, "��᫠� ᮮ�饭��: ", msg
            end if
            
            if ((myid .eq. (j + 1)) .and. ((j + 1) .ne. numprocs)) then
                call mpi_irecv(msg, n, mpi_character, j, msgTag, mpi_comm_world, request, ierr)
                call mpi_wait(request, status, ierr)
                !write(*, *), "����� �", j + 1, "�ਭ� ᮮ�饭��: ", msg
            else if ((myid .eq. 0) .and. ((j + 1) .eq. numprocs)) then
                call mpi_irecv(msg, n, mpi_character, j, msgTag, mpi_comm_world, request, ierr)
                call mpi_wait(request, status, ierr)
                !write(*, *), "����� �", 0, "�ਭ� ᮮ�饭��: ", msg
            end if
        end do
    end do
end subroutine subSendUnblock
