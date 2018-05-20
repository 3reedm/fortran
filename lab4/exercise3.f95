! ��ࠫ���쭠� MPI-����� ���᫥��� ��।��񭭮�� ��⥣ࠫ�
program exercise3
    implicit none
    include 'mpif.h'
    
    real(kind = 8) :: resi, res, h, sum, x, f, a, b
    integer :: n, myid, numprocs, i, ierr
    
    ! ������� ����⥣ࠫ쭮� �㭪樨
    f(x) = sin(x)

    call mpi_init(ierr)
    
    ! �뢮� id ����� � ��饥 �᫮ ����ᮢ (��� �������)
    call mpi_comm_rank(mpi_comm_world, myid, ierr)
    call mpi_comm_size(mpi_comm_world, numprocs, ierr)
    print *, "����� ", myid, " �� ", numprocs
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    ! ����� � id = 0 ����訢��� ��室�� �����
    if (myid .eq. 0) then
        print *, ''
        print *, '������ �筮��� (�᫮ ࠧ�������) :'
        read *, n
        print *, ''
        
        print *, '������ �࠭��� ��⥣ࠫ�: [a, b]'
        read *, a, b
        print *, ''
        
        ! ���᫥��� 蠣� ���ࢠ��
        h = (b - a) / n
    endif
    
    ! ����뫪� ����ᮬ id = 0 ������, ����祭��� �� �।��饬 蠣�, �ᥬ ��⠫�� ����ᠬ ��㯯�
    call mpi_bcast(n, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(a, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
    call mpi_bcast(h, 1, mpi_double_precision, 0, mpi_comm_world, ierr)

    ! ���᫥��� ��।��񭭮�� ��⥣ࠫ� �ᥬ� ����ᠬ� (��� ������� ᢮� 㧫�) 
    sum  = 0.0
    do i = myid + 1, n, numprocs
        x = a + h * (dble(i) - 0.5)
        sum = sum + f(x)
    end do
    resi = h * sum 

    ! �㬬�஢���� १���⮢ ��� ����ᮢ
    call mpi_reduce(resi, res, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)

    ! ����� � id = 0 �뢮��� १����
    if (myid .eq. 0) then
        write(*, '(A, f12.5)'), ' ��⥣ࠫ = ', res
    endif

    call mpi_finalize(ierr)
end program exercise3