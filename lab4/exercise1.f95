! ��ࠫ���쭠� MPI-����� 㬭������ ����� � ������� �⠭���⭮�� ��⮤� � ��⮤� �����ࠤ�
program exercise1
    implicit none
    include "mpif.h" 
    
    ! �����
    real, allocatable :: A(:, :), B(:, :), C(:, :)
    integer :: m, n, k, myid, numprocs, ierr
    real(kind = 8) :: t1, t2, tl, tp
    interface
        subroutine printMatrix(B)
            real, dimension(:, :), intent(inout) :: B 
        end subroutine printMatrix
    end interface
    
    ! ����⢨�
    call mpi_init(ierr)
    
    call mpi_comm_rank(mpi_comm_world, myid, ierr)
    
    ! ���뢠��� ������
    if (myid .eq. 0) then
        call mpi_comm_size(mpi_comm_world, numprocs, ierr)
        
        write(*,  *) "                          ��������� �����"
        write(*,  *) ""
        write(*,  *) "        ������ ࠧ��୮�� ����� (A = m * n, B = n * k, C = m * k)"
        write(*,  *) ""
        
        write(*,  *) "������ m ="   
        read(*, *) m
    
        write(*,  *) "������ n ="
        read(*, *) n
    
        write(*,  *) "������ k ="
        read(*, *) k
    end if
    
    ! ����뫪� ������
    call mpi_bcast(m, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(n, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(k, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(numprocs, 1, mpi_integer, 0, mpi_comm_world, ierr)
    
    ! �뤥����� ����� ��� ���ᨢ A
    allocate(A(1 : m, 1 : n))  
    
    if (myid .eq. 0) then    
        ! �뤥����� ����� ��� ���ᨢ� B, C
        allocate(B(1 : n, 1 : k), C(1 : m, 1 : k))
    
        ! ���������� ���ᨢ�� A, B (�᫠�� �� 0 �� 1)
        call random_number(A)
        call random_number(B)
    end if

    ! ����뫪� ������ A
    call mpi_bcast(A, m * n, mpi_real, 0, mpi_comm_world, ierr) 
    
    if (myid .eq. 0) then   
        ! ����� ����� (A, B)
        !print *, ""
        !print *, "����� A"
        !call printMatrix(A)
        !print *, ""
        !print *, "����� B"
        !call printMatrix(B)
        print *, ""
        print *, ""
    
        print *, "                  �஢�ઠ �ந�����⥫쭮�� �����⬮� "
        print *, ""
    
        print *, "�⠭����� ������ 㬭������ �����"
        
        ! �ਬ������ �⠭���⭮� ��᫥����⥫쭮� ��楤��� 0 ����ᮬ
        t1 = mpi_wtime(ierr)
        call stdMM(A, B, C, m, k, n)
        t2 = mpi_wtime(ierr)
        tl = t2 - t1
        
        !call printMatrix(C)
        write(*, '(A, f8.3)') " ��᫥����⥫쭮� �६� (ᥪ㭤�): ", tl
    end if
    
    ! �ਬ������ �⠭���⭮� ��ࠫ���쭮� ��楤��� �ᥬ� ����ᠬ�
    t1 = mpi_wtime(ierr)
    call stdPMM(A, B, C, m, k, n, myid, numprocs)
    t2 = mpi_wtime(ierr)
    tp = t2 - t1
    
    if (myid .eq. 0) then
        !call printMatrix(C)
        write(*, '(A, f8.3)') " ��ࠫ���쭮� �६� (ᥪ㭤�): ", tp
        write(*, '(A, f8.3)') " �먣��� ��ࠫ���쭮�� ��ਠ�� (ᥪ㭤�): ", tl - tp
        print *, ""
    
        print *, "������ �����ࠤ� 㬭������ �����"
        
        ! �ਬ������ ��᫥����⥫쭮� ��楤��� �����ࠤ� 0 ����ᮬ
        t1 = mpi_wtime(ierr)
        call winogradMM(A, B, C, m, k, n)
        t2 = mpi_wtime(ierr)
        tl = t2 - t1
        
        !call printMatrix(C)
        write(*, '(A, f8.3)') " ��᫥����⥫쭮� �६� (ᥪ㭤�): ", tl
    end if
    
    ! �ਬ������ ��ࠫ���쭮� ��楤��� �����ࠤ� �ᥬ� ����ᠬ�
    t1 = mpi_wtime(ierr)
    call winogradPMM(A, B, C, m, k, n, myid, numprocs)
    t2 = mpi_wtime(ierr)
    tp = t2 - t1
    
    if (myid .eq. 0) then
        !call printMatrix(C)
        write(*, '(A, f8.3)') " ��ࠫ���쭮� �६� (ᥪ㭤�): ", tp
        write(*, '(A, f8.3)') " �먣��� ��ࠫ���쭮�� ��ਠ�� (ᥪ㭤�): ", tl - tp
        print *, ""
    end if
    
    call mpi_finalize(ierr)
end program exercise1

subroutine stdMM(A, B, C, nRows, nColumns, nAB)
    ! ���ᠭ�� ��६�����
    integer, intent(in) :: nRows, nColumns, nAB
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: A
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: B
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: C
    integer :: l, i, j
    
    ! ����⢨� 
    do i = 1, nRows
        do j = 1, nColumns
            C(i, j) = A(i, 1) * B(1, j)
            do l = 2, nAB
                C(i, j) = C(i, j) + A(i, l) * B(l, j)
            end do
        end do
    end do
end subroutine stdMM

subroutine winogradMM(A, B, C, nRows, nColumns, nAB)
    ! ���ᠭ�� ��६�����
    integer, intent(in) :: nRows, nColumns, nAB
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: A
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: B
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: C
    real, dimension(1 : nRows) :: rowFactor
    real, dimension(1 : nColumns) :: columnFactor
    integer :: l, i, j, d 
    
    ! ����⢨�   
    d = nAB / 2
    
    ! ���᫥��� rowFactor ��� ������ A
    do i = 1, nRows
        rowFactor(i) = A(i, 1) * A(i, 2)
        do j = 2, d
            rowFactor(i) = rowFactor(i) + A(i, 2 * j - 1) * A(i, 2 * j)
        end do
    end do
    
    ! ���᫥��� columnFactor ��� ������ B
    do i = 1, nColumns
        columnFactor(i) = B(1, i) * B(2, i)
        do j = 2, d
            columnFactor(i) = columnFactor(i) + B(2 * j - 1, i) * B(2 * j, i)
        end do
    end do
    
    ! ���᫥��� ������ C
    do i = 1, nRows
        do j = 1, nColumns
            C(i, j) = - rowFactor(i) - columnFactor(j)
            do l = 1, d
                C(i, j) = C(i, j) + (A(i, 2 * l - 1) + B(2 * l, j)) * (A(i, 2 * l) + B(2 * l - 1, j))
            end do
        end do
    end do
          
    ! �ਡ������� 童��� � ��砥 ����⭮� ࠧ��⪨
    if (mod(nAB, 2) .eq. 1) then  
        do i = 1, nRows
            do j = 1, nColumns
                C(i, j) = C(i, j) + A(i, nAB) * B(nAB, j)
            end do
        end do
    end if
end subroutine winogradMM

subroutine stdPMM(Ares, Bres, Cres, nRows, nColumns, nAB, myid, numprocs)
    include 'mpif.h'
    
    ! ���ᠭ�� ��६�����
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: Ares
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: Bres
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: Cres
    integer, intent(in) :: nRows, nColumns, nAB, myid, numprocs
    real, allocatable :: B(:, :), C(:, :)
    integer, dimension(0 : numprocs - 1) :: scounts, displs, offsets
    integer :: i, ierr, offset
   
    ! ����⢨�  
    ! ���������� �� ����ᠬ (� ��।�祩 ࠧ��� ���権 �ᥬ ����ᠬ)
    if (myid .eq. 0) call getStats(scounts, displs, offsets, nColumns, nAB, numprocs)
    call mpi_scatter(offsets, 1, mpi_integer, offset, 1, mpi_integer, 0, mpi_comm_world, ierr)
    
    ! �뤥����� ����� ��� ���ᨢ� B, C
    allocate(B(1 : nAB, 1 : offset), C(1 : nRows, 1 : offset))
    
    ! ����뫪� ���権 ����� ������ (Bres -> B) �ᥬ ����ᠬ
    call mpi_scatterv(Bres, scounts, displs, mpi_real, B, nAB * offset, mpi_real, 0, mpi_comm_world, ierr)
    
    ! �ਬ������ ��᫥����⥫쭮� ��楤��� �ᥬ� ����ᠬ�
    call stdMM(Ares, B, C, nRows, offset, nAB)
    
    ! ��������� ��ࠬ��஢ ��� ᡮન ���権 (� -< Cres)
    if (myid .eq. 0) call getStats(scounts, displs, offsets, nColumns, nRows, numprocs)
    call mpi_scatter(offsets, 1, mpi_integer, offset, 1, mpi_integer, 0, mpi_comm_world, ierr)
 
    ! ���ઠ ���権
    call mpi_gatherv(C, nRows * offset, mpi_real, Cres, scounts, displs, mpi_real, 0, mpi_comm_world, ierr)   
end subroutine stdPMM

subroutine winogradPMM(Ares, Bres, Cres, nRows, nColumns, nAB, myid, numprocs)
    include 'mpif.h'
    
    ! ���ᠭ�� ��६�����
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: Ares
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: Bres
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: Cres
    integer, intent(in) :: nRows, nColumns, nAB, myid, numprocs
    real, allocatable :: B(:, :), C(:, :)
    integer, dimension(0 : numprocs - 1) :: scounts, displs, offsets
    integer :: i, ierr, offset
    
    ! ����⢨� 
    ! ���������� �� ����ᠬ (� ��।�祩 ࠧ��� ���権 �ᥬ ����ᠬ)
    if (myid .eq. 0) call getStats(scounts, displs, offsets, nColumns, nAB, numprocs)
    call mpi_scatter(offsets, 1, mpi_integer, offset, 1, mpi_integer, 0, mpi_comm_world, ierr)
    
    ! �뤥����� ����� ��� ���ᨢ� B, C
    allocate(B(1 : nAB, 1 : offset), C(1 : nRows, 1 : offset))
    
    ! ����뫪� ���権 ����� ������ (Bres -> B) �ᥬ ����ᠬ
    call mpi_scatterv(Bres, scounts, displs, mpi_real, B, nAB * offset, mpi_real, 0, mpi_comm_world, ierr)
    
    ! �ਬ������ ��᫥����⥫쭮� ��楤��� �����ࠤ� �ᥬ� ����ᠬ�
    call winogradMM(Ares, B, C, nRows, offset, nAB)
    
    ! ��������� ��ࠬ��஢ ��� ᡮન ���権 (� -< Cres)
    if (myid .eq. 0) call getStats(scounts, displs, offsets, nColumns, nRows, numprocs)
    call mpi_scatter(offsets, 1, mpi_integer, offset, 1, mpi_integer, 0, mpi_comm_world, ierr)
 
    ! ���ઠ ���権
    call mpi_gatherv(C, nRows * offset, mpi_real, Cres, scounts, displs, mpi_real, 0, mpi_comm_world, ierr)
end subroutine winogradPMM

! ����ணࠬ��, ��� ����祭�� ������ ࠧ������� �� ����ᠬ
subroutine getStats(scounts, displs, offsets, nColumns, nRows, numprocs)  

    ! ���ᠭ�� ��६�����
    integer, intent(in) :: nColumns, nRows, numprocs
    integer, dimension(0 : numprocs - 1), intent(out) :: scounts, displs, offsets
    integer :: tmpOffset
    
    if (numprocs .ge. nColumns) then
        do i = 0, nColumns - 1
            scounts(i) = nRows
            displs(i) = i * nRows
            offsets(i) = 1
        end do
            
        do i = nColumns, numprocs - 1
            scounts(i) = 0
            displs(i) = 0
            offsets(i) = 0
        end do   
    else 
        tmpOffset = nint(real(nColumns) / real(numprocs))
            
        do i = 0, numprocs - 2
            scounts(i) = nRows * tmpOffset
            displs(i) = i * nRows * tmpOffset
            offsets(i) = tmpOffset
        end do
            
        scounts(numprocs - 1) = nRows * (nColumns - (numprocs - 1) * tmpOffset)
        displs(numprocs - 1) = (numprocs - 1) * nRows * tmpOffset
        offsets(numprocs - 1) = nColumns - (numprocs - 1) * tmpOffset
    end if
end subroutine getStats

! ����ணࠬ�� ��� ���� �� ��㬥୮� ������
subroutine printMatrix(A) 
    integer :: nRows, nColumns, i, j
    real, dimension(:, :), intent(inout) :: A 
    character(len = 20) :: fmt
       
    nRows = size(A, 1)
    nColumns = size(A, 2)
    
    do i = 1, nRows
        do j = 1, nColumns
            if (A(i, j) .eq. 0) then
                A(i, j) = abs(A(i, j))
            end if
        end do
    end do 
    
    do i = 1, nRows
        write(fmt, *) nColumns
        write(*, "("//adjustl(fmt)//"(f8.3), t1)") (A(i, j), j = 1, nColumns)
    end do  
end subroutine printMatrix