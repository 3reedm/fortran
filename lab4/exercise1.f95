! Параллельная MPI-версия умножения матриц с помощью стандартного метода и метода Винограда
program exercise1
    implicit none
    include "mpif.h" 
    
    ! Данные
    real, allocatable :: A(:, :), B(:, :), C(:, :)
    integer :: m, n, k, myid, numprocs, ierr
    real(kind = 8) :: t1, t2, tl, tp
    interface
        subroutine printMatrix(B)
            real, dimension(:, :), intent(inout) :: B 
        end subroutine printMatrix
    end interface
    
    ! Действия
    call mpi_init(ierr)
    
    call mpi_comm_rank(mpi_comm_world, myid, ierr)
    
    ! Считывание данных
    if (myid .eq. 0) then
        call mpi_comm_size(mpi_comm_world, numprocs, ierr)
        
        write(*,  *) "                          Умножение матриц"
        write(*,  *) ""
        write(*,  *) "        Введите размерности матриц (A = m * n, B = n * k, C = m * k)"
        write(*,  *) ""
        
        write(*,  *) "Введите m ="   
        read(*, *) m
    
        write(*,  *) "Введите n ="
        read(*, *) n
    
        write(*,  *) "Введите k ="
        read(*, *) k
    end if
    
    ! Рассылка данных
    call mpi_bcast(m, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(n, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(k, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(numprocs, 1, mpi_integer, 0, mpi_comm_world, ierr)
    
    ! Выделение памяти под массив A
    allocate(A(1 : m, 1 : n))  
    
    if (myid .eq. 0) then    
        ! Выделение памяти под массивы B, C
        allocate(B(1 : n, 1 : k), C(1 : m, 1 : k))
    
        ! Заполнение массивов A, B (числами от 0 до 1)
        call random_number(A)
        call random_number(B)
    end if

    ! Рассылка матрицы A
    call mpi_bcast(A, m * n, mpi_real, 0, mpi_comm_world, ierr) 
    
    if (myid .eq. 0) then   
        ! Печать матриц (A, B)
        !print *, ""
        !print *, "Матрица A"
        !call printMatrix(A)
        !print *, ""
        !print *, "Матрица B"
        !call printMatrix(B)
        print *, ""
        print *, ""
    
        print *, "                  Проверка производительности алгоритмов "
        print *, ""
    
        print *, "Стандартный алгоритм умножения матриц"
        
        ! Применение стандартной последовательной процедуры 0 процессом
        t1 = mpi_wtime(ierr)
        call stdMM(A, B, C, m, k, n)
        t2 = mpi_wtime(ierr)
        tl = t2 - t1
        
        !call printMatrix(C)
        write(*, '(A, f8.3)') " последовательное время (секунды): ", tl
    end if
    
    ! Применение стандартной параллельной процедуры всеми процессами
    t1 = mpi_wtime(ierr)
    call stdPMM(A, B, C, m, k, n, myid, numprocs)
    t2 = mpi_wtime(ierr)
    tp = t2 - t1
    
    if (myid .eq. 0) then
        !call printMatrix(C)
        write(*, '(A, f8.3)') " параллельное время (секунды): ", tp
        write(*, '(A, f8.3)') " выигрыш параллельного варианта (секунды): ", tl - tp
        print *, ""
    
        print *, "Алгоритм Винограда умножения матриц"
        
        ! Применение последовательной процедуры Винограда 0 процессом
        t1 = mpi_wtime(ierr)
        call winogradMM(A, B, C, m, k, n)
        t2 = mpi_wtime(ierr)
        tl = t2 - t1
        
        !call printMatrix(C)
        write(*, '(A, f8.3)') " последовательное время (секунды): ", tl
    end if
    
    ! Применение параллельной процедуры Винограда всеми процессами
    t1 = mpi_wtime(ierr)
    call winogradPMM(A, B, C, m, k, n, myid, numprocs)
    t2 = mpi_wtime(ierr)
    tp = t2 - t1
    
    if (myid .eq. 0) then
        !call printMatrix(C)
        write(*, '(A, f8.3)') " параллельное время (секунды): ", tp
        write(*, '(A, f8.3)') " выигрыш параллельного варианта (секунды): ", tl - tp
        print *, ""
    end if
    
    call mpi_finalize(ierr)
end program exercise1

subroutine stdMM(A, B, C, nRows, nColumns, nAB)
    ! Описание переменных
    integer, intent(in) :: nRows, nColumns, nAB
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: A
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: B
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: C
    integer :: l, i, j
    
    ! Действия 
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
    ! Описание переменных
    integer, intent(in) :: nRows, nColumns, nAB
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: A
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: B
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: C
    real, dimension(1 : nRows) :: rowFactor
    real, dimension(1 : nColumns) :: columnFactor
    integer :: l, i, j, d 
    
    ! Действия   
    d = nAB / 2
    
    ! Вычисление rowFactor для матрицы A
    do i = 1, nRows
        rowFactor(i) = A(i, 1) * A(i, 2)
        do j = 2, d
            rowFactor(i) = rowFactor(i) + A(i, 2 * j - 1) * A(i, 2 * j)
        end do
    end do
    
    ! Вычисление columnFactor для матрицы B
    do i = 1, nColumns
        columnFactor(i) = B(1, i) * B(2, i)
        do j = 2, d
            columnFactor(i) = columnFactor(i) + B(2 * j - 1, i) * B(2 * j, i)
        end do
    end do
    
    ! Вычисление матрицы C
    do i = 1, nRows
        do j = 1, nColumns
            C(i, j) = - rowFactor(i) - columnFactor(j)
            do l = 1, d
                C(i, j) = C(i, j) + (A(i, 2 * l - 1) + B(2 * l, j)) * (A(i, 2 * l) + B(2 * l - 1, j))
            end do
        end do
    end do
          
    ! Прибавление членов в случае нечётной разметки
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
    
    ! Описание переменных
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: Ares
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: Bres
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: Cres
    integer, intent(in) :: nRows, nColumns, nAB, myid, numprocs
    real, allocatable :: B(:, :), C(:, :)
    integer, dimension(0 : numprocs - 1) :: scounts, displs, offsets
    integer :: i, ierr, offset
   
    ! Действия  
    ! Разделение по процессам (с передачей размера порций всем процессам)
    if (myid .eq. 0) call getStats(scounts, displs, offsets, nColumns, nAB, numprocs)
    call mpi_scatter(offsets, 1, mpi_integer, offset, 1, mpi_integer, 0, mpi_comm_world, ierr)
    
    ! Выделение памяти под массивы B, C
    allocate(B(1 : nAB, 1 : offset), C(1 : nRows, 1 : offset))
    
    ! Рассылка порций левой матрицы (Bres -> B) всем процессам
    call mpi_scatterv(Bres, scounts, displs, mpi_real, B, nAB * offset, mpi_real, 0, mpi_comm_world, ierr)
    
    ! Применение последовательной процедуры всеми процессами
    call stdMM(Ares, B, C, nRows, offset, nAB)
    
    ! Изменение параметров для сборки порций (С -< Cres)
    if (myid .eq. 0) call getStats(scounts, displs, offsets, nColumns, nRows, numprocs)
    call mpi_scatter(offsets, 1, mpi_integer, offset, 1, mpi_integer, 0, mpi_comm_world, ierr)
 
    ! Сборка порций
    call mpi_gatherv(C, nRows * offset, mpi_real, Cres, scounts, displs, mpi_real, 0, mpi_comm_world, ierr)   
end subroutine stdPMM

subroutine winogradPMM(Ares, Bres, Cres, nRows, nColumns, nAB, myid, numprocs)
    include 'mpif.h'
    
    ! Описание переменных
    real, dimension(1 : nRows, 1 : nAB), intent(in) :: Ares
    real, dimension(1 : nAB, 1 : nColumns), intent(in) :: Bres
    real, dimension(1 : nRows, 1 : nColumns), intent(out) :: Cres
    integer, intent(in) :: nRows, nColumns, nAB, myid, numprocs
    real, allocatable :: B(:, :), C(:, :)
    integer, dimension(0 : numprocs - 1) :: scounts, displs, offsets
    integer :: i, ierr, offset
    
    ! Действия 
    ! Разделение по процессам (с передачей размера порций всем процессам)
    if (myid .eq. 0) call getStats(scounts, displs, offsets, nColumns, nAB, numprocs)
    call mpi_scatter(offsets, 1, mpi_integer, offset, 1, mpi_integer, 0, mpi_comm_world, ierr)
    
    ! Выделение памяти под массивы B, C
    allocate(B(1 : nAB, 1 : offset), C(1 : nRows, 1 : offset))
    
    ! Рассылка порций левой матрицы (Bres -> B) всем процессам
    call mpi_scatterv(Bres, scounts, displs, mpi_real, B, nAB * offset, mpi_real, 0, mpi_comm_world, ierr)
    
    ! Применение последовательной процедуры Винограда всеми процессами
    call winogradMM(Ares, B, C, nRows, offset, nAB)
    
    ! Изменение параметров для сборки порций (С -< Cres)
    if (myid .eq. 0) call getStats(scounts, displs, offsets, nColumns, nRows, numprocs)
    call mpi_scatter(offsets, 1, mpi_integer, offset, 1, mpi_integer, 0, mpi_comm_world, ierr)
 
    ! Сборка порций
    call mpi_gatherv(C, nRows * offset, mpi_real, Cres, scounts, displs, mpi_real, 0, mpi_comm_world, ierr)
end subroutine winogradPMM

! Подпрограмма, для получения данных разделения по процессам
subroutine getStats(scounts, displs, offsets, nColumns, nRows, numprocs)  

    ! Описание переменных
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

! Подпрограмма для печати любой двумерной матрицы
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