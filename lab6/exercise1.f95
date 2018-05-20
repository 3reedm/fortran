! Определение времени, необходимого на создание и уничтожение нити с помощью OpenMP
program exercise1
    implicit none 
    integer :: t1, t2, th_num
    
    write(*,  *) "   --Определение времени, необходимого для создания и уничтожения потока--" 
    write(*,  *) "" 

    th_num = 2
    call system_clock(t1) 
    !$omp parallel num_threads(th_num)
    !$omp end parallel
    call system_clock(t2) 
    
    write(*, '(A, f5.3, A)') " Время, потраченное на создание и уничтожение потока: ", real(t2-t1) / (th_num*1000), " секунд"
end program exercise1