! ��।������ �६���, ����室����� �� ᮧ����� � 㭨�⮦���� ��� � ������� OpenMP
program exercise1
    implicit none 
    integer :: t1, t2, th_num
    
    write(*,  *) "   --��।������ �६���, ����室����� ��� ᮧ����� � 㭨�⮦���� ��⮪�--" 
    write(*,  *) "" 

    th_num = 2
    call system_clock(t1) 
    !$omp parallel num_threads(th_num)
    !$omp end parallel
    call system_clock(t2) 
    
    write(*, '(A, f5.3, A)') " �६�, ����祭��� �� ᮧ����� � 㭨�⮦���� ��⮪�: ", real(t2-t1) / (th_num*1000), " ᥪ㭤"
end program exercise1