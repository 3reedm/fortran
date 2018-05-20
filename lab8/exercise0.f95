! Вычисление корней квадратного уравнения
program exercise0
    implicit none
    
    ! Раздел описания переменных
    real :: a, b, c
    complex :: x1, x2
     
    ! Раздел операторов
    print *, 'У вас имеется уравнение вида: ax^2 + bx + c = 0,  где a - ненулевой' 
    do while (abs(a - 0) .lt. 1e-10)
        print *, 'Введите коэффициент [a] при [x^2] ([a] .ne. 0)'
        read *, a 
    end do
    print *, 'Введите коэффициент [b] при [x]'
    read *, b
    print *, 'Введите свободный член (коэффициент) [c]'
    read *, c
      
    call sqrsolution(a, b, c, x1, x2) 
    print *, 'Корни уравнения в формате комплексных чисел (real, image):' 
    print *, 'x1 =', x1, 'и x2 =', x2
    pause
end program exercise0

subroutine sqrsolution(a, b, c, x1, x2)
    real, intent(in) :: a, b, c
    complex, intent(out) :: x1, x2
    complex :: D
        
    D = b * b - 4 * a * c
    x1 = (- b + sqrt(D)) / (2 * a)
    x2 = (- b - sqrt(D)) / (2 * a)
end subroutine 