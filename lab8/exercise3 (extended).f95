! Получение массива с произвольным числом элементов и нахождение его евклидовой нормы
program loop
    implicit none
    
    ! Раздел описания переменных
    integer :: i
    real :: nv
    real, dimension(1:4) :: v
    interface 
        real function norm(v) result(res)
            integer :: nRows, i
            real, dimension(:) :: v 
        end function norm
    end interface
    
    ! Раздел операторов
    data v/4,3,2,1/
    
    print *, 'Массив: ', (v(i), i = 1, 4)
    
    nv = norm(v)
    
    print *, ''
    print *, 'Евклидова норма массива:', nv
end program loop

real function norm(v) result(res)
    ! Входной параметр
    integer :: nRows, i
    real, dimension(:) :: v 
    
    nRows = size(v)
    forall(i=1:nRows) v(i) = v(i)**2
    res = sqrt(sum(v))
end function norm