! Функция для нахождения евклидовой нормы массива произвольной длины
real function norm(v) result(res)
    ! Входной параметр
    integer :: nRows, i
    real, dimension(:) :: v 
    
    nRows = size(v)
    forall(i=1:nRows) v(i) = v(i)**2
    res = sqrt(sum(v))
end function norm