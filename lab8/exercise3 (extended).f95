! ��������� ������� � ������������ ������ ��������� � ���������� ��� ���������� �����
program loop
    implicit none
    
    ! ������ �������� ����������
    integer :: i
    real :: nv
    real, dimension(1:4) :: v
    interface 
        real function norm(v) result(res)
            integer :: nRows, i
            real, dimension(:) :: v 
        end function norm
    end interface
    
    ! ������ ����������
    data v/4,3,2,1/
    
    print *, '������: ', (v(i), i = 1, 4)
    
    nv = norm(v)
    
    print *, ''
    print *, '��������� ����� �������:', nv
end program loop

real function norm(v) result(res)
    ! ������� ��������
    integer :: nRows, i
    real, dimension(:) :: v 
    
    nRows = size(v)
    forall(i=1:nRows) v(i) = v(i)**2
    res = sqrt(sum(v))
end function norm