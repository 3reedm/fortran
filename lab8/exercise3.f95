! ������� ��� ���������� ���������� ����� ������� ������������ �����
real function norm(v) result(res)
    ! ������� ��������
    integer :: nRows, i
    real, dimension(:) :: v 
    
    nRows = size(v)
    forall(i=1:nRows) v(i) = v(i)**2
    res = sqrt(sum(v))
end function norm