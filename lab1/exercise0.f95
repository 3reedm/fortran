! ���������� ������ ����������� ���������
program exercise0
    implicit none
    
    ! ������ �������� ����������
    real :: a, b, c
    complex :: x1, x2
     
    ! ������ ����������
    print *, '� ��� ������� ��������� ����: ax^2 + bx + c = 0,  ��� a - ���������' 
    do while (abs(a - 0) .lt. 1e-10)
        print *, '������� ����������� [a] ��� [x^2] ([a] .ne. 0)'
        read *, a 
    end do
    print *, '������� ����������� [b] ��� [x]'
    read *, b
    print *, '������� ��������� ���� (�����������) [c]'
    read *, c
      
    call sqrsolution(a, b, c, x1, x2) 
    print *, '����� ��������� � ������� ����������� ����� (real, image):' 
    print *, 'x1 =', x1, '� x2 =', x2
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