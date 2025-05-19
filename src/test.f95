pure function add(a, b) result(c)
    use variables, only : sp
    implicit none
    real(sp), intent(in) :: a, b
    real(sp) :: c
    c = a + b
end function

pure function subtract(a, b) result(c)
    use variables, only : sp
    implicit none
    real(sp), intent(in) :: a, b
    real(sp) :: c
    c = a - b
end function

