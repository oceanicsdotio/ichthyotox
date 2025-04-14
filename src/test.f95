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

function random() result(c)
    use variables, only : sp
    use random, only : random_number_generator
    implicit none
    real(sp) :: c
    allocate(random_number_generator)
    c = random_number_generator%init()
end function
