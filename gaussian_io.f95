program gaussian

    use MOD_RAND, only : LAG_RAND, random
    implicit none

    integer :: ii

    allocate(random)
    call random%init()

    open(unit=101, file='./random.dat', status='replace')

    do ii = 1, 5000
        write(101, "(1F20.6)") random%get()
    end do

    call random%stats()

end program