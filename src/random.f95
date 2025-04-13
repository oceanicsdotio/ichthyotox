module random
    ! Module With Random Number Generator
    use variables, only : ZERO, sp
    implicit none
    save
    private

    ! Simple uniform using system methods
    type, public :: Uniform_Random_Number
        real(sp), private :: ru(2)
        logical, private  :: current
    contains
        procedure, public :: uniform => random_uniform
    end type
  
    ! Creating gaussian distributions
    type, public, extends(Uniform_Random_Number) :: Gaussian_Random_Number
        real(sp), private :: rn(2), sumMeanDiffSq, mean
        logical, private  :: statistics
        integer, private  :: samples
    contains
        ! Initialize random number generator
        procedure, public :: init => random_initialize
        ! Generates two gaussian randoms, stored in the struct, calculate statistics on the fly
        procedure, private :: normal => random_normal 
        ! Returns value of one stored gaussian, and trigger new calculation if needed
        procedure, public :: get => random_get 
        ! Prints summary statistics to logging output
        procedure, public :: stats => random_displayStatistics 
        ! Single random Gaussian value
        procedure, public :: gaussian => random_get
        ! Array of random values truncated to -1, 1
        procedure, public :: clipped => random_clipped_normal
    end type
    
    class(Gaussian_Random_Number), public, allocatable :: random_number_generator
  
  contains
    subroutine random_initialize(self)
        ! init random seed from computer clock
        class(Gaussian_Random_Number), intent(inout) :: self
        integer :: uu, seed_size
        integer(kind = 4) :: clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        call system_clock(clock)
        seed = clock + 37 * (/ (uu - 1, uu = 1, seed_size) /)
        call random_seed(put=seed)
        deallocate(seed)

        self%current = .true.
        self%samples = 0
        self%rn = (/ zero, zero /)
        self%ru = (/ zero, zero /)
        self%mean = zero
        self%sumMeanDiffSq = zero

        call self%normal()
    end subroutine
  
    subroutine random_normal(self)
        ! generate two random normal numbers using Box-Muller method
        class(Gaussian_Random_Number), intent(inout) :: self
        real(sp) :: sq, values(0:1), residual
    
        do
            values(:) = (/ self%uniform(), self%uniform() /)
            sq = sum(values*values)
            if (sq <= 1.0_sp) exit
        end do
    
        self%rn(:) = sqrt(-2.0_SP*log(sq)/sq)*values(:) ! generate normal rands from uniform
        self%samples = self%samples + 2
        residual = sum(self%rn - self%mean)
        self%mean = self%mean + residual / self%samples  ! TODO: make sure this is still correct
        self%sumMeanDiffSq = self%sumMeanDiffSq + residual*(sum(self%rn - self%mean))
    end subroutine

    real(sp) function random_uniform(self) result(rand)
        ! Use system level uniform `random_number()` transformed in -1 to 1 range
        class(Uniform_Random_Number), intent(inout) :: self
        call random_number(self%ru(1))
        rand = 2.0_SP*self%ru(1) - 1.0_SP
    end function  

    real(sp) function random_get(self)
        ! returns one of stored gaussian random numbers and generates new ones when used
        class(Gaussian_Random_Number), intent(inout) :: self
        random_get = merge(self%rn(1), self%rn(2), self%current)
        if (self%current) call self%normal()
        self%current = .not. self%current
    end function
  
    real(sp) function random_clipped_normal(self)
        ! get new random normal and reassign if out of range
        class(Gaussian_Random_Number), intent(inout) :: self
        random_clipped_normal = self%get() 
        do while (abs(random_clipped_normal) > 1.0_SP)
            random_clipped_normal = self%get() 
        end do
    end function

    subroutine random_displayStatistics(self)
        ! calculate and display distribution statistics for the random number system
        class(Gaussian_Random_Number), intent(in) :: self
        print *, "Statistics for random gaussian numbers.";
        print *, "    Samples:       ", self%samples
        print *, "    Mean:          ", self%mean
        print *, "    Variance:      ", self%sumMeanDiffSq / float(self%samples - 1)
        print *, "    Std Deviation: ", sqrt(self%sumMeanDiffSq / float(self%samples - 1))
    end subroutine
end module
