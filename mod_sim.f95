module MOD_SIM

  use MOD_PREC, only : sp
  use ALL_VARS, only : ZERO
  implicit none
  save ! State is saved in the compiled object

  private
  public LAG_SIM, domain ! only type information and 

  type LAG_SIM
    character(len = 100), public :: simID
    integer, public :: nnodes=0, nelements=0, nlayers=0, lines_read=0 ! simulation id, size and state
    real(sp), public :: globalIrradiance=ZERO, meshArea=ZERO, layerDepth=ZERO, layerSigma=ZERO, time=ZERO, daytime=ZERO, clocktime=ZERO ! domain variables and clocks: elapsed, divided days, and twenty four hour periodic
    real(sp), allocatable, dimension(:), private :: elementSigmaVolume, elementArea ! mesh stats
    real(sp), allocatable, dimension(:), public :: verticaltox, verticaldiff, verticaltemp, verticalrho ! uniform horizontal fields

  contains
    ! call in this order
    procedure, public :: init => simulation_initialize
    procedure, public :: load => simulation_read
    procedure, public :: geo => simulation_geometry
    procedure, public :: vdiff => simulation_diffusion

  end type LAG_SIM; class(LAG_SIM), allocatable :: domain ! domain structure imported from this module
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine simulation_initialize(self, exp_type)
    ! initialize variables for current simulation taking id number as input
    use MOD_PREC, only : sp
    use ALL_VARS, only : zero
    class(LAG_SIM), intent(inout) :: self
    integer, intent(in) :: exp_type


    allocate(self%verticaltox(0:self%nlayers+1));
    if (exp_type == 4) then
      self%verticaltox(:) = 6324.0_SP ! experiment D
    else
      self%verticaltox(:) = zero ! experiments A-C
    end if

    allocate(self%elementArea(0:self%nelements)); self%elementArea(:)=zero
    allocate(self%elementSigmaVolume(0:self%nelements)); self%elementSigmaVolume(:)=zero
    allocate(self%verticaldiff(0:self%nlayers+1)); self%verticaldiff=zero
    allocate(self%verticaltemp(0:self%nlayers+1)); self%verticaltemp=zero
    allocate(self%verticalrho(0:self%nlayers+1)); self%verticalrho=zero

  end subroutine simulation_initialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine simulation_geometry(self, vertx, verty, node_indices, bathymetry)
    ! called once during simulation setup, calculates area of any triangular mesh or subregion, 
    use MOD_PREC, only : sp
    use LIMS, only : KBM1
    class(LAG_SIM), intent(inout) :: self
    real(sp), dimension(0:self%nnodes), intent(in) :: vertx, verty
    integer, dimension(0:self%nelements, 4), intent(in) :: node_indices
    real(sp), dimension(0:self%nnodes) :: bathymetry
    real(sp), dimension(0:self%nelements) :: length1, length2, length3, ss, average_thickness
    integer :: ii

    do ii = 1, self%nelements
      length1(ii) = sqrt(  (vertx(node_indices(ii, 1)) - vertx(node_indices(ii, 2)) )**(2.0_sp) + ( verty(node_indices(ii, 1)) - verty(node_indices(ii, 2)) )**(2.0_sp)  )
      length2(ii) = sqrt(  (vertx(node_indices(ii, 2)) - vertx(node_indices(ii, 3)) )**(2.0_sp) + ( verty(node_indices(ii, 2)) - verty(node_indices(ii, 3)) )**(2.0_sp)  )
      length3(ii) = sqrt(  (vertx(node_indices(ii, 3)) - vertx(node_indices(ii, 1)) )**(2.0_sp) + ( verty(node_indices(ii, 3)) - verty(node_indices(ii, 1)) )**(2.0_sp)  )
      average_thickness = abs(sum( bathymetry(node_indices(ii,1:3)) )/3.0_sp) ! must be positive
    end do

    ss(:) = 0.5_sp*(length1(:) + length2(:)+ length3(:))
    self%elementArea(:) = sqrt( ss * (ss(:) - length1(:)) * (ss(:) - length1(:)) * (ss(:) - length3(:)) ) ! herons's formula for area
    self%meshArea = sum(self%elementArea(:))
    self%elementSigmaVolume(:) = self%elementArea(:)*average_thickness(:)
    self%layerDepth = average_thickness(1)
    self%layerSigma = float(KBM1)**(-1.0_SP)

  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine simulation_read(self, u_vel, v_vel, w_vel, diffusivity, elevation, salinity, temperature, density)

    use MOD_PREC, only : sp
    use ALL_VARS, only : ZERO
    use LIMS, only : KB, M, N
    use parameters, only : iophys

    class(LAG_SIM), intent(inout) :: self
    real(sp), dimension(0:N, KB), intent(inout) :: u_vel, v_vel, w_vel
    real(sp), dimension(0:M, KB), intent(inout) :: diffusivity, salinity, temperature, density
    real(sp), dimension(0:M), intent(inout) :: elevation

    character(len = 100) :: vert_format
    real(sp) :: time
    integer :: ii

    write(vert_format, "(A7,I6,A7)") "(F10.3,", 3*KB, "F20.10)"
    read(iophys, vert_format) time, self%verticaltemp(1:KB), self%verticalrho(1:KB), self%verticaldiff(1:KB)

    u_vel(:,:)=ZERO; v_vel(:,:)=ZERO; w_vel(:,:)=ZERO
    elevation(:)=-abs(ZERO); salinity(:,:)=ZERO

    node_based: do ii = 1, self%nnodes
      temperature(ii, 1:KB) = self%verticaltemp(1:KB)
      diffusivity(ii, 1:KB) = self%verticaldiff(1:KB)
      density(ii, 1:KB) = self%verticalrho(1:KB)
    end do node_based

    if (self%lines_read == 0) rewind(unit=iophys)
    self%lines_read = self%lines_read + 1

  end subroutine simulation_read

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine simulation_diffusion(self) ! NK 2/11/16
    ! one dimensional vertical diffusion
    use MOD_PREC, only : sp ! for real precision
    use LIMS, only : KB, KBM1
    use ALL_VARS, only : dti, ZERO ! time interpolation step

    class(LAG_SIM), intent(inout) :: self ! domain structure
    integer :: ii, jj, nsteps ! iteration parameters
    real(sp) :: dt_substep, dt_stable, secondDerivative, wdiff ! finite element variables
    real(sp), dimension(0:self%nlayers+1) :: TP1; TP1=ZERO ! dissolved profile at time plus one

    wdiff = 60.0_SP*60.0_SP*10.0_SP**(-5.0_SP) ! includes hour and meter conversion from 0.1 cm^2/s
    dt_stable = 0.5_SP*abs(self%layerDepth)**(2.0_SP) ! longest step for conditional stability
    nsteps = ceiling(dti/dt_stable) ! min number of steps to achieve stability, cannot be less than one
    dt_substep = dti/float(nsteps) ! automatically substeps nsteps times, only valid for dK/dt=0

    time_loop: do ii = 1, nsteps
      self%verticaltox(1) = self%verticaltox(2); self%verticaltox(KB) = self%verticaltox(KBM1)
      do jj = 2, KBM1
        secondDerivative = (self%verticaltox(jj-1) + self%verticaltox(jj+1) - 2.0_SP*self%verticaltox(jj)) / self%layerDepth**(2.0_SP) ! central in space, forward in time second derivative with depth
        TP1(jj) = self%verticaltox(jj) + wdiff*dt_substep*secondDerivative
      end do
      TP1(1) = TP1(2); TP1(KB) = TP1(KBM1) ! copy in domain value to boundary nodes
      self%verticaltox(:) = TP1(:)
    end do time_loop

  end subroutine simulation_diffusion
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module MOD_SIM


module MOD_RAND

  ! random number class, used for economically creating gaussian distributions, etc.
  use MOD_PREC, only : sp ! for real precision
  implicit none
  save
  private
  public LAG_RAND, random

  type :: LAG_RAND

    real(sp), private :: rn1, rn2, ru1, ru2
    logical, private  :: current, statistics
    real(sp), private :: sumMeanDiffSq, mean
    integer, private  :: samples

  contains

    procedure, public :: init => random_initialize ! subroutine initializes generator from system clock
    procedure, private :: normal => random_normal ! subroutine generates two gaussian randoms which are stored in the object, also calculates statisitcs on the fly
    procedure, public :: get => random_get ! function returns the value of one stored gaussian, and will trigger a new calculatlion when necessary
    procedure, public :: gaussian => random_get
    procedure, public :: array => random_array
    procedure, public :: uniform => random_uniform
    procedure, public :: clipped => random_clipped_normal
    procedure, public :: stats => random_displayStatistics ! subroutine prints summary statistics to command line
    procedure, private :: test => random_test ! subroutine iterates over N calls, then prints summary statistics

  end type LAG_RAND

  class(LAG_RAND), allocatable :: random

contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine random_initialize(self)

    ! init random seed from computer clock
    class(LAG_RAND), intent(inout) :: self
    integer :: uu, seed_size
    integer(kind = 4) :: clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = seed_size)
    allocate(seed(seed_size))
    call system_clock(clock)
    seed = clock + 37 * (/ (uu - 1, uu = 1, seed_size) /)
    call random_seed(put = seed)
    deallocate(seed)

    self%current = .true.
    self%statistics = .true.
    self%samples = 0
    self%rn1 = 0.0_SP
    self%rn2 = 0.0_SP
    self%ru1 = 0.0_SP
    self%ru2 = 0.0_SP
    self%mean = 0.0_SP
    self%sumMeanDiffSq = 0.0_SP

    call self%normal()

  end subroutine random_initialize
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine random_normal(self) ! OK
    ! generate two random normal numbers using Box-Muller method
    use MOD_PREC, only : sp
    class(LAG_RAND), intent(inout) :: self
    real(sp) :: V1, V2, S, meanDiff

    do
      V1 = self%uniform() ! uniform random centered on zero
      V2 = self%uniform()
      S = V1**2.0_SP + V2**2.0_SP ! check the sum of the squares and reject if outside range
      if (S <= 1.0_SP) exit
    end do

    self%rn1 = sqrt(-2.0_SP*log(S)/S)*V1 ! generate normal rands from uniform
    self%rn2 = sqrt(-2.0_SP*log(S)/S)*V2

    if (self%statistics) then ! calculate mean and stddev on the fly if desired
      self%samples = self%samples + 1
      meanDiff = self%rn1 - self%mean
      self%mean = self%mean + meanDiff/self%samples
      self%sumMeanDiffSq = self%sumMeanDiffSq + meanDiff*(self%rn1 - self%mean)

      self%samples = self%samples + 1
      meanDiff = self%rn2 - self%mean
      self%mean = self%mean + meanDiff/self%samples
      self%sumMeanDiffSq = self%sumMeanDiffSq + meanDiff*(self%rn2 - self%mean)
    end if

  end subroutine random_normal
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function random_array(self, nn)

    use MOD_PREC, only : sp
    class(LAG_RAND), intent(inout) :: self
    integer, intent(in) :: nn

    real(sp), dimension(nn) :: random_array
    real(sp), dimension(nn) :: temp_array
    integer :: ii

    do ii = 1, nn
      temp_array(ii) = self%get()
    end do

    random_array(:) = temp_array(:)

  end function random_array
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(sp) function random_uniform(self) result(rand)

    use MOD_PREC, only : sp
    class(LAG_RAND), intent(inout) :: self

    call random_number(self%ru1)

    rand = 2.0_SP*self%ru1 - 1.0_SP ! returns uniform pseudorandom in -1 to 1 range

  end function random_uniform
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(sp) function random_get(self) result(rand)! OK
    ! returns one of stored gaussian random numbers and generates new ones when used
    use MOD_PREC, only : sp
    class(LAG_RAND), intent(inout) :: self

    if (self%current) then ! if stored randoms haven't been used,
      rand = self%rn1 ! return the first
      self%current = .false.
    else ! if one stored random has been used,
      rand = self%rn2 ! return the second
      self%current = .true.
      call self%normal() ! generate new numbers
    end if

  end function random_get
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(sp) function random_clipped_normal(self) result(rand)! OK
    ! returns one of stored gaussian random numbers and generates new ones when used
    use MOD_PREC, only : sp
    class(LAG_RAND), intent(inout) :: self
    real(sp) :: deviate

    deviate = self%get() ! get new random normal
    do while (abs(deviate) > 1.0_SP) ! check if -1,1
      deviate = self%get() ! reassign if out of range
    end do
    rand = deviate

  end function random_clipped_normal


  subroutine random_displayStatistics(self) ! OK
    ! calculate and display distribution statistics for the random number system
    use MOD_PREC ! for real precision
    implicit none

    class(LAG_RAND), intent(in) :: self
    write(*, *); write(*, *) "Statistics for random gaussian numbers generated so far..."; write(*, *)
    write(*, *) "    Samples:       ", self%samples
    write(*, *) "    Mean:          ", self%mean
    write(*, *) "    Variance:      ", self%sumMeanDiffSq / float(self%samples - 1)
    write(*, *) "    Std Deviation: ", sqrt(self%sumMeanDiffSq / float(self%samples - 1))

  end subroutine random_displayStatistics


  subroutine random_test(self) ! OK
    ! iteratively generate random gaussian numbers and calculate statistics
    use MOD_PREC, only : sp ! for real precision

    class(LAG_RAND), intent(inout) :: self
    integer :: ii, nSample = 1000
    real(sp) :: randomNumber

    do ii = 1, nSample-1
      randomNumber = self%get() ! return next gaussian random
    end do
    call self%stats() ! display distribution statistics

  end subroutine

end module MOD_RAND