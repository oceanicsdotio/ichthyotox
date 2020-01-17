module behavior

  use ALL_VARS
  use MOD_LAG, only : LAG_OBJ
  use ALL_VARS, only : zero, pi, pi2
  implicit none
  save

  real(sp), dimension(0:4), parameter :: speedtable = (/0.50_SP, 1.00_SP, 0.50_SP, 0.25_SP, 0.33_SP/)
  real(sp), dimension(0:4), parameter :: angletable = (/2.00_SP, 0.25_SP, 0.25_SP, 1.00_SP, 0.50_SP/)
  real(sp), dimension(0:1), parameter :: memory = (/0.5_SP, 0.96_SP/) ! unitless, memory coefficients
  real(sp), dimension(1:2), parameter :: threshold = (/0.005_SP*10.0_SP**(-6.0_SP), 0.5_SP/) ! detection thresholds
  real(sp), dimension(1:2), parameter :: weight = (/0.7_SP, 1.0_SP/) ! sensitivity analyis @ (/0.1_SP, 1.0_SP/)

  real(sp), parameter :: initBodylength = 0.1_SP ! meters
  real(sp), parameter :: growthMax = 0.0025_SP * 12.0_SP * 0.001_SP ! conversion to meters per hour from mm per 5min
  real(sp), parameter :: util_cutoff = 0.01 ! level at which default behavior is chosen
  real(sp), parameter :: absorptionRate = 0.01_SP*10.0_SP*0.046748_SP ! grams of toxin / m^2 / hour / [toxin]
  real(sp), parameter :: depurationRate = 0.01_SP ! sensitivity analysis @ 0.005
  real(sp), parameter :: ingestionRate = 0.001_SP*0.02_SP
  real(sp), parameter :: toxfrac = 0.015_SP*10.0_SP**(-6.0_SP)
  real(sp), parameter :: speedimpair = 0.9_SP ! sensitivity analysis @ 0.5

  logical, parameter :: enforce_default = .false.
  logical, parameter :: no_flight = .false.
  logical, parameter :: ingestion_multiplier = .true.

  public LAG_FISH, FISH

  type, extends(LAG_OBJ) :: LAG_FISH

    logical, allocatable, dimension(:), private :: impaired
    integer, allocatable, dimension(:), private :: last_rule

    real(sp), allocatable, dimension(:), private :: &
            & reverse, &
            & suitability, & ! spatial varying growth rate
            & length, &
            & effective_length, & ! impairment scalar
            & mass, &
            & microcystin, & ! body toxin
            & dissolved, & ! in situ toxin concentration
            & angle, &
            & pathway

    real(sp), allocatable, dimension(:, :), private :: &
            & event, &        ! fish x agents
            & probability, &  ! fish x (agents x timescales)
            & utility         ! fish x (agents x timescales)

  contains
    ! Call in the order: dynamics, toxicity, movement
    procedure, public :: init => fish_initialize
    procedure, public :: writeState => fish_writeState
    procedure, public :: movement => fish_movement ! behavior selection and movement

  end type LAG_FISH; type(LAG_FISH), allocatable :: FISH

contains

  subroutine fish_initialize(self)

    use MOD_RAND, only : random
    use ALL_VARS, only : pi2, zero, sp
    class(LAG_FISH), intent(inout) :: self

    self%species = 'fish'
    call self%readPosition() ! read particles counts and allocates position variables

    allocate( self%suitability(self%ndrft), &
            & self%impaired(self%ndrft), &
            & self%microcystin(self%ndrft), &
            & self%dissolved(self%ndrft), &
            & self%angle(self%ndrft), &
            & self%pathway(self%ndrft), &
            & self%last_rule(self%ndrft), &
            & self%reverse(self%ndrft), &
            & self%mass(self%ndrft), &
            & self%length(self%ndrft), &
            & self%utility(self%ndrft, 0:4), &
            & self%probability(self%ndrft, 1:4), &
            & self%event(self%ndrft, 2), &
            & self%effective_length(self%ndrft))

    self%impaired = .false.
    self%event = 0
    self%last_rule = 0

    self%angle = pi2*random%uniform()
    self%length = initBodylength
    self%effective_length = 1.0_sp
    self%mass = 2.0_SP * 10.0_SP**(-6.0_SP) * (1000.0*self%length(:))**(3.38_SP)

    self%microcystin = zero
    self%dissolved = zero
    self%pathway = zero
    self%probability = zero
    self%utility = zero
    self%reverse = zero
    self%suitability = zero

  end subroutine

  subroutine fish_writeState(self, fid)
    use MOD_SIM, only : domain ! domain structure for elapsed time
    use ALL_VARS, only : zero
    class(LAG_FISH), intent(in) :: self ! cyanobacteria extended type
    integer, intent(in) :: fid ! persistent file unit number
    integer :: ii

    write(fid, "(1F10.2,9000(I6,3F20.6))") domain%time, &
            & (self%itag(ii), self%mass(ii), 1000.0_SP*self%microcystin(ii), zero, ii=1,self%ndrft)

  end subroutine


  subroutine fish_movement(self)

    use ALL_VARS, only : sp
    use MOD_RAND, only : random
    use MOD_SIM, only : domain
    use MOD_TOX, only : cyano
    use ALL_VARS, only : pi, pi2, dti, vxmin, vxmax, vymin, vymax, zero, M, KBM1

    class(LAG_FISH), intent(inout) :: self
    integer :: ii, jj, rule_index
    real(sp) :: maxutil, speed
    real(sp), dimension(self%ndrft) :: absorption, depuration, ingestion, mass

    maxutil = zero
    speed = zero
    self%suitability(:) = 0.5_SP*(1.0_sp + sin(pi2*(self%xp(:) - vxmin - (vxmax/4.0_SP))/(vxmax - vxmin)))

    where (self%microcystin(:) / self%mass(:) > toxfrac) ! induce impairment if toxin level above some threshold
      self%impaired(:) = .true.
      self%effective_length(:) = speedimpair
    elsewhere
      self%impaired(:) = .false.
      self%effective_length(:) = 1.0_SP
    end where

    self%event(:, :) = zero
    where (self%microcystin(:) > threshold(1)) self%event(:, 1) = 1.0_SP ! intoxication/mortality
    where (self%suitability(:) > threshold(2)) self%event(:, 2) = 1.0_SP ! current suitability

    do ii = 1, self%ndrft

      self%probability(ii, 1:2) = (1.0_SP - memory(0:1))*self%event(ii, 1) + memory(0:1)*self%probability(ii, 1:2)
      self%probability(ii, 3:4) = (1.0_SP - memory(0:1))*self%event(ii, 2) + memory(0:1)*self%probability(ii, 3:4)

      self%utility(ii, 0) = 0.01
      self%utility(ii, 1:2) = weight(1)*self%probability(ii, 1:2)
      self%utility(ii, 3:4) = weight(2)*self%probability(ii, 3:4)
      
      if (no_flight) then
        self%utility(ii, 1) = zero
        self%utility(ii, 2) = zero
      end if

      rule_index = 0 ! default behavior
      maxutil = self%utility(ii, rule_index)
      if (.not. enforce_default) then
        do jj = 1, 4
          if (self%utility(ii, jj) > maxutil) then
            maxutil = self%utility(ii, jj) ! select highest util, or default behaviors
            rule_index = jj
          end if
        end do
      end if

      self%last_rule(ii) = rule_index ! store last behavior
      self%reverse(ii) = merge(1.0_SP, zero,((rule_index == 1) .and. (self%reverse(ii) < 0.5_SP))) ! reverse direction for avoidance
      self%angle(ii) = self%angle(ii) + self%reverse(ii)*pi + &
              &merge(random%uniform(), random%clipped(), self%impaired(ii)) * pi * angletable(rule_index)

      if (self%angle(ii) < -pi) then
        self%angle(ii) = self%angle(ii) + pi2 ! normalize angles to -pi, pi]
      elseif (self%angle(ii) > pi) then
        self%angle(ii) = self%angle(ii) - pi2
      end if
      
      speed = dti * 3600.0_sp * speedtable(rule_index) * self%length(ii) * self%effective_length(ii)
      self%xp(ii) = self%xp(ii) + cos(self%angle(ii)) * speed
      self%yp(ii) = self%yp(ii) + sin(self%angle(ii)) * speed
      
    end do

    ! wrap positions
    where (self%xp(:) > vxmax)
      self%xp(:) = vxmin + (self%xp(:) - vxmax)
    elsewhere (self%xp(:) < vxmin)
      self%xp(:) = vxmax - (vxmin - self%xp(:))
    end where

    where (self%yp(:) > vymax)
      self%yp(:) = vymin + (self%yp(:) - vymax)
    elsewhere (self%yp(:) < vymin)
      self%yp(:) = vymax - (vymin - self%yp(:))
    end where

    ! add consumed biomass to gut, and a portion of that to fish mass
    ! length growth of individuals in meters per hour based on small pelagic fish
    self%length(:) = self%length(:) + growthMax * self%suitability(:) * dti
    mass(:) = 2.0_SP * 10.0_SP**(-6.0_SP) * (1000.0_SP*self%length(:))**(3.38_SP)
    ingestion(:) = (self%mass(:) - mass(:)) * ingestionRate * merge(10.0_SP, 1.0_SP, ingestion_multiplier) * &
            & sum(cyano%microcystin(:)) / sum(cyano%carbohydrate(:) + cyano%protein(:))

    depuration(:) = depurationRate * self%microcystin(:)
    absorption(:) = self%zinterp(domain%verticaltox(:)) / domain%meshArea / 500.0_SP * 500.0_SP * &
            & absorptionRate * self%length(:)

    self%microcystin(:) = self%microcystin(:) + ingestion(:) + (absorption(:) - depuration(:))*dti
    self%mass(:) = mass(:)
    self%pathway = zero

  end subroutine

end module