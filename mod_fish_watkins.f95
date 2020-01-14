module MOD_FISH

  use MOD_PREC
  use MOD_LAG, only : LAG_OBJ
  use ALL_VARS, only : zero, pi, pi2
  implicit none
  save

  ! watkins and rose parameters
  real(sp), dimension(0:4), parameter :: speedtable = (/0.50_SP, 1.00_SP, 0.50_SP, 0.25_SP, 0.33_SP/)
  real(sp), dimension(0:4), parameter :: angletable = (/2.00_SP, 0.25_SP, 0.25_SP, 1.00_SP, 0.50_SP/)
  real(sp), dimension(0:1), parameter :: memory = (/0.5_SP, 0.96_SP/) ! unitless, memory coefficients, m0, m1
  real(sp), dimension(1:2), parameter :: threshold = (/0.005_SP*10.0_SP**(-6.0_SP), 0.5_SP/) ! detection thresholds, r1 and r2
  !real(sp), dimension(1:2), parameter :: weight = (/0.1_SP, 1.0_SP/)
  real(sp), dimension(1:2), parameter :: weight = (/0.7_SP, 1.0_SP/)
    
  real(sp), parameter :: initBodylength = 0.1_SP ! meters
  real(sp), parameter :: growthMax = 0.0025_SP * 12.0_SP * 0.001_SP ! conversion to meters per hour from mm per 5min
  real(sp), parameter :: util_cutoff = 0.01 ! level at which default behavior is chosen
  real(sp), parameter :: absorptionRate = 0.01_SP*10.0_SP*0.046748_SP ! grams of toxin / m^2 / hour / [toxin]
  real(sp), parameter :: depurationRate = 0.01_SP
  !real(sp), parameter :: depurationRate = 0.005_SP ! lower depuration rate by 50%, for sensitivity
  real(sp), parameter :: ingestionRate = 0.001_SP*0.02_SP
  real(sp), parameter :: toxfrac = 0.015_SP*10.0_SP**(-6.0_SP)
  real(sp), parameter :: speedimpair = 0.9_SP
  !real(sp), parameter :: speedimpair = 0.5_SP ! decrease speed when impaired, for sensitivity
  logical, parameter :: enforce_default = .false.
  logical, parameter :: no_flight = .false.
  logical, parameter :: ingestion_multiplier = .true.
  
  public LAG_FISH, FISH

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type, extends(LAG_OBJ) :: LAG_FISH

    ! fish behavioral data
    logical, allocatable, dimension(:), private :: impaired
    integer, allocatable, dimension(:), private :: last_rule
    real(sp), allocatable, dimension(:), private :: reverse
    real(sp), allocatable, dimension(:, :), private :: event ! dim: fish x agents
    real(sp), allocatable, dimension(:, :), private :: probability ! dim: fish x (agents x timescales)
    real(sp), allocatable, dimension(:, :), private :: utility ! dim: fish x (agents x timescales)

    ! fish  state variables
    real(sp), allocatable, dimension(:), private :: suitability ! spatial varying growth rate
    real(sp), allocatable, dimension(:), private :: length ! length of fish
    real(sp), allocatable, dimension(:), private :: effective_length ! impairment scalar
    real(sp), allocatable, dimension(:), private :: mass ! mass of fish
    real(sp), allocatable, dimension(:), private :: microcystin ! biocative microcystin
    real(sp), allocatable, dimension(:), private :: dissolved ! in situ toxin concentration
    real(sp), allocatable, dimension(:), private :: angle ! current facing
    real(sp), allocatable, dimension(:), private :: pathway

    contains
      ! When using procedures in main program, include subroutine calls in the order: dynamics, toxicity, movement
      ! Initialization aliases
      procedure, public :: init => fish_initialize ! allocate additioanl varaibles for lagrangian particles taking the fish type
      procedure, public :: writeState => fish_writeState

      ! Fish physiology procedure aliases
      procedure, private :: growth => fish_growth
      procedure, public :: movement => fish_movement ! subroutine collects all procedures necessary for behavior selection and movement under one subroutine call

      
  end type LAG_FISH; type(LAG_FISH), allocatable :: FISH
  contains

  ! INITIALIZATION / IO
  subroutine fish_initialize(self) ! OK
    ! Allocate fish specific variables
    use MOD_PREC, only : sp
    use MOD_RAND, only : random
    use ALL_VARS, only : pi2, zero
    class(LAG_FISH), intent(inout) :: self

    self%species = 'fish'
    call self%readPosition() ! finds ndrft and allocates position variables

    allocate( self%suitability(self%ndrft) ); self%suitability(:) = zero
    allocate( self%impaired(self%ndrft) ); self%impaired(:) = .false.
    allocate( self%event(self%ndrft, 2) ); self%event = 0
    allocate( self%probability(self%ndrft, 1:4) ); self%probability = zero
    allocate( self%utility(self%ndrft, 0:4) ); self%utility = zero
    allocate( self%length(self%ndrft) ); self%length = initBodylength
    allocate( self%effective_length(self%ndrft) ); self%effective_length = 1.0_sp
    allocate( self%mass(self%ndrft) ); self%mass = 2.0_SP * 10.0_SP**(-6.0_SP) * (1000.0*self%length(:))**(3.38_SP)
    allocate( self%microcystin(self%ndrft) ); self%microcystin = zero
    allocate( self%dissolved(self%ndrft) ); self%dissolved = zero
    allocate( self%angle(self%ndrft) ); self%angle = pi2*random%uniform()
    allocate( self%pathway(self%ndrft) ); self%pathway = zero
    allocate( self%last_rule(self%ndrft) ); self%last_rule = 0
    allocate( self%reverse(self%ndrft) ); self%reverse = zero
    
    

  end subroutine fish_initialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fish_writeState(self, fid) ! OK
    use MOD_SIM, only : domain ! domain structure for elapsed time
    use ALL_VARS, only : zero
    class(LAG_FISH), intent(in) :: self ! cyanobacteria extended type
    integer, intent(in) :: fid ! persistent file unit number
    integer :: ii
    
    write(fid, "(1F10.2,9000(I6,3F20.6))") domain%time, (self%itag(ii), self%mass(ii), 1000.0_SP*self%microcystin(ii), zero, ii=1,self%ndrft)
  
  end subroutine fish_writeState
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fish_growth(self) ! OK
    ! add consumed biomass to gut, and a portion of that to fish mass
    use MOD_PREC, only : sp ! for real precision
    use ALL_VARS, only : zero, dti
    use MOD_SIM, only : domain
    use MOD_TOX, only : cyano
    class(LAG_FISH), intent(inout) :: self ! uses only data contained in and inhereted by LAG_FISH
    real(sp), dimension(self%ndrft) :: growth, toxinLoad, old_mass, delta_mass ! array of carbon consumed by each fish particle
    real(sp), dimension(self%ndrft) :: absorption, depuration, ingestion, ptox
    integer :: ii

    toxinLoad(:) = sum(cyano%microcystin(:)) / sum(cyano%carbohydrate(:)+cyano%protein(:))
    old_mass(:) = 2.0_SP * 10.0_SP**(-6.0_SP) * (1000.0_SP*self%length(:))**(3.38_SP)
    growth(:) = growthMax*self%suitability(:) ! length growth of individuals in meters per hour, based on small pelagic fish
    self%length(:) = self%length(:) + growth(:)*dti ! update length
    self%mass(:) = 2.0_SP * 10.0_SP**(-6.0_SP) * (1000.0_SP*self%length(:))**(3.38_SP) ! update mass, based on small pelagic fish
    delta_mass(:) = self%mass(:) - old_mass(:)
    
    
    ingestion(:) =  delta_mass(:)*toxinLoad(:)*ingestionRate
    if (ingestion_multiplier) then
    	ingestion(:) = ingestion(:)*10.0_SP
    end if
    ptox(:) = self%zinterp(domain%verticaltox(:))/domain%meshArea
    absorption(:) = ptox(:)/(500.0_SP*500.0_SP)*absorptionRate*self%length(:)
    depuration(:) = depurationRate*self%microcystin(:)
    
    self%pathway = zero
    
    self%microcystin(:) = self%microcystin(:) + ingestion(:) + absorption(:)*dti - depuration(:)*dti 
    
  end subroutine fish_growth

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fish_movement(self)

      use MOD_PREC, only : sp
      use LIMS, only : M, KBM1
      use MOD_RAND, only : random
      use MOD_SIM, only : domain
      use ALL_VARS, only : pi, pi2, dti, VXMIN, VXMAX, VYMIN, VYMAX, zero

      class(LAG_FISH), intent(inout) :: self
      integer :: ii, jj, rule_index
      real(sp) :: maxutil = zero
      
      self%suitability(:) = 0.5_SP*(1.0_sp + sin( pi2*(self%xp(:)-VXMIN-(VXMAX/4.0_SP))/(VXMAX-VXMIN) ))
      
      do ii = 1, self%ndrft
        
        self%event(ii,:) = zero ; self%impaired(ii) = .false. ; self%effective_length(ii) = 1.0_SP
        if ( self%microcystin(ii) > threshold(1) ) self%event(ii, 1) = 1.0_SP ! intoxication/mortality
        if ( self%suitability(ii) > threshold(2) ) self%event(ii, 2) = 1.0_SP ! current suitability
        if ( self%microcystin(ii)/self%mass(ii) > toxfrac ) then ! induce impairment if toxin level above some threshold
          self%impaired(ii) = .true.
          self%effective_length(ii) = speedimpair
          
        end if
          
        self%probability(ii, 1) = (1.0_SP - memory(0))*self%event(ii, 1) + memory(0)*self%probability(ii, 1)
        self%probability(ii, 2) = (1.0_SP - memory(1))*self%event(ii, 1) + memory(1)*self%probability(ii, 2)
        self%probability(ii, 3) = (1.0_SP - memory(0))*self%event(ii, 2) + memory(0)*self%probability(ii, 3)
        self%probability(ii, 4) = (1.0_SP - memory(1))*self%event(ii, 2) + memory(1)*self%probability(ii, 4)
        self%utility(ii, 0) = 0.01
        self%utility(ii, 1) = weight(1)*self%probability(ii, 1)
        self%utility(ii, 2) = weight(1)*self%probability(ii, 2)
        self%utility(ii, 3) = weight(2)*self%probability(ii, 3)
        self%utility(ii, 4) = weight(2)*self%probability(ii, 4)
        
        maxutil = self%utility(ii,0)
        
        if (no_flight) then
            self%utility(ii,1) = zero
            self%utility(ii,2) = zero
        end if
        
        rule_index = 0
        do jj = 1,4
          if (self%utility(ii,jj) > maxutil) then
            maxutil = self%utility(ii,jj) ! select highest util, or default behaviors
            rule_index = jj
          end if
          
        end do
        
        
        if (enforce_default) then
            rule_index = 0
        end if
        
        
        if ((rule_index == 1).and.(self%reverse(ii) < 0.5_SP)) then
          self%reverse(ii) = 1.0_SP ! reverse direction for avoidance
        else
          self%reverse(ii) = zero
        end if
        
        self%last_rule(ii) = rule_index ! store last behavior
        if (self%impaired(ii)) then
          self%angle(ii) = self%angle(ii) + self%reverse(ii)*pi + random%uniform()*pi*angletable(rule_index) ! uniform random angle if impaired
        else
          self%angle(ii) = self%angle(ii) + self%reverse(ii)*pi + random%clipped()*pi*angletable(rule_index) ! clipped gaussian random angle if sober, more directed
        end if
      
        if (self%angle(ii) < -pi) self%angle(ii) = self%angle(ii) + pi2 ! normalize angles to -pi, pi]
        if (self%angle(ii) > pi) self%angle(ii) = self%angle(ii) - pi2
        self%xp(ii) = self%xp(ii) + cos(self%angle(ii)) * dti * 3600.0_sp * speedtable(rule_index) * self%length(ii) * self%effective_length(ii)
        self%yp(ii) = self%yp(ii) + sin(self%angle(ii)) * dti * 3600.0_sp * speedtable(rule_index) * self%length(ii) * self%effective_length(ii)
      end do
      
      ! when particles leave domain have them enter from the opposite boundary with no change to trajectory
      where (self%xp(:) > VXMAX)
        self%xp(:) = VXMIN + (self%xp(:) - VXMAX)
      elsewhere (self%xp(:) < VXMIN)
        self%xp(:) = VXMAX - (VXMIN - self%xp(:))
      end where
      
      where (self%yp(:) > VYMAX)
        self%yp(:) = VYMIN + (self%yp(:) - VYMAX)
      elsewhere (self%yp(:) < VYMIN)
        self%yp(:) = VYMAX - (VYMIN - self%yp(:))
      end where
      
      call self%growth()
      
  end subroutine fish_movement

end module MOD_FISH
