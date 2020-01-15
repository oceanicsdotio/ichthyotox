module MOD_TOX
  use parameters, only : sp ! real precision
  use MOD_LAG, only : LAG_OBJ ! module extends the lagrangian particle type
  use ALL_VARS, only : zero
  implicit none
  save ! object file will contain persistent data

  real(sp), parameter :: colonyBaseRadius = 75.0_SP * 10.0_SP**(-6.0_SP) ! meters
  real(sp), parameter :: tempRef = 25.0_SP ! reference temperature for limit fcn
  real(sp), parameter :: tempOpt = 28.0_SP ! optimal growth temperature
  real(sp), parameter :: tempLethal = 35.0_SP ! lethal temperature
  real(sp), parameter :: excretionFrac = 0.1_SP ! unitless
  real(sp), parameter :: fixationMax = 11.4_SP ! maximum carbon fixation, per hour rate
  real(sp), parameter :: fixationBeta = 0.02_SP ! shape factor in fixation calculation
  real(sp), parameter :: respirationBasic = 0.004_SP ! basic respiration rate, per hour
  real(sp), parameter :: respirationActive = 0.2_SP ! active respiration rate, unitless
  real(sp), parameter :: densityMax = 1150.0_SP ! maximum empirical density of algal cells, kg/m^3
  real(sp), parameter :: densityMin = 1037.0_SP ! minimum empirical density of algal cells, kg/m^3
  real(sp), parameter :: vesicleDensity = 150.0_SP ! density of gas filled vesicles, kg/m^3
  real(sp), parameter :: cellFrac = 0.25_SP ! fraction of colony volume composed of cell material, unitless
  real(sp), parameter :: carbonRatioMax = 4.0_SP ! maximum empirical ratio of carbon reservoirs in algal cells, unitless
  real(sp), parameter :: vesicleFrac = 0.08_SP ! fraction of cell volume occupied by vesicles, unitless
  real(sp), parameter :: irradOpt = 250.0_SP ! optimal irradiance W/M^-2
  real(sp), parameter :: synthesisMax = 0.05_SP ! per hour rate
  real(sp), parameter :: tempFcnAlpha = 0.286_SP ! shape factor, unitless
  real(sp), parameter :: tempFcnBeta = 0.05_SP   ! shape coefficient, unitless
  real(sp), parameter :: cellDensityCoefficient = 0.7_SP ! shape coefficient, unitless
  real(sp), parameter :: lightExtinctionBiomass = 14.0_SP ! light extinction due to overlying biomass
  real(sp), parameter :: lightAttenuationWater = 0.15_SP ! light extinction coefficient due to coastal waters
  real(sp), parameter :: shading_upscale = 1.0_SP

  !private ! all private by default
  public LAG_TOX, CYANO ! type and instance used by other modules and main program

  type, extends(LAG_OBJ) :: LAG_TOX ! derived from lagrangian particle class
    ! algal data
    real(sp) :: mclrProductionRate = zero
    real(sp) :: mclrExcretionRate = zero

    real(sp), allocatable, dimension(:), private :: radius ! radius for vertical movement of spherical colonies, meters
    real(sp), allocatable, dimension(:), private :: irradiance ! irradiance at particle position, watts per sq meter
    real(sp), allocatable, dimension(:), private :: biomass ! overlying biomass used for light attenuation calculations, gram
    real(sp), allocatable, dimension(:), public :: carbohydrate ! mass of carbohydration ballast, grams
    real(sp), allocatable, dimension(:), public :: protein ! mass of buoyant cellular material, grams
    real(sp), allocatable, dimension(:), public :: microcystin ! mass of cellular toxins, grams
    real(sp), allocatable, dimension(:), private :: delta_rho ! difference in density between water and colony

  contains
    ! user interface subroutines
    procedure, public :: init => colony_initialize ! read position and count, allocate state variables
    procedure, public :: writeState => colony_writeState ! write internal state variables
    procedure, public :: movement => colony_verticalMovement ! subroutine that updates position of colony particles
    procedure, public :: random => colony_random_walk ! vertical random walk for cyanobacteria

    ! private utility procedures
    procedure, private :: readState => colony_readState ! read state variables from initial value file
    procedure, private :: tempLimit => colony_temperatureLimit ! function that returns array scaling coefficients
    procedure, private :: tempFunction => colony_temperatureFunction ! function that returns array scaling coefficients

    ! movement function call aliases
    procedure, private :: velocity => colony_stokesVelocity ! function that returns array of particle vertical velocities
    procedure, private :: density => colony_algaeDensity ! function that returns array of colony densities
    procedure, private :: viscosity => colony_dynamicViscosity ! function that returns array of viscosity at particle locations

    ! mass transfer function call aliases
    procedure, private :: fixation => colony_carbonFixation ! subroutine adds mass to carbohydrate ballast
    procedure, private :: synthesis => colony_carbonSynthesis ! subroutine moves mass from carbohydrate to protein
    procedure, private :: excretion => colony_carbonExcretion ! subroutine subtracts mass from protein
    procedure, private :: respiration => colony_carbonRespiration ! subroutine subtracts mass from carbohydrate
    procedure, private :: production => colony_microcystinProduction ! constant production
    procedure, private :: release => colony_microcystinExcretion ! constant excretion

  end type LAG_TOX; type(LAG_TOX), allocatable :: CYANO

contains

  subroutine colony_initialize(self, exp_type) ! OK
    ! read position and state, and allocate internal variables
    use ALL_VARS, only : zero
    use parameters, only : sp
    class(LAG_TOX), intent(inout) :: self
    integer :: exp_type

    self%species = 'cyanobacteria' ! file name prefix
    if (exp_type == 1) then
      self%mclrProductionRate = zero
      self%mclrExcretionRate = zero
    else if (exp_type == 2) then
      self%mclrProductionRate = 0.0001_sp
      self%mclrExcretionRate = zero
    else if (exp_type == 3) then
      self%mclrProductionRate = 0.0001_sp
      self%mclrExcretionRate = 0.0001_sp
    else if (exp_type == 4) then
      self%mclrProductionRate = zero
      self%mclrExcretionRate = 0.0001_sp
    else
      write(*,*) "Problem with experiment type, stopping..."; stop
    end if

    call self%readPosition() ! finds ndrft and allocates position arrays
    allocate( self%radius(self%ndrft) ); self%radius = ZERO ! read from file
    allocate( self%irradiance(self%ndrft) ); self%irradiance = ZERO ! calculated at runtime
    allocate( self%biomass(self%ndrft) ); self%biomass = ZERO ! calculated at runtime
    allocate( self%carbohydrate(self%ndrft) ); self%carbohydrate = ZERO ! read from file
    allocate( self%protein(self%ndrft) ); self%protein = ZERO ! read from file
    allocate( self%microcystin(self%ndrft) ); self%microcystin = ZERO ! read from file
    allocate (self%delta_rho(self%ndrft) ); self%delta_rho = ZERO ! calculated at runtime
    call self%readState() ! read initial values from file

  end subroutine colony_initialize


  subroutine colony_readState(self) ! OK
    use parameters, only : iovar
    use ALL_VARS, only : folderprefix
    class(LAG_TOX), intent(inout) :: self
    integer :: ii, indexMatch
    character(len = 50) :: filename
    logical :: fexist

    ! Check for state variable file, allocate storage and read values
    write(filename, "(A)") "./"//trim(folderprefix)//"/"//trim(self%species)//'_var.dat' ! variables read from "cyanobacteria_var.dat"
    inquire(file=trim(filename), exist=fexist) ! check for file
    if (.not. fexist) then
      write(*, *) 'State variable file: ', filename,' does not exist, halting...'
      stop
    end if

    open(unit=iovar, file=filename, form='formatted')
    read(iovar, "(I6)") indexMatch
    if (indexMatch /= self%ndrft) then
      write(*, *) 'Dimensions of position and state variable files are not equal, halting...'
      stop
    end if

    do ii = 1, self%ndrft
      read(iovar, "(4F20.6)") self%radius(ii), self%carbohydrate(ii), self%protein(ii), self%microcystin(ii)
    end do
    close(iovar)
  end subroutine colony_readState


  subroutine colony_writeState(self, fid)
    use MOD_SIM, only : domain ! domain structure for elapsed time
    class(LAG_TOX), intent(in) :: self ! cyanobacteria extended type
    integer, intent(in) :: fid ! persistent file unit number
    integer :: ii

    write(fid, "(1F10.2,9000(I6,3F20.3))") domain%time, (self%itag(ii), self%carbohydrate(ii), self%protein(ii), self%microcystin(ii), ii=1,self%ndrft)

  end subroutine colony_writeState


  recursive subroutine colony_QsortC(absdepth, order) ! OK
    ! Recursive Fortran 95 quicksort routine sorts real numbers into ascending numerical order
    ! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
    ! Based on algorithm from Cormen et al., Introduction to Algorithms, 1997 printing
    ! Made F conformant by Walt Brainerd http://www.fortran.com/qsort_c.f95
    use parameters, only : sp
    real(sp), intent(inout), dimension(:) :: absdepth
    integer, intent(inout), dimension(:) :: order
    integer :: iq

    if (size(absdepth) > 1) then
      call colony_Partition(absdepth, order, iq)
      call colony_QsortC(absdepth(:iq-1), order(:iq-1))
      call colony_QsortC(absdepth(iq:), order(iq:))
    end if
  end subroutine colony_QsortC


  subroutine colony_Partition(A, B, marker) ! OK
    use parameters, only : sp

    real(sp), intent(inout), dimension(:) :: A
    integer, intent(inout), dimension(:) :: B
    integer, intent(out) :: marker
    integer :: ii, jj, itemp ! iterators and temporary index 
    real(sp) :: temp ! temporary sorting value
    real(sp) :: x ! pivot point
    x = A(1)
    ii= 0
    jj= size(A) + 1

    do
      jj = jj-1
      do
        if (A(jj) <= x) exit
        jj = jj-1
      end do
      ii = ii+1
      do
        if (A(ii) >= x) exit
        ii = ii+1
      end do
      if (ii < jj) then
        ! exchange A(ii) and A(jj)
        temp = A(ii)
        A(ii) = A(jj)
        A(jj) = temp
        itemp = B(ii)
        B(ii) = B(jj)
        B(jj) = itemp
      elseif (ii == jj) then
        marker = ii+1
        return
      else
        marker = ii
        return
      endif
    end do
  end subroutine colony_Partition


  function colony_carbonFixation(self) ! OK
    ! Updates carbohydrate ballast state variable due to fixation (alias is "fixation")
    use parameters, only : SP ! for single or double precision
    use MOD_SIM, only : domain

    class(LAG_TOX), intent(inout) :: self
    real(sp), dimension(self%ndrft) :: colony_carbonFixation ! result array
    real(sp), dimension(self%ndrft) :: irradRatio ! array of particle specific actual:optimal light ratios
    real(sp), dimension(self%ndrft) :: fixationCoef ! array transport for each particle
    integer :: ii
    real(sp), dimension(self%ndrft) :: proxy_depth, avg_self_shade ! copy array for updating depth by recursive binary partitioning
    integer, dimension(self%ndrft) :: indices ! copy array for sorting indices by recursive binary partitioning

    indices(:) = self%ITAG(:) ! copy of particle indices
    proxy_depth(:) = abs(self%ZP(:)) ! use sigma to avoid positive ZPT values in sorting
    call colony_QsortC(proxy_depth(:), indices(:)) ! sort indices by position from shallowest to deepest

    self%biomass(:) = (self%carbohydrate(:) + self%protein(:))/domain%meshArea ! contribution to shading by particle
    avg_self_shade(:) = (exp(-lightExtinctionBiomass*self%biomass(:)) - 1.0_SP) / (-lightExtinctionBiomass*self%biomass(:)) ! average area under self shading curve from zero to self biomass
    self%irradiance(indices(1)) = domain%globalIrradiance ! first particle unshaded
    do ii = 2, self%ndrft
      self%irradiance(indices(ii)) = self%irradiance(indices(ii-1)) * exp(-lightExtinctionBiomass*shading_upscale*self%biomass(indices(ii-1))) ! attenuate by next biomass
    end do
    self%irradiance(:) = self%irradiance(:) * avg_self_shade(:) ! multiply by particle self shading component
    self%irradiance(:) = self%irradiance(:) * exp(self%zpt(:)*lightAttenuationWater) ! find irradiance at particle position after attenuation
    irradRatio(:) = self%irradiance(:) / irradOpt ! substitution function
    fixationCoef(:) = (2.0_SP + fixationBeta) * irradRatio(:) / (irradRatio(:)**2.0_SP + fixationBeta*irradRatio(:) + 1.0_SP) ! scaling coefficient for transfer
    colony_carbonFixation(:) = fixationMax*fixationCoef(:)*self%protein(:)*(1.0_SP - vesicleFrac)*(carbonRatioMax - self%carbohydrate(:)/self%protein(:))/carbonRatioMax ! actual mass transfer

  end function colony_carbonFixation


  function colony_carbonSynthesis(self) ! OK
    ! updates carbohydrate and protein state variables due to synthesis transport (alias is "synthesis")
    ! calls tempLimit()
    use parameters, only : sp ! real precision
    class(LAG_TOX), intent(inout) :: self
    real(sp), dimension(self%ndrft) :: colony_carbonSynthesis

    colony_carbonSynthesis(:) =  self%carbohydrate(:) * synthesisMax * self%tempLimit()

  end function colony_carbonSynthesis


  function colony_carbonExcretion(self) ! OK
    ! update protein and dissolved pools due to excretion transport (alias is "excretion")
    ! temperature is tracked for all particles, so function uses algae array subset
    use parameters, only : sp ! for single or double precision
    class(LAG_TOX), intent(inout) :: self
    real(sp), dimension(self%ndrft) :: colony_carbonExcretion

    colony_carbonExcretion(:) = excretionFrac * self%tempFunction() * (respirationBasic*self%carbohydrate(:) + synthesisMax*self%protein(:))

  end function colony_carbonExcretion


  function colony_carbonRespiration(self)
    ! update carbohydrate and dissolved pools due to respiration transport (alias is "respiration")
    ! temperature is tracked for all particles, so function calls use algae array subset
    use parameters, only : sp ! for single or double precision
    class(LAG_TOX), intent(inout) :: self
    real(sp), dimension(self%ndrft) :: colony_carbonRespiration

    colony_carbonRespiration(:) = respirationBasic*self%tempFunction()*self%protein(:) + respirationActive*synthesisMax*self%tempLimit()*self%carbohydrate(:)

  end function


  function colony_temperatureLimit(self)
    ! Returns array of temperature limitation coefficents (0,1) for C synthesis
    use parameters, only : sp ! for single or double precision
    class(LAG_TOX), intent(in) :: self
    real(sp), dimension(self%ndrft) :: colony_temperatureLimit ! array of output coefficients for each particle

    colony_temperatureLimit(:) = (self%TEMP(:) / tempOpt * (  ((self%TEMP(:) - tempLethal)/(tempOpt - tempLethal))**( (tempRef - tempOpt) / tempOpt )  ))**(4.0_SP)

  end function


  function colony_temperatureFunction(self)
    ! returns array of scaling coefficents for biometric fcns
    use parameters, only : sp! for single and double precision

    class(LAG_TOX), intent(in) :: self
    real(sp), dimension(self%ndrft) :: colony_temperatureFunction ! array of output coefficients for each particles

    colony_temperatureFunction(:) = tempFcnAlpha * exp(tempFcnBeta * (self%temp(:) - tempOpt + tempRef))

  end function


  function colony_microcystinProduction(self)
    ! calculates toxin production per time step
    use parameters, only : sp ! for precision

    class(LAG_TOX), intent(inout) :: self
    real(sp), dimension(self%ndrft) :: colony_microcystinProduction ! array of microcystin production for colony particles
    colony_microcystinProduction(:) = self%mclrProductionRate * self%protein(:)

  end function colony_microcystinProduction


  function colony_microcystinExcretion(self)
    ! calculates temperature dependent toxin loss and moves mass to host element
    use parameters, only : sp

    class(LAG_TOX), intent(inout) :: self
    real(sp), dimension(self%ndrft) :: colony_microcystinExcretion
    colony_microcystinExcretion(:) = self%mclrExcretionRate * self%protein(:)

  end function colony_microcystinExcretion


  subroutine colony_verticalMovement(self)
    ! update position due to buoyant movement: calls velocity(), zinterp(), zlocate(), sigma()
    use parameters, only : sp ! single precision
    use ALL_VARS, only : KB, KBM1
    use MOD_SIM, only : domain
    use ALL_VARS, only : dti, zero ! integration step
    use parameters, only : A_RK, B_RK, MSTAGE, strict_integration ! runge-kutta integration parameters

    class(LAG_TOX), intent(inout) :: self
    integer :: ii, jj
    real(sp) :: mcoef
    real(sp), dimension(self%ndrft) :: idz, ini_position, ini_sigma, ini_microcystin, ini_carbohydrate, ini_protein ! initial valuess 
    real(sp), dimension(self%ndrft, 0:MSTAGE) :: chi_dissolved, chi_position, chi_carbohydrate, chi_protein, chi_microcystin ! stage fcn evaluations
    real(sp), dimension(self%ndrft) :: synthesis, excretion, respiration, fixation, mc_excretion, mc_production, calc_array
    chi_position=ZERO; chi_carbohydrate=ZERO; chi_protein=ZERO; chi_microcystin=abs(ZERO)


    ! save initial values
    ini_carbohydrate(:) = self%carbohydrate(:)
    ini_protein(:) = self%protein(:)
    ini_microcystin(:) = self%microcystin(:)
    ini_sigma(:) = self%zp(:)
    ini_position(:) = self%zpt(:)

    do ii = 1, MSTAGE
      ! get new stage values
      self%zpt(:)           = ini_position(:)     + A_RK(ii)*dti*chi_position(:, ii-1)      ! new stage position
      self%carbohydrate(:)  = ini_carbohydrate(:) + A_RK(ii)*dti*chi_carbohydrate(:, ii-1)  ! new carb stage value
      self%protein(:)       = ini_protein(:)      + A_RK(ii)*dti*chi_protein(:, ii-1)       ! new protein stage value
      self%microcystin(:)   = ini_microcystin(:)  + A_RK(ii)*dti*chi_microcystin(:, ii-1)   ! new toxin stage value

      ! enforce strict limits on state variables 
      self%carbohydrate(:)  = max(self%carbohydrate(:), zero) ! keep carbs non-negative
      self%protein(:)       = max(self%protein(:), zero) ! keep protein non-negative
      self%microcystin(:)   = max(self%microcystin(:), zero) ! keep toxin non-negative
      self%zpt(:)           = min(self%zpt(:), self%ep(:)) ! stop at surface
      self%zpt(:)           = max(self%zpt(:), self%hp(:)) ! stop at sediment

      ! update particle values at new stage for next calculation
      self%zp(:) = self%sigma(self%zpt(:)) ! update sigma position
      self%layer(:) = self%zlocate(self%zp) ! update layers
      self%rho(:) = self%zinterp(domain%verticalrho(:)) ! interp of density at stage position
      self%temp(:) = self%zinterp(domain%verticaltemp(:)) ! interp of temperature at stage position

      ! calculate mass transfer and update differential equations at stage starting position
      fixation(:) = self%fixation()
      synthesis(:) = self%synthesis()
      respiration(:) = self%respiration()
      excretion(:) = self%excretion()
      mc_excretion(:) = self%release()
      mc_production(:) = self%production()


      fixation(:) = self%fixation()
      synthesis(:) = self%synthesis()
      respiration(:) = self%respiration()
      excretion(:) = self%excretion()
      mc_excretion(:) = self%release()
      mc_production(:) = self%production()

      if (ii < MSTAGE) then
        mcoef = A_RK(ii+1)
      else
        mcoef = 1.0_SP
      end if

      calc_array(:) = synthesis(:)/respiration(:)
      where ( ((synthesis(:)+respiration(:))*dti*mcoef) > (self%carbohydrate(:) + fixation(:)*dti*mcoef) )
        synthesis(:) = (self%carbohydrate(:)/dti/mcoef + fixation(:))/(1.0_SP + 1.0_SP/calc_array(:))
        respiration(:) = (self%carbohydrate(:)/dti/mcoef + fixation(:))/(1.0_SP + calc_array(:))
      end where

      where ( excretion(:)*dti*mcoef > (self%protein(:) + synthesis(:)*dti*mcoef) )
        excretion(:) = self%protein(:)/dti + synthesis(:)
      end where

      where ( mc_excretion(:)*dti*mcoef > (self%microcystin(:) + mc_production(:)*dti*mcoef) )
        mc_excretion(:) = self%microcystin(:)/dti/mcoef + mc_production(:)
      end where

      chi_position(:, ii) = self%velocity()
      chi_carbohydrate(:,ii) = fixation(:) - respiration(:) - synthesis(:)
      chi_protein(:,ii) = synthesis(:) - excretion(:)
      chi_microcystin(:,ii) = mc_production(:) - mc_excretion(:)
      chi_dissolved(:,ii) = mc_excretion(:)

    end do


    ! restore initial values
    self%carbohydrate(:) = ini_carbohydrate(:)
    self%protein(:) = ini_protein(:)
    self%microcystin(:) = ini_microcystin(:)
    self%zp(:) = ini_sigma(:)
    self%zpt(:) = ini_position(:)

    do ii = 1, MSTAGE ! add weighted stages to initial values

      where (chi_dissolved(:, ii)*dti*B_RK(ii) > self%microcystin(:)) chi_dissolved(:,ii) = self%microcystin(:)/dti/B_RK(ii)
      where (-chi_carbohydrate(:, ii)*dti*B_RK(ii) > self%carbohydrate(:)) chi_carbohydrate(:,ii) = self%carbohydrate(:)/dti/B_RK(ii)
      where (-chi_protein(:, ii)*dti*B_RK(ii) > self%protein(:)) chi_protein(:,ii) = self%protein(:)/dti/B_RK(ii)

      self%zpt(:) = self%zpt(:) + chi_position(:,ii)*B_RK(ii)*dti ! update depth
      self%carbohydrate(:) = self%carbohydrate(:) + chi_carbohydrate(:, ii)*dti*B_RK(ii) ! update carbs
      self%protein(:) = self%protein(:) + chi_protein(:, ii)*dti*B_RK(ii) ! update protein
      self%microcystin(:) = self%microcystin(:) + chi_microcystin(:, ii)*dti*B_RK(ii) ! update toxin

      where (self%zpt(:) > self%ep(:))
        self%zpt(:) = self%ep(:)
      else where (self%zpt(:) < self%hp(:))
        self%zpt(:) = self%hp(:)
      end where


      ! update dissolved toxin  at particle position for each stage
      self%zp(:) = self%sigma(self%zpt(:)) ! update sigma
      self%layer(:) = self%zlocate(self%zp(:)) ! update layer
      idz(:) = float(KBM1)*abs((1.0_SP/float(KBM1) * (self%layer(:)-1)) - self%zp(:)) ! relative distance from layer above
      where (idz(:) > 1.0_SP)
        idz(:) = 1.0_SP
      elsewhere (idz(:) < zero)
        idz(:) = zero
      end where

      where (self%layer(:) == 1)
        domain%verticaltox(self%layer(:)) = domain%verticaltox(self%layer(:)) + 2.0_SP*B_RK(ii)*dti*chi_dissolved(:,ii)*(1.0_SP - idz(:))/domain%layerDepth ! add to sigma above
        domain%verticaltox(self%layer(:)+1) = domain%verticaltox(self%layer(:)+1) + B_RK(ii)*dti*chi_dissolved(:,ii)*idz(:)/domain%layerDepth ! add to sigma below
      elsewhere (self%layer(:) == KBM1)
        domain%verticaltox(self%layer(:)) = domain%verticaltox(self%layer(:)) + B_RK(ii)*dti*chi_dissolved(:,ii)*(1.0_SP - idz(:))/domain%layerDepth ! add to sigma above
        domain%verticaltox(self%layer(:)+1) = domain%verticaltox(self%layer(:)+1) + 2.0_SP*B_RK(ii)*dti*chi_dissolved(:,ii)*idz(:)/domain%layerDepth ! add to sigma below
      elsewhere
        domain%verticaltox(self%layer(:)) = domain%verticaltox(self%layer(:)) + B_RK(ii)*dti*chi_dissolved(:,ii)*(1.0_SP - idz(:))/domain%layerDepth ! add to sigma above
        domain%verticaltox(self%layer(:)+1) = domain%verticaltox(self%layer(:)+1) + B_RK(ii)*dti*chi_dissolved(:,ii)*idz(:)/domain%layerDepth ! add to sigma below
      end where

    end do

    self%carbohydrate(:) = max(self%carbohydrate(:), zero)
    self%protein(:) = max(self%protein(:), zero)
    where (self%microcystin(:) < zero) self%microcystin(:) = zero
    self%zpt(:) = min(self%zpt(:), self%ep(:)) ! collect at surface, mixing during random walk
    self%zpt(:) = max(self%zpt(:), self%hp(:)) ! collect at sediment, mixing during random walk
    self%zp(:) = self%sigma(self%zpt(:)) ! convert to sigma coordinate for finding current layer
    self%layer(:) = self%zlocate(self%zp(:)) ! update layer
    self%rho(:) = self%zinterp(domain%verticalrho(:)) ! vertical interp of density at final position
    self%temp(:) = self%zinterp(domain%verticaltemp(:)) ! vertical interp of density at final position

  end subroutine


  subroutine colony_random_walk(self)
    ! vertical and horizontal random walk
    use ALL_VARS, only : dti, dtrw, z, KB, KBM1, KBM2
    use MOD_RAND, only : random
    use MOD_SIM, only : domain
    class(LAG_TOX), intent(inout) :: self

    real(sp), dimension(self%ndrft) :: wdiff, kzp, dkzp ! diffusivity and first derivative at particle positions
    real(sp), dimension(KB) :: kspline_in, kspline_out, zspline, dkspline, smooth1, smooth2
    integer :: substeps, ii
    real(sp), parameter :: variance = 1.0_SP
    real(sp), parameter :: AC = 1.0_SP/6.0_SP
    real(sp), parameter :: BIG = 1.0E30

    ! Spline by Ross and Sharples (2004) to creates continuous and differentiable diffusivity profile, to meet time step criterion DT<<MIN(1/K")
    !smooth2(:) = domain%verticaldiff(1:KB) ! smoothed array of first sigma layer from initial values
    !    smooth2(2:KBM1) = (domain%verticaldiff(1:KBM2) + domain%verticaldiff(2:KBM1)) / 2.0_SP ! two point filter
    !    smooth1(:) = smooth2(:) ! copy first two point filter
    !    smooth1(2:KBM1) = AC*smooth2(1:KBM2) + (1.0_SP - 2.0_SP*AC)*smooth2(2:KBM1) + AC*smooth2(3:KB) ! three point smoothing filter [1/6 2/3 1/6]
    !
    !    do ii = 1, KB
    !       kspline_in(KB-ii+1) = smooth1(ii) ! reverse smoothed array for spline subroutine
    !       zspline(KB-ii+1) = z(ii) ! reverse order of sigma layer depths for spline subroutine
    !    end do
    !    call spline(zspline, kspline_in, KB, BIG, BIG, kspline_out) ! create spline, natural (K'=0) at boundaries
    !    do ii = 1, KB
    !      dkspline(ii) = kspline_out(KB-ii+1) ! revert output to original order
    !    end do

    !vertical random walk
    do substeps = 1, int(dti/dtrw)
      !kzp(:) = self%zinterp(smooth1) ! interpolate diffusivity at particle depth
      !dkzp(:) = self%zinterp(dkspline) ! interpolate first derivative at particle depth
      !kzp(:) = max(kzp(:), 0.0_SP) ! remove negative values
      kzp(:) = 60.0_SP*60.0_SP*10.0_SP**(-4.0_SP)*0.1
      dkzp(:) = 0.0_SP
      wdiff(:) = dkzp(:)*dtrw + random%array(self%ndrft)*sqrt((2.0_sp*kzp(:) + dkzp(:)**(2.0_SP))*dtrw/variance) ! Ross and Sharples 2004 eqn 1
      self%zpt(:) = self%zpt(:) + wdiff(:)
      self%zpt(:) = min(self%zpt(:), self%ep(:)) ! stop at surface during integration
      self%zpt(:) = max(self%zpt(:), self%hp(:)) ! stop at sediment during integration
      self%zp(:) = self%sigma(self%zpt(:))
      self%layer(:) = self%zlocate(self%zp(:))
    end do
  end subroutine


  function colony_stokesVelocity(self)
    ! returns stokes velocity of particle in m/hr, if lighter than water result is positive
    ! calls density() and viscosity()
    use parameters, only : sp ! precision
    use ALL_VARS, only : grav ! grav is positive, m/s2
    class(LAG_TOX), intent(inout) :: self
    real(sp), dimension(self%ndrft) :: colony_stokesVelocity ! output array of particle vertical velocities

    self%delta_rho(:) = self%rho(:) - self%density() ! water density array - colony density fcn
    colony_stokesVelocity(:) = 60.0_SP*60.0_SP * (2.0_SP/9.0_SP) * grav * self%radius(:)**(2.0_SP) * self%delta_rho(:) * self%viscosity()**(-1.0_SP)
  end function


  function colony_algaeDensity(self) ! OK
    ! returns actual colony density including contibutions of mucus and gas vacuoles
    use parameters, only : sp
    !use parameters, only : densityMin, densitymax, cellFrac, vesicleFrac, vesicleDensity, cellDensityCoefficient
    class(LAG_TOX), intent(in) :: self
    real(sp), dimension(self%ndrft) :: colony_algaeDensity ! output array of overall colony density
    real(sp), dimension(self%ndrft) :: cellDensity ! array of density without mucus and vacuoles

    cellDensity(:) = densityMin + (densityMax - densityMin)*(1.0_SP - exp(-cellDensityCoefficient * self%carbohydrate(:)/self%protein(:)))
    colony_algaeDensity(:) = (1.0_SP - cellFrac)*(self%rho(:) + 0.7_SP) + cellFrac*((1.0_SP - vesicleFrac)*cellDensity(:) + vesicleFrac*vesicleDensity)

  end function


  function colony_dynamicViscosity(self)
    ! returns array of dynamic viscosity values at particle locations
    use parameters, only : sp ! for precision
    class(LAG_TOX), intent(in) :: self
    real(sp), dimension(self%ndrft) :: colony_dynamicViscosity ! output array of viscosity values
    ! real(sp), dimension(self%ndrft) :: A, B, visc_pure
    ! Sharqway et al 2010
    !    A = 1.541_SP + 19.998_SP*10.0_SP**(-2.0_SP)*self%temp - 9.52_SP*10.0_SP**(-5.0_SP)*self%temp**(2.0_SP)
    !    B = 7.974_SP - 7.561_SP*10.0_SP**(-2.0_SP) + 4.724_SP*10.0_SP**(-4.0_SP)*self%temp**(2.0_SP)
    !    visc_pure = 4.2844_SP*10.0_SP**(-5.0_SP) + (0.157_SP*(self%temp+64.993_SP)**(2.0_SP)-91.296_SP)**(-1.0_SP)
    !    colony_dynamicViscosity = visc_pure*(1.0_SP + A*self%sal + B*self%sal**(2.0_SP))
    colony_dynamicViscosity(:) = 10.0_SP**(-3.0_SP) * 10.0_SP**(-1.65_SP + 262.0_SP/(self%temp(:) + 169.0_SP))
  end function

end module MOD_TOX
