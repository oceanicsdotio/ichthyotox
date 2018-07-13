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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine colony_initialize(self, exp_type) ! OK
    ! read position and state, and allocate internal variables
    use ALL_VARS, only : zero
    use MOD_PREC, only : sp
    class(LAG_TOX), intent(inout) :: self
    integer :: exp_type
    
    self%species = 'cyanobacteria' ! file name prefix
    if (exp_type .eq. 1) then
      self%mclrProductionRate = zero
      self%mclrExcretionRate = zero
    else if (exp_type .eq. 2) then
      self%mclrProductionRate = 0.0001_sp
      self%mclrExcretionRate = zero
    else if (exp_type .eq. 3) then
      self%mclrProductionRate = 0.0001_sp
      self%mclrExcretionRate = 0.0001_sp
    else if (exp_type .eq. 4) then
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function colony_carbonFixation(self) ! OK
    ! Updates carbohydrate ballast state variable due to fixation (alias is "fixation")
    use MOD_PREC, only : SP ! for single or double precision
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function colony_carbonSynthesis(self) ! OK
    ! updates carbohydrate and protein state variables due to synthesis transport (alias is "synthesis")
    ! calls tempLimit()
    use MOD_PREC, only : sp ! real precision
    class(LAG_TOX), intent(inout) :: self
    real(sp), dimension(self%ndrft) :: colony_carbonSynthesis
    
    colony_carbonSynthesis(:) =  self%carbohydrate(:) * synthesisMax * self%tempLimit()

  end function colony_carbonSynthesis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function colony_carbonExcretion(self) ! OK
    ! update protein and dissolved pools due to excretion transport (alias is "excretion")
    ! temperature is tracked for all particles, so function uses algae array subset
    use MOD_PREC, only : sp ! for single or double precision
    class(LAG_TOX), intent(inout) :: self
    real(sp), dimension(self%ndrft) :: colony_carbonExcretion
    
    colony_carbonExcretion(:) = excretionFrac * self%tempFunction() * (respirationBasic*self%carbohydrate(:) + synthesisMax*self%protein(:))
    
  end function colony_carbonExcretion
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function colony_carbonRespiration(self) ! OK
    ! update carbohydrate and dissolved pools due to respiration transport (alias is "respiration")
    ! temperature is tracked for all particles, so function calls use algae array subset
    use MOD_PREC, only : sp ! for single or double precision
    class(LAG_TOX), intent(inout) :: self
    real(sp), dimension(self%ndrft) :: colony_carbonRespiration

    colony_carbonRespiration(:) = respirationBasic*self%tempFunction()*self%protein(:) + respirationActive*synthesisMax*self%tempLimit()*self%carbohydrate(:)
    
  end function colony_carbonRespiration
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function colony_temperatureLimit(self) ! OKAY
    ! Returns array of temperature limitation coefficents (0,1) for C synthesis
    use MOD_PREC, only : sp ! for single or double precision
    class(LAG_TOX), intent(in) :: self
    real(sp), dimension(self%ndrft) :: colony_temperatureLimit ! array of output coefficients for each particle
    
    colony_temperatureLimit(:) = (self%TEMP(:) / tempOpt * (  ((self%TEMP(:) - tempLethal)/(tempOpt - tempLethal))**( (tempRef - tempOpt) / tempOpt )  ))**(4.0_SP)
    
  end function colony_temperatureLimit
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function colony_temperatureFunction(self) ! OK
    ! returns array of scaling coefficents for biometric fcns
    use MOD_PREC, only : sp! for single and double precision
    
    class(LAG_TOX), intent(in) :: self
    real(sp), dimension(self%ndrft) :: colony_temperatureFunction ! array of output coefficients for each particles

    colony_temperatureFunction(:) = tempFcnAlpha * exp(tempFcnBeta * (self%temp(:) - tempOpt + tempRef))

  end function colony_temperatureFunction
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function colony_microcystinProduction(self) ! OK
    ! calculates toxin production per time step
    use MOD_PREC, only : sp ! for precision
    
    class(LAG_TOX), intent(inout) :: self
    real(sp), dimension(self%ndrft) :: colony_microcystinProduction ! array of microcystin production for colony particles
    colony_microcystinProduction(:) = self%mclrProductionRate * self%protein(:)
    
  end function colony_microcystinProduction
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function colony_microcystinExcretion(self) ! OK
    ! calculates temperature dependent toxin loss and moves mass to host element
    use MOD_PREC, only : sp
    
    class(LAG_TOX), intent(inout) :: self
    real(sp), dimension(self%ndrft) :: colony_microcystinExcretion
    colony_microcystinExcretion(:) = self%mclrExcretionRate * self%protein(:)
    
    
  end function colony_microcystinExcretion
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine colony_verticalMovement(self) ! OK
    ! update position due to buoyant movement: calls velocity(), zinterp(), zlocate(), sigma()
    use MOD_PREC, only : sp ! single precision
    use LIMS, only : KB, KBM1
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
    
    runge_kutta_integration: do ii = 1, MSTAGE
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
      
      if (ii .lt. MSTAGE) then 
        mcoef = A_RK(ii+1)
      else
        mcoef = 1.0_SP
      end if
      
      calc_array(:) = synthesis(:)/respiration(:)
      where ( ((synthesis(:)+respiration(:))*dti*mcoef) .gt. (self%carbohydrate(:) + fixation(:)*dti*mcoef) )
        synthesis(:) = (self%carbohydrate(:)/dti/mcoef + fixation(:))/(1.0_SP + 1.0_SP/calc_array(:))
        respiration(:) = (self%carbohydrate(:)/dti/mcoef + fixation(:))/(1.0_SP + calc_array(:))
      end where
      
      where ( excretion(:)*dti*mcoef .gt. (self%protein(:) + synthesis(:)*dti*mcoef) )
        excretion(:) = self%protein(:)/dti + synthesis(:)
      end where

      where ( mc_excretion(:)*dti*mcoef .gt. (self%microcystin(:) + mc_production(:)*dti*mcoef) )
        mc_excretion(:) = self%microcystin(:)/dti/mcoef + mc_production(:)
      end where

      chi_position(:, ii) = self%velocity()
      chi_carbohydrate(:,ii) = fixation(:) - respiration(:) - synthesis(:)
      chi_protein(:,ii) = synthesis(:) - excretion(:)
      chi_microcystin(:,ii) = mc_production(:) - mc_excretion(:)
      chi_dissolved(:,ii) = mc_excretion(:)
      
    end do runge_kutta_integration
    
    
     ! restore initial values
    self%carbohydrate(:) = ini_carbohydrate(:)
    self%protein(:) = ini_protein(:)
    self%microcystin(:) = ini_microcystin(:)
    self%zp(:) = ini_sigma(:)
    self%zpt(:) = ini_position(:)
    
    stage_summation: do ii = 1, MSTAGE ! add weighted stages to initial values
    
    
      where (chi_dissolved(:, ii)*dti*B_RK(ii) .gt. self%microcystin(:)) chi_dissolved(:,ii) = self%microcystin(:)/dti/B_RK(ii)
      where (-chi_carbohydrate(:, ii)*dti*B_RK(ii) .gt. self%carbohydrate(:)) chi_carbohydrate(:,ii) = self%carbohydrate(:)/dti/B_RK(ii)
      where (-chi_protein(:, ii)*dti*B_RK(ii) .gt. self%protein(:)) chi_protein(:,ii) = self%protein(:)/dti/B_RK(ii)
    
      self%zpt(:) = self%zpt(:) + chi_position(:,ii)*B_RK(ii)*dti ! update depth
      self%carbohydrate(:) = self%carbohydrate(:) + chi_carbohydrate(:, ii)*dti*B_RK(ii) ! update carbs
      self%protein(:) = self%protein(:) + chi_protein(:, ii)*dti*B_RK(ii) ! update protein
      self%microcystin(:) = self%microcystin(:) + chi_microcystin(:, ii)*dti*B_RK(ii) ! update toxin
      
      where (self%zpt(:) .gt. self%ep(:)) 
        self%zpt(:) = self%ep(:)
      else where (self%zpt(:) .lt. self%hp(:)) 
        self%zpt(:) = self%hp(:)
      end where
      
      
      ! update dissolved toxin  at particle position for each stage
      self%zp(:) = self%sigma(self%zpt(:)) ! update sigma
      self%layer(:) = self%zlocate(self%zp(:)) ! update layer
      idz(:) = float(KBM1)*abs((1.0_SP/float(KBM1) * (self%layer(:)-1)) - self%zp(:)) ! relative distance from layer above
      where (idz(:) .gt. 1.0_SP) 
        idz(:) = 1.0_SP
      elsewhere (idz(:) .lt. zero) 
        idz(:) = zero
      end where
      
      where (self%layer(:) .eq. 1) 
        domain%verticaltox(self%layer(:)) = domain%verticaltox(self%layer(:)) + 2.0_SP*B_RK(ii)*dti*chi_dissolved(:,ii)*(1.0_SP - idz(:))/domain%layerDepth ! add to sigma above
        domain%verticaltox(self%layer(:)+1) = domain%verticaltox(self%layer(:)+1) + B_RK(ii)*dti*chi_dissolved(:,ii)*idz(:)/domain%layerDepth ! add to sigma below
      elsewhere (self%layer(:) .eq. KBM1)
        domain%verticaltox(self%layer(:)) = domain%verticaltox(self%layer(:)) + B_RK(ii)*dti*chi_dissolved(:,ii)*(1.0_SP - idz(:))/domain%layerDepth ! add to sigma above
        domain%verticaltox(self%layer(:)+1) = domain%verticaltox(self%layer(:)+1) + 2.0_SP*B_RK(ii)*dti*chi_dissolved(:,ii)*idz(:)/domain%layerDepth ! add to sigma below
      elsewhere 
        domain%verticaltox(self%layer(:)) = domain%verticaltox(self%layer(:)) + B_RK(ii)*dti*chi_dissolved(:,ii)*(1.0_SP - idz(:))/domain%layerDepth ! add to sigma above
        domain%verticaltox(self%layer(:)+1) = domain%verticaltox(self%layer(:)+1) + B_RK(ii)*dti*chi_dissolved(:,ii)*idz(:)/domain%layerDepth ! add to sigma below
      end where
      
    end do stage_summation


    self%carbohydrate(:) = max(self%carbohydrate(:), zero)
    self%protein(:) = max(self%protein(:), zero)
    where (self%microcystin(:) .lt. zero) self%microcystin(:) = zero
    self%zpt(:) = min(self%zpt(:), self%ep(:)) ! collect at surface, mixing during random walk
    self%zpt(:) = max(self%zpt(:), self%hp(:)) ! collect at sediment, mixing during random walk
    self%zp(:) = self%sigma(self%zpt(:)) ! convert to sigma coordinate for finding current layer
    self%layer(:) = self%zlocate(self%zp(:)) ! update layer
    self%rho(:) = self%zinterp(domain%verticalrho(:)) ! vertical interp of density at final position
    self%temp(:) = self%zinterp(domain%verticaltemp(:)) ! vertical interp of density at final position
    

    
  end subroutine colony_verticalMovement
 
