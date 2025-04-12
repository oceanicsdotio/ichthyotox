module Cyanobacteria
    ! The Cyanobacteria module handles state and calculations for the growth,
    ! buoyancy, and toxicity at the individual colony level for organisms
    ! like cyanobacteria of Anabaena and Microcystis genera. State is saved 
    ! between executions, allowing seperate invokations to do different things (init/run)

    use lagrangian, only: Agent ! module extends the lagrangian particle type
    use variables, only: zero, sp
    implicit none
    save

    real(sp), parameter :: &
            & colonyBaseRadius = 75.0*10.0**(-6.0), & ! meters
            & tempRef = 25.0, & ! reference temperature for limit fcn
            & tempOpt = 28.0, & ! optimal growth temperature
            & tempLethal = 35.0, & ! lethal temperature
            & excretionFrac = 0.1, & ! unitless
            & fixationMax = 11.4, & ! maximum carbon fixation, per hour rate
            & fixationBeta = 0.02, & ! shape factor in fixation calculation
            & respirationBasic = 0.004, & ! basic respiration rate, per hour
            & respirationActive = 0.2, & ! active respiration rate, unitless
            & densityMax = 1150.0, & ! maximum empirical density of algal cells, kg/m^3
            & densityMin = 1037.0, & ! minimum empirical density of algal cells, kg/m^3
            & vesicleDensity = 150.0, & ! density of gas filled vesicles, kg/m^3
            & cellFrac = 0.25, & ! fraction of colony volume composed of cell material, unitless
            & carbonRatioMax = 4.0, & ! maximum empirical ratio of carbon reservoirs in algal cells, unitless
            & vesicleFrac = 0.08, & ! fraction of cell volume occupied by vesicles, unitless
            & irradOpt = 250.0, & ! optimal irradiance W/M^-2
            & synthesisMax = 0.05, & ! per hour rate
            & tempFcnAlpha = 0.286, & ! shape factor, unitless
            & tempFcnBeta = 0.05, &   ! shape coefficient, unitless
            & cellDensityCoefficient = 0.7, & ! shape coefficient, unitless
            & lightExtinctionBiomass = 14.0, & ! light extinction due to overlying biomass
            & lightAttenuationWater = 0.15, & ! light extinction coefficient due to coastal waters
            & shading_upscale = 1.0

    type, public, extends(Agent) :: CyanobacteriaAgent
        ! algal state inherited from lagrangian particle class
        real(sp) :: &
                & mclrProductionRate = zero, &
                & mclrExcretionRate = zero

        real(sp), allocatable, dimension(:), private :: &
                & radius, & ! radius for vertical movement of spherical colonies, meters
                & irradiance, & ! irradiance at particle position, watts per sq meter
                & biomass, & ! overlying biomass used for light attenuation calculations, gram
                & delta_rho ! difference in density between water and colony

        real(sp), allocatable, dimension(:), public :: &
                & carbohydrate, & ! mass of carbohydration ballast, grams
                & protein, & ! mass of buoyant cellular material, grams
                & microcystin ! mass of cellular toxins, grams

    contains
        ! user interface subroutines
        procedure, public :: init => colony_initialize ! read position and count, allocate state variables
        procedure, public :: writeState => colony_writeState ! write internal state variables
        procedure, public :: movement => colony_verticalMovement ! subroutine that updates position of colony particles
        procedure, public :: random => colony_random_walk ! vertical random walk for cyanobacteria

        ! private utility procedures
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

    end type CyanobacteriaAgent; 

contains

    subroutine colony_initialize(self, experimentType)
        ! read position and state, and allocate internal variables
        use variables, only: iovar, folderprefix
        class(CyanobacteriaAgent), intent(inout) :: self
        integer :: experimentType, ii, indexMatch
        character(len=50) :: filename
        logical :: fexist

        self%species = 'cyanobacteria' ! file name prefix
        if (experimentType == 1) then
            self%mclrProductionRate = zero
            self%mclrExcretionRate = zero
        else if (experimentType == 2) then
            self%mclrProductionRate = 0.0001
            self%mclrExcretionRate = zero
        else if (experimentType == 3) then
            self%mclrProductionRate = 0.0001
            self%mclrExcretionRate = 0.0001
        else if (experimentType == 4) then
            self%mclrProductionRate = zero
            self%mclrExcretionRate = 0.0001
        else
            write (*, *) "Problem with experiment type, stopping..."; stop
        end if

        call self%readPosition() ! finds ndrft and allocates position arrays

        allocate (self%radius(self%ndrft), &
                & self%irradiance(self%ndrft), &
                & self%biomass(self%ndrft), &
                & self%carbohydrate(self%ndrft), &
                & self%protein(self%ndrft), &
                & self%microcystin(self%ndrft), &
                & self%delta_rho(self%ndrft)); 
        self%radius = zero ! read from file
        self%irradiance = zero ! calculated at runtime
        self%biomass = zero ! calculated at runtime
        self%carbohydrate = zero ! read from file
        self%protein = zero ! read from file
        self%microcystin = zero ! read from file
        self%delta_rho = zero ! calculated at runtime

        ! Check for state variable file, allocate storage and read values
        write (filename, "(A)") "./"//trim(folderprefix)//"/"//trim(self%species)//'_var.dat' ! variables read from "cyanobacteria_var.dat"
        inquire (file=trim(filename), exist=fexist) ! check for file
        if (.not. fexist) then
            write (*, *) 'State variable file: ', filename, ' does not exist.'
            stop
        end if

        open (unit=iovar, file=filename, form='formatted')
        read (iovar, "(I6)") indexMatch
        if (indexMatch /= self%ndrft) then
            write (*, *) 'Dimensions of position and state variable files are not equal.'
            stop
        end if

        do ii = 1, self%ndrft
            read (iovar, "(4F20.6)") self%radius(ii), self%carbohydrate(ii), self%protein(ii), self%microcystin(ii)
        end do
        close (iovar)

    end subroutine

    subroutine colony_writeState(self, fid, time)
        use simulation, only: domain ! domain structure for elapsed time
        class(CyanobacteriaAgent), intent(in) :: self ! cyanobacteria extended type
        integer, intent(in) :: fid ! persistent file unit number
        real(sp), intent(in) :: time ! elapsed time in seconds
        integer :: ii

        write(fid, "(1F10.2,9000(I6,3F20.3))") time, & 
            & (self%itag(ii), self%carbohydrate(ii), self%protein(ii), &
            & self%microcystin(ii), ii=1,self%ndrft)

    end subroutine

    recursive subroutine colony_QsortC(absdepth, order) ! OK
        ! Recursive Fortran 95 quicksort routine sorts real numbers into ascending numerical order
        ! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
        ! Based on algorithm from Cormen et al., Introduction to Algorithms, 1997 printing
        ! Made F conformant by Walt Brainerd http://www.fortran.com/qsort_c.f95
        use variables, only: sp
        real(sp), intent(inout), dimension(:) :: absdepth
        integer, intent(inout), dimension(:) :: order
        integer :: iq

        if (size(absdepth) > 1) then
            call colony_Partition(absdepth, order, iq)
            call colony_QsortC(absdepth(:iq - 1), order(:iq - 1))
            call colony_QsortC(absdepth(iq:), order(iq:))
        end if
    end subroutine

    subroutine colony_Partition(A, B, marker) ! OK
        use variables, only: sp

        real(sp), intent(inout), dimension(:) :: A
        integer, intent(inout), dimension(:) :: B
        integer, intent(out) :: marker
        integer :: ii, jj, itemp ! iterators and temporary index
        real(sp) :: temp ! temporary sorting value
        real(sp) :: x ! pivot point
        x = A(1)
        ii = 0
        jj = size(A) + 1

        do
            jj = jj - 1
            do
                if (A(jj) <= x) exit
                jj = jj - 1
            end do
            ii = ii + 1
            do
                if (A(ii) >= x) exit
                ii = ii + 1
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
                marker = ii + 1
                return
            else
                marker = ii
                return
            end if
        end do
    end subroutine

    function colony_carbonFixation(self)
        ! Updates carbohydrate ballast state variable due to fixation (alias is "fixation")
        use variables, only: SP ! for single or double precision
        use simulation, only: domain

        class(CyanobacteriaAgent), intent(inout) :: self
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
        avg_self_shade(:) = (exp(-lightExtinctionBiomass*self%biomass(:)) - 1.0)/(-lightExtinctionBiomass*self%biomass(:)) ! average area under self shading curve from zero to self biomass
        self%irradiance(indices(1)) = domain%globalIrradiance ! first particle unshaded

        do ii = 2, self%ndrft
            ! attenuate by next biomass
            self%irradiance(indices(ii)) = self%irradiance(indices(ii-1)) * &
                & exp(-lightExtinctionBiomass*shading_upscale*self%biomass(indices(ii-1))) 
        end do

        self%irradiance(:) = self%irradiance(:)*avg_self_shade(:) ! multiply by particle self shading component
        self%irradiance(:) = self%irradiance(:)*exp(self%zpt(:)*lightAttenuationWater) ! find irradiance at particle position after attenuation
        irradRatio(:) = self%irradiance(:)/irradOpt ! substitution function
        fixationCoef(:) = (2.0 + fixationBeta)*irradRatio(:)/(irradRatio(:)**2.0 + &
            & fixationBeta*irradRatio(:) + 1.0) ! scaling coefficient for transfer
        ! actual mass transfer
        colony_carbonFixation(:) = fixationMax*fixationCoef(:)*self%protein(:)*(1.0 - vesicleFrac) * &
            & (carbonRatioMax - self%carbohydrate(:)/self%protein(:))/carbonRatioMax 

    end function

    function colony_carbonSynthesis(self)
        ! updates carbohydrate and protein state variables due to synthesis transport (alias is "synthesis")
        ! calls tempLimit()
        use variables, only: sp ! real precision
        class(CyanobacteriaAgent), intent(inout) :: self
        real(sp), dimension(self%ndrft) :: colony_carbonSynthesis

        colony_carbonSynthesis(:) = self%carbohydrate(:)*synthesisMax*self%tempLimit()

    end function

    function colony_carbonExcretion(self)
        ! update protein and dissolved pools due to excretion transport (alias is "excretion")
        ! temperature is tracked for all particles, so function uses algae array subset
        use variables, only: sp ! for single or double precision
        class(CyanobacteriaAgent), intent(inout) :: self
        real(sp), dimension(self%ndrft) :: colony_carbonExcretion

colony_carbonExcretion(:) = excretionFrac*self%tempFunction()*(respirationBasic*self%carbohydrate(:) + synthesisMax*self%protein(:))

    end function

    function colony_carbonRespiration(self)
        ! update carbohydrate and dissolved pools due to respiration transport (alias is "respiration")
        ! temperature is tracked for all particles, so function calls use algae array subset
        class(CyanobacteriaAgent), intent(inout) :: self
        real(sp), dimension(self%ndrft) :: colony_carbonRespiration

    colony_carbonRespiration(:) = respirationBasic*self%tempFunction()*self%protein(:) + &
            & respirationActive*synthesisMax*self%tempLimit()*self%carbohydrate(:)

    end function

    function colony_temperatureLimit(self) result(limit)
        ! Returns array of temperature limitation coefficents (0,1) for C synthesis
        class(CyanobacteriaAgent), intent(in) :: self
        real(sp), dimension(self%ndrft) :: limit ! array of output coefficients for each particle

        limit(:) = (self%TEMP(:) / tempOpt * (  ((self%TEMP(:) - tempLethal)/(tempOpt - tempLethal))** &
            & ( (tempRef - tempOpt) / tempOpt )  ))**(4.0)

    end function

    function colony_temperatureFunction(self)
        ! returns array of scaling coefficents for biometric fcns
        use variables, only: sp! for single and double precision

        class(CyanobacteriaAgent), intent(in) :: self
        real(sp), dimension(self%ndrft) :: colony_temperatureFunction ! array of output coefficients for each particles

        colony_temperatureFunction(:) = tempFcnAlpha*exp(tempFcnBeta*(self%temp(:) - tempOpt + tempRef))

    end function

    function colony_microcystinProduction(self)
        ! calculates toxin production per time step
        use variables, only: sp ! for precision

        class(CyanobacteriaAgent), intent(inout) :: self
        real(sp), dimension(self%ndrft) :: colony_microcystinProduction ! array of microcystin production for colony particles
        colony_microcystinProduction(:) = self%mclrProductionRate*self%protein(:)

    end function

    function colony_microcystinExcretion(self)
        ! calculates temperature dependent toxin loss and moves mass to host element
        use variables, only: sp

        class(CyanobacteriaAgent), intent(inout) :: self
        real(sp), dimension(self%ndrft) :: colony_microcystinExcretion
        colony_microcystinExcretion(:) = self%mclrExcretionRate*self%protein(:)

    end function

    subroutine colony_verticalMovement(self)
        ! update position due to buoyant movement: calls velocity(), zinterp(), zlocate(), sigma()
        use simulation, only: domain
        use variables, only: dti, zero, KBM1 ! integration step
        use variables, only: A_RK, B_RK, MSTAGE, strict_integration ! runge-kutta integration parameters

        class(CyanobacteriaAgent), intent(inout) :: self
        integer :: ii
        real(sp) :: mcoef
        real(sp), dimension(self%ndrft) :: idz, ini_position, ini_sigma, ini_microcystin, ini_carbohydrate, ini_protein ! initial valuess
        real(sp), dimension(self%ndrft, 0:MSTAGE) :: chi_dissolved, chi_position, chi_carbohydrate, chi_protein, chi_microcystin ! stage fcn evaluations
        real(sp), dimension(self%ndrft) :: synthesis, excretion, respiration, fixation, mc_excretion, mc_production, calc_array
        chi_position = ZERO; chi_carbohydrate = ZERO; chi_protein = ZERO; chi_microcystin = abs(ZERO)

        ! save initial values
        ini_carbohydrate(:) = self%carbohydrate(:)
        ini_protein(:) = self%protein(:)
        ini_microcystin(:) = self%microcystin(:)
        ini_sigma(:) = self%zp(:)
        ini_position(:) = self%zpt(:)

        do ii = 1, MSTAGE
            ! get new stage values
            self%zpt(:) = ini_position(:) + A_RK(ii)*dti*chi_position(:, ii - 1)      ! new stage position
            self%carbohydrate(:) = ini_carbohydrate(:) + A_RK(ii)*dti*chi_carbohydrate(:, ii - 1)  ! new carb stage value
            self%protein(:) = ini_protein(:) + A_RK(ii)*dti*chi_protein(:, ii - 1)       ! new protein stage value
            self%microcystin(:) = ini_microcystin(:) + A_RK(ii)*dti*chi_microcystin(:, ii - 1)   ! new toxin stage value

            ! enforce strict limits on state variables
            self%carbohydrate(:) = max(self%carbohydrate(:), zero) ! keep carbs non-negative
            self%protein(:) = max(self%protein(:), zero) ! keep protein non-negative
            self%microcystin(:) = max(self%microcystin(:), zero) ! keep toxin non-negative
            self%zpt(:) = min(self%zpt(:), self%ep(:)) ! stop at surface
            self%zpt(:) = max(self%zpt(:), self%hp(:)) ! stop at sediment

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

            mcoef = merge(A_RK(ii + 1), 1.0, ii < mstage)

            calc_array(:) = synthesis(:)/respiration(:)
            where (((synthesis(:) + respiration(:))*dti*mcoef) > (self%carbohydrate(:) + fixation(:)*dti*mcoef))
                synthesis(:) = (self%carbohydrate(:)/dti/mcoef + fixation(:))/(1.0 + 1.0/calc_array(:))
                respiration(:) = (self%carbohydrate(:)/dti/mcoef + fixation(:))/(1.0 + calc_array(:))
            end where

            where (excretion(:)*dti*mcoef > (self%protein(:) + synthesis(:)*dti*mcoef))
                excretion(:) = self%protein(:)/dti + synthesis(:)
            end where

            where (mc_excretion(:)*dti*mcoef > (self%microcystin(:) + mc_production(:)*dti*mcoef))
                mc_excretion(:) = self%microcystin(:)/dti/mcoef + mc_production(:)
            end where

            chi_position(:, ii) = self%velocity()
            chi_carbohydrate(:, ii) = fixation(:) - respiration(:) - synthesis(:)
            chi_protein(:, ii) = synthesis(:) - excretion(:)
            chi_microcystin(:, ii) = mc_production(:) - mc_excretion(:)
            chi_dissolved(:, ii) = mc_excretion(:)

        end do

        ! restore initial values
        self%carbohydrate(:) = ini_carbohydrate(:)
        self%protein(:) = ini_protein(:)
        self%microcystin(:) = ini_microcystin(:)
        self%zp(:) = ini_sigma(:)
        self%zpt(:) = ini_position(:)

        do ii = 1, MSTAGE ! add weighted stages to initial values

            where (chi_dissolved(:, ii)*dti*B_RK(ii) > self%microcystin(:)) chi_dissolved(:, ii) = self%microcystin(:)/dti/B_RK(ii)
    where (-chi_carbohydrate(:, ii)*dti*B_RK(ii) > self%carbohydrate(:)) chi_carbohydrate(:, ii) = self%carbohydrate(:)/dti/B_RK(ii)
            where (-chi_protein(:, ii)*dti*B_RK(ii) > self%protein(:)) chi_protein(:, ii) = self%protein(:)/dti/B_RK(ii)

            self%zpt(:) = self%zpt(:) + chi_position(:, ii)*B_RK(ii)*dti ! update depth
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
            idz(:) = float(KBM1)*abs((1.0/float(KBM1)*(self%layer(:) - 1)) - self%zp(:)) ! relative distance from layer above
            where (idz(:) > 1.0)
                idz(:) = 1.0
            elsewhere(idz(:) < zero)
                idz(:) = zero
            end where

            where (self%layer(:) == 1)
                domain%verticaltox(self%layer(:)) = domain%verticaltox(self%layer(:)) + &
                    & 2.0*B_RK(ii)*dti*chi_dissolved(:,ii)*(1.0 - idz(:))/domain%layerDepth ! add to sigma above
                domain%verticaltox(self%layer(:)+1) = domain%verticaltox(self%layer(:)+1) + &
                    & B_RK(ii)*dti*chi_dissolved(:,ii)*idz(:)/domain%layerDepth ! add to sigma below
            elsewhere(self%layer(:) == KBM1)
                domain%verticaltox(self%layer(:)) = domain%verticaltox(self%layer(:)) + &
                    & B_RK(ii)*dti*chi_dissolved(:,ii)*(1.0 - idz(:))/domain%layerDepth ! add to sigma above
                domain%verticaltox(self%layer(:)+1) = domain%verticaltox(self%layer(:)+1) + &
                    & 2.0*B_RK(ii)*dti*chi_dissolved(:,ii)*idz(:)/domain%layerDepth ! add to sigma below
            elsewhere
                domain%verticaltox(self%layer(:)) = domain%verticaltox(self%layer(:)) + &
                    & B_RK(ii)*dti*chi_dissolved(:,ii)*(1.0 - idz(:))/domain%layerDepth ! add to sigma above
                domain%verticaltox(self%layer(:)+1) = domain%verticaltox(self%layer(:)+1) + &
                    & B_RK(ii)*dti*chi_dissolved(:,ii)*idz(:)/domain%layerDepth ! add to sigma below
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

    subroutine colony_random_walk(self, noise)
        ! vertical and horizontal random walk
        use variables, only: dti, dtrw
        class(CyanobacteriaAgent), intent(inout) :: self

        real(sp), dimension(self%ndrft) :: noise
        real(sp), dimension(self%ndrft) :: wdiff, kzp, dkzp ! diffusivity and first derivative at particle positions
        integer :: substeps
        real(sp), parameter :: variance = 1.0, AC = 1.0/6.0, BIG = 1.0E30

        ! vertical random walk
        do substeps = 1, int(dti/dtrw)
           
            kzp(:) = 60.0*60.0*10.0**(-4.0)*0.1
            dkzp(:) = 0.0
            wdiff(:) = dkzp(:)*dtrw + noise * &
                & sqrt((2.0*kzp(:) + dkzp(:)**2.0)*dtrw/variance) ! Ross and Sharples 2004 eqn 1
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
        use variables, only: grav ! grav is positive, m/s2
        class(CyanobacteriaAgent), intent(inout) :: self
        real(sp), dimension(self%ndrft) :: colony_stokesVelocity ! output array of particle vertical velocities

        self%delta_rho(:) = self%rho(:) - self%density() ! water density array - colony density fcn
        colony_stokesVelocity(:) = &
                & (60.0**2)*(2.0/9.0)*grav*self%radius(:)**2.0* &
                & self%delta_rho(:)*self%viscosity()**(-1.0)

    end function

    elemental function algae_density(carbohydrate, protein, water_density)
        ! colony density including contibutions of mucus and gas vacuoles
        real(sp) :: algae_density, carbohydrate, protein, water_density, cell_density
        
        ! array of density without mucus and vacuoles
        cell_density = densityMin + (densityMax - densityMin) * &
            & (1.0 - exp(-cellDensityCoefficient * carbohydrate/protein))
        algae_density = (1.0 - cellFrac) * (water_density + 0.7) + &
            & cellFrac * ((1.0 - vesicleFrac) * cell_density + vesicleFrac*vesicleDensity)
    end function

    function colony_algaeDensity(self)
        ! colony density including contibutions of mucus and gas vacuoles
        class(CyanobacteriaAgent), intent(in) :: self
        real(sp), dimension(self%ndrft) :: colony_algaeDensity ! output array of overall colony density

        colony_algaeDensity(:) = algae_density(self%carbohydrate, self%protein, self%rho)
    end function

    function colony_dynamicViscosity(self)
        ! returns array of dynamic viscosity values at particle locations
        class(CyanobacteriaAgent), intent(in) :: self
        logical :: simple = .true.
        real(sp), dimension(self%ndrft) :: A, B, visc_pure, colony_dynamicViscosity

        if (simple) then
            colony_dynamicViscosity(:) = 10.0**(-3.0)*10.0**(-1.65 + 262.0/(self%temp(:) + 169.0))
        else ! Sharqway et al 2010
            A(:) = 1.541 + 19.998*10.0**(-2.0)*self%temp - 9.52*10.0**(-5.0)*self%temp**(2.0)
            B(:) = 7.974 - 7.561*10.0**(-2.0) + 4.724*10.0**(-4.0)*self%temp**(2.0)
            visc_pure(:) = 4.2844*10.0**(-5.0) + (0.157*(self%temp + 64.993)**(2.0) - 91.296)**(-1.0)
            colony_dynamicViscosity(:) = visc_pure*(1.0 + A*self%sal + B*self%sal**(2.0))
        end if
    end function

end module