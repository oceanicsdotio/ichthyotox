module Fish

    use lagrangian, only: agent
    use variables, only: SP, M, KB
    
    implicit none
    save

    real(SP), parameter :: &
        & PI = 3.141592653, &
        & travel_distance = 0.5787, & ! m/s = 50 km/day, test @ 20 and 2
        & epsx = sqrt((travel_distance**2.0)*0.5), &
        & epsx_sigma = 0.5*travel_distance, &
        & salinity_optimal = 30.0, & ! test @ 2.0
        & salinity_sigma_coef = 5.0, &
        & speed_table(0:4) = (/0.50, 1.00, 0.50, 0.25, 0.33/), &
        & angle_table(0:4) = (/2.00, 0.25, 0.25, 1.00, 0.50/), &
        & memory(0:1) = (/0.5, 0.96/), & ! unitless, memory coefficients
        & threshold(1:2) = (/0.005*10.0**(-6.0), 0.5/), & ! detection thresholds
        & weight(1:2) = (/0.7, 1.0/), & ! sensitivity analyis @ (/0.1, 1.0/)
        & initBodylength = 0.1, & ! meters
        & growthMax = 0.0025*12.0*0.001, & ! conversion to meters per hour from mm per 5min
        & util_cutoff = 0.01, & ! level at which default behavior is chosen
        & absorptionRate = 0.01*10.0*0.046748, & ! grams of toxin / m^2 / hour / [toxin]
        & depurationRate = 0.01, & ! sensitivity analysis @ 0.005
        & ingestionRate = 0.001*0.02, &
        & toxfrac = 0.015*10.0**(-6.0), &
        & speedimpair = 0.9 ! sensitivity analysis @ 0.5
        & h1h1 = 0.75_sp, &
        & h2h2 = 0.9_sp, &

    logical, parameter :: &
        & enforce_default = .false., &
        & no_flight = .false., &
        & ingestion_multiplier = .true.

    type, public, extends(Agent) :: FishAgent

        logical, allocatable, dimension(:), private :: &
            & impaired

        integer, allocatable, dimension(:), private :: &
            & last_rule

        real(SP), allocatable, dimension(:), private :: &
            & reverse, &
            & suitability, & ! spatial varying growth rate
            & length, &
            & effective_length, & ! impairment scalar
            & mass, &
            & microcystin, & ! body toxin
            & dissolved, & ! in situ toxin concentration
            & angle

        real(SP), allocatable, dimension(:, :), private :: &
            & event, & ! fish x agents
            & probability, & ! fish x (agents x timescales)
            & utility ! fish x (agents x timescales)

    contains
        ! Call in the order: dynamics, toxicity, movement
        procedure, public :: init => fish_initialize
        procedure, public :: writeState => fish_writeState
        procedure, public :: movement => fish_movement ! behavior selection and movement
    end type; 

contains

    subroutine fish_initialize(self, random_angle)

        class(FishAgent), intent(inout) :: self
        real(SP), intent(in) :: random_angle

        self%species = 'fish'
        call self%readPosition() ! read particles counts and allocates position variables

        allocate (&
            & self%suitability(self%ndrft), &
            & self%impaired(self%ndrft), &
            & self%microcystin(self%ndrft), &
            & self%dissolved(self%ndrft), &
            & self%angle(self%ndrft), &
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

        self%angle = 2.0*PI*random_angle
        self%length = initBodylength
        self%effective_length = 1.0
        self%mass = 2.0*10.0**(-6.0)*(1000.0*self%length)**(3.38)

        self%microcystin = 0.0_sp
        self%dissolved = 0.0_sp
        self%probability = 0.0_sp
        self%utility = 0.0_sp
        self%reverse = 0.0_sp
        self%suitability = 0.0_sp

    end subroutine

    subroutine fish_writeState(self, fid, time)
        real(SP), intent(in) :: time
        class(FishAgent), intent(in) :: self ! cyanobacteria extended type
        integer, intent(in) :: fid ! persistent file unit number
        integer :: ii
        write (fid, "(1F10.2,9000(3F20.6))") time, &
            & (self%mass(ii), 1000.0*self%microcystin(ii), 0.0_sp, ii=1, self%ndrft)
    end subroutine

    subroutine kinesis(self, noise, deltat, salinity, temperature, density, HIN, EIN)
        class(FishAgent), intent(inout) :: self
        real(SP), intent(in) :: noise(self%ndrft, 2)
        real(SP), intent(in) :: deltat ! time step, usually DTI
        real(SP), dimension(0:M, KB), intent(in) :: salinity, temperature, density ! grid based field for kinesis (usually salinity)
        real(SP), dimension(0:M), intent(in) :: HIN, EIN ! grid based field for kinesis (usually salinity)

        real(SP), dimension(self%ndrft) :: PDXT, PDYT
        logical, dimension(self%ndrft) :: inwater
        integer :: ii
        real(SP) :: pp1, p1

        do ii = 1, self%ndrft
            pp1 = (self%sal(ii) - salinity_optimal) / salinity_sigma_coef
            p1 = exp(-0.5 * (pp1 * pp1))
            ! Update U and V velocities
            self%up(ii) = self%UP(ii) * h1h1 * p1 + (noise(ii,1)*epsx_sigma + epsx) * (1.0 - h2h2 * p1) 
            self%up(ii) = self%VP(ii) * h1h1 * p1 + (noise(ii,2)*epsx_sigma + epsx) * (1.0 - h2h2 * p1)
            ! Update position
            pdxt(ii) = self%xp(ii) + self%up(ii) * deltat
            pdyt(ii) = self%yp(ii) + self%vp(ii) * deltat
        end do

        ! Evaluate Temporary Location
        inwater = .true.

        ! Update only particles still in water
        call self%find_host_element(pdxt, pdyt, inwater)
        where (inwater)
            self%xp = PDXT
            self%yp = PDYT
        end where
       
        ! interpolate bathymetry, elevation, and fields at new position
        call self%INTERP_ELH(self%xp, self%yp, HIN, EIN, 1) 
        call self%INTERP_FIELDS(self%xp, self%yp, self%zp, salinity, temperature, density, 0) 

    end subroutine

    subroutine fish_movement(self, noise, carbon_ratio)

        use simulation, only: domain
        use variables, only: dti, vxmin, vxmax, vymin, vymax

        class(FishAgent), intent(inout) :: self
        integer :: ii, jj, rule_index
        real(SP), intent(in) :: noise(self%ndrft), carbon_ratio(self%ndrft)
        real(SP) :: maxutil, speed
        real(SP), dimension(self%ndrft) :: absorption, depuration, ingestion, mass

        maxutil = 0.0_sp
        speed = 0.0_sp
        self%suitability = 0.5*(1.0 + sin(2.0*PI*(self%xp - vxmin - (vxmax/4.0))/(vxmax - vxmin)))

        where (self%microcystin/self%mass > toxfrac) ! induce impairment if toxin level above some threshold
            self%impaired = .true.
            self%effective_length = speedimpair
        elsewhere
            self%impaired = .false.
            self%effective_length = 1.0
        end where

        self%event(:, :) = 0.0_sp
        where (self%microcystin > threshold(1)) self%event(:, 1) = 1.0 ! intoxication/mortality
        where (self%suitability > threshold(2)) self%event(:, 2) = 1.0 ! current suitability

        do ii = 1, self%ndrft

            self%probability(ii, 1:2) = (1.0 - memory(0:1))*self%event(ii, 1) + memory(0:1)*self%probability(ii, 1:2)
            self%probability(ii, 3:4) = (1.0 - memory(0:1))*self%event(ii, 2) + memory(0:1)*self%probability(ii, 3:4)

            self%utility(ii, 0) = 0.01
            self%utility(ii, 1:2) = weight(1)*self%probability(ii, 1:2)
            self%utility(ii, 3:4) = weight(2)*self%probability(ii, 3:4)

            if (no_flight) then
                self%utility(ii, 1) = 0.0_sp
                self%utility(ii, 2) = 0.0_sp
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
            self%reverse(ii) = merge(1.0_sp, 0.0_sp, ((rule_index == 1) .and. (self%reverse(ii) < 0.5))) ! reverse direction for avoidance
            self%angle(ii) = self%angle(ii) + self%reverse(ii)*pi + &
                    & noise(ii)*pi*angle_table(rule_index)

            if (self%angle(ii) < -pi) then
                self%angle(ii) = self%angle(ii) + 2*PI ! normalize angles to -pi, pi]
            elseif (self%angle(ii) > pi) then
                self%angle(ii) = self%angle(ii) - 2*PI
            end if

            speed = dti*3600.0*speed_table(rule_index)*self%length(ii)*self%effective_length(ii)
            self%xp(ii) = self%xp(ii) + cos(self%angle(ii))*speed
            self%yp(ii) = self%yp(ii) + sin(self%angle(ii))*speed

        end do

        ! wrap positions
        where (self%xp > vxmax)
            self%xp = vxmin + (self%xp - vxmax)
        elsewhere(self%xp < vxmin)
            self%xp = vxmax - (vxmin - self%xp)
        end where

        where (self%yp > vymax)
            self%yp = vymin + (self%yp - vymax)
        elsewhere(self%yp < vymin)
            self%yp = vymax - (vymin - self%yp)
        end where

        ! length growth of individuals in meters per hour based on small pelagic fish
        self%length = self%length + growthMax * self%suitability * dti

        ! add consumed biomass to gut, and a portion of that to fish mass
        mass = 2.0*10.0**(-6.0)*(1000.0*self%length)**(3.38)
        ingestion = (self%mass - mass)*ingestionRate*merge(10.0, 1.0, ingestion_multiplier) * carbon_ratio
        depuration = depurationRate * self%microcystin
        absorption = self%zinterp(domain%verticaltox)/domain%meshArea/500.0*500.0* &
            & absorptionRate*self%length
        self%microcystin = self%microcystin + ingestion + (absorption - depuration)*dti
        self%mass = mass
    end subroutine
end module
