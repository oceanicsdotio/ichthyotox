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
        & speed_table(0:4) = [0.50, 1.00, 0.50, 0.25, 0.33], &
        & angle_table(0:4) = [2.00, 0.25, 0.25, 1.00, 0.50], &
        & memory(0:1) = [0.5, 0.96], & ! unitless, memory coefficients
        & threshold(1:2) = [0.005*10.0**(-6.0), 0.5], & ! detection thresholds
        & weight(1:2) = [0.7, 1.0], & ! sensitivity analyis @ (/0.1, 1.0/)
        & initBodylength = 0.1, & ! meters
        & growthMax = 0.0025*12.0*0.001, & ! conversion to meters per hour from mm per 5min
        & util_cutoff = 0.01, & ! level at which default behavior is chosen
        & absorptionRate = 0.01*10.0*0.046748, & ! grams of toxin / m^2 / hour / [toxin]
        & depurationRate = 0.01, & ! sensitivity analysis @ 0.005
        & ingestionRate = 0.01*0.02, &
        & toxfrac = 0.015*10.0**(-6.0), &
        & speedimpair = 0.9, & ! sensitivity analysis @ 0.5
        & h1h1 = 0.75_sp, & ! kinesis parameters
        & h2h2 = 0.9_sp

    integer, parameter :: &
        & INTOXICATION_CUE = 1, &
        & SUITABILITY_CUE = 2
    type, public, extends(Agent) :: FishAgent

        integer, allocatable, dimension(:), private :: &
            & last_rule

        real(SP), allocatable, dimension(:), private :: &
            & reverse, & ! flag for direction switching
            & suitability, & ! spatial varying growth rate
            & length, & ! body length proportional to swimming speed
            & effective_length, & ! impairment scalar
            & mass, & ! individual mass
            & toxin, & ! body toxin
            & dissolved, & ! in situ toxin concentration
            & angle ! orientation

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
            & self%suitability(self%count), &
            & self%toxin(self%count), &
            & self%dissolved(self%count), &
            & self%angle(self%count), &
            & self%last_rule(self%count), &
            & self%reverse(self%count), &
            & self%mass(self%count), &
            & self%length(self%count), &
            & self%utility(self%count, 0:4), &
            & self%probability(self%count, 1:4), &
            & self%event(self%count, 2), &
            & self%effective_length(self%count))

        self%event = 0
        self%last_rule = 0

        self%angle = 2.0*PI*random_angle
        self%length = initBodylength
        self%effective_length = 1.0
        self%mass = 2.0*10.0**(-6.0)*(1000.0*self%length)**(3.38)

        self%toxin = 0.0_sp
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
            & (self%mass(ii), 1000.0*self%toxin(ii), 0.0_sp, ii=1, self%count)
    end subroutine

    subroutine kinesis(self, salinity, dt, noise)
        class(FishAgent), intent(inout) :: self
        real(SP), intent(in) :: noise(self%count, 2) ! random components of U and V, allows reproducibility
        real(SP), intent(in) :: dt ! time step
        real(SP), dimension(0:M, KB), intent(in) :: salinity
        integer :: ii
        real(SP) :: suitability, signal
        do ii = 1, self%count
            suitability = (self%sal(ii) - salinity_optimal) / salinity_sigma_coef
            signal = exp(-0.5 * (suitability * suitability))
            ! Update U and V velocities
            self%up(ii) = self%UP(ii) * h1h1 * signal + (noise(ii,1)*epsx_sigma + epsx) * (1.0 - h2h2 * signal) 
            self%vp(ii) = self%VP(ii) * h1h1 * signal + (noise(ii,2)*epsx_sigma + epsx) * (1.0 - h2h2 * signal)
            ! Update position
            self%xp = self%xp(ii) + self%up(ii) * dt
            self%xp = self%yp(ii) + self%vp(ii) * dt
        end do
    end subroutine

    subroutine intoxication(self, noise, carbon_ratio)
        use simulation, only: domain
        use variables, only: dti, vxmin, vxmax, vymin, vymax

        class(FishAgent), intent(inout) :: self
        integer :: ii
        real(SP), intent(in) :: noise(self%count), carbon_ratio(self%count)
        real(SP), dimension(self%count) :: absorption, depuration, ingestion

    end subroutine

    function periodic_suitability(self) result(suitability)
        use variables, only: dti, vxmin, vxmax, vymin, vymax
        class(FishAgent), intent(inout) :: self
        real(SP), dimension(self%count) :: suitability
        suitability = 0.5*(1.0 + sin(2.0*PI*(self%xp - vxmin - (vxmax/4.0))/(vxmax - vxmin)))
    end function

    subroutine fish_movement(self, noise, carbon_ratio)

        use simulation, only: domain
        use variables, only: dti, vxmin, vxmax, vymin, vymax

        class(FishAgent), intent(inout) :: self
        integer :: ii, jj, rule_index
        real(SP), intent(in) :: noise(self%count), carbon_ratio(self%count)
        real(SP) :: maxutil, speed
        real(SP), dimension(self%count) :: absorption, depuration, ingestion, mass

        maxutil = 0.0_sp
        speed = 0.0_sp
        self%suitability = self%periodic_suitability()

        where (self%toxin/self%mass > toxfrac) ! induce impairment if toxin level above some threshold
            self%effective_length = speedimpair
        elsewhere
            self%effective_length = 1.0
        end where

        self%event(:, :) = 0.0_sp
        where (self%toxin > threshold(INTOXICATION_CUE)) self%event(:,INTOXICATION_CUE) = 1.0 ! intoxication/mortality
        where (self%suitability > threshold(SUITABILITY_CUE)) self%event(:, SUITABILITY_CUE) = 1.0 ! current suitability

        do ii = 1, self%count

            self%probability(ii, 1:2) = (1.0 - memory(0:1))*self%event(ii, 1) + memory(0:1)*self%probability(ii, 1:2)
            self%probability(ii, 3:4) = (1.0 - memory(0:1))*self%event(ii, 2) + memory(0:1)*self%probability(ii, 3:4)

            self%utility(ii, 0) = 0.01
            self%utility(ii, 1:2) = weight(1)*self%probability(ii, 1:2)
            self%utility(ii, 3:4) = weight(2)*self%probability(ii, 3:4)

            rule_index = 0 ! default behavior
            maxutil = self%utility(ii, rule_index)
            do jj = 1, 4
                if (self%utility(ii, jj) > maxutil) then
                    maxutil = self%utility(ii, jj) ! select highest util, or default behaviors
                    rule_index = jj
                end if
            end do

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
        ingestion = (self%mass - mass) * ingestionRate * carbon_ratio
        depuration = depurationRate * self%toxin
        absorption = self%zinterp(domain%verticaltox)/domain%meshArea/500.0*500.0* &
            & absorptionRate*self%length
        self%toxin = self%toxin + ingestion + (absorption - depuration)*dti
        self%mass = mass
    end subroutine
end module
