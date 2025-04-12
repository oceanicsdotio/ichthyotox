module Lagrangian
  ! type and variables for lagrangian particle system
  use variables, only : sp, ZERO
  implicit none
  save

  private

  type, public :: Agent

    character(len = 20) :: species ! string for identification
    logical :: fixed_depth = .false. ! fixed particle depth option in cartesian
    integer :: ndrft  ! total particles

    logical, allocatable, dimension(:) :: &
        & indomain ! Particle is in the domain

    integer, allocatable, dimension(:) :: &
        & ITAG, &     ! Label for the particle
        & host, &     ! Element containing particle
        & layer, &    ! sigma layer
        & found, &    ! Host element is found
        & SBOUND      ! Host element has a solid boundary node

    real(SP), allocatable, dimension(:) :: &
        & XP, YP, ZP, &       ! position of particle
        & XPT, YPT, ZPT, &    ! absolute position of particle
        & HP, &               ! Bathymetry at particle position
        & EP, &               ! Free surface height at particle
        & UP, VP, WP, &       ! velocity of particle
        & TEMP, &             ! temperature at particle position
        & SAL, &              ! salinity at particle position
        & RHO                 ! density at particle position

  contains

    procedure, public :: lag_alloc => lag_alloc
    procedure, public :: stats => lag_printStatistics ! print particle information if desired
    procedure, public :: readPosition => readLocations ! read position and ndrft from file
    procedure, public :: writePosition => lag_writePosition ! write position to file
    procedure, public :: sigma => lag_cart2sig ! convert cartesian to sigma
    procedure, public :: cartesian => lag_sig2cart ! convert sigma to cartesian
    procedure, public :: zlocate => lag_getlayers ! update the sigma layer that particles are currently below
    procedure, public :: zinterp => lag_zinterp
    procedure, public :: find_host_element => find_host_element
    procedure, public :: traject => traject
    procedure, public :: interp_elh => interp_elh
    procedure, public :: interp_fields => interp_fields
    procedure, public :: interp_v => interp_v
    procedure, public :: interp_kh => interp_kh

  end type

contains

  subroutine find_host_element(lag, x, y, inwater)
    !  Find host elements of particles by searching progressively further elements
    use variables, only: N, NV, ISONB, VX, NV, VY, XC, YC, NTVE, NBVE
  

    class(Agent), intent(inout) :: lag
    real(sp), dimension(lag%ndrft), intent(in) :: x, y
    logical, dimension(lag%ndrft), intent(out) :: INWATER

    real(sp), dimension(1:N, 1) :: distance
    real(sp) :: previous
    integer :: ii, jj, kk, ind
    integer, dimension(2) :: nearby

    do ii = 1, lag%ndrft

      if (.not. lag%indomain(ii)) cycle
      if (isintriangle(VX(NV(lag%HOST(ii), 1:3)), VY(NV(lag%HOST(ii), 1:3)), x(ii), y(ii))) then
        lag%FOUND(ii) = 1
        cycle
      end if

      do jj = 1, 3
        do kk = 1, NTVE(NV(lag%HOST(ii), jj))
          ind = NBVE(NV(lag%HOST(ii), jj), kk)
          if (isintriangle(VX(NV(ind, 1:3)), VY(NV(ind, 1:3)), x(ii), y(ii))) then
            lag%FOUND(ii) = 1
            lag%HOST(ii) = ind
            lag%SBOUND(ii) = merge(1, 0, &
                    & (ISONB(NV(ind, 1)) == 1) .or. (ISONB(NV(ind, 2)) == 1) .or. (ISONB(NV(ind, 3)) == 1))
            exit
          end if
        end do
      end do

      if (lag%found(ii) == 1) cycle

      distance(1:n, 1) = sqrt((xc(1:n) - x(ii))**2 + (yc(1:n) - y(ii))**2)
      previous = zero

      do jj = 1, 16
        nearby(:) = minloc(distance, distance > previous)
        if (nearby(1) == 0) exit
        previous = distance(nearby(1), 1)
        if (ISINTRIANGLE(VX(NV(nearby(1), 1:3)), VY(NV(nearby(1), 1:3)) , x(ii), y(ii))) then
          lag%found(ii) = 1
          lag%host(ii) = nearby(1)
          lag%sbound(ii) = merge(1, 0, &
                  & (ISONB(NV(lag%host(ii), 1)) == 1) .or. &
                  & (ISONB(NV(lag%host(ii), 2)) == 1) .or. &
                  & (ISONB(NV(lag%host(ii), 3)) == 1))
          exit
        end if
      end do
    end do
  end subroutine


  logical function isintriangle(XT, YT, X0, Y0)
    ! Determine if point is in triangle defined by nodes (XT(3),YT(3))
    real(sp), intent(in) :: X0, Y0
    real(sp), intent(in) :: XT(3), YT(3)
    real(sp) :: F1, F2, F3

    F1 = (Y0-YT(1))*(XT(2)-XT(1)) - (X0-XT(1))*(YT(2)-YT(1))
    F2 = (Y0-YT(3))*(XT(1)-XT(3)) - (X0-XT(3))*(YT(1)-YT(3))
    F3 = (Y0-YT(2))*(XT(3)-XT(2)) - (X0-XT(2))*(YT(3)-YT(2))

    isintriangle = ((F1*F3 >= 0.0_sp) .and. (F3*F2 >= 0.0_sp))

  end function


  subroutine lag_alloc(self)

    class(Agent), intent(inout) :: self

    allocate( self%XP( self%ndrft ), &
            & self%YP( self%ndrft ), &
            & self%ZP( self%ndrft ), &
            & self%HP( self%ndrft ), &
            & self%EP( self%ndrft ), &
            & self%FOUND( self%ndrft ), &
            & self%HOST( self%ndrft ), &
            & self%LAYER( self%ndrft ), &
            & self%SBOUND( self%ndrft ), &
            & self%INDOMAIN( self%ndrft ), &
            & self%TEMP( self%ndrft ), &
            & self%SAL( self%ndrft ), &
            & self%RHO( self%ndrft ), &
            & self%UP( self%ndrft ), &
            & self%VP( self%ndrft ), &
            & self%WP( self%ndrft ))

    self%XP(:) = zero
    self%YP(:) = zero
    self%ZP(:) = zero
    self%HP(:) = zero
    self%EP(:) = zero
    self%FOUND(:) = 0
    self%HOST(:) = 1
    self%LAYER(:) = 1
    self%SBOUND(:) = 0
    self%INDOMAIN(:) = .true.
    self%TEMP(:) = zero
    self%SAL(:) = zero
    self%RHO(:) = zero
    self%UP(:) = zero
    self%VP(:) = zero
    self%WP(:) = zero

  end subroutine


  subroutine lag_printStatistics(self)

    implicit none
    class(Agent), intent(inout) :: self

    ! Report status of particles
    write(*, *)
    write(*, *) '    Particle Class'
    write(*, *) '        Species       : ', self%species
    write(*, *) '        Quantity      : ', self%ndrft
    write(*, *) '        Out of bounds : ', self%ndrft - count(self%INDOMAIN)
    write(*, *) '        Stopped       : ', self%ndrft - sum(self%FOUND)
    write(*, *)

  end subroutine


  subroutine readLocations(self)
    use variables, only : folderprefix

    class(Agent), intent(inout) :: self

    integer :: fid = 1, ii
    character(len = 50) :: filename
    logical :: fexist

    write(filename, "(A)") "./"//trim(folderprefix)//"/"//trim(self%species)//"_ini.dat"
    inquire(file=trim(filename), exist=fexist)
    if (.not. fexist) then
      print *, 'Initial position file ', filename, ' does not exist.'
      stop
    end if

    open(unit=fid, file=filename, form='formatted')
    read(fid, "(I6)") self%ndrft ! read number of particles
    allocate( self%itag(self%ndrft), &
            & self%xpt(self%ndrft), &
            & self%ypt(self%ndrft), &
            & self%zpt(self%ndrft))

    do ii = 1, self%ndrft
      read(fid, "(I6, 3F20.6)") self%itag(ii), self%XPT(ii), self%YPT(ii), self%ZPT(ii) ! read identifier and position for each particle
    end do

    close(fid)

  end subroutine


  subroutine lag_writePosition(self, fid, time)
    ! write time and particle id/position to already open file
    use simulation, only: domain ! domain structure for current time only
    class(Agent), intent(inout) :: self ! lagrangian particle structure
    integer, intent(in) :: fid ! unit number of open output file
    real(sp), intent(in) :: time ! time to write
    integer :: ii
    write(fid, "(1F10.2,9000(I6,3F20.3))") time, (self%ITAG(ii), self%XPT(ii), self%YPT(ii), self%ZPT(ii), ii=1,self%ndrft)
  end subroutine

  function lag_cart2sig(self, cartesian)
    ! Calculate sigma vertical position from cartesian
    class(Agent), intent(in) :: self
    real(sp), dimension(self%ndrft), intent(in) :: cartesian
    real(sp), dimension(self%ndrft) :: lag_cart2sig
    lag_cart2sig(:) = -1.0_SP * abs(cartesian(:)-self%EP(:)) / abs(self%EP(:)-self%HP(:))
  end function

  function lag_sig2cart(self, sigma)
    ! calculate cartesian vertical from sigma coordinate
    class(Agent), intent(in) :: self
    real(sp), dimension(self%ndrft), intent(in) :: sigma
    real(sp), dimension(self%ndrft) :: lag_sig2cart
    lag_sig2cart(:) = sigma(:)*(self%EP(:) - self%HP(:)) + self%EP(:)
  end function

  function lag_getlayers(self, sigma) ! OK
    ! update current sigma layer of particles between moves
    use variables, only : KBM1
    class(Agent), intent(in) :: self
    real(sp), dimension(self%ndrft), intent(in) :: sigma
    integer, dimension(self%ndrft) :: lag_getlayers, layers
    layers(:) = ceiling(-KBM1*sigma(:))
    where (layers(:) < 1) layers(:) = 1
    lag_getlayers(:) = layers(:)
  end function

  function lag_zinterp(self, verticalvar)
    use variables, only : KB ! sigma layers
    use simulation, only : domain ! domain structure
    class(Agent), intent(inout) :: self ! lagrangian particle swarm object
    real(sp), dimension(0:KB+1), intent(in) :: verticalvar
    real(sp), dimension(self%ndrft) :: idz, lag_zinterp

    idz(:) = (domain%layerSigma*(self%layer(:)-1) - self%zp(:))/domain%layerSigma ! relative distance from layer above
    lag_zinterp(:) = verticalvar(self%layer(:))*(1.0_SP - idz(:)) + verticalvar(self%layer(:)+1)*idz(:)

  end function


  subroutine TRAJECT(self, dt, U1, U2, V1, V2, W1, W2, HL, EL1, EL2, salinity, temperature, density)

    ! integrate particle position from x0 to xn using velocity fields at time t0 and time tn
    use variables, only : N, M, KB, MSTAGE, A_RK, B_RK, C_RK ! runge-kutta parameters

    class(Agent), intent(inout) :: self ! lagrangian particle object
    real(SP), intent(in) :: dt ! time step, usually DTI
    real(SP), dimension(0:N, KB), intent(in) :: U1, U2, V1, V2, W1, W2 ! velocity fields at start and end of time step
    real(SP), dimension(0:M), intent(in) :: HL, EL1, EL2 ! bathymetry and free surface height
    real(SP), dimension(0:M, KB), intent(in) :: temperature, salinity, density

    real(SP), dimension(self%ndrft) :: PDXT, PDYT, PDZT, PDX, PDY, PDZ ! RK stage positions
    logical, dimension(self%ndrft) :: inwater
    integer :: stage
    real(SP), dimension(0:N, KB, 0:2) :: velocity ! ERK stage velocity field
    real(SP), dimension(0:M) :: ELL ! ERK stage freesurface height
    real(SP), dimension(self%ndrft, 0:MSTAGE, 3) :: CHI ! ERK stage function evaluation for velocities
    real(SP), parameter :: EPS  = 10.0 ** (-5.0) ! depth of dry element

    CHI = zero ! Initialize Stage Functional Evaluations

    PDXT(:) = self%xp(:) ! Assign position at previous time to current position
    PDYT(:) = self%yp(:)
    PDZT(:) = self%zp(:)

    do stage = 1, MSTAGE ! Runge-Kutta integration stages
      ! New stage position updated from time zero position
      PDX(:) = self%yp(:) + A_RK(stage) * dt * CHI(:, stage - 1, 1)
      PDY(:) = self%yp(:) + A_RK(stage) * dt * CHI(:, stage - 1, 2)
      PDZ(:) = self%zp(:) + A_RK(stage) * dt * CHI(:, stage - 1, 3)
      PDZ(:) = max(PDZ(:), -(2.0 + PDZ(:))) ! reflect sigma depth off bottom
      PDZ(:) = min(PDZ(:), zero) ! keep sigma depth below free surface

      ! Calculate velocity field for stage using c_rk coefficients
      velocity(:, :, 0) = (1.0 - C_RK(stage)) * U1 + C_RK(stage) * U2
      velocity(:, :, 1) = (1.0 - C_RK(stage)) * V1 + C_RK(stage) * V2
      velocity(:, :, 2) = (1.0 - C_RK(stage)) * W1 + C_RK(stage) * W2
      ELL = (1.0 - C_RK(stage)) * EL1 + C_RK(stage) * EL2

      call self%INTERP_V(PDX, PDY, PDZ, velocity)
      ! interpolate particle velocity, automatically updates host elements
      call self%INTERP_ELH(PDX, PDY, HL, ELL, 0) ! interpolate elevation and bathymetry at stage particle position, zero denotes not to search for host elements

      CHI(:, stage, 1) = self%UP ! Update CHI values for next time step
      CHI(:, stage, 2) = self%VP

      where ((self%EP - self%HP) < EPS) 
        CHI(:, stage, 3) = zero ! Limit vertical motion in very shallow water
      elsewhere
        CHI(:, stage, 3) = self%WP / (self%HP - self%EP)    ! delta_sigma/deltaT = ww/D
      end where
    end do

    do stage = 1, MSTAGE
      where (self%indomain)
        PDXT(:) = PDXT(:) + dt * CHI(:, stage, 1) * B_RK(stage) ! Update current position if particle is in domain
        PDYT(:) = PDYT(:) + dt * CHI(:, stage, 2) * B_RK(stage)
        PDZT(:) = PDZT(:) + dt * CHI(:, stage, 3) * B_RK(stage)
      end where
    end do

    self%FOUND = 0
    inwater(:) = .true.

    ! Perform robust progressive-topology search
    call self%find_host_element(PDXT, PDYT, inwater) 
    where (inwater)
      self%xp(:) = PDXT(:) ! Update only particles still in water
      self%yp(:) = PDYT(:)
      self%zp(:) = PDZT(:)
    end where

    self%zp = max(self%zp, -(2.0_SP + self%zp)) ! reflect off bottom, sigma
    self%zp = min(self%zp, -self%zp) ! reflect off free surface, sigma
    self%zpt(:) = self%cartesian(self%zp(:)) ! Calculate particle location in cartesian vertical coordinate
    self%layer(:) = self%zlocate(self%zp(:)) ! only valid for sigma layers of equal separation

    call self%INTERP_ELH(self%xp, self%yp, HL, ELL, 1) ! interpolate bathymetry and elevation
    call self%INTERP_FIELDS(self%xp, self%yp, self%zp, salinity, temperature, density, 0) ! interpolate salinity and temperature
  end subroutine

  subroutine INTERP_V(lag, XP, YP, ZP, velocity)
    ! linear interpolation of velocity field at particle positions
    use variables, only : A1U, A2U, NBE, YC, XC, DZ, ZZ, KBM1, N, KB

    real(sp), intent(in), dimension(0:N, 1:KB, 0:2) :: velocity

    real(SP) :: delta(3, 2), interpolated(2, 3)


    class(Agent), intent(inout) :: lag
    real(sp), intent(in), dimension(lag%ndrft) :: XP, YP, ZP ! ZP is sigma depth
    real(sp), dimension(0:N, 1:KB) :: UIN, VIN, WIN

    logical, dimension(lag%ndrft) ::  INWATER
    integer :: ii, host, E1, E2, E3, K1, K2, K
    real(SP) :: DUDX, DUDY, DVDX, DVDY, DWDX, DWDY, UE01, UE02, VE01, VE02, WE01, WE02
    real(SP) :: ZF1, ZF2, X0C, Y0C

    UIN = velocity(:, :, 0)
    VIN = velocity(:, :, 1)
    WIN = velocity(:, :, 2)
    INWATER(:) = .true.
    lag%FOUND(:) = 0
    call lag%find_host_element(XP, YP, INWATER) ! determine host element

    particle_loop: do ii = 1, lag%ndrft
      if ( (.not. lag%INDOMAIN(ii)) .or. (.not. INWATER(ii))) cycle ! skip particles outside domain
      host = lag%HOST(ii)
      E1  = NBE(host,1)
      E2  = NBE(host,2)
      E3  = NBE(host,3)
      X0C = XP(ii) - XC(host)
      Y0C = YP(ii) - YC(host)

      ! Determine sigma layers above and below particle
      if (ZP(ii) .gt. ZZ(1)) then ! Particle near surface
         K1  = 1
         K2  = 1
         ZF1 = 1.0_SP
         ZF2 = 0.0_SP
      else if (ZP(ii) .lt. ZZ(KBM1)) then ! Particle near bottom
         K1 = KBM1
         K2 = KBM1
         ZF1 = 1.0_SP
         ZF2 = 0.0_SP
      else
         K1 = int( (ZZ(1) - ZP(ii)) / DZ(1) ) + 1;
         K2 = K1 + 1
         ZF1 = (ZP(ii) - ZZ(K2)) / DZ(1)
         ZF2 = (ZZ(K1) - ZP(ii)) / DZ(1)
      end if


      ! Linear interpolation of velocity in sigma level above particle
      K = K1
      DUDX = A1U(host,1)*UIN(host,K)+A1U(host,2)*UIN(E1,K)+A1U(host,3)*UIN(E2,K)+A1U(host,4)*UIN(E3,K)
      DUDY = A2U(host,1)*UIN(host,K)+A2U(host,2)*UIN(E1,K)+A2U(host,3)*UIN(E2,K)+A2U(host,4)*UIN(E3,K)
      DVDX = A1U(host,1)*VIN(host,K)+A1U(host,2)*VIN(E1,K)+A1U(host,3)*VIN(E2,K)+A1U(host,4)*VIN(E3,K)
      DVDY = A2U(host,1)*VIN(host,K)+A2U(host,2)*VIN(E1,K)+A2U(host,3)*VIN(E2,K)+A2U(host,4)*VIN(E3,K)
      DWDX = A1U(host,1)*WIN(host,K)+A1U(host,2)*WIN(E1,K)+A1U(host,3)*WIN(E2,K)+A1U(host,4)*WIN(E3,K)
      DWDY = A2U(host,1)*WIN(host,K)+A2U(host,2)*WIN(E1,K)+A2U(host,3)*WIN(E2,K)+A2U(host,4)*WIN(E3,K)
      UE01 = UIN(host,K) + DUDX*X0C + DUDY*Y0C
      VE01 = VIN(host,K) + DVDX*X0C + DVDY*Y0C
      WE01 = WIN(host,K) + DWDX*X0C + DWDY*Y0C

      ! Linear interpolation of velocity in sigma level below particle
      K = K2
      DUDX = A1U(host,1)*UIN(host,K)+A1U(host,2)*UIN(E1,K)+A1U(host,3)*UIN(E2,K)+A1U(host,4)*UIN(E3,K)
      DUDY = A2U(host,1)*UIN(host,K)+A2U(host,2)*UIN(E1,K)+A2U(host,3)*UIN(E2,K)+A2U(host,4)*UIN(E3,K)
      DVDX = A1U(host,1)*VIN(host,K)+A1U(host,2)*VIN(E1,K)+A1U(host,3)*VIN(E2,K)+A1U(host,4)*VIN(E3,K)
      DVDY = A2U(host,1)*VIN(host,K)+A2U(host,2)*VIN(E1,K)+A2U(host,3)*VIN(E2,K)+A2U(host,4)*VIN(E3,K)
      DWDX = A1U(host,1)*WIN(host,K)+A1U(host,2)*WIN(E1,K)+A1U(host,3)*WIN(E2,K)+A1U(host,4)*WIN(E3,K)
      DWDY = A2U(host,1)*WIN(host,K)+A2U(host,2)*WIN(E1,K)+A2U(host,3)*WIN(E2,K)+A2U(host,4)*WIN(E3,K)
      UE02 = UIN(host,K) + DUDX*X0C + DUDY*Y0C
      VE02 = VIN(host,K) + DVDX*X0C + DVDY*Y0C
      WE02 = WIN(host,K) + DWDX*X0C + DWDY*Y0C

      ! Interpolate particle velocity between two sigma layers
      lag%UP(ii) = UE01*ZF1 + UE02*ZF2
      lag%VP(ii) = VE01*ZF1 + VE02*ZF2
      lag%WP(ii) = WE01*ZF1 + WE02*ZF2

    end do particle_loop

  end subroutine INTERP_V


  subroutine INTERP_ELH(self, XP, YP, HIN, EIN, FHE) ! OK
    ! Linearly interpolate elevation and bathymetry at a set of particle positions
    use variables, only : AW0, AWX, AWY, NV, XC, YC
    use variables, only : M

    class(Agent), intent(inout) :: self
    real(sp), intent(in), dimension(self%ndrft) :: XP, YP ! position arrays, needed because subroutine is used between updates of LAG structure
    real(sp), intent(in), dimension(0:M) :: HIN, EIN ! baythmetry and elevation inputs
    integer, intent(in) :: FHE ! Find host elements: 0 if host has correct elements; 1 if host should be updated

    integer :: ii, host, N1, N2, N3
    logical, dimension(self%ndrft) :: inwater
    real(sp) :: H0, HX, HY, E0, EX, EY, offset(2)

    if (FHE == 1) then
      inwater(:) = .true.
      self%FOUND(:) = 0
      call self%find_host_element(XP, YP, inwater)
    end if

    do ii = 1, self%ndrft
      if (.not. self%INDOMAIN(ii)) cycle ! skip particle outside domain
      host = self%host(ii) ! element containing particle
      N1 = NV(self%host(ii), 1); 
      N2 = NV(self%host(ii), 2); 
      N3 = NV(self%host(ii), 3) ! node indices
      offset(1) = XP(ii) - XC(host); 
      offset(2) = YP(ii) - YC(host) ! distance from element center

      H0 = AW0(host, 1)*HIN(N1) + AW0(host, 2)*HIN(N2) + AW0(host, 3)*HIN(N3)
      HX = AWX(host, 1)*HIN(N1) + AWX(host, 2)*HIN(N2) + AWX(host, 3)*HIN(N3)
      HY = AWY(host, 1)*HIN(N1) + AWY(host, 2)*HIN(N2) + AWY(host, 3)*HIN(N3)
      self%HP(ii) = -1.0*(H0 + HX*offset(1) + HY*offset(2)) ! Linear interpolation of bathymetry, forced to be negtaive?

      E0 = AW0(host, 1)*EIN(N1) + AW0(host, 2)*EIN(N2) + AW0(host, 3)*EIN(N3)
      EX = AWX(host, 1)*EIN(N1) + AWX(host, 2)*EIN(N2) + AWX(host, 3)*EIN(N3)
      EY = AWY(host, 1)*EIN(N1) + AWY(host, 2)*EIN(N2) + AWY(host, 3)*EIN(N3)
      self%EP(ii) = -1.0*(E0 + EX*offset(1) + EY*offset(2)) ! Linear interpolation of free surface height, forced to be positive
    end do
  end subroutine


  subroutine INTERP_FIELDS(self, XP, YP, ZP, SAL, TEMP, RHO, FHE) ! OK
    ! Linearly interpolates salinity, temperature and density
    use variables, only : AW0, AWX, AWY, NV, XC, YC, KB, M

    class(Agent), intent(inout) :: self
    real(sp), intent(in), dimension(self%ndrft) :: XP, YP, ZP
    real(sp), intent(in), dimension(0:M, KB) :: SAL, TEMP, RHO
    integer, intent(in) :: FHE ! Find host elements: 0 if host has correct elements; 1 if host should be updated

    logical, dimension(self%ndrft) :: inwater
    integer :: ii, host, N1, N2, N3
    real(sp) :: S0, SX, SY, T0, TX, TY, D0, DX, DY, offset(2), ZTMP(self%ndrft)

    ZTMP(:) = ZP(:)
    find_host: if (FHE == 1) then
      inwater(:) = .true.
      self%FOUND(:) = 0
      call self%find_host_element(XP, YP, inwater)
    end if find_host

    host = self%host(ii) ! element containing particle
    N1 = NV(host, 1); 
    N2 = NV(host, 2); 
    N3 = NV(host, 3) ! node indices
    offset(1) = XP(ii) - XC(host); 
    offset(2) = YP(ii) - YC(host) ! distance from element center

    S0 = AW0(host, 1)*SAL(N1, 1) + AW0(host, 2)*SAL(N2, 1) + AW0(host, 3)*SAL(N3, 1)
    SX = AWX(host, 1)*SAL(N1, 1) + AWX(host, 2)*SAL(N2, 1) + AWX(host, 3)*SAL(N3, 1)
    SY = AWY(host, 1)*SAL(N1, 1) + AWY(host, 2)*SAL(N2, 1) + AWY(host, 3)*SAL(N3, 1)
    self%SAL(ii) = abs(S0 + SX*offset(1) + SY*offset(2)) ! Linear interpolation of salinity field

    T0 = AW0(host, 1)*TEMP(N1, 1) + AW0(host, 2)*TEMP(N2, 1) + AW0(host, 3)*TEMP(N3, 1)
    TX = AWX(host, 1)*TEMP(N1, 1) + AWX(host, 2)*TEMP(N2, 1) + AWX(host, 3)*TEMP(N3, 1)
    TY = AWY(host, 1)*TEMP(N1, 1) + AWY(host, 2)*TEMP(N2, 1) + AWY(host, 3)*TEMP(N3, 1)
    self%TEMP(ii) = abs(T0 + TX*offset(1) + TY*offset(2)) ! Linear interpolation of temperature field

    D0 = AW0(host, 1)*RHO(N1, 1) + AW0(host, 2)*RHO(N2, 1) + AW0(host, 3)*RHO(N3, 1)
    DX = AWX(host, 1)*RHO(N1, 1) + AWX(host, 2)*RHO(N2, 1) + AWX(host, 3)*RHO(N3, 1)
    DY = AWY(host, 1)*RHO(N1, 1) + AWY(host, 2)*RHO(N2, 1) + AWY(host, 3)*RHO(N3, 1)
    self%RHO(ii) = abs(D0 + DX*offset(1) + DY*offset(2)) ! Linear interpolation of temperature field


  end subroutine INTERP_FIELDS


  subroutine INTERP_KH(self, ZKH, DZKH, KHOUT, DKHOUT)
    ! Obtain a spline interpolation on the vertical with the provided eddy diffusivity (ZKH) and its derivative (DZKH) at grid point points
    ! then linear interpolation on the horizontal
    !  RETURNS:   Both dkh/dz (DKHOUT) and kh (KHOUT)

    use variables, only: M, KB, NV, XC, Z, KBM1, YC, DTRW, AW0, AWX, AWY
    class(Agent), intent(inout) :: self
    real(SP), intent(out), dimension(self%ndrft) :: DKHOUT, KHOUT
    real(SP), intent(in), dimension(0:M, KB) :: DZKH, ZKH

    real(SP) :: offset(2), COF1, COF2, COF3
    integer :: N1, N2, N3, ii
    real(SP) :: DKHR1, DKHR2, DKHR3 !, DKHR4
    real(SP) :: KHR1, KHR2, KHR3 !, KHR4
    !real(SP) :: DDKHR1, DDKHR2, DDKHR3, DDKHR4
    real(SP) :: DKHTMP, DZP
    integer, dimension(0:KB+1) :: NZRINDX
    real(SP) :: HK, HK2, AK, AK2, AK3, BK, BK2, BK3
    integer :: KLO, KHI, NZR, host

    ! Interpolate eddy diffusivity and its derivative
    KHOUT  = zero
    DKHOUT = zero
    NZRINDX(0) = 1
    do ii = 1, KB
      NZRINDX(ii) = ii
    end do
    NZRINDX(KB+1) = KB

    do ii = 1, self%ndrft
      if (.not. self%indomain(ii)) cycle
      host = self%host(ii) ! element containing particle
      N1 = NV(host, 1); N2 = NV(host, 2); N3 = NV(host, 3) ! get node indices of host element
      offset(1) = self%XP(ii) - XC(host); offset(2) = self%YP(ii) - YC(host) ! distance from element center

      ! DERIVATIVE OF THE DIFFUSION
      ! find vertical location
      NZR = floor( float(KBM1)*abs(self%zp(ii)) ) ! guess value for hunt
      !call HUNT(Z, KB, self%ZP(ii), NZR) ! Z is sigma coordinate value from variables, only necessary if sigma layers are different thicknesses

      KHI = NZRINDX(NZR)
      KLO = NZRINDX(NZR + 1)

      ! as k in z(k) increases, sigma decreases.
      HK=Z(KHI)-Z(KLO)
      AK=(Z(KHI)-self%ZP(ii))/HK
      AK2=AK**2
      BK=(self%ZP(ii)-Z(KLO))/HK
      BK2=BK**2

      DKHR1=(-ZKH(N1,KLO)+ZKH(N1,KHI))/HK+((-3*AK2+1)*DZKH(N1,KLO)+(3*BK2-1)*DZKH(N1,KHI))*HK/6
      DKHR2=(-ZKH(N2,KLO)+ZKH(N2,KHI))/HK+((-3*AK2+1)*DZKH(N2,KLO)+(3*BK2-1)*DZKH(N2,KHI))*HK/6
      DKHR3=(-ZKH(N3,KLO)+ZKH(N3,KHI))/HK+((-3*AK2+1)*DZKH(N3,KLO)+(3*BK2-1)*DZKH(N3,KHI))*HK/6

      COF1=AW0(self%HOST(ii),1)*DKHR1 +AW0(self%HOST(ii),2)*DKHR2+AW0(self%HOST(ii),3)*DKHR3
      COF2=AWX(self%HOST(ii),1)*DKHR1 +AWX(self%HOST(ii),2)*DKHR2+AWX(self%HOST(ii),3)*DKHR3
      COF3=AWY(self%HOST(ii),1)*DKHR1 +AWY(self%HOST(ii),2)*DKHR2+AWY(self%HOST(ii),3)*DKHR3

      DKHTMP = COF1 + COF2*offset(1) + COF3*offset(2)
      DKHOUT(ii) = DKHTMP / (self%HP(ii) - self%EP(ii))     !--want the answer to be dkh/dz not dkh/dsigma, changed sign, keeney 1/7

      ! DIFFUSION ITSELF
      ! find z in grid again, as per visser, but in sigma
      DZP = self%ZP(ii) + 0.5*DKHOUT(ii)*DTRW/(self%HP(ii) - self%EP(ii)) ! changed sign, keeney 1/7
      ! adding 0.5*dkhtmp*dtrw, new z can be out of [0;-1]
      DZP = min(DZP, ZERO)
      DZP = max(DZP, -1.0_SP)

      ! find vertical location
      !call HUNT(Z, KB, DZP, NZR)
      NZR = floor( float(KBM1)*abs(DZP ) ) ! guess value for hunt
      KHI = NZRINDX(NZR)
      KLO = NZRINDX(NZR+1)

      HK = Z(KHI) - Z(KLO)
      HK2 = HK**2
      AK = (Z(KHI)-DZP)/HK
      AK3 = AK**3
      BK = (DZP-Z(KLO))/HK
      BK3 = BK**3

      KHR1 = AK*(ZKH(N1,KLO)) + BK*(ZKH(N1,KHI)) + ((AK3-AK)*(DZKH(N1,KLO))+(BK3-BK)*(DZKH(N1,KHI)))*HK2/6
      KHR2 = AK*(ZKH(N2,KLO)) + BK*(ZKH(N2,KHI)) + ((AK3-AK)*(DZKH(N2,KLO))+(BK3-BK)*(DZKH(N2,KHI)))*HK2/6
      KHR3 = AK*(ZKH(N3,KLO)) + BK*(ZKH(N3,KHI)) + ((AK3-AK)*(DZKH(N3,KLO))+(BK3-BK)*(DZKH(N3,KHI)))*HK2/6

      COF1 = AW0(self%HOST(ii),1)*KHR1 + AW0(self%HOST(ii),2)*KHR2+AW0(self%HOST(ii),3)*KHR3
      COF2 = AWX(self%HOST(ii),1)*KHR1 + AWX(self%HOST(ii),2)*KHR2+AWX(self%HOST(ii),3)*KHR3
      COF3 = AWY(self%HOST(ii),1)*KHR1 + AWY(self%HOST(ii),2)*KHR2+AWY(self%HOST(ii),3)*KHR3

      KHOUT(ii) = COF1 + COF2*offset(1) + COF3*offset(2)

    end do
  end subroutine
end module
