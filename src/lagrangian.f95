module Lagrangian
  ! type and variables for lagrangian particle system
  use variables, only : sp
  implicit none
  save
  private

  type, public :: Agent
    character(len = 20) :: species    ! string for identification
    integer :: ndrft                  ! total particles
    integer, allocatable, dimension(:) :: &
        & host, &                     ! Element containing particle
        & layer                    ! sigma layer
    real(SP), allocatable, dimension(:) :: &
        & XP, YP, ZP, &               ! position of particle
        & HP, &                       ! Bathymetry at particle position
        & EP, &                       ! Free surface height at particle
        & UP, VP, WP, &               ! velocity of particle
        & TEMP, &                     ! temperature at particle position
        & SAL, &                      ! salinity at particle position
        & RHO                         ! density at particle position

  contains

    procedure, public :: lag_alloc => lag_alloc
    procedure, public :: stats => lag_printStatistics ! print particle information if desired
    procedure, public :: readPosition => readPositions ! read position and ndrft from file
    procedure, public :: writePosition => lag_writePosition ! write position to file
    procedure, public :: sigma => lag_cart2sig ! convert cartesian to sigma
    procedure, public :: cartesian => lag_sig2cart ! convert sigma to cartesian
    procedure, public :: zlocate => lag_getlayers ! update the sigma layer that particles are currently below
    procedure, public :: zinterp => lag_zinterp
    procedure, public :: find_host_element => find_host_element
    procedure, public :: traject => traject
    procedure, public :: interp_fields => interp_fields
    procedure, public :: interp_v => interp_v
    procedure, public :: interp_kh => interp_kh

  end type

contains

  subroutine find_host_element(lag, x, y)
    !  Find host elements of particles by searching progressively further elements
    use variables, only: N, NV, VX, NV, VY, XC, YC, NTVE, NBVE
    class(Agent), intent(inout) :: lag
    real(sp), dimension(lag%ndrft), intent(in) :: x, y

    real(sp), dimension(1:N, 1) :: distance
    real(sp) :: previous
    integer :: ii, jj, kk, ind
    integer, dimension(2) :: nearby

    do ii = 1, lag%ndrft
      if (isintriangle(VX(NV(lag%HOST(ii), 1:3)), VY(NV(lag%HOST(ii), 1:3)), x(ii), y(ii))) then
        cycle
      end if

      do jj = 1, 3
        do kk = 1, NTVE(NV(lag%HOST(ii), jj))
          ind = NBVE(NV(lag%HOST(ii), jj), kk)
          if (isintriangle(VX(NV(ind, 1:3)), VY(NV(ind, 1:3)), x(ii), y(ii))) then
            lag%HOST(ii) = ind
            exit
          end if
        end do
      end do

      distance(1:n, 1) = sqrt((xc(1:n) - x(ii))**2 + (yc(1:n) - y(ii))**2)
      previous = 0.0_sp

      do jj = 1, 16
        nearby(:) = minloc(distance, distance > previous)
        if (nearby(1) == 0) exit
        previous = distance(nearby(1), 1)
        if (ISINTRIANGLE(VX(NV(nearby(1), 1:3)), VY(NV(nearby(1), 1:3)) , x(ii), y(ii))) then
          lag%host(ii) = nearby(1)
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
            & self%HOST( self%ndrft ), &
            & self%LAYER( self%ndrft ), &
            & self%TEMP( self%ndrft ), &
            & self%SAL( self%ndrft ), &
            & self%RHO( self%ndrft ), &
            & self%UP( self%ndrft ), &
            & self%VP( self%ndrft ), &
            & self%WP( self%ndrft ))

    self%XP(:) = 0.0_sp
    self%YP(:) = 0.0_sp
    self%ZP(:) = 0.0_sp
    self%HP(:) = 0.0_sp
    self%EP(:) = 0.0_sp
    self%HOST(:) = 1
    self%LAYER(:) = 1
    self%TEMP(:) = 0.0_sp
    self%SAL(:) = 0.0_sp
    self%RHO(:) = 0.0_sp
    self%UP(:) = 0.0_sp
    self%VP(:) = 0.0_sp
    self%WP(:) = 0.0_sp

  end subroutine

  subroutine lag_printStatistics(self)
    class(Agent), intent(inout) :: self
    write(*, *)
    write(*, *) '    Particle Class'
    write(*, *) '        Species       : ', self%species
    write(*, *) '        Quantity      : ', self%ndrft
    write(*, *)
  end subroutine

  subroutine readPositions(self)
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
    allocate(self%xp(self%ndrft), &
            & self%yp(self%ndrft), &
            & self%zp(self%ndrft))
    do ii = 1, self%ndrft
      read(fid, "(3F20.6)") self%xp(ii), self%yp(ii), self%zp(ii)
    end do

    close(fid)
  end subroutine

  subroutine lag_writePosition(self, fid, time)
    ! write time and particle id/position to already open file
    class(Agent), intent(inout) :: self ! lagrangian particle structure
    integer, intent(in) :: fid ! unit number of open output file
    real(sp), intent(in) :: time ! time to write
    integer :: ii
    write(fid, "(1F10.2,9000(3F20.3))") time, (self%XP(ii), self%YP(ii), self%ZP(ii), ii=1,self%ndrft)
  end subroutine

  pure function lag_cart2sig(self, cartesian) result(sigma)
    ! Calculate sigma vertical position from cartesian
    class(Agent), intent(in) :: self
    real(sp), dimension(self%ndrft), intent(in) :: cartesian
    real(sp), dimension(self%ndrft) :: sigma
    sigma(:) = -1.0_SP * abs(cartesian(:) - self%EP(:)) / abs(self%EP(:) - self%HP(:))
  end function

  pure function lag_sig2cart(self, sigma) result(depth)
    ! calculate cartesian vertical from sigma coordinate
    class(Agent), intent(in) :: self
    real(sp), dimension(self%ndrft), intent(in) :: sigma
    real(sp), dimension(self%ndrft) :: depth
    depth(:) = sigma(:)*(self%EP(:) - self%HP(:)) + self%EP(:)
  end function

  pure function lag_getlayers(self, sigma) result(layers)
    ! update current sigma layer of particles between moves
    use variables, only : KBM1
    class(Agent), intent(in) :: self
    real(sp), dimension(self%ndrft), intent(in) :: sigma
    integer, dimension(self%ndrft) :: layers
    layers(:) = ceiling(-KBM1*sigma(:))
    where (layers(:) < 1) layers(:) = 1
  end function

  pure function lag_zinterp(self, verticalvar) result(interpolated)
    use variables, only : KB ! sigma layers
    use simulation, only : domain ! domain structure
    class(Agent), intent(in) :: self ! lagrangian particle swarm object
    real(sp), dimension(0:KB+1), intent(in) :: verticalvar
    real(sp), dimension(self%ndrft) :: idz, interpolated

    idz(:) = (domain%layerSigma*(self%layer(:)-1) - self%zp(:))/domain%layerSigma ! relative distance from layer above
    interpolated(:) = verticalvar(self%layer(:))*(1.0_SP - idz(:)) + verticalvar(self%layer(:)+1)*idz(:)

  end function


  subroutine TRAJECT(self, dt, U1, U2, V1, V2, W1, W2, HL, EL1, EL2, salinity, temperature, density)

    ! integrate particle position from x0 to xn using velocity fields at time t0 and time tn
    use variables, only : N, M, KB, MSTAGE, A_RK, B_RK, C_RK ! runge-kutta parameters

    class(Agent), intent(inout) :: self ! lagrangian particle object
    real(SP), intent(in) :: dt ! time step, usually DTI
    real(SP), dimension(0:N, KB), intent(in) :: U1, U2, V1, V2, W1, W2 ! velocity fields at start and end of time step
    real(SP), dimension(0:M), intent(in) :: HL, EL1, EL2 ! bathymetry and free surface height
    real(SP), dimension(0:M, KB), intent(in) :: temperature, salinity, density

    real(SP), dimension(self%ndrft) :: PDX, PDY, PDZ ! RK stage positions
    integer :: stage
    real(SP), dimension(0:N, KB, 0:2) :: velocity ! ERK stage velocity field
    real(SP), dimension(0:M) :: ELL ! ERK stage freesurface height
    real(SP), dimension(self%ndrft, 0:MSTAGE, 3) :: CHI ! ERK stage function evaluation for velocities
    real(SP), parameter :: EPS  = 10.0 ** (-5.0) ! depth of dry element

    CHI = 0.0_sp ! Initialize Stage Functional Evaluations

    do stage = 1, MSTAGE ! Runge-Kutta integration stages
      ! New stage position updated from time zero position
      PDX(:) = self%yp(:) + A_RK(stage) * dt * CHI(:, stage - 1, 1)
      PDY(:) = self%yp(:) + A_RK(stage) * dt * CHI(:, stage - 1, 2)
      PDZ(:) = self%zp(:) + A_RK(stage) * dt * CHI(:, stage - 1, 3)
      PDZ(:) = max(PDZ(:), -(2.0 + PDZ(:))) ! reflect sigma depth off bottom
      PDZ(:) = min(PDZ(:), 0.0_sp) ! keep sigma depth below free surface

      ! Calculate velocity field for stage using c_rk coefficients
      velocity(:, :, 0) = (1.0 - C_RK(stage)) * U1 + C_RK(stage) * U2
      velocity(:, :, 1) = (1.0 - C_RK(stage)) * V1 + C_RK(stage) * V2
      velocity(:, :, 2) = (1.0 - C_RK(stage)) * W1 + C_RK(stage) * W2
      ELL = (1.0 - C_RK(stage)) * EL1 + C_RK(stage) * EL2

      call self%INTERP_V(PDX, PDY, PDZ, velocity)
      call self%INTERP_FIELDS(PDX, PDY, PDZ, salinity, temperature, density, HL, ELL)

      CHI(:, stage, 1) = self%UP ! Update CHI values for next time step
      CHI(:, stage, 2) = self%VP

      where ((self%EP - self%HP) < EPS) 
        CHI(:, stage, 3) = 0.0_sp ! Limit vertical motion in very shallow water
      elsewhere
        CHI(:, stage, 3) = self%WP / (self%HP - self%EP)    ! delta_sigma/deltaT = ww/D
      end where
    end do

    do stage = 1, MSTAGE
      self%xp(:) = self%xp(:) + dt * CHI(:, stage, 1) * B_RK(stage)
      self%yp(:) = self%yp(:) + dt * CHI(:, stage, 2) * B_RK(stage)
      self%zp(:) = self%zp(:) + dt * CHI(:, stage, 3) * B_RK(stage)
    end do

    self%zp = max(self%zp, -(2.0_SP + self%zp)) ! reflect off bottom
    self%zp = min(self%zp, -self%zp) ! reflect off free surface
    self%layer(:) = self%zlocate(self%zp(:)) ! only valid for sigma layers of equal separation
    call self%INTERP_FIELDS(self%xp, self%yp, self%zp, salinity, temperature, density, HL, ELL)
  end subroutine

  subroutine INTERP_V(lag, XP, YP, ZP, velocity)
    ! linear interpolation of velocity field at particle positions
    use variables, only : A1U, A2U, NBE, YC, XC, DZ, ZZ, KBM1, N, KB

    real(sp), intent(in), dimension(0:N, 1:KB, 0:2) :: velocity
    class(Agent), intent(inout) :: lag
    real(sp), intent(in), dimension(lag%ndrft) :: XP, YP, ZP ! ZP is sigma depth
    real(sp), dimension(0:N, 1:KB) :: UIN, VIN, WIN
    integer :: ii, host, E1, E2, E3, K1, K2, K
    real(SP) :: DUDX, DUDY, DVDX, DVDY, DWDX, DWDY, UE01, UE02, VE01, VE02, WE01, WE02
    real(SP) :: ZF1, ZF2, X0C, Y0C

    UIN = velocity(:, :, 0)
    VIN = velocity(:, :, 1)
    WIN = velocity(:, :, 2)
    call lag%find_host_element(XP, YP) ! determine host element

    particle_loop: do ii = 1, lag%ndrft
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

  pure function linear_planar_interpolation(offset, host, stencil) result(sample)
    ! Linearly interpolate a 2D planar field at a point, base for 3D
    use variables, only : AW0, AWX, AWY
    real(sp), intent(in) :: offset(2), stencil(3)
    integer, intent(in) :: host ! element and layer containing particle
    real(sp) :: f, fx, fy
    real(sp) :: sample
    f = sum(AW0(host, 1:3) * stencil)
    fx = sum(AWX(host, 1:3) * stencil) * offset(1)
    fy = sum(AWY(host, 1:3) * stencil) * offset(2)
    sample = abs(f + fx + fy)
  end function

  pure function linear_interpolation_at_layer(offset, host, field, layer) result(sample)
    ! Linearly interpolates a field at a point in discrete layer
    use variables, only : NV, KB, M
    real(sp), intent(in) :: offset(2)
    real(sp), intent(in), dimension(0:M, KB) :: field
    integer, intent(in) :: host, layer ! element and layer containing particle
    real(sp) :: sample, stencil(3)
    stencil = field(NV(host, 1:3), layer)
    sample = linear_planar_interpolation(offset, host, stencil)
  end function

  pure function linear_interpolation_of_surface(offset, host, field) result(sample)
    ! Linearly interpolate a 2D surface at a point
    use variables, only : NV, M, XC, YC
    real(sp), intent(in) :: offset(2)
    real(sp), intent(in), dimension(0:M) :: field
    integer, intent(in) :: host ! element and layer containing particle
    real(sp) :: sample, stencil(3)
    stencil = field(NV(host, 1:3))
    sample = linear_planar_interpolation(offset, host, stencil)
  end function

  subroutine INTERP_FIELDS(self, XP, YP, ZP, SAL, TEMP, RHO, HIN, EIN) ! OK
    ! Linearly interpolates salinity, temperature and density
    use variables, only : KB, M, XC, YC

    class(Agent), intent(inout) :: self
    real(sp), intent(in), dimension(self%ndrft) :: XP, YP, ZP
    real(sp), intent(in), dimension(0:M, KB) :: SAL, TEMP, RHO
    real(sp), intent(in), dimension(0:M) :: HIN, EIN
    real(sp), dimension(2) :: offset

    integer :: ii, host
    call self%find_host_element(XP, YP)
    do ii = 1, self%ndrft
      host = self%host(ii) ! element containing particle
      offset = [XP(ii) - XC(host), YP(ii) - YC(host)]
      self%SAL(ii) = linear_interpolation_at_layer(offset, host, SAL, 1)
      self%TEMP(ii) = linear_interpolation_at_layer(offset, host, TEMP, 1)
      self%RHO(ii) = linear_interpolation_at_layer(offset, host, RHO, 1)
      self%HP(ii) = -linear_interpolation_of_surface(offset, host, HIN)
      self%EP(ii) = -linear_interpolation_of_surface(offset, host, EIN)
    end do
  end subroutine

  pure subroutine INTERP_KH(self, ZKH, DZKH, KHOUT, DKHOUT)
    ! Obtain a spline interpolation on the vertical with eddy diffusivity (ZKH) and 
    ! its derivative (DZKH) at grid points, then linear interpolation on the horizontal
    !  RETURNS:   Both dkh/dz (DKHOUT) and kh (KHOUT)

    use variables, only: M, KB, NV, XC, Z, KBM1, YC, DTRW, AW0, AWX, AWY
    class(Agent), intent(in) :: self
    real(SP), intent(out), dimension(self%ndrft) :: DKHOUT, KHOUT
    real(SP), intent(in), dimension(0:M, KB) :: DZKH, ZKH
    real(SP) :: offset(2)
    integer :: N1, N2, N3, ii
    real(SP) :: DKHR1, DKHR2, DKHR3
    real(SP) :: KHR1, KHR2, KHR3
    real(SP) :: DKHTMP, DZP
    integer, dimension(0:KB+1) :: NZRINDX
    real(SP) :: HK, HK2, AK, AK2, AK3, BK, BK2, BK3, stencil(3)
    integer :: KLO, KHI, NZR, host

    ! Interpolate eddy diffusivity and its derivative
    KHOUT  = 0.0_sp
    DKHOUT = 0.0_sp
    NZRINDX(0) = 1
    do ii = 1, KB
      NZRINDX(ii) = ii
    end do
    NZRINDX(KB+1) = KB

    do ii = 1, self%ndrft

      host = self%host(ii)
      N1 = NV(host, 1)
      N2 = NV(host, 2)
      N3 = NV(host, 3)
      offset(1) = self%XP(ii) - XC(host)
      offset(2) = self%YP(ii) - YC(host) ! distance from element center

      ! DERIVATIVE OF THE DIFFUSION
      ! find vertical location
      NZR = floor( float(KBM1) * abs(self%zp(ii)) )
      KHI = NZRINDX(NZR)
      KLO = NZRINDX(NZR + 1)

      ! as k in z(k) increases, sigma decreases.
      HK= Z(KHI) - Z(KLO)
      HK2 = HK**2
      AK = (Z(KHI) - self%ZP(ii))/HK
      AK2 = AK**2
      BK = (self%ZP(ii) - Z(KLO))/HK
      BK2=BK**2

      DKHR1=(-ZKH(N1,KLO)+ZKH(N1,KHI))/HK+((-3*AK2+1)*DZKH(N1,KLO)+(3*BK2-1)*DZKH(N1,KHI))*HK/6
      DKHR2=(-ZKH(N2,KLO)+ZKH(N2,KHI))/HK+((-3*AK2+1)*DZKH(N2,KLO)+(3*BK2-1)*DZKH(N2,KHI))*HK/6
      DKHR3=(-ZKH(N3,KLO)+ZKH(N3,KHI))/HK+((-3*AK2+1)*DZKH(N3,KLO)+(3*BK2-1)*DZKH(N3,KHI))*HK/6
      stencil = [DKHR1, DKHR2, DKHR3]
      DKHTMP = linear_planar_interpolation(offset, host, stencil)
      DKHOUT(ii) = DKHTMP / (self%HP(ii) - self%EP(ii))

      ! DIFFUSION ITSELF
      ! find z in grid again, as per visser, but in sigma
      DZP = self%ZP(ii) + 0.5*DKHOUT(ii)*DTRW/(self%HP(ii) - self%EP(ii))
      DZP = min(DZP, 0.0_sp)
      DZP = max(DZP, -1.0_SP)

      ! find vertical location
      NZR = floor( float(KBM1)*abs(DZP) ) ! guess value for hunt
      KHI = NZRINDX(NZR)
      KLO = NZRINDX(NZR+1)

    
      AK = (Z(KHI)-DZP)/HK
      AK3 = AK**3
      BK = (DZP-Z(KLO))/HK
      BK3 = BK**3

      KHR1 = AK*(ZKH(N1,KLO)) + BK*(ZKH(N1,KHI)) + ((AK3-AK)*(DZKH(N1,KLO))+(BK3-BK)*(DZKH(N1,KHI)))*HK2/6
      KHR2 = AK*(ZKH(N2,KLO)) + BK*(ZKH(N2,KHI)) + ((AK3-AK)*(DZKH(N2,KLO))+(BK3-BK)*(DZKH(N2,KHI)))*HK2/6
      KHR3 = AK*(ZKH(N3,KLO)) + BK*(ZKH(N3,KHI)) + ((AK3-AK)*(DZKH(N3,KLO))+(BK3-BK)*(DZKH(N3,KHI)))*HK2/6
      stencil = [KHR1, KHR2, KHR3]
      
      KHOUT(ii) = linear_planar_interpolation(offset, host, stencil)

    end do
  end subroutine
end module
