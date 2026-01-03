module Lagrangian
  ! type and variables for lagrangian particle system
  use variables, only : sp
  implicit none
  save
  private

  type, public :: Agent
    character(len = 20) :: species    ! string for identification
    integer :: count                  ! total particles
    integer, allocatable, dimension(:) :: &
        & host, &                     ! Element containing particle
        & layer                    ! sigma layer
    real(SP), allocatable, dimension(:) :: &
        & XP, YP, ZP, &               ! position of particle
        & HP, &                       ! bathymetry at particle position
        & EP, &                       ! surface height at particle
        & UP, VP, &                   ! velocity of particle
        & TEMP, &                     ! temperature at particle position
        & SAL, &                      ! salinity at particle position
        & RHO                         ! density at particle position

  contains

    procedure, public :: lag_alloc => lag_alloc
    procedure, public :: readPosition => readPositions ! read position and ndrft from file
    procedure, public :: writePosition => lag_writePosition ! write position to file
    procedure, public :: sigma => lag_cart2sig ! convert cartesian to sigma
    procedure, public :: zlocate => layer ! the layer the particle is below
    procedure, public :: zinterp => zinterp
    procedure, public :: find_host_element => find_host_element
    procedure, public :: interp_fields => interp_fields

  end type

contains

  subroutine find_host_element(self, x, y)
    !  Find host elements of particles by searching progressively further elements
    use variables, only: N, NV, VX, NV, VY, XC, YC, NTVE, NBVE
    class(Agent), intent(inout) :: self
    real(sp), dimension(self%count), intent(in) :: x, y
    real(sp), dimension(1:N, 1) :: distance
    real(sp) :: previous
    integer :: ii, jj, kk, ind, vertices(3), triangle, vertex
    integer, dimension(2) :: nearby

    do ii = 1, self%count
      triangle = self%host(ii)
      vertices = NV(triangle, 1:3)
      if (isintriangle(VX(vertices), VY(vertices), x(ii), y(ii))) then
        ! is in same triangle
        cycle
      end if
      ! search neighbors
      do jj = 1, 3
        do kk = 1, NTVE(vertices(jj))
          triangle = NBVE(vertices(jj), kk)
          vertices = NV(triangle, 1:3)
          if (isintriangle(VX(vertices), VY(vertices), x(ii), y(ii))) then
            self%host(ii) = ind
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
        triangle = nearby(1)
        vertices = NV(triangle, 1:3)
        if (ISINTRIANGLE(VX(vertices), VY(vertices) , x(ii), y(ii))) then
          self%host(ii) = nearby(1)
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
    allocate( self%XP( self%count ), &
            & self%YP( self%count ), &
            & self%ZP( self%count ), &
            & self%HP( self%count ), &
            & self%EP( self%count ), &
            & self%HOST( self%count ), &
            & self%LAYER( self%count ), &
            & self%TEMP( self%count ), &
            & self%SAL( self%count ), &
            & self%RHO( self%count ), &
            & self%UP( self%count ), &
            & self%VP( self%count ))
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
    read(fid, "(I6)") self%count ! read number of particles
    allocate(self%xp(self%count), &
            & self%yp(self%count), &
            & self%zp(self%count))
    do ii = 1, self%count
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
    write(fid, "(1F10.2,9000(3F20.3))") time, (self%XP(ii), self%YP(ii), self%ZP(ii), ii=1,self%count)
  end subroutine

  pure function lag_cart2sig(self, cartesian) result(sigma)
    ! Calculate sigma vertical position from cartesian
    class(Agent), intent(in) :: self
    real(sp), dimension(self%count), intent(in) :: cartesian
    real(sp), dimension(self%count) :: sigma
    sigma(:) = -1.0_SP * abs(cartesian(:) - self%EP(:)) / abs(self%EP(:) - self%HP(:))
  end function

  pure function layer(self, sigma) result(layers)
    ! update current sigma layer of particles between moves
    use variables, only : KBM1
    class(Agent), intent(in) :: self
    real(sp), dimension(self%count), intent(in) :: sigma
    integer, dimension(self%count) :: layers
    layers(:) = ceiling(-KBM1*sigma(:))
    where (layers(:) < 1) layers(:) = 1
  end function

  pure function zinterp(self, verticalvar) result(interpolated)
    use variables, only : layers ! sigma layers
    use simulation, only : domain ! domain structure
    class(Agent), intent(in) :: self ! lagrangian particle swarm object
    real(sp), dimension(0:layers+1), intent(in) :: verticalvar
    real(sp), dimension(self%count) :: idz, interpolated
    idz(:) = (domain%layerSigma*(self%layer(:)-1) - self%zp(:))/domain%layerSigma ! relative distance from layer above
    interpolated(:) = verticalvar(self%layer(:))*(1.0_SP - idz(:)) + verticalvar(self%layer(:)+1)*idz(:)
  end function

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

  pure function linear_interpolation_at_layer(offset, triangle, field, layer) result(sample)
    ! Linearly interpolates a field at a point in discrete layer
    use variables, only : NV, layers, M
    real(sp), intent(in) :: offset(2)
    real(sp), intent(in), dimension(0:M, layers) :: field
    integer, intent(in) :: triangle, layer
    real(sp) :: sample, stencil(3)
    stencil = field(NV(triangle, 1:3), layer)
    sample = linear_planar_interpolation(offset, triangle, stencil)
  end function

  pure function linear_interpolation_of_surface(offset, triangle, field) result(sample)
    ! Linearly interpolate a 2D surface at a point
    use variables, only : NV, M, XC, YC
    real(sp), intent(in) :: offset(2)
    real(sp), intent(in), dimension(0:M) :: field
    integer, intent(in) :: triangle ! element and layer containing particle
    real(sp) :: sample, stencil(3)
    stencil = field(NV(triangle, 1:3))
    sample = linear_planar_interpolation(offset, triangle, stencil)
  end function

  subroutine INTERP_FIELDS(self, XP, YP, ZP, salinity, temperature, density, depth, surface)
    ! Linearly interpolates salinity, temperature and density
    use variables, only : layers, M, XC, YC

    class(Agent), intent(inout) :: self
    real(sp), intent(in), dimension(self%count) :: XP, YP, ZP
    real(sp), intent(in), dimension(0:M, layers) :: salinity, temperature, density
    real(sp), intent(in), dimension(0:M) :: depth, surface
    real(sp), dimension(2) :: offset

    integer :: ii, host
    call self%find_host_element(XP, YP)
    do ii = 1, self%count
      host = self%host(ii) ! element containing particle
      offset = [XP(ii) - XC(host), YP(ii) - YC(host)]
      self%SAL(ii) = linear_interpolation_at_layer(offset, host, salinity, 1)
      self%TEMP(ii) = linear_interpolation_at_layer(offset, host, temperature, 1)
      self%RHO(ii) = linear_interpolation_at_layer(offset, host, density, 1)
      self%HP(ii) = -linear_interpolation_of_surface(offset, host, depth)
      self%EP(ii) = -linear_interpolation_of_surface(offset, host, surface)
    end do
  end subroutine
end module
