module variables
    implicit none
    save

    integer, parameter :: &
        & iocp = 101, &
        & iocs = 102, &
        & iotox = 103, &
        & iofp = 201, &
        & iofs = 202, &
        & iophys = 301, &
        & sp = SELECTED_REAL_KIND(12,300), & ! double precision, single -> (6,30)
        & MSTAGE = 4 ! number of Runge-Kutta integration stages
  
    ! physical and mathematical constants
    real(sp), parameter :: &
        & ZERO = 0.0_sp, &
        & irradSurf = 650.0_sp ! W/M^2

    ! Runge-Kutta integration coefficients
    real(sp), parameter, dimension(4) :: &
        & A_RK = (/ 0.0_sp, 0.5_sp, 0.5_sp, 1.0_sp/), &
        & B_RK = (/ 1.0_sp/6.0_sp, 1.0_sp/3.0_sp, 1.0_sp/3.0_sp, 1.0_sp/6.0_sp /), &
        & C_RK = (/ 0.0_sp, 0.5_sp, 0.5_sp, 1.0_sp /)

    character(LEN = 80) :: CASENAME, FOLDERPREFIX ! File Specifiers

    integer :: &
        & HOURLAG, &
        & IELAG, &
        & ISLAG, &
        & TDRIFT, &
        & N, &        ! Number of elements
        & M, &        ! Number of nodes
        & KB, &       ! Number of sigma levels
        & KBM1, &     ! Number of sigma levels-1
        & KBM2, &     ! Number of sigma levels-2
        & NE, &       ! Number of unique edges
        & MX_NBR_ELEM ! Max number of elements surrounding a node

    real(sp) :: &
        & DTOUT, &
        & INSTP, &
        & DHOR, &
        & DTRW, &
        & DTI, &
        & VXMIN, &
        & VYMIN, &
        & VXMAX, &
        & VYMAX ! Grid metrics
  
    ! Shape coefficient arrays and control volume metrics
    real(sp), allocatable, dimension(:, :) :: &
        & A1U, &
        & A2U, &
        & AWX, &
        & AWY, &
        & AW0

    ! Node, boundary condition, and control volume
    integer, allocatable :: &
        & NV(:, :), &   ! Node numbering for elements
        & NBE(:, :), &  ! Indices of element neighbors
        & NTVE(:), &
        & ISONB(:), &   ! Node marker = 0,1,2
        & ISBCE(:), &
        & NBVE(:, :), &
        & NBVT(:, :)

    ! 1-d arrays for the sigma coordinate and surfaces
    real(sp), allocatable, dimension(:) :: &
        & Z, &    ! Sigma coordinate value
        & ZZ, &   ! Intra level sigma value
        & DZ, &   ! Delta-sigma value
        & H, &    ! bathymetry
        & D, &    ! depth
        & EL, &   ! surface elevation
        & ET, &   ! surface elevation, previous time step
        & XC, &   ! cell center
        & YC, &
        & VX, &   ! vertices
        & VY

    ! internal mode arrays (cell-based)
    real(sp), allocatable, dimension(:, :) :: &
        & U, &  ! x-velocity
        & V, &  ! y-velocity
        & WW, & ! z-velocity
        & UT, & ! x-velocity, previous
        & VT, & ! y-velocity, previous
        & WWT   ! z-velocity, previous

    ! vertex-based arrays
    real(sp), allocatable, dimension(:, :) :: &
        & T1, &   ! temperature at nodes
        & S1, &   ! salinity at nodes
        & R1, &   ! density at nodes
        & TT1, &  ! temperature, previous
        & ST1, &  ! salinity, previous
        & RT1, &  ! density, previous
        & KH      ! turbulent diffusivity
end module
