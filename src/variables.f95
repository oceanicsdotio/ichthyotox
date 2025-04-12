module variables
    implicit none
    save
  
    logical, parameter :: strict_integration =  .false. ! set mass transfer

    integer, parameter :: &
        & iocp = 101, iocs = 102, iotox = 103, iofp = 201, &
        & iofs = 202, iophys = 301, iorun = 302, iovar=303, &
        & sp = SELECTED_REAL_KIND(12,300), & ! double precision, single -> (6,30)
        & MSTAGE = 4 ! number of Runge-Kutta integration stages
  
    ! physical and mathematical constants
    real(sp), parameter :: &
        & GRAV = 9.81_sp, & ! note that this is positive
        & ZERO = 0.0_sp, &
        & ONE_THIRD = 1.0_sp/3.0_sp, &
        & h1h1 = 0.75_sp, &
        & h2h2 = 0.9_sp, &
        & irradSurf = 650.0_sp ! W/M^2

    ! Runge-Kutta integration coefficients
    real(sp), parameter, dimension(4) :: &
        & A_RK = (/ 0.0_sp, 0.5_sp, 0.5_sp, 1.0_sp/), &
        & B_RK = (/ 1.0_sp/6.0_sp, ONE_THIRD, ONE_THIRD, 1.0_sp/6.0_sp /), &
        & C_RK = (/ 0.0_sp, 0.5_sp, 0.5_sp, 1.0_sp /)

    logical :: P_SIGMA, OUT_SIGMA, F_DEPTH ! I/O position coordinate systems
 
    character(LEN = 80) :: CASENAME, GEOAREA, OUTDIR, INPDIR, INFOFILE, LAGINI, FOLDERPREFIX ! File Specifiers

    integer :: &
        & IOPAR, IPT, INLAG, &  ! File Unit Specifiers
        & YEARLAG, MONTHLAG, DAYLAG, HOURLAG, IELAG, ISLAG, TDRIFT, ITOUT, IRW, &
        & N, &        ! Number of elements
        & M, &        ! Number of nodes
        & KB, &       ! Number of sigma levels
        & KBM1, &     ! Number of sigma levels-1
        & KBM2, &     ! Number of sigma levels-2
        & NE, &       ! Number of unique edges
        & MX_NBR_ELEM ! Max number of elements surrounding a node
  
    real(sp) :: &
        & DTOUT, INSTP, DHOR, DTRW, DTI, &
        & VXMIN, VYMIN, VXMAX, VYMAX ! Grid metrics
  
    ! Shape coefficient arrays and control volume metrics
    real(sp), allocatable, dimension(:, :) :: A1U, A2U, AWX, AWY, AW0

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
        & DZZ, &  ! Delta of intra level sigma
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
        & W, &  ! sigma-velocity
        & WW, & ! z-velocity
        & UT, & ! x-velocity, previous timestep
        & VT, & ! y-velocity, previous timestep
        & WT, & ! sigma-velocity, previous timestep
        & WWT   ! z-velocity, previous timesteps

    ! vertex-based arrays
    real(sp), allocatable, dimension(:, :) :: &
        & T1, &   ! temperature at nodes
        & S1, &   ! salinity at nodes
        & R1, &   ! density at nodes
        & TT1, &  ! temperature, previous time
        & ST1, &  ! salinity, previous time
        & RT1, &  ! density, previous time
        & WTS, &  ! vertical velocity in sigma system
        & KH      ! turbulent diffusivity
end module