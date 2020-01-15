module ALL_VARS
  ! Variables
  use parameters
  implicit none
  save

  real(sp) :: DTOUT, INSTP, DHOR, DTRW, DTI
  logical :: P_SIGMA, OUT_SIGMA, F_DEPTH ! I/O position coordinate systems
  character(LEN = 80) :: CASENAME, GEOAREA, OUTDIR, INPDIR, INFOFILE, LAGINI, FOLDERPREFIX ! File Specifiers
  integer :: IOPAR, IPT, INLAG  ! File Unit Specifiers
  integer :: YEARLAG, MONTHLAG, DAYLAG, HOURLAG ! Starting date of the tracking
  integer :: IELAG, ISLAG, TDRIFT, ITOUT, IRW

  integer :: N ! Number of elements
  integer :: M ! Number of nodes
  integer :: KB ! Number of sigma levels
  integer :: KBM1 ! Number of sigma levels-1
  integer :: KBM2 ! Number of sigma levels-2
  integer :: NE ! Number of unique edges
  integer :: MX_NBR_ELEM ! Max number of elements surrounding a node

  ! Constants
  real(SP), parameter :: GRAV = 9.81_SP ! note that this is positive
  real(SP), parameter :: PI = 3.141592653_SP
  real(SP), parameter :: PI2 = 6.283185307_SP
  real(SP), parameter :: ZERO = 0.0_SP
  real(SP), parameter :: ONE_THIRD = 1.0_SP/3.0_SP
  !real(SP), parameter :: traveld = 0.023148148148_SP  ! m/s = 2000 m/day
  !real(SP), parameter :: traveld = 0.23148148148_SP  ! m/s = 20 km/day
  real(SP), parameter :: traveld = 0.5787_SP  ! m/s = 50 km/day
  real(SP), parameter :: epsx = SQRT((traveld**2.0)*0.5_SP)
  real(SP), parameter :: epsx_sigma= 0.5_SP*traveld
  !real(SP), parameter :: sal_opt = 2.0_SP
  real(SP), parameter :: sal_opt = 30.0_SP
  real(SP), parameter :: sal_sigma = 5.0_SP
  real(SP), parameter :: w1w1 = 0.5_SP
  real(SP), parameter :: h1h1 = 0.75_SP
  real(SP), parameter :: h2h2 = 0.9_SP

  ! Grid Metrics
  real(SP) :: VXMIN, VYMIN, VXMAX, VYMAX
  real(SP), allocatable :: XC(:) ! X-coord at face center
  real(SP), allocatable :: YC(:) ! Y-coord at face center
  real(SP), allocatable :: VX(:) ! X-coord at grid point
  real(SP), allocatable :: VY(:) ! Y-coord at grid point

  ! Node, Boundary Condition, and Control Volume
  integer, allocatable :: NV(:,:) ! Node numbering for elements
  integer, allocatable :: NBE(:,:) ! Indices of elmnt neighbors
  integer, allocatable :: NTVE(:)
  integer, allocatable :: ISONB(:) ! Node marker = 0,1,2
  integer, allocatable :: ISBCE(:)
  integer, allocatable :: NBVE(:,:)
  integer, allocatable :: NBVT(:,:)

  ! 1-d arrays for the sigma coordinate
  real(SP), allocatable :: Z(:) ! Sigma coordinate value
  real(SP), allocatable :: ZZ(:) ! Intra level sigma value
  real(SP), allocatable :: DZ(:) ! Delta-sigma value
  real(SP), allocatable :: DZZ(:) ! Delta of intra level sigma

  ! 2-d flow variable arrays at nodes
  real(SP), allocatable :: H(:) ! Bathymetric depth
  real(SP), allocatable :: D(:) ! Current depth
  real(SP), allocatable :: EL(:) ! Current surface elevation
  real(SP), allocatable :: ET(:) ! Surface elevation at previous time step

  ! internal mode arrays (element based)
  real(SP), allocatable :: U(:,:) ! x-velocity
  real(SP), allocatable :: V(:,:) ! Y-velocity
  real(SP), allocatable :: W(:,:) ! Vertical velocity in sigma system
  real(SP), allocatable :: WW(:,:) ! Z-velocity
  real(SP), allocatable :: UT(:,:) ! X-velocity from previous timestep
  real(SP), allocatable :: VT(:,:) ! Y-velocity from previous timestep
  real(SP), allocatable :: WT(:,:) ! Velocity sigma from previous timestep
  real(SP), allocatable :: WWT(:,:) ! Z-velocity from previous timesteps

  ! 3d variable arrays-(node based)
  real(SP), allocatable :: T1(:,:) ! Temperature at nodes
  real(SP), allocatable :: S1(:,:) ! Salinity at nodes
  real(SP), allocatable :: R1(:,:) ! density at nodes
  real(SP), allocatable :: TT1(:,:) ! Temperature from previous time
  real(SP), allocatable :: ST1(:,:) ! Salinity from previous time
  real(SP), allocatable :: RT1(:,:) ! density from previous time
  real(SP), allocatable :: WTS(:,:) ! Vertical velocity in sigma system
  real(SP), allocatable :: KH(:,:) ! Turbulent diffusivity

  ! shape coefficient arrays and control volume metrics
  real(SP), allocatable :: A1U(:,:)
  real(SP), allocatable :: A2U(:,:)
  real(SP), allocatable :: AWX(:,:)
  real(SP), allocatable :: AWY(:,:)
  real(SP), allocatable :: AW0(:,:)

end module ALL_VARS