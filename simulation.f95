module simulation

  use variables, only : ZERO, sp
  implicit none
  save ! State is saved in the compiled object

  private

  type, public :: LAG_SIM
    character(len = 100), public :: simID
    integer, public :: nnodes=0, nelements=0, nlayers=0, lines_read=0 ! simulation id, size and state
    real(sp), public :: globalIrradiance=ZERO, meshArea=ZERO, layerDepth=ZERO, layerSigma=ZERO, time=ZERO, daytime=ZERO, clocktime=ZERO ! domain variables and clocks: elapsed, divided days, and twenty four hour periodic
    real(sp), allocatable, dimension(:), private :: elementSigmaVolume, elementArea ! mesh stats
    real(sp), allocatable, dimension(:), public :: verticaltox, verticaldiff, verticaltemp, verticalrho ! uniform horizontal fields

  contains
    ! call in this order
    procedure, public :: init => simulation_initialize
    procedure, public :: load => simulation_read
    procedure, public :: geo => simulation_geometry
    procedure, public :: vdiff => simulation_diffusion

  end type LAG_SIM;
  class(LAG_SIM), allocatable, public :: domain ! domain structure imported from this module


  ! creating gaussian distributions, etc.
  type, public :: LAG_RAND

    real(sp), private :: rn1, rn2, ru1, ru2, sumMeanDiffSq, mean
    logical, private  :: current, statistics
    integer, private  :: samples

  contains

    procedure, public :: init => random_initialize ! subroutine initializes generator from system clock
    procedure, private :: normal => random_normal ! subroutine generates two gaussian randoms which are stored in the object, also calculates statisitcs on the fly
    procedure, public :: get => random_get ! function returns the value of one stored gaussian, and will trigger a new calculatlion when necessary
    procedure, public :: gaussian => random_get
    procedure, public :: array => random_array
    procedure, public :: uniform => random_uniform
    procedure, public :: clipped => random_clipped_normal
    procedure, public :: stats => random_displayStatistics ! subroutine prints summary statistics to command line
    procedure, private :: test => random_test ! subroutine iterates over N calls, then prints summary statistics

  end type LAG_RAND

  class(LAG_RAND), allocatable, public :: random


contains

  subroutine simulation_initialize(self, exp_type)

    use variables, only : zero, sp
    class(LAG_SIM), intent(inout) :: self
    integer, intent(in) :: exp_type

    allocate( self%verticaltox(0:self%nlayers+1), &
            & self%elementArea(0:self%nelements), &
            & self%elementSigmaVolume(0:self%nelements), &
            & self%verticaldiff(0:self%nlayers+1), &
            & self%verticaltemp(0:self%nlayers+1), &
            & self%verticalrho(0:self%nlayers+1))

    self%verticaltox(:) = merge(6324.0_sp, zero, exp_type == 4) ! experiment D or A-C, TODO: move this outside
    self%elementArea = zero
    self%elementSigmaVolume = zero
    self%verticaldiff = zero
    self%verticaltemp = zero
    self%verticalrho = zero

  end subroutine


  subroutine simulation_geometry(self, vertx, verty, node_indices, bathymetry)
    ! called once during simulation setup, calculates area of any triangular mesh or subregion,
    use variables, only : KBM1, sp
    class(LAG_SIM), intent(inout) :: self
    real(sp), dimension(0:self%nnodes), intent(in) :: vertx, verty
    integer, dimension(0:self%nelements, 4), intent(in) :: node_indices
    real(sp), dimension(0:self%nnodes) :: bathymetry
    real(sp), dimension(0:self%nelements) :: length1, length2, length3, ss, average_thickness
    integer :: ii

    do ii = 1, self%nelements
      length1(ii) = sqrt(  (vertx(node_indices(ii, 1)) - vertx(node_indices(ii, 2)) )**(2.0_sp) + ( verty(node_indices(ii, 1)) - verty(node_indices(ii, 2)) )**(2.0_sp)  )
      length2(ii) = sqrt(  (vertx(node_indices(ii, 2)) - vertx(node_indices(ii, 3)) )**(2.0_sp) + ( verty(node_indices(ii, 2)) - verty(node_indices(ii, 3)) )**(2.0_sp)  )
      length3(ii) = sqrt(  (vertx(node_indices(ii, 3)) - vertx(node_indices(ii, 1)) )**(2.0_sp) + ( verty(node_indices(ii, 3)) - verty(node_indices(ii, 1)) )**(2.0_sp)  )
      average_thickness = abs(sum( bathymetry(node_indices(ii,1:3)) )/3.0_sp) ! must be positive
    end do

    ss(:) = 0.5_sp*(length1(:) + length2(:)+ length3(:))
    self%elementArea(:) = sqrt( ss * (ss(:) - length1(:)) * (ss(:) - length1(:)) * (ss(:) - length3(:)) ) ! herons's formula for area
    self%meshArea = sum(self%elementArea(:))
    self%elementSigmaVolume(:) = self%elementArea(:)*average_thickness(:)
    self%layerDepth = average_thickness(1)
    self%layerSigma = float(KBM1)**(-1.0_SP)

  end subroutine


  subroutine simulation_read(self, u_vel, v_vel, w_vel, diffusivity, elevation, salinity, temperature, density)

    use variables, only : ZERO, sp, KB, M, N, iophys

    class(LAG_SIM), intent(inout) :: self
    real(sp), dimension(0:N, KB), intent(inout) :: u_vel, v_vel, w_vel
    real(sp), dimension(0:M, KB), intent(inout) :: diffusivity, salinity, temperature, density
    real(sp), dimension(0:M), intent(inout) :: elevation

    character(len = 100) :: vert_format
    real(sp) :: time
    integer :: ii

    write(vert_format, "(A7,I6,A7)") "(F10.3,", 3*KB, "F20.10)"
    read(iophys, vert_format) time, self%verticaltemp(1:KB), self%verticalrho(1:KB), self%verticaldiff(1:KB)

    u_vel(:,:)=ZERO; v_vel(:,:)=ZERO; w_vel(:, :)=ZERO
    elevation(:)=-abs(ZERO); salinity(:,:)=ZERO

    do ii = 1, self%nnodes
      temperature(ii, 1:KB) = self%verticaltemp(1:KB)
      diffusivity(ii, 1:KB) = self%verticaldiff(1:KB)
      density(ii, 1:KB) = self%verticalrho(1:KB)
    end do

    if (self%lines_read == 0) rewind(unit=iophys)
    self%lines_read = self%lines_read + 1

  end subroutine


  subroutine simulation_diffusion(self)
    ! one dimensional vertical diffusion
    use variables, only : KB, KBM1, dti, ZERO, sp ! time interpolation step

    class(LAG_SIM), intent(inout) :: self ! domain structure
    integer :: ii, jj, nsteps ! iteration parameters
    real(sp) :: dt_substep, dt_stable, secondDerivative, wdiff ! finite element variables
    real(sp), dimension(0:self%nlayers+1) :: TP1; TP1=ZERO ! dissolved profile at time plus one

    wdiff = 60.0_SP*60.0_SP*10.0_SP**(-5.0_SP) ! includes hour and meter conversion from 0.1 cm^2/s
    dt_stable = 0.5_SP*abs(self%layerDepth)**(2.0_SP) ! longest step for conditional stability
    nsteps = ceiling(dti/dt_stable) ! min number of steps to achieve stability, cannot be less than one
    dt_substep = dti/float(nsteps) ! automatically substeps nsteps times, only valid for dK/dt=0

    do ii = 1, nsteps
      self%verticaltox(1) = self%verticaltox(2); self%verticaltox(KB) = self%verticaltox(KBM1)
      do jj = 2, KBM1
        secondDerivative = (self%verticaltox(jj-1) + self%verticaltox(jj+1) - 2.0_SP*self%verticaltox(jj)) / self%layerDepth**(2.0_SP) ! central in space, forward in time second derivative with depth
        TP1(jj) = self%verticaltox(jj) + wdiff*dt_substep*secondDerivative
      end do
      TP1(1) = TP1(2); TP1(KB) = TP1(KBM1) ! copy in domain value to boundary nodes
      self%verticaltox(:) = TP1(:)
    end do

  end subroutine

  subroutine TRIANGLE_GRID_EDGE

    !  Define triangular mesh used for flux computations.

    !     variable list:
    !  vx(m)    :: vx(i) = x-coordinate of node i (input from mesh)
    !  vy(m)    :: vy(i) = y-coordinate of node i (input from mesh)
    !  nv(n,3)  :: nv(i:1-3) = 3 node numbers of element i
    !  xc(n)    :: xc(i) = x-coordinate of element i (calculated from vx)          !
    !  yc(n)    :: yc(i) = y-coordinate of element i (calculated from vy)          !
    !                                                                              !
    !  nbe(n,3) :: nbe(i,1->3) = element index of 1->3 neighbors of element i      !
    !  isbce(n) :: flag if element is on the boundary, see below for values        !
    !  isonb(m) :: flag is node is on the boundary, see below for values           !
    !                                                                              !
    !  ntve(m)  :: the number of neighboring elements of node m                    !
    !  nbve(m,ntve(m)) :: nbve(i,1->ntve(i)) = ntve elements containing node i     !
    !  nbvt(m,ntve(m)) :: nbvt(i,j) = the node number of node i in element         !
    !                     nbve(i,j) (has a value of 1,2,or 3)                      !
    !                                                                              !
    !classification of the triangles nodes, and edges                         !
    !                                                                              !
    !     isonb(i)=0:  node in the interior computational domain                   !
    !     isonb(i)=1:  node on the solid boundary                                  !
    !     isonb(i)=2:  node on the open boundary                                   !
    !                                                                              !
    !     isbce(i)=0:  element in the interior computational domain                !
    !     isbce(i)=1:  element on the solid boundary                               !
    !     isbce(i)=2:  element on the open boundary                                !
    !     isbce(i)=3:  element with 2 solid boundary edges                         !

    use variables
    implicit none

    integer, allocatable, dimension(:, :) :: NB_TMP, CELLS, NBET
    integer, allocatable, dimension(:) :: CELLCNT
    integer :: ii, jj, kk, ll, NTMP, NCNT, NFLAG, JJB, N1, N2, N3, J1, J2, J3
    real(sp) :: X1, X2, X3, Y1, Y2, Y3, DELT, AI1, AI2, AI3, BI1, BI2, BI3, CI1, CI2, CI3, DELTX, DELTY, B1, B2, ART(N)

    ! SET UP MESH (HORIZONTAL COORDINATES)
    ! CALCULATE GLOBAL MINIMUMS AND MAXIMUMS
    VXMIN = MINVAL(VX(1:M)) ; VXMAX = MAXVAL(VX(1:M))
    VYMIN = MINVAL(VY(1:M)) ; VYMAX = MAXVAL(VY(1:M))

    ! SHIFT GRID TO UPPER RIGHT CARTESIAN
    VX = VX - VXMIN
    VY = VY - VYMIN

    ! CALCULATE GLOBAL ELEMENT CENTER GRID COORDINATES
    do ii = 1, N
      XC(ii) = (VX(NV(ii, 1)) + VX(NV(ii, 2)) + VX(NV(ii, 3)))/3.0_SP
      YC(ii) = (VY(NV(ii, 1)) + VY(NV(ii, 2)) + VY(NV(ii, 3)))/3.0_SP
    end do

    XC(0) = zero
    YC(0) = zero
    ART(:)  = zero
    do ii = 1, N
      ART(ii) = (VX(NV(ii, 2)) - VX(NV(ii, 1))) * (VY(NV(ii, 3)) - VY(NV(ii, 1))) - (VX(NV(ii, 3)) - VX(NV(ii, 1))) * (VY(NV(ii, 2)) - VY(NV(ii, 1)))
    end do
    ART = ABS(0.5_SP*ART)

    ! INITIALIZE
    ISBCE = 0
    ISONB = 0
    NBE   = 0

    ! DETERMINE NBE(i=1:n,j=1:3): INDEX OF 1 to 3 NEIGHBORING ELEMENTS
    allocate(NBET(N, 3)) ; NBET = 0
    allocate(CELLS(M, 50)) ; CELLS = 0
    allocate(CELLCNT(M))  ; CELLCNT = 0

    do ii = 1, N
      N1 = NV(ii, 1) ; CELLCNT(N1) = CELLCNT(N1)+1
      N2 = NV(ii, 2) ; CELLCNT(N2) = CELLCNT(N2)+1
      N3 = NV(ii, 3) ; CELLCNT(N3) = CELLCNT(N3)+1
      CELLS(NV(ii, 1), CELLCNT(N1)) = ii
      CELLS(NV(ii, 2), CELLCNT(N2)) = ii
      CELLS(NV(ii, 3), CELLCNT(N3)) = ii
    end do

    if (maxval(cellcnt) > 50) write(*, *) 'bad', maxval(cellcnt)

    DO ii = 1, N
      N1 = NV(ii,1)
      N2 = NV(ii,2)
      N3 = NV(ii,3)
      DO J1 = 1, CELLCNT(N1)
        DO J2 = 1, CELLCNT(N2)
          IF ((CELLS(N1, J1) == CELLS(N2, J2)) .AND. CELLS(N1,J1) /= ii) NBE(ii, 3) = CELLS(N1, J1)
        END DO
      END DO
      DO J2 = 1, CELLCNT(N2)
        DO J3 = 1, CELLCNT(N3)
          IF ((CELLS(N2, J2) == CELLS(N3, J3)) .AND. CELLS(N2, J2) /= ii) NBE(ii,1) = CELLS(N2, J2)
        END DO
      END DO
      DO J1 = 1, CELLCNT(N1)
        DO J3 = 1, CELLCNT(N3)
          IF((CELLS(N1, J1) == CELLS(N3,J3)) .AND. CELLS(N1,J1) /= ii) NBE(ii, 2) = CELLS(N3, J3)
        END DO
      END DO
    END DO
    DEALLOCATE(CELLS,CELLCNT)

    !   IF(MSR)WRITE(IPT,*)  '!  NEIGHBOR FINDING      :    COMPLETE'
    !
    ! Ensure all elements have at least one neighbor
    NFLAG = 0
    do ii = 1, N
      if (sum(NBE(ii, 1:3)) == 0) then
        NFLAG = 1
        write(*, *) 'ELEMENT ', ii, ' AT ', XC(ii), YC(ii), ' HAS NO NEIGHBORS'
        stop
      end if
    end do
    if (NFLAG == 1) stop


    ! if element on boundary set isbce(i)=1 and isonb(j)=1 for boundary nodes j
    do ii = 1, N
      if ( MIN(NBE(ii, 1), NBE(ii, 2), NBE(ii, 3)) == 0 ) then
        ISBCE(ii) = 1  ! element on boundary
        if (NBE(ii, 1) == 0) then
          ISONB(NV(ii, 2)) = 1 ; ISONB(NV(ii, 3)) = 1
        end if
        if (NBE(ii,2) == 0) then
          ISONB(NV(ii, 1)) = 1 ; ISONB(NV(ii, 3)) = 1
        end if
        if (NBE(ii,3) == 0) then
          ISONB(NV(ii, 1)) = 1 ; ISONB(NV(ii, 2)) = 1
        end if
      end if
    end do

    ! DEFINE NTVE, NBVE, NBVT
    ! ntve(1:m): total number of the surrounding triangles connected to the given node
    ! nbve(1:m, 1:ntve+1): the identification number of surrounding triangles with a common node (counted clockwise)
    ! nbvt(1:m,ntve(1:m)): the idenfication number of a given node over each individual surrounding triangle(counted clockwise)                                              !

    ! Determine max number of surrounding elements
    MX_NBR_ELEM = 0
    element_loop: do ii = 1, M
      NCNT = 0
      node_loop: do jj = 1, N
        if ( float(NV(jj, 1) - ii) * float(NV(jj, 2) - ii) * float(NV(jj, 3) - ii) == 0.0_SP ) NCNT = NCNT + 1
      end do node_loop
      MX_NBR_ELEM = MAX(MX_NBR_ELEM, NCNT)
    end do element_loop


    ! allocate arrays based on mx_nbr_elem
    ALLOCATE( NBVE(M, MX_NBR_ELEM + 1))
    ALLOCATE( NBVT(M, MX_NBR_ELEM + 1))


    ! Determine number of surrounding elements for node i = ntve(i)
    ! determine nbve - indices of neighboring elements of node i
    ! determine nbvt - index (1,2, or 3) of node i in neighboring element

    do ii = 1, M
      NCNT = 0
      do jj = 1, N
        if ( float(NV(jj, 1) - ii) * float(NV(jj, 2) - ii) * float(NV(jj, 3) - ii) == 0.0_SP) then
          NCNT = NCNT+1
          NBVE(ii, NCNT) = jj
          if ((NV(jj, 1) - ii) == 0) NBVT(ii, NCNT) = 1
          if ((NV(jj, 2) - ii) == 0) NBVT(ii, NCNT) = 2
          if ((NV(jj, 3) - ii) == 0) NBVT(ii, NCNT) = 3
        end if
      end do
      NTVE(ii) = NCNT
    end do

    !
    !--Reorder Order Elements Surrounding a Node to Go in a Cyclical Procession----!
    !--Determine NTSN  = Number of Nodes Surrounding a Node (+1)-------------------!
    !--Determine NBSN  = Node Numbers of Nodes Surrounding a Node------------------!

    allocate(NB_TMP(M,MX_NBR_ELEM+1))
    do ii = 1, M
      if (ISONB(ii) == 0) then
        NB_TMP(1, 1) = NBVE(ii, 1)
        NB_TMP(1, 2) = NBVT(ii, 1)
        do jj = 2, NTVE(ii) + 1
          kk = NB_TMP(jj - 1, 1)
          ll = NB_TMP(jj - 1, 2)
          NB_TMP(jj, 1) = NBE(kk, ll + 1 - INT((ll + 1)/4)*3)
          ll=NB_TMP(jj, 1)
          IF ((NV(ll, 1) - ii) == 0) NB_TMP(jj, 2) = 1
          IF ((NV(ll, 2) - ii) == 0) NB_TMP(jj, 2) = 2
          IF ((NV(ll, 3) - ii) == 0) NB_TMP(jj, 2) = 3
        end do

        do jj = 2, NTVE(ii) + 1
          NBVE(ii, jj) = NB_TMP(jj, 1)
        end do

        do jj = 2, NTVE(ii) + 1
          NBVT(ii, jj) = NB_TMP(jj, 2)
        end do

        NTMP = NTVE(ii) + 1
        if (NBVE(ii, 1) /= NBVE(ii, NTMP)) then
          print*, ii,'nbve(ii) not correct!!'
          stop
        end if
        if (NBVT(ii,1) /= NBVT(ii, NTMP)) then
          print*, ii,'NBVT(ii) NOT CORRECT!!'
          stop
        end if

      else
        JJB = 0

        do jj = 1, NTVE(ii)
          ll = NBVT(ii, jj)
          if (NBE(NBVE(ii, jj), ll + 2 - int((ll + 2)/4)*3) == 0) then
            JJB = JJB + 1
            NB_TMP(JJB, 1) = NBVE(ii, jj)
            NB_TMP(JJB, 2) = NBVT(ii, jj)
          end if
        end do

        if (JJB /= 1) then
          print*, 'ERROR IN ISONB !, I, J', ii, jj
          stop
        end if

        do jj = 2, NTVE(ii)
          kk = NB_TMP(jj - 1, 1)
          ll = NB_TMP(jj - 1, 2)
          NB_TMP(jj, 1) = NBE(kk, ll + 1 - int((ll + 1)/4)*3)
          ll = NB_TMP(jj, 1)
          if ((NV(ll, 1) - ii) == 0) NB_TMP(jj, 2) = 1
          if ((NV(ll, 2) - ii) == 0) NB_TMP(jj, 2) = 2
          if ((NV(ll, 3) - ii) == 0) NB_TMP(jj, 2) = 3
        end do

        do jj = 1, NTVE(ii)
          NBVE(ii, jj) = NB_TMP(jj, 1)
          NBVT(ii, jj) = NB_TMP(jj, 2)
        enddo

        NBVE(ii, NTVE(ii) + 1) = 0
      end if
    end do
    deallocate(NB_TMP)

    !--Huang change 911
    ! This is a special part. Since in original FVCOM subroutine tge.F
    !       ISONB(I_OBC_N(I)) = 2
    ! I_OBC_N is read from input file xxx_obc.dat
    ! here Martin deleted that part
    ! This line of code needs to be changed for each new case
    !   do i=1,135
    !      ISONB(i) = 2
    !   enddo


    do ii = 1, N
      if (ISBCE(ii) == 0) then
        Y1 = YC(NBE(ii, 1)) - YC(ii)
        Y2 = YC(NBE(ii, 2)) - YC(ii)
        Y3 = YC(NBE(ii, 3)) - YC(ii)
        X1 = XC(NBE(ii, 1)) - XC(ii)
        X2 = XC(NBE(ii, 2)) - XC(ii)
        X3 = XC(NBE(ii, 3)) - XC(ii)
        X1 = X1/1000.0_SP
        X2 = X2/1000.0_SP
        X3 = X3/1000.0_SP
        Y1 = Y1/1000.0_SP
        Y2 = Y2/1000.0_SP
        Y3 = Y3/1000.0_SP

        delt=(x1*y2-x2*y1)**2+(x1*y3-x3*y1)**2+(x2*y3-x3*y2)**2
        delt=delt*1000.0

        a1u(ii, 1) = (y1+y2+y3)*(x1*y1+x2*y2+x3*y3)- (x1+x2+x3)*(y1**2+y2**2+y3**2)
        a1u(ii, 1) = a1u(ii, 1)/delt
        a1u(ii, 2) = (y1**2+y2**2+y3**2)*x1 - (x1*y1+x2*y2+x3*y3)*y1
        a1u(ii, 2) = a1u(ii, 2)/delt
        a1u(ii, 3) = (y1**2+y2**2+y3**2)*x2 - (x1*y1+x2*y2+x3*y3)*y2
        a1u(ii, 3) = a1u(ii, 3)/delt
        a1u(ii, 4) = (y1**2+y2**2+y3**2)*x3 - (x1*y1+x2*y2+x3*y3)*y3
        a1u(ii, 4) = a1u(ii, 4)/delt

        a2u(ii, 1) = (x1+x2+x3)*(x1*y1+x2*y2+x3*y3) - (y1+y2+y3)*(x1**2+x2**2+x3**2)
        a2u(ii, 1) = a2u(ii, 1)/delt
        a2u(ii, 2) = (x1**2+x2**2+x3**2)*y1-(x1*y1+x2*y2+x3*y3)*x1
        a2u(ii, 2) = a2u(ii, 2)/delt
        a2u(ii, 3) = (x1**2+x2**2+x3**2)*y2-(x1*y1+x2*y2+x3*y3)*x2
        a2u(ii, 3) = a2u(ii, 3)/delt
        a2u(ii, 4) = (x1**2+x2**2+x3**2)*y3-(x1*y1+x2*y2+x3*y3)*x3
        a2u(ii, 4) = a2u(ii, 4)/delt
      end if

      x1 = vx(nv(ii, 1)) - xc(ii)
      x2 = vx(nv(ii, 2)) - xc(ii)
      x3 = vx(nv(ii, 3)) - xc(ii)
      y1 = vy(nv(ii, 1)) - yc(ii)
      y2 = vy(nv(ii, 2)) - yc(ii)
      y3 = vy(nv(ii, 3)) - yc(ii)


      ai1 = y2 - y3
      ai2 = y3 - y1
      ai3 = y1 - y2
      bi1 = x3 - x2
      bi2 = x1 - x3
      bi3 = x2 - x1
      ci1 = x2*y3 - x3*y2
      ci2 = x3*y1 - x1*y3
      ci3 = x1*y2 - x2*y1

      aw0(ii, 1) = -ci1 / 2.0 / art(ii)
      aw0(ii, 2) = -ci2 / 2.0 / art(ii)
      aw0(ii, 3) = -ci3 / 2.0 / art(ii)
      awx(ii, 1) = -ai1 /2.0 / art(ii)
      awx(ii, 2) = -ai2 / 2.0 / art(ii)
      awx(ii, 3) = -ai3 / 2.0 / art(ii)
      awy(ii, 1) = -bi1 / 2.0 / art(ii)
      awy(ii, 2) = -bi2 / 2.0 / art(ii)
      awy(ii, 3) = -bi3 / 2.0 / art(ii)

    end do

    do ii = 1, n
      if (isbce(ii) > 1) then
        do jj = 1, 4
          a1u(ii, jj) = 0.0_SP
          a2u(ii, jj) = 0.0_SP
        end do
      else if (isbce(ii) == 1) then
        do jj = 1, 3
          if (nbe(ii, jj) == 0) ll = jj
        end do
        j1 = ll + 1 - int((ll + 1)/4)*3
        j2 = ll + 2 - int((ll + 2)/4)*3
        x1 = vx(nv(ii, j1)) - xc(ii)
        x2 = vx(nv(ii, j2)) - xc(ii)
        y1 = vy(nv(ii, j1)) - yc(ii)
        y2 = vy(nv(ii, j2)) - yc(ii)

        delt = x1*y2 - x2*y1
        b1 = (y2 - y1)/delt
        b2 = (x1 - x2)/delt
        deltx = vx(nv(ii, j1)) - vx(nv(ii, j2))
        delty = vy(nv(ii, j1)) - vy(nv(ii, j2))


        a1u(ii, 1) = 0.0_SP
        a1u(ii, ll + 1) = 0.0_SP
        a1u(ii, j1 + 1) = 0.0_SP
        a1u(ii, j2 + 1) = 0.0_SP

        a2u(ii, 1) = 0.0_SP
        a2u(ii, ll + 1) = 0.0_SP
        a2u(ii, j1 + 1) = 0.0_SP
        a2u(ii, j2 + 1) = 0.0_SP
      end if
    end do

  end subroutine

  subroutine hunt(sigma_nodes, KB, sigma_particle, jlo) ! Z, KB, self%ZP(ii), NZR
    ! from numerical recipies vol 2
    use variables

    integer, intent(inout) :: jlo ! sigma layer below particle?
    integer, intent(in) :: KB ! number of sigma layers
    real(sp), dimension(KB), intent(in) :: sigma_nodes ! sigma coordinate value
    real(sp), intent(in) :: sigma_particle ! depth of particle

    integer :: inc, jhi, jm
    logical :: ascnd

    ascnd = (sigma_nodes(KB) > sigma_nodes(1)) ! bottom sigma layer greater than first layer
    if ((jlo <= 0) .or. (jlo > KB)) then
      jlo = 0
      jhi = KB + 1
      goto 3
    endif
    inc = 1

    if ((sigma_particle >= sigma_nodes(jlo)) .eqv. ascnd) then
      1    jhi = jlo + inc
      if(jhi > KB)then
        jhi=n+1
      else if ((sigma_particle >= sigma_nodes(jhi)) .eqv. ascnd) then
        jlo = jhi
        inc = inc + inc
        goto 1
      end if
    else
      jhi = jlo
      2    jlo = jhi - inc
      if (jlo < 1) then
        jlo=0
      else if ((sigma_particle < sigma_nodes(jlo)) .eqv. ascnd) then
        jhi = jlo
        inc = inc + inc
        go to 2
      end if
    end if
    3 if (jhi-jlo == 1) then
      if (sigma_particle == sigma_nodes(KB)) jlo = KB - 1
      if (sigma_particle == sigma_nodes(1)) jlo = 1
      return
    end if
    jm = (jhi + jlo)/2
    if (sigma_particle >= sigma_nodes(jm) .eqv. ascnd) then
      jlo = jm
    else
      jhi = jm
    end if
    go to 3

  end subroutine hunt


  subroutine spline(x, y, n2, yp1, ypn, y2)
    ! from numerical recipies vol 2, but modfied so that nmax=50
    use variables, only : sp
    implicit none
    integer  :: n2
    real(sp), intent(in) :: x(n2), y(n2), yp1, ypn
    real(sp), intent(out), dimension(n2) :: y2

    integer  :: i, k
    integer, parameter :: nmax=50
    real(sp) :: p, qn, sig, un
    real(sp), dimension(nmax) :: u

    if (yp1 > 0.99e30) then ! force natural lower boundary
      y2(1) = 0.0_SP
      u(1) = 0.0_SP
    else ! or set specific values
      y2(1) = -0.5_SP
      u(1) = (3.0_SP/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    end if
    do i = 2, n2-1 ! tridiagonal algorithm decomp
      sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
      p = sig * y2(i-1) + 2.0_SP
      y2(i) = (sig-1.) / p
      u(i) = (6.0*((y(i+1) - y(i)) / (x(i+1) - x(i)) - (y(i) - y(i-1))/(x(i) - x(i-1)))/(x(i+1) - x(i-1)) - (- sig*u(i-1)))/p
    end do
    if (ypn > 0.99e30) then ! force natural upper boundary
      qn = 0.0_SP
      un = 0.0_SP
    else
      qn = 0.5_SP
      un = (3.0_SP/(x(n2)-x(n2-1)))*(ypn-(y(n2)-y(n2-1))/(x(n2)-x(n2-1)))
    end if
    y2(n2) = (un - qn * u(n2-1)) / (qn * y2(n2-1) + 1.0_SP)
    do k = n2-1, 1, -1
      y2(k) = y2(k) * y2(k+1) + u(k)
    end do

  end subroutine spline


  subroutine random_initialize(self)

    ! init random seed from computer clock
    class(LAG_RAND), intent(inout) :: self
    integer :: uu, seed_size
    integer(kind = 4) :: clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = seed_size)
    allocate(seed(seed_size))
    call system_clock(clock)
    seed = clock + 37 * (/ (uu - 1, uu = 1, seed_size) /)
    call random_seed(put = seed)
    deallocate(seed)

    self%current = .true.
    self%statistics = .true.
    self%samples = 0
    self%rn1 = 0.0_SP
    self%rn2 = 0.0_SP
    self%ru1 = 0.0_SP
    self%ru2 = 0.0_SP
    self%mean = 0.0_SP
    self%sumMeanDiffSq = 0.0_SP

    call self%normal()

  end subroutine


  subroutine random_normal(self) ! OK
    ! generate two random normal numbers using Box-Muller method
    use variables, only : sp
    class(LAG_RAND), intent(inout) :: self
    real(sp) :: V1, V2, S, meanDiff

    do
      V1 = self%uniform() ! uniform random centered on zero
      V2 = self%uniform()
      S = V1**2.0_SP + V2**2.0_SP ! check the sum of the squares and reject if outside range
      if (S <= 1.0_SP) exit
    end do

    self%rn1 = sqrt(-2.0_SP*log(S)/S)*V1 ! generate normal rands from uniform
    self%rn2 = sqrt(-2.0_SP*log(S)/S)*V2

    if (self%statistics) then ! calculate mean and stddev on the fly if desired
      self%samples = self%samples + 1
      meanDiff = self%rn1 - self%mean
      self%mean = self%mean + meanDiff/self%samples
      self%sumMeanDiffSq = self%sumMeanDiffSq + meanDiff*(self%rn1 - self%mean)

      self%samples = self%samples + 1
      meanDiff = self%rn2 - self%mean
      self%mean = self%mean + meanDiff/self%samples
      self%sumMeanDiffSq = self%sumMeanDiffSq + meanDiff*(self%rn2 - self%mean)
    end if

  end subroutine


  function random_array(self, nn)

    use variables, only : sp
    class(LAG_RAND), intent(inout) :: self
    integer, intent(in) :: nn

    real(sp), dimension(nn) :: random_array
    real(sp), dimension(nn) :: temp_array
    integer :: ii

    do ii = 1, nn
      temp_array(ii) = self%get()
    end do

    random_array(:) = temp_array(:)

  end function


  real(sp) function random_uniform(self) result(rand)

    use variables, only : sp
    class(LAG_RAND), intent(inout) :: self

    call random_number(self%ru1)

    rand = 2.0_SP*self%ru1 - 1.0_SP ! returns uniform pseudorandom in -1 to 1 range

  end function


  real(sp) function random_get(self) result(rand)! OK
    ! returns one of stored gaussian random numbers and generates new ones when used
    use variables, only : sp
    class(LAG_RAND), intent(inout) :: self

    if (self%current) then ! if stored randoms haven't been used,
      rand = self%rn1 ! return the first
      self%current = .false.
    else ! if one stored random has been used,
      rand = self%rn2 ! return the second
      self%current = .true.
      call self%normal() ! generate new numbers
    end if

  end function


  real(sp) function random_clipped_normal(self) result(rand)! OK
    ! returns one of stored gaussian random numbers and generates new ones when used
    use variables, only : sp
    class(LAG_RAND), intent(inout) :: self
    real(sp) :: deviate

    deviate = self%get() ! get new random normal
    do while (abs(deviate) > 1.0_SP) ! check if -1,1
      deviate = self%get() ! reassign if out of range
    end do
    rand = deviate

  end function random_clipped_normal


  subroutine random_displayStatistics(self) ! OK
    ! calculate and display distribution statistics for the random number system
    use variables ! for real precision
    implicit none

    class(LAG_RAND), intent(in) :: self
    write(*, *); write(*, *) "Statistics for random gaussian numbers generated so far..."; write(*, *)
    write(*, *) "    Samples:       ", self%samples
    write(*, *) "    Mean:          ", self%mean
    write(*, *) "    Variance:      ", self%sumMeanDiffSq / float(self%samples - 1)
    write(*, *) "    Std Deviation: ", sqrt(self%sumMeanDiffSq / float(self%samples - 1))

  end subroutine random_displayStatistics


  subroutine random_test(self)
    ! iteratively generate random gaussian numbers and calculate statistics
    use variables, only : sp

    class(LAG_RAND), intent(inout) :: self
    integer :: ii, nSample = 1000
    real(sp) :: randomNumber

    do ii = 1, nSample-1
      randomNumber = self%get() ! return next gaussian random
    end do
    call self%stats() ! display distribution statistics

  end subroutine



  subroutine getInteger(file, key, target, fid)

    character(len = *), intent(in) :: key, file
    integer, intent(out) :: target
    integer, intent(in), optional :: fid
    integer :: errorCode

    errorCode = find_key(file, key, iscal=target)
    if (errorCode /= 0) then
      write(fid, *) 'Error reading '//key//': ', errorCode
      stop
    end if

  end subroutine


  subroutine getString(file, key, target, fid)

    character(len = *), intent(in) :: key, file
    character(len = 80), intent(out) :: target
    integer, intent(in) :: fid
    integer :: errorCode, ii

    errorCode = find_key(file, key, cval = target)
    if (errorCode /= 0) then
      write(fid, *) 'Error reading '//key//': ', errorCode
      stop
    else
      ! remove trailing directory '/'
      ii = len_trim(target)
      if (target(ii:ii) == "/") target(ii:ii) = " "
    end if

  end subroutine


  integer function find_key(FNAME, VNAME, ISCAL, FSCAL, IVEC, FVEC, CVEC, NSZE, CVAL, LVAL)
    !   Scan an Input File for a Variable
    !   RETURN VALUE:
    !        0 = FILE FOUND, VARIABLE VALUE FOUND
    !       -1 = FILE DOES NOT EXIST OR PERMISSIONS ARE INCORRECT
    !       -2 = VARIABLE NOT FOUND OR IMPROPERLY SET
    !       -3 = VARIABLE IS OF DIFFERENT TYPE, CHECK INPUT FILE
    !       -4 = VECTOR PROVIDED BUT DATA IS SCALAR TYPE
    !       -5 = NO DATATYPE DESIRED, EXITING

    !   REQUIRED INPUT:
    !        FNAME = File Name
    !        FSIZE = Length of Filename

    !   optional (MUST PROVIDE ONE)
    !        ISCAL = integer SCALAR
    !        FSCAL = FLOAT SCALAR
    !        CVAL = character VARIABLE
    !        LVAL = LOGICAL VARIABLE
    !        IVEC = integer VECTOR **
    !        FVEC = FLOAT VECTOR **
    !        CVEC = STRING VECTOR **
    !      **NSZE = ARRAY SIZE (MUST BE PROVIDED WITH IVEC/FVEC)

    use variables
    implicit none
    character(LEN = *) :: FNAME, VNAME
    integer, intent(inout), optional :: ISCAL, IVEC(*)
    REAL(SP), intent(inout), optional :: FSCAL, FVEC(*)
    character(LEN=80), optional :: CVAL, CVEC(*)
    LOGICAL, intent(inout), optional :: LVAL
    integer, intent(inout), optional :: NSZE

    REAL(SP) REALVAL(150)
    integer  INTVAL(150)
    character(LEN=20 ) :: key
    character(LEN=80 ) :: STRINGVAL(150),TITLE
    character(LEN=80 ) :: line
    character(LEN=400) :: buffer
    character(LEN=7  ) :: type_of
    character(LEN=20 ), DIMENSION(200) :: SET
    integer :: last, NVAL, lines, NREP
    logical :: SETYES, ALLSET, CHECK, LOGVAL
    character(len=*), parameter :: continue_line = "////"
    character(len = len_trim(copy)) :: text
    character(len = len_trim(copy)) :: value, TEMP, fragments(200)
    character(len = 80) :: TSTRING
    character(len = 6) :: ERRSTRING
    character(len = 16) :: NUMCHARS = "0123456789+-Ee. "
    integer :: EQLOC, length, ii, LOCEX, NP
    logical :: flag

    find_key = 0

    ! OPEN THE INPUT FILE
    inquire(file=TRIM(FNAME), exist=CHECK)
    if (.not. CHECK) then
      find_key = -1
    end if

    open(10, file=trim(FNAME))
    rewind(10)

    lines = 0
    do while (.true.)

      buffer(1:len(buffer)) = ' '
      NREP  = 0
      lines = lines + 1
      read(10,'(a)', end=20) line
      buffer(1:80) = line(1:80)

      ! PROCESS LINE CONTINUATIONS
      110 CONTINUE
      last = len_trim(line)
      if (last /= 0) then
        if ( line(last-1:last) == '\\\\') then

          NREP = NREP + 1
          read(10, '(a)', end=20) line
          lines = lines + 1
          buffer( NREP*80 + 1 : NREP*80 +80) = line(1:80)
          GOTO 110
        end if
      end if

      ! REMOVE LINE CONTINUATION character \\
      if (NREP > 0) then
        do last = 2, LEN_TRIM(buffer)
          if ( buffer(last-1:last) == '\\\\') buffer(last-1:last) = '  '
        end do
      end if

      fragments = " "
      type_of = "error"
      LOGVAL = .false.

      write(ERRSTRING, "(I6)") lines

      LOCEX = index(text, "!")
      if (LOCEX /= 0) text = text(1:LOCEX-1)
      length = len_trim(text)

      if (length == 0) then
        type_of = "none"
        key = "none"
        return
      end if

      ! Commas to spaces
      where (text(:) == ",")
        text(:) = " "
      end where

      ! Find assignment "="
      EQLOC = index(text, "=")
      if (EQLOC == 0) call raise(6,'DATA LINE '//ERRSTRING//' MUST CONTAIN "=" ')

      ! split name and value substrings
      key = text(1:EQLOC-1)
      value  = adjustl(text(EQLOC+1:LENGTH))
      length = len_trim(value)

      if (length == 0) call raise(6,'IN DATA PARAMETER FILE', 'VARIABLE LINE'//ERRSTRING//' HAS NO ASSOCIATED VALUE')

      ! check for logical
      if ((value(1:1) == "T" .or. value(1:1) == "F") .and. length == 1) then
        type_of = "logical"
        if (value(1:1) == "T") LOGVAL = .true.
        return
      end if

      ! is string if contains non-numeric characters
      do ii = 1, length
        if (index(NUMCHARS, value(ii:ii)) == 0) then

          type_of = "string"
          TSTRING = value
          stringval(1) = TSTRING
          NVAL = 1
          flag = .true.

          do ii = 1, length
            if (value(ii:ii) /= " ") then
              fragments(NVAL) = trim(fragments(NVAL)) // value(ii:ii)
              flag = .true.
            else
              if (flag) NVAL = NVAL + 1
              flag = .false.
            end if
          end do

          do ii = 1, NVAL
            stringval(ii + 1) = trim(fragments(ii))
          end do
          return

        end if
      end do

      type_of = merge("float  ", "integer", index(value, ".") /= 0)

      ! Split lines
      NP = 1
      flag = .true.
      do ii = 1, length
        if (value(ii:ii) /= " ") then
          fragments(NP) = trim(fragments(NP)) // value(ii:ii)
          flag = .true.
        else
          if (flag) NP = NP + 1
          flag = .false.
        end if
      end do

      ! numerical
      NVAL = NP
      do ii = 1, NP
        if (type_of == "float") then
          read(trim(fragments(ii)), *) realval(ii)
        else
          read(trim(fragments(ii)), *) intval(ii)
        end if
      end do


      if (trim(name) == trim(VNAME)) then

        if (PRESENT(ISCAL)) then
          if (type_of == 'integer') then
            ISCAL = INTVAL(1)
            return
          else
            find_key = -3
          end if
        elseif(present(FSCAL)) then
          if (type_of == 'float') then
            FSCAL = REALVAL(1)
            return
          else
            find_key = -3
          end if
        elseif(present(CVAL))THEN
          if (type_of == 'string') then
            CVAL = STRINGVAL(1)
            return
          else
            find_key = -3
          end if
        elseif (present(LVAL)) THEN
          if (type_of == 'logical') then
            LVAL = LOGVAL
            return
          else
            find_key = -3
          end if
        else if (present(IVEC)) then
          if (NVAL > 1) then
            if (type_of == 'integer') then
              IVEC(1:NVAL) = INTVAL(1:NVAL)
              NSZE = NVAL
              return
            else
              find_key = -3
            end if
          else
            find_key = -4
          end if
        elseif (present(FVEC)) then
          if (NVAL > 1) then
            IF (type_of == 'float') then
              FVEC(1:NVAL) = REALVAL(1:NVAL)
              NSZE = NVAL
              return
            else
              find_key = -3
            end if
          else
            find_key = -4
          end if
        elseif (present(CVEC)) then
          if (NVAL > 0) then
            if (type_of == 'string') then
              CVEC(1:NVAL) = STRINGVAL(2:NVAL+1)
              NSZE = NVAL
              return
            else
              find_key = -3
            end if
          else
            find_key = -4
          end if
        else
          find_key = -5
        end if
      end if
    end do
    20 close(10)
    find_key = -2
  end function

end module
