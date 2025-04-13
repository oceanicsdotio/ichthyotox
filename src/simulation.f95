module simulation
    use variables, only : ZERO, sp
    use random, only : random_number_generator
    
    implicit none
    save
    private
  
    type, public :: Mesh
    contains
        procedure, public :: init => TRIANGLE_GRID_EDGE
    end type
  
    type, public :: Experiment
        character(len = 100), public :: simID
        integer, public :: &
            & nnodes=0, &
            & nelements=0, &
            & nlayers=0, &
            & lines_read=0
        real(sp), public :: &
            & globalIrradiance=ZERO, &
            & meshArea=ZERO, &
            & layerDepth=ZERO, &
            & layerSigma=ZERO, &
            & time=ZERO, & ! total time
            & daytime=ZERO, & ! time in days
            & clocktime=ZERO ! twenty four hour periodic
        real(sp), allocatable, dimension(:), private :: elementSigmaVolume, elementArea ! mesh stats
        real(sp), allocatable, dimension(:), public :: verticaltox, verticaldiff, verticaltemp, verticalrho ! uniform horizontal fields
    contains
        ! call in this order
        procedure, public :: init => initSimulation
        procedure, public :: read => readSimulation
        procedure, public :: diffuse => verticalDiffusion
    end type
  
    type(Experiment), allocatable, public :: domain ! domain structure imported from this module

contains
    subroutine initSimulation(self, toxin)
        ! allocate and initialize data structures
        class(Experiment), intent(inout) :: self
        real(sp), intent(in) :: toxin ! 6324 for D

        allocate( &
            & self%verticaltox(0:self%nlayers+1), &
            & self%elementArea(0:self%nelements), &
            & self%elementSigmaVolume(0:self%nelements), &
            & self%verticaldiff(0:self%nlayers+1), &
            & self%verticaltemp(0:self%nlayers+1), &
            & self%verticalrho(0:self%nlayers+1))

        self%verticaltox(:) = toxin
        self%elementArea = zero
        self%elementSigmaVolume = zero
        self%verticaldiff = zero
        self%verticaltemp = zero
        self%verticalrho = zero
    end subroutine
  
    subroutine readSimulation(self, u_vel, v_vel, w_vel, diffusivity, elevation, salinity, temperature, density)
        ! read physical drivers from file
        use variables, only : KB, M, N, iophys
    
        class(Experiment), intent(inout) :: self
        real(sp), dimension(0:N, KB), intent(inout) :: u_vel, v_vel, w_vel
        real(sp), dimension(0:M, KB), intent(inout) :: diffusivity, salinity, temperature, density
        real(sp), dimension(0:M), intent(inout) :: elevation
    
        character(len = 100) :: vert_format
        real(sp) :: time
        integer :: ii
    
        write(vert_format, "(A7,I6,A7)") "(F10.3,", 3*KB, "F20.10)"
        read(iophys, vert_format) time, self%verticaltemp(1:KB), self%verticalrho(1:KB), self%verticaldiff(1:KB)
    
        u_vel(:, :) = zero
        v_vel(:, :) = zero 
        w_vel(:, :) = zero
        elevation(:) = -abs(zero)
        salinity(:,:) = zero
        do ii = 1, self%nnodes
            temperature(ii, 1:KB) = self%verticaltemp(1:KB)
            diffusivity(ii, 1:KB) = self%verticaldiff(1:KB)
            density(ii, 1:KB) = self%verticalrho(1:KB)
        end do
        if (self%lines_read == 0) rewind(unit=iophys)
        self%lines_read = self%lines_read + 1
    end subroutine

    subroutine verticalDiffusion(self)
        ! One-dimensional vertical diffusion
        use variables, only : KB, KBM1, dti
        class(Experiment), intent(inout) :: self
        integer :: ii, jj, steps_to_stability
        real(sp) :: &
            & substep, & ! automatic time substep, only valid for dK/dt=0
            & stable, &  ! longest time step for conditional stability
            & diffusivity = 60.0*60.0*10.0**(-5.0), &  ! 0.1 cm^2/s
            & profile(0:self%nlayers+1); 
        
        profile = zero
        stable = 0.5 * self%layerDepth**2.0
        steps_to_stability = ceiling(dti/stable) ! min number of steps to achieve stability, cannot be less than one
        substep = dti / float(steps_to_stability)
    
        do ii = 1, steps_to_stability
            self%verticaltox(1) = self%verticaltox(2)
            self%verticaltox(KB) = self%verticaltox(KBM1)
    
            ! central-in-space (z), forward-in-time, second derivative
            do jj = 2, KBM1
            profile(jj) = self%verticaltox(jj) + diffusivity * substep * &
                    & (self%verticaltox(jj-1) + self%verticaltox(jj+1) - 2.0*self%verticaltox(jj)) / &
                    & self%layerDepth * self%layerDepth
            end do
    
            profile(1) = profile(2);
            profile(KB) = profile(KBM1) ! copy in domain value to boundary nodes
            self%verticaltox(:) = profile(:)
        end do
    end subroutine

    subroutine TRIANGLE_GRID_EDGE(self)
        !  Define triangular mesh used for flux computations.
    
        !     variable list:
        !  vx(m)    :: vx(i) = x-coordinate of node i (input from mesh)
        !  vy(m)    :: vy(i) = y-coordinate of node i (input from mesh)
        !  nv(n,3)  :: nv(i:1-3) = 3 node numbers of element i
        !  xc(n)    :: xc(i) = x-coordinate of element i (calculated from vx)          !
        !  yc(n)    :: yc(i) = y-coordinate of element i (calculated from vy)          !
        !                                                                              !
        !  nbe(n,3) :: nbe(i,1:3) = element index of 1->3 neighbors of element i      !
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
        !     isonb
        !       0: interior computational domain                   !
        !       1: solid boundary                                  !
        !       2: open boundary                                   !
        !                                                                              !
        !     isbce
        !       0:  interior computational domain                !
        !       1:  1 solid boundary                               !
        !       2:  open boundary                                !
        !       3:  2 solid boundary edges                         !
    
        use variables, only : xc, yc, vx, vy, nv, M, isbce, nbe, isonb, awx, awy
        use variables, only : n, a1u, NBVE, NBVT, vxmin, vxmax, vymin, VYMAX, aw0, a2u, MX_NBR_ELEM, NTVE
        class(Mesh), intent(inout) :: self
        integer, allocatable, dimension(:, :) :: NB_TMP, CELLS, NBET
        integer, allocatable, dimension(:) :: CELLCNT
        integer :: ii, jj, kk, ll, NTMP, NCNT, JJB, N1, N2, N3, J1, J2, J3, tri(3)
        real(sp) :: X1, X2, X3, Y1, Y2, Y3, DELT, AI(3), BI(3), CI(3), DELTX, DELTY, B1, B2, ART(N)
    
        ! SET UP MESH (HORIZONTAL COORDINATES)
        ! CALCULATE GLOBAL MINIMUMS AND MAXIMUMS
        VXMIN = MINVAL(VX(1:M))
        VXMAX = MAXVAL(VX(1:M))
        VYMIN = MINVAL(VY(1:M))
        VYMAX = MAXVAL(VY(1:M))
    
        ! CALCULATE GLOBAL ELEMENT CENTER GRID COORDINATES
        xc = zero
        yc = zero
        do ii = 1, N
            xc(ii) = sum(VX(NV(ii, :))) / 3.0_SP
            yc(ii) = sum(VY(NV(ii, :))) / 3.0_SP
        end do
    
        ART(:) = abs(((VX(NV(:, 2)) - VX(NV(:, 1))) * (VY(NV(:, 3)) - VY(NV(:, 1))) - &
                & (VX(NV(:, 3)) - VX(NV(:, 1))) * (VY(NV(:, 2)) - VY(NV(:, 1))))*0.5)
    
        ISBCE = 0
        ISONB = 0
        NBE   = 0
    
        allocate( &
            &NBET(N, 3), & ! index of neighbors
            & CELLS(M, 50), &
            & CELLCNT(M))
    
        NBET = 0
        cellcnt = 0
        cells = 0
    
        do ii = 1, n
            tri = nv(ii, 1:3)
            cellcnt(tri) = cellcnt(tri) + 1
            cells(tri, cellcnt(tri)) = ii
        end do
    
        do ii = 1, N
            N1 = NV(ii,1)
            N2 = NV(ii,2)
            N3 = NV(ii,3)
            do J1 = 1, CELLCNT(N1)
            do J2 = 1, CELLCNT(N2)
                if ((CELLS(N1, J1) == CELLS(N2, J2)) .and. CELLS(N1,J1) /= ii) NBE(ii, 3) = CELLS(N1, J1)
            end do
            end do
            do J2 = 1, CELLCNT(N2)
            do J3 = 1, CELLCNT(N3)
                IF ((CELLS(N2, J2) == CELLS(N3, J3)) .and. CELLS(N2, J2) /= ii) NBE(ii,1) = CELLS(N2, J2)
            end do
            end do
            do J1 = 1, CELLCNT(N1)
            do J3 = 1, CELLCNT(N3)
                IF((CELLS(N1, J1) == CELLS(N3,J3)) .and. CELLS(N1,J1) /= ii) NBE(ii, 2) = CELLS(N3, J3)
            end do
            end do
        end do
        deallocate(CELLS, CELLCNT)
    
        ! Ensure all elements have at least one neighbor
        do ii = 1, N
            if (sum(NBE(ii, 1:3)) == 0) then
            print *, 'cell ', ii, ' @ ', xc(ii), yc(ii), ' has no neighbors'
            stop
            end if
        end do
    
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
        do ii = 1, M
            NCNT = 0
            do jj = 1, N
            if ( float(NV(jj, 1) - ii) * float(NV(jj, 2) - ii) * float(NV(jj, 3) - ii) == 0.0_SP ) NCNT = NCNT + 1
            end do
            MX_NBR_ELEM = MAX(MX_NBR_ELEM, NCNT)
        end do
    
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

        !--Reorder Order Elements Surrounding a Node to Go in a Cyclical Procession----!
        !--Determine NTSN  = Number of Nodes Surrounding a Node (+1)-------------------!
        !--Determine NBSN  = Node Numbers of Nodes Surrounding a Node------------------!
        allocate(NB_TMP(M, MX_NBR_ELEM+1))
        do ii = 1, M
            if (ISONB(ii) == 0) then
            NB_TMP(1, 1) = NBVE(ii, 1)
            NB_TMP(1, 2) = NBVT(ii, 1)
            do jj = 2, NTVE(ii) + 1
                kk = NB_TMP(jj - 1, 1)
                ll = NB_TMP(jj - 1, 2)
                NB_TMP(jj, 1) = NBE(kk, ll + 1 - INT((ll + 1)/4)*3)
                ll=NB_TMP(jj, 1)
                if ((NV(ll, 1) - ii) == 0) NB_TMP(jj, 2) = 1
                if ((NV(ll, 2) - ii) == 0) NB_TMP(jj, 2) = 2
                if ((NV(ll, 3) - ii) == 0) NB_TMP(jj, 2) = 3
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
    
            ai = (/ y2 - y3, y3 - y1, y1 - y2 /)
            bi = (/ x3 - x2, x1 - x3, x2 - x1 /)
            ci = (/ x2*y3 - x3*y2, x3*y1 - x1*y3, x1*y2 - x2*y1 /)
    
            aw0(ii, 1:3) = -1.0 * ci / 2.0 / art(ii)
            awx(ii, 1:3) = -1.0 * ai / 2.0 / art(ii)
            awy(ii, 1:3) = -1.0 * bi / 2.0 / art(ii)
        end do
    
        do ii = 1, n
            if (isbce(ii) > 1) then
    
            a1u(ii, 1:4) = zero
            a2u(ii, 1:4) = zero
    
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
    
            a1u(ii, (/ 0, ll, j1, j2 /) + 1) = zero
            a2u(ii, (/ 0, ll, j1, j2 /) + 1) = zero
            end if
        end do
    end subroutine
end module
