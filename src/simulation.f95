module simulation
    use variables, only : sp
    use random, only : random_number_generator
    
    implicit none
    save
    private

    type, public :: Experiment
        character(len = 100), public :: simID
        integer, public :: &
            & nnodes=0, &
            & nelements=0, &
            & nlayers=0, &
            & lines_read=0
        real(sp), public :: &
            & globalIrradiance=0.0_sp, &
            & meshArea=0.0_sp, &
            & layerDepth=0.0_sp, &
            & layerSigma=0.0_sp, &
            & time=0.0_sp, & ! total time
            & daytime=0.0_sp, & ! time in days
            & clocktime=0.0_sp ! twenty four hour periodic
        real(sp), allocatable, dimension(:), public :: elementSigmaVolume, elementArea ! mesh stats
        real(sp), allocatable, dimension(:), public :: verticaltox, verticaldiff, verticaltemp, verticalrho ! uniform horizontal fields
    contains
        ! call in this order
        procedure, public :: init => initSimulation
        procedure, public :: read => readSimulation
        procedure, public :: diffuse => verticalDiffusion
    end type
  
    type(Experiment), allocatable, public :: domain ! domain structure imported from this module
    public :: topology, globalIrradiance, write_rectangular_mesh_files
contains
    pure elemental function globalIrradiance(time, irradiance) result(light)
        ! calculate light based on time and maximum value for region
        real(sp) :: light, daytime, clock
        real(sp), intent(in) :: time, irradiance
        daytime = time / 24.0_sp
        clock = time - 24.0_sp*floor(daytime)
        light = 0.0_sp
        if ((clock > 6.0_sp) .and. (clock < 18.0_sp)) then
            light = 0.5_SP*irradiance*(1.0_SP + cos(4.0_SP*3.14159*daytime))
        end if
    end function

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
        self%elementArea = 0.0_sp
        self%elementSigmaVolume = 0.0_sp
        self%verticaldiff = 0.0_sp
        self%verticaltemp = 0.0_sp
        self%verticalrho = 0.0_sp
    end subroutine
  
    subroutine readSimulation(self, u_vel, v_vel, w_vel, diffusivity, elevation, salinity, temperature, density)
        ! read physical drivers from file
        use variables, only : layers, M, N, iophys
    
        class(Experiment), intent(inout) :: self
        real(sp), dimension(0:N, layers), intent(inout) :: u_vel, v_vel, w_vel
        real(sp), dimension(0:M, layers), intent(inout) :: diffusivity, salinity, temperature, density
        real(sp), dimension(0:M), intent(inout) :: elevation
    
        character(len = 100) :: vert_format
        real(sp) :: time
        integer :: ii
    
        write(vert_format, "(A7,I6,A7)") "(F10.3,", 3*layers, "F20.10)"
        read(iophys, vert_format) time, self%verticaltemp(1:layers), self%verticalrho(1:layers), self%verticaldiff(1:layers)
    
        u_vel(:, :) = 0.0_sp
        v_vel(:, :) = 0.0_sp 
        w_vel(:, :) = 0.0_sp
        elevation(:) = -abs(0.0_sp)
        salinity(:,:) = 0.0_sp
        do ii = 1, self%nnodes
            temperature(ii, 1:layers) = self%verticaltemp(1:layers)
            diffusivity(ii, 1:layers) = self%verticaldiff(1:layers)
            density(ii, 1:layers) = self%verticalrho(1:layers)
        end do
        if (self%lines_read == 0) rewind(unit=iophys)
        self%lines_read = self%lines_read + 1
    end subroutine

    subroutine verticalDiffusion(self)
        ! One-dimensional vertical diffusion
        use variables, only : layers, KBM1, dti
        class(Experiment), intent(inout) :: self
        integer :: ii, jj, steps_to_stability
        real(sp) :: &
            & substep, & ! automatic time substep, only valid for dK/dt=0
            & stable, &  ! longest time step for conditional stability
            & diffusivity = 60.0*60.0*10.0**(-5.0), &  ! 0.1 cm^2/s
            & profile(0:self%nlayers+1); 
        
        profile = 0.0_sp
        stable = 0.5 * self%layerDepth**2.0
        steps_to_stability = ceiling(dti/stable) ! min number of steps to achieve stability, cannot be less than one
        substep = dti / float(steps_to_stability)
    
        do ii = 1, steps_to_stability
            self%verticaltox(1) = self%verticaltox(2)
            self%verticaltox(layers) = self%verticaltox(KBM1)
    
            ! central-in-space (z), forward-in-time, second derivative
            do jj = 2, KBM1
            profile(jj) = self%verticaltox(jj) + diffusivity * substep * &
                    & (self%verticaltox(jj-1) + self%verticaltox(jj+1) - 2.0*self%verticaltox(jj)) / &
                    & self%layerDepth * self%layerDepth
            end do
    
            profile(1) = profile(2);
            profile(layers) = profile(KBM1) ! copy in domain value to boundary nodes
            self%verticaltox(:) = profile(:)
        end do
    end subroutine

    subroutine topology
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
        use variables, only : n, a1u, NBVE, NBVT, vxmin, vxmax, vymin, VYMAX, aw0, a2u, NTVE

        integer, allocatable, dimension(:, :) :: NB_TMP, CELLS
        integer :: ii, jj, kk, ll, count, JJB, N1, N2, N3, J1, J2, J3, tri(3), MX_NBR_ELEM
        real(sp) :: X1, X2, X3, Y1, Y2, Y3, DELT, AI(3), BI(3), CI(3), area(N)

        ! calculate global coordinate ranges
        VXMIN = minval(VX(1:M))
        VXMAX = maxval(VX(1:M))
        VYMIN = minval(VY(1:M))
        VYMAX = maxval(VY(1:M))
    
        ! calculate global element center grid coordinates
        xc = 0.0_sp
        yc = 0.0_sp
        do ii = 1, N
            xc(ii) = sum(VX(NV(ii, :))) / 3.0_SP
            yc(ii) = sum(VY(NV(ii, :))) / 3.0_SP
        end do
    
        area(:) = abs(((VX(NV(:, 2)) - VX(NV(:, 1))) * (VY(NV(:, 3)) - VY(NV(:, 1))) - &
                & (VX(NV(:, 3)) - VX(NV(:, 1))) * (VY(NV(:, 2)) - VY(NV(:, 1))))*0.5)

        ISBCE = 0
        ISONB = 0
        NBE   = 0

        allocate(CELLS(M, 50))

        cells = 0
    
        ! For each element, get the node indices, and increment the count of
        ! triangles in which the node appears.
        do ii = 1, n
            tri = nv(ii, 1:3)
            ntve(tri) = ntve(tri) + 1
            cells(tri, ntve(tri)) = ii
        end do
    
        do ii = 1, N
            N1 = NV(ii,1)
            N2 = NV(ii,2)
            N3 = NV(ii,3)
            do J1 = 1, ntve(N1)
                do J2 = 1, ntve(N2)
                    if ((CELLS(N1, J1) == CELLS(N2, J2)) .and. CELLS(N1,J1) /= ii) NBE(ii, 3) = CELLS(N1, J1)
                end do
            end do
            do J2 = 1, ntve(N2)
                do J3 = 1, ntve(N3)
                    IF ((CELLS(N2, J2) == CELLS(N3, J3)) .and. CELLS(N2, J2) /= ii) NBE(ii,1) = CELLS(N2, J2)
                end do
            end do
            do J1 = 1, ntve(N1)
                do J3 = 1, ntve(N3)
                    IF((CELLS(N1, J1) == CELLS(N3,J3)) .and. CELLS(N1,J1) /= ii) NBE(ii, 2) = CELLS(N3, J3)
                end do
            end do
        end do
        deallocate(CELLS)

        ! Flag boundary elements and nodes
        do ii = 1, N
            if (NBE(ii, 1) == 0) then
                ISONB(NV(ii, 2)) = 1
                ISONB(NV(ii, 3)) = 1
            end if
            if (NBE(ii, 2) == 0) then
                ISONB(NV(ii, 1)) = 1
                ISONB(NV(ii, 3)) = 1
            end if
            if (NBE(ii, 3) == 0) then
                ISONB(NV(ii, 1)) = 1
                ISONB(NV(ii, 2)) = 1
            end if
            if (any(NBE(ii, :) == 1)) then
                ISBCE(ii) = 1
            end if
        end do

        MX_NBR_ELEM = maxval(ntve)
        allocate(NBVE(M, MX_NBR_ELEM + 1))
        allocate(NBVT, NB_TMP, mold=NBVE)
        ! For each node, scan through all elements, to determine number of times the node
        ! appears in the element list.  If it appears, then store the bi-directional mapping.
        do ii = 1, M
            count = 0
            do jj = 1, N
                if ( (NV(jj, 1) - ii) * (NV(jj, 2) - ii) * (NV(jj, 3) - ii) == 0) then
                    count = count + 1
                    NBVE(ii, count) = jj
                    if ((NV(jj, 1) - ii) == 0) NBVT(ii, count) = 1
                    if ((NV(jj, 2) - ii) == 0) NBVT(ii, count) = 2
                    if ((NV(jj, 3) - ii) == 0) NBVT(ii, count) = 3
                end if
            end do
        end do

        ! Reorder Order Elements Surrounding a Node to Go in a Cyclical Procession
        ! Determine NTSN = Number of Nodes Surrounding a Node (+1)
        ! Determine NBSN = Indices of neighboring nodes

        ! Loop through each node
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
                    NBVT(ii, jj) = NB_TMP(jj, 2)
                end do
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

        ! Calculate shape coefficients for each element
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
            else if (isbce(ii) == 1) then
                do jj = 1, 3
                    if (nbe(ii, jj) == 0) ll = jj
                end do
                j1 = ll + 1 - int((ll + 1)/4)*3
                j2 = ll + 2 - int((ll + 2)/4)*3
                a1u(ii, [0, ll, j1, j2] + 1) = 0.0_sp
                a2u(ii, [0, ll, j1, j2] + 1) = 0.0_sp
            else if (isbce(ii) > 1) then
                a1u(ii, 1:4) = 0.0_sp
                a2u(ii, 1:4) = 0.0_sp
            end if

            x1 = vx(nv(ii, 1)) - xc(ii)
            x2 = vx(nv(ii, 2)) - xc(ii)
            x3 = vx(nv(ii, 3)) - xc(ii)
            y1 = vy(nv(ii, 1)) - yc(ii)
            y2 = vy(nv(ii, 2)) - yc(ii)
            y3 = vy(nv(ii, 3)) - yc(ii)

            ai = [y2 - y3, y3 - y1, y1 - y2]
            bi = [x3 - x2, x1 - x3, x2 - x1]
            ci = [x2*y3 - x3*y2, x3*y1 - x1*y3, x1*y2 - x2*y1]

            aw0(ii, 1:3) = -1.0 * ci / 2.0 / area(ii)
            awx(ii, 1:3) = -1.0 * ai / 2.0 / area(ii)
            awy(ii, 1:3) = -1.0 * bi / 2.0 / area(ii)
        end do
    end subroutine


    subroutine write_rectangular_mesh_files(prefix, bottom_depth, div_inc, x_range, y_range, layers)
        ! This program generates a simple triangular grid for the simulation
        character(len=50), intent(in) :: prefix
        real(sp), intent(in) :: &
            & bottom_depth, div_inc, x_range, y_range
        integer, intent(in) :: layers
        integer :: &
            & nodes, elements, node_index=1, element_index=1, start_pattern=1, &
            & alternate_pattern, x_div, y_div, ii, jj, fid=500, fid2=501, &
            & x_nodes, y_nodes
        integer, dimension(3) :: vertices

        x_div = floor( x_range / div_inc )
        y_div = floor( y_range / div_inc )
        x_nodes = x_div + 1
        y_nodes = y_div + 1
        nodes = (x_div + 1) * (y_div + 1)
        elements = 2 * x_div * y_div

        print *, "Saving mesh data... "
        write(*, '(A,I5)') "    Nodes:    ", nodes
        write(*, '(A,I5)') "    Elements: ", elements
        open(unit=fid, file=trim(prefix)//"/nodes.txt", status='replace')
        open(unit=fid2, file=trim(prefix)//"/elements.txt", status='replace')
        write(fid, *) nodes
        write(fid2, *) elements, layers
        do ii = 1, y_nodes
            alternate_pattern = start_pattern
            do jj = 1, x_nodes
                write(fid, *) float(jj-1) * div_inc, float(ii-1) * div_inc, -abs(bottom_depth)
                if ( (jj < x_nodes) .and. (ii < y_nodes) ) then
                    ! bottom row of triangles
                    vertices(1) = node_index
                    vertices(2) = node_index + x_nodes ! note index
                    vertices(3) = node_index + 1
                    if (alternate_pattern < 0) vertices(2) = vertices(2) + 1
                    write(fid2, *) vertices(1:3)
                    element_index = element_index + 1
                
                    ! upper row of triangles
                    vertices(1) = node_index + x_nodes
                    vertices(2) = node_index + x_nodes + 1
                    vertices(3) = node_index + 1
                    if (alternate_pattern < 0) vertices(3) = vertices(3) - 1
                    write(fid2, *) vertices(1:3)
                    element_index = element_index + 1
                    alternate_pattern = -alternate_pattern ! switch triangle pattern as you move x=0 to x=x_range
                end if
                node_index = node_index + 1 
            end do
            start_pattern = -start_pattern ! switch triangle pattern between rows as grid is built
        end do
        close(fid)
        close(fid2)
    end subroutine

end module
