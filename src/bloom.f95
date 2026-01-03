program bloom

    use variables, only: sp, layers, M, N, KBM1, KBM2, irradSurf, folderprefix, duration, DTI, DTOUT, INSTP, DTRW, xc, yc, vx, vy, nv, isbce, nbe, isonb, awx, awy, u, v, ww, ut, vt, wwt, kh, z, zz, dz, a1u, NBVE, NBVT, vxmin, vxmax, vymin, VYMAX, aw0, a2u, NTVE, h, d, el, et, t1, tt1, iophys, iotox, iocs, iocp
    use random, only : random_number_generator
    use simulation, only : domain, topology
    use cyanobacteria, only : CyanobacteriaAgent
    use io, only : scanInteger, scanReal, scanString, scanLogical

    implicit none

    type(CyanobacteriaAgent) :: particles

    real(sp), allocatable, dimension(:) :: length1, length2, length3, ss, average_thickness
    ! Element-based velocity fields start and end of hour
    real(sp), allocatable, dimension(:, :) :: UNC, UNC2, VNC, VNC2, WNC, WNC2 
    ! Node-based diffusion and physical fields start and end of hour
    real(sp), allocatable, dimension(:, :) :: KHNC, KHNC1, KHNC2, salinity, SALNC2, temperature, TEMPNC2, density, RHONC2
    ! Node-based free surface height field read, start and end of hour
    real(sp), allocatable, dimension(:) :: ELNC, ELNC1, ELNC2
    integer :: NH, IT, HOUR, IINT, ii, index, ionode=100, ioelem=101, exp_type, start=1, end
    character(len = 80) :: input_file, foldername, state_format, file

    call get_command_argument(1, foldername) ! Import case name from command line
    if (len_trim(foldername) == 0) then
        print *, 'Provide simulation ID as argument'; stop
    end if
    folderprefix = "data/"//adjustl(foldername)


    ! Read in variables and set values
    write(*, *);
    write(*, "(A)", advance='no') 'Importing simulation parameters... '
    file = trim(folderprefix)//"/parameters.dat"
    exp_type = scanInteger(file, "EXPTYPE") ! Experiment type
    HOUR = scanInteger(file, "HOUR") ! day
    duration = scanInteger(file, "TDRIFT") ! Simulation duration
    DTI = scanReal(file, "DTI") ! Input time step
    DTOUT = scanReal(file, "DTOUT") ! Output time step
    INSTP = scanReal(file, "INSTP") ! Input time step of flow fields
    DTRW = scanReal(file, "DTRW") ! Random walk time step
    write(*, *) "Finished"

    ! Determine number of elements and nodes in the model
    write(*, "(A)", advance='no') "Getting element data... "
    allocate(domain)
    file = trim(folderprefix)//"/elements.txt"
    open(ioelem, file=file)
    read(ioelem, *) N, layers ! get number of nodes, elements, and sigma layers

    domain%nelements = N
    domain%nlayers = layers
    KBM1 = layers - 1
    KBM2 = layers - 2
    allocate(NV(0:N, 4)); NV = 0  ! child nodes
    allocate(NBE(0:N, 3)); NBE = 0  ! element neighbors
    allocate(ISBCE(0:N)); ISBCE = 0
    ! Element based 2D fields
    allocate(length1(0:N)); length1 = 0.0_sp
    allocate(&
        & length2, &
        & length3, &
        & ss, &
        & average_thickness, &
        & XC, & ! center x-coordinate
        & YC, & ! center y-coordinate
        & source=length1)
    ! Element based 3D fields
    allocate(U(0:N, layers)); U = 0.0_sp  ! x-velocity
    allocate(&
        & V, & ! y-velocity
        & WW, & ! z-velocity
        & UT, & ! previous x-velocity
        & VT, & ! previous y-velocity
        & WWT, & ! previous z-velocity
        & KH, & ! turbulent quantity
        & UNC, &
        & VNC, &
        & WNC, &
        & UNC2, &
        & VNC2, &
        & WNC2, &
        & source=U)

    ! Shape coefficient arrays and control volume metrics
    allocate(A1U(0:N, 4)); A1U = 0.0_sp
    allocate(A2U, source=A1U)
    allocate(AWX(0:N, 3)); AWX = 0.0_sp
    allocate(AWY, AW0, source=AWX)
    do ii = 1, N
        read(ioelem, *) index, NV(ii, 1), NV(ii, 3), NV(ii, 2)
    end do
    NV(:, 4) = NV(:, 1) ! duplicate node for computation
    close(ioelem)
    ! Sigma based 1D fields and coordinates
    allocate(Z(layers)); Z = 0.0_sp ! s-position
    allocate(&
        & ZZ, & ! INTRA LEVEL SIGMA VALUE
        & DZ, & ! DELTA-SIGMA VALUE
        & source=Z)

    write(*, *) "Finished"
    write(*, "(A)", advance='no') "Getting node data... "
    file = trim(folderprefix)//"/nodes.txt"
    open(ionode, file=file)
    read(ionode, *) M
    domain%nnodes = M

    ! Node, Boundary Condition, and Control Volume
    allocate(NTVE(0:M)); NTVE = 0
    allocate(ISONB, source=NTVE);  ! NODE MARKER = 0,1,2
    ! Node based 2D floating point fields and coordinates
    allocate(H(0:M)); H = 0.0_sp ! bathymetry
    allocate(&
        & D, & ! depth
        & EL, & ! surface elevation
        & ET, & ! previous surface elevation
        & VX, & ! x-position
        & VY, & ! y-position
        & ELNC, &
        & ELNC1, &
        & ELNC2, &
        & source=H)
    ! Node based floating point 3D fields
    allocate(T1(0:M, layers)); T1 = 0.0_sp ! temperature
    allocate(&
        & TT1, & ! previous temperature
        & KHNC, &
        & KHNC1, &
        & KHNC2, &
        & salinity, &
        & SALNC2, &
        & temperature, &
        & TEMPNC2, &
        & density, &
        & RHONC2, &
        & source=T1)

    do ii = 1, M
        read(ionode, *) index, VX(ii), VY(ii), H(ii) ! get node position and depth
    end do
    close(ionode)

    write(*, *) "Finished"

    call domain%init(toxin=0.0_sp)
    
    write(*, "(A)", advance='no') "Computing mesh topology... "
    call topology
    write(*, *) "Finished"

    state_format="(1F12.6, 9000(I6,3F12.6))"

    allocate(random_number_generator)
    call random_number_generator%init()

    write(*, "(A)", advance='no') "Loading physical field data... "
    open(unit=iophys, file=input_file, status='old', position='rewind')
    call domain%read(UNC, VNC, WNC, KHNC, ELNC, salinity, temperature, density)
    write(*, *) "Finished"

    start = HOUR ! tracking begin iteration is 0
    end = start + duration - 1 ! tracking end iteration

    write(*, "(A)", advance='no') "Initializing particle structures... "
    call particles%init(exp_type) ! initialize and allocate type specific structures
    write(*, *) "Finished"
    write(*, "(A)", advance='no') "Allocating common variables... "
    call particles%lag_alloc() ! allocate common variables other than position

    write(*, *) "Finished"
    write(*, "(A)", advance='no') "Finding host elements... "
    call particles%find_host_element(particles%XP, particles%YP)
    write(*, *) "Finished"

    write(*, "(A)", advance='no') "Interpolating physical fields... "
    call particles%INTERP_FIELDS(particles%XP, particles%YP, particles%ZP, salinity, temperature, density, H, ELNC)
    write(*, *) "Finished"

    write(*, "(A)", advance='no') "Adjusting vertical domain... "
    particles%ZP(:) = -1.0_sp*abs(particles%ZP(:)) ! make depth negative
    particles%LAYER(:) = particles%zlocate(particles%ZP(:)) ! valid when sigma layers are equal thickness
    write(*, *) "Finished"

    write(*, "(A)", advance='no') "Writing position and state variables to file... "
    open(unit=iocp, file="./"//trim(folderprefix)//"/"//"cyanobacteria_position.dat", status='replace') ! create new position file
    open(unit=iocs, file="./"//trim(folderprefix)//"/"//"cyanobacteria_state.dat", status='replace') ! create new state variable file
    domain%time = float(HOUR) ! decimal simulation time for output timestamp
    call particles%writeState(iocs, domain%time)
    call particles%writePosition(iocp, domain%time) ! write particle positions to output file
    write(*, *) "Finished"
    write(*, *) "Starting simulation loop... "
    call domain%read(UT, VT, WWT, KHNC1, ET, ST1, TT1, RT1) ! read physical fields
    HOUR = HOUR + 1
    do ii = 1, domain%nelements
        length1(ii) = sqrt(  (vx(nv(ii, 1)) - vx(nv(ii, 2)) )**(2.0_sp) + ( vy(nv(ii, 1)) - vy(nv(ii, 2)) )**(2.0_sp)  )
        length2(ii) = sqrt(  (vx(nv(ii, 2)) - vx(nv(ii, 3)) )**(2.0_sp) + ( vy(nv(ii, 2)) - vy(nv(ii, 3)) )**(2.0_sp)  )
        length3(ii) = sqrt(  (vx(nv(ii, 3)) - vx(nv(ii, 1)) )**(2.0_sp) + ( vy(nv(ii, 3)) - vy(nv(ii, 1)) )**(2.0_sp)  )
        average_thickness = abs(sum( H(nv(ii,1:3)) )/3.0_sp) ! must be positive
    end do

    ss(:) = 0.5_sp*(length1(:) + length2(:)+ length3(:)) ! herons's formula for area
    domain%elementArea = sqrt( ss * (ss(:) - length1(:)) * (ss(:) - length1(:)) * (ss(:) - length3(:)) ) 
    domain%meshArea = sum(domain%elementArea)
    domain%elementSigmaVolume = domain%elementArea * average_thickness
    domain%layerDepth = average_thickness(1)
    domain%layerSigma = float(KBM1)**(-1.0_SP)

    open(unit=iotox, file=trim(folderprefix)//"/dissolved_toxin.txt", status='replace')
    IINT = 0
    do NH = start, end
        write(*,*);
        write(*, "(I4,A,I4,A)", advance='no') NH-start+1, ' / ', end-start+1, ' steps: '
        call domain%read(UNC2, VNC2, WNC2, KHNC2, ELNC2, SALNC2, TEMPNC2, RHONC2)
        HOUR = HOUR + 1
        ! length (float) of flow field interpolation divided by length (float) of inner time step
        do IT = 1, int(INSTP/DTI)

            write(*, "(A)", advance='no') "|"
            IINT = IINT + 1 ! increment total step count
            domain%time = float(IINT)*DTI + float(start) ! decimal simulation time for output timestamp
            domain%daytime = domain%time / 24.0_sp
            domain%clocktime = domain%time - 24.0_sp*floor(domain%daytime)

            if ((domain%clocktime > 6.0_sp) .and. (domain%clocktime < 18.0_sp)) then
                domain%globalIrradiance = 0.5_SP*irradSurf*(1.0_SP + cos(4.0_SP*3.14159*domain%daytime))
            else
                domain%globalIrradiance = 0.0_sp
            end if

            call domain%diffuse() ! vertical diffusion of dissolved toxin
            call particles%movement()
            call particles%random(noise=particles%noise) ! random walk
            if ( mod(IINT, int(DTOUT/DTI) ) == 0) then
                call particles%writePosition(iocp, time=domain%time) ! output position, same for all particle types
                call particles%writeState(iocs, time=domain%time)
                write(iotox,"(1f20.3,51F20.6)") domain%time, domain%verticaltox(1:layers)
            end if
        end do
    end do

    write(*,*); write(*,*); write(*,*) "Simulation finished."; write(*,*)
    ! call random%stats()

end program