program bloom

    use variables
    use random, only : random_number_generator
    use simulation, only : domain, TRIANGLE_GRID_EDGE
    use cyanobacteria, only : CyanobacteriaAgent
    use io, only : scanInteger, scanReal, scanString, scanLogical
    
    implicit none

    type(CyanobacteriaAgent) :: particles

    real(sp), allocatable, dimension(:) :: length1, length2, length3, ss, average_thickness
    real(sp), dimension(0:N, KB) :: UNC, UNC2, VNC, VNC2, WNC, WNC2 ! velocity fields start and end of hour
    real(sp), dimension(0:M, KB) :: KHNC, KHNC1, KHNC2, SALNC, SALNC2, TEMPNC, TEMPNC2, RHONC, RHONC2 ! diffusion and physical fields start and end of hour
    real(sp), dimension(0:M) :: ELNC, ELNC1, ELNC2 ! free surface height field read, start and end of hour
    real(sp) :: LAG_TIME

    integer , allocatable, dimension(:) :: INWATER
    integer :: NH, IT, HOUR, IINT, ii, index, ionode=100, ioelem=101, exp_type, NCT, ISCAN
    character(len = 100) :: input_file, input, meshfile, foldername, state_format
    character(len = 80) :: file

    call get_command_argument(1, foldername) ! Import case name from command line
    if (len_trim(foldername) == 0) then
        print *, 'Provide simulation ID as argument'; stop
    end if
    folderprefix = "data/"//trim(adjustl(foldername))


    ! Read in variables and set values
    write(*, *);
    write(*, "(A)", advance='no') 'Importing simulation parameters... '
    file = trim(folderprefix)//"/parameters.dat"
    exp_type = scanInteger(file, "EXPTYPE") ! Experiment type
    INFOFILE = scanString(file, "INFOFILE") ! Info file
    HOURLAG = scanInteger(file, "HOURLAG") ! day
    TDRIFT = scanInteger(file, "TDRIFT") ! Total time to move drifters (TDRIFT)
    DTI = scanReal(file, "DTI") ! External time step (DTI)
    DTOUT = scanReal(file, "DTOUT") ! External time step (DTOUT)
    INSTP = scanReal(file, "INSTP") ! Input time step of flow fields (instp)
    DTRW = scanReal(file, "DTRW") ! Random walk time step
    DHOR = scanReal(file, "DHOR") ! Horizontal diffusion coefficient (DHOR)
    P_SIGMA = scanLogical(file, "P_SIGMA") ! Sigma coordinate system
    OUT_SIGMA = scanLogical(file, "OUT_SIGMA") ! Sigma coordinate system
    F_DEPTH = scanLogical(file, "F_DEPTH") ! fixed depth


    ! Determine number of elements and nodes in the model
    write(*, "(A)", advance='no') "Reading mesh files... "
    allocate(domain)

    write(*, *) "Finished"
    write(*, "(A)", advance='no') "Getting element data... "
    meshfile = trim(folderprefix)//"/elements.txt"
    open(ioelem, file = meshfile)
    read(ioelem, *) N, KB ! get number of nodes, elements, and sigma layers
    domain%nelements = N
    domain%nlayers = KB
    NCT = N*3
    allocate(length1(0:N)); length1 = zero
    allocate(length2(0:N)); length2 = zero
    allocate(length3(0:N)); length3 = zero
    allocate(ss(0:N)); ss = zero
    allocate(average_thickness(0:N)); average_thickness = zero
    allocate(XC(0:N)); XC = zero ! X-COORD AT FACE CENTER
    allocate(YC(0:N)); YC = zero ! Y-COORD AT FACE CENTER
    allocate(NV(0:N,4))           ;NV = 0  ! NODE NUMBERING FOR ELEMENTS
    allocate(NBE(0:N,3))          ;NBE = 0  ! INDICES OF ELEMENT NEIGHBORS
    allocate(ISBCE(0:N))          ;ISBCE = 0
    ! internal mode arrays-(element based)
    allocate(U(0:N, KB))       ;U = zero   ! X-VELOCITY
    allocate(V(0:N, KB))       ;V = zero   ! Y-VELOCITY
    allocate(W(0:N, KB))       ;W = zero   ! VERTICAL VELOCITY IN SIGMA SYSTEM
    allocate(WW(0:N, KB))      ;WW = zero   ! Z-VELOCITY
    allocate(UT(0:N, KB))      ;UT = zero   ! X-VELOCITY FROM PREVIOUS TIMESTEP
    allocate(VT(0:N, KB))      ;VT = zero   ! Y-VELOCITY FROM PREVIOUS TIMESTEP
    allocate(WT(0:N, KB))      ;WT = zero   ! VERTICAL VELOCITY FROM PREVIOUS TIMESTEP
    allocate(WWT(0:N, KB))     ;WWT = zero   ! Z-VELOCITY FROM PREVIOUS TIMESTEP
    allocate(KH(0:N, KB))     ;KH = zero   ! TURBULENT QUANTITY
    ! Shape coefficient arrays and control volume metrics
    allocate(A1U(0:N, 4))         ;A1U = zero
    allocate(A2U(0:N, 4))         ;A2U = zero
    allocate(AWX(0:N, 3))         ;AWX = zero
    allocate(AWY(0:N, 3))         ;AWY = zero
    allocate(AW0(0:N, 3))         ;AW0 = zero
    do ii = 1, N
        read(ioelem, *) index, NV(ii, 1), NV(ii, 3), NV(ii, 2)
    end do
    NV(:, 4) = NV(:, 1) ! duplicate node for computation
    close(ioelem)
    ! 1-d arrays for the sigma coordinate
    allocate(Z(KB))               ; Z = zero    ! SIGMA COORDINATE VALUE
    allocate(ZZ(KB))              ; ZZ = zero    ! INTRA LEVEL SIGMA VALUE
    allocate(DZ(KB))              ; DZ = zero    ! DELTA-SIGMA VALUE
    allocate(DZZ(KB))             ; DZZ = zero    ! DELTA OF INTRA LEVEL SIGMA
    KBM1 = KB - 1
    KBM2 = KB - 2
    write(*, *) "Finished"


    write(*, "(A)", advance='no') "Getting node data... "
    open(ionode, file = trim(folderprefix)//"/nodes.txt")

    read(ionode, *) M

    domain%nnodes = M


  
    
    ! Grid Metrics
    allocate(VX(0:M)); VX = zero ! X-COORD AT GRID POINT
    allocate(VY(0:M)); VY = zero ! Y-COORD AT GRID POINT
    ! Node, Boundary Condition, and Control Volume
    allocate(NTVE(0:M))           ;NTVE = 0
    allocate(ISONB(0:M))          ;ISONB = 0  ! NODE MARKER = 0,1,2
    ! 2-d flow variable arrays at nodes
    allocate(H(0:M))       ;H = zero       ! BATHYMETRIC DEPTH
    allocate(D(0:M))       ;D = zero       ! DEPTH
    allocate(EL(0:M))      ;EL = zero       ! SURFACE ELEVATION
    allocate(ET(0:M))      ;ET = zero       ! SURFACE ELEVATION PREVIOUS TIMESTEP
    ! 3d variable arrays-(node based)
    allocate(T1(0:M, KB))       ;T1 = zero  ! TEMPERATURE AT NODES
    allocate(S1(0:M, KB))       ;S1 = zero  ! SALINITY AT NODES
    allocate(R1(0:M, KB))       ;R1 = zero  ! DENSITY AT NODES
    allocate(TT1(0:M, KB))      ;TT1 = zero  ! TEMPERATURE FROM PREVIOUS TIME
    allocate(ST1(0:M, KB))      ;ST1 = zero  ! SALINITY FROM PREVIOUS TIME
    allocate(RT1(0:M, KB))      ;RT1 = zero
    allocate(WTS(0:M, KB))      ;WTS = zero  ! VERTICAL VELOCITY IN SIGMA SYSTEM

    do ii = 1, M
        read(ionode, *) index, VX(ii), VY(ii), H(ii) ! get node position and depth
    end do
    close(ionode)

    write(*, *) "Finished"
    write(*, *)
    write(*, *) '    Nodes        :', M
    write(*, *) '    Elements     :', N
    write(*, *) '    Sigma layers :', KB
    write(*, *)

    call domain%init(exp_type) ! allocate and initialize global environmental variables, and additional mesh-based variables
    
    write(*, "(A)", advance='no') "Computing mesh topology... "
    call TRIANGLE_GRID_EDGE
    write(*, *) "Finished"

    state_format="(1F12.6, 9000(I6,3F12.6))"

    write(*, "(A)", advance='no') "Allocating cyanobacteria structures... "
    allocate(particles)
    allocate(random_number_generator)
    call random_number_generator%init()
    write(*, *) "Finished"

    write(*, "(A)", advance='no') "Loading physical field data... "
    input_file = trim(folderprefix)//"/forcing.txt"
    open(unit=iophys, file=input_file, status='old', position='rewind')
    call domain%read(UNC, VNC, WNC, KHNC, ELNC, SALNC, TEMPNC, RHONC)
    write(*, *) "Finished"

    HOUR = HOURLAG ! start time (int) is 0
    ISLAG = HOURLAG ! tracking begin iteration is 0
    IELAG = ISLAG + TDRIFT - 1 ! tracking end iteration

    write(*, "(A)", advance='no') "Initializing particle structures... "
    call particles%init(exp_type) ! initialize and allocate type specific structures, also reads position and state variables from file
    write(*, *) "Finished"
    write(*, "(A)", advance='no') "Allocating common variables... "
    call particles%lag_alloc() ! allocate common variables other than position and itag

    allocate( INWATER(particles%ndrft) ); INWATER(:) = 1
    write(*, *) "Finished"

    particles%XP(:) = particles%XPT(:) ! Shift x to model coordinate system
    particles%YP(:) = particles%YPT(:) ! Shift y to model coordinate system

    write(*, "(A)", advance='no') "Finding host elements... "
    call particles%find_host_element(particles%XP, particles%YP, INWATER)
    where (particles%FOUND == 0) particles%indomain(:) = .false.
    write(*, *) "Finished"

    write(*, "(A)", advance='no') "Interpolating physical fields... "
    call particles%INTERP_ELH(particles%XP, particles%YP, H, ELNC, 1) ! interpolate elevation and bathymetry
    call particles%INTERP_FIELDS(particles%XP, particles%YP, particles%ZP, SALNC, TEMPNC, RHONC, 0) ! interpolate salinity and temperature
    write(*, *) "Finished"

    write(*, "(A)", advance='no') "Adjusting vertical domain... "
    particles%ZPT(:) = -1.0_sp*abs(particles%ZPT(:)) ! make depth negative
    particles%ZP(:) = particles%sigma(particles%ZPT(:)) ! convert to sigma coordinate
    particles%LAYER(:) = particles%zlocate(particles%ZP(:)) ! valid when sigma layers are equal thickness

    write(*, "(A)", advance='no') "Writing position and state variables to file... "
    open(unit=iocp, file="./"//trim(folderprefix)//"/"//"cyanobacteria_position.dat", status='replace') ! create new position file
    open(unit=iocs, file="./"//trim(folderprefix)//"/"//"cyanobacteria_state.dat", status='replace') ! create new state variable file
    call particles%writeState(iocs)
    call particles%writePosition(iocp) ! write particle positions to output file
    write(*, *) "Finished"

    write(*, *) ! Print particle statistics
    write(*, *) '    Tracking Info'
    write(*, *) '        Start iteration :', ISLAG
    write(*, *) '        Final iteration :', IELAG
    write(*, *)
    call particles%stats

    write(*, *) "Starting simulation loop... "

    HOUR = HOURLAG ! First reading of velocity fields from netcdf file
    call domain%read(UT, VT, WWT, KHNC1, ET, ST1, TT1, RT1) ! read physical fields
    call domain%geo(VX, VY, NV, H)

    do ii = 1, domain%nelements
        length1(ii) = sqrt(  (vx(nv(ii, 1)) - vx(nv(ii, 2)) )**(2.0_sp) + ( vy(nv(ii, 1)) - vy(nv(ii, 2)) )**(2.0_sp)  )
        length2(ii) = sqrt(  (vx(nv(ii, 2)) - vx(nv(ii, 3)) )**(2.0_sp) + ( vy(nv(ii, 2)) - vy(nv(ii, 3)) )**(2.0_sp)  )
        length3(ii) = sqrt(  (vx(nv(ii, 3)) - vx(nv(ii, 1)) )**(2.0_sp) + ( vy(nv(ii, 3)) - vy(nv(ii, 1)) )**(2.0_sp)  )
        average_thickness = abs(sum( H(nv(ii,1:3)) )/3.0_sp) ! must be positive
    end do

    ss(:) = 0.5_sp*(length1(:) + length2(:)+ length3(:))
    domain%elementArea(:) = sqrt( ss * (ss(:) - length1(:)) * (ss(:) - length1(:)) * (ss(:) - length3(:)) ) ! herons's formula for area
    domain%meshArea = sum(domain%elementArea(:))
    domain%elementSigmaVolume(:) = domain%elementArea(:)*average_thickness(:)
    domain%layerDepth = average_thickness(1)
    domain%layerSigma = float(KBM1)**(-1.0_SP)

    HOUR = HOUR + 1

    open(unit=iotox, file=trim(folderprefix)//"/dissolved_toxin.dat", status='replace')
    IINT = 0
    do NH = ISLAG, IELAG ! timestep units are hours, but not necessarily whole numbers
        write(*,*);
        write(*, "(I4,A,I4,A)", advance='no') NH-ISLAG+1, ' / ', IELAG-ISLAG+1, ' steps: '
        call domain%read(UNC2, VNC2, WNC2, KHNC2, ELNC2, SALNC2, TEMPNC2, RHONC2)

        HOUR = HOUR + 1
        ! length (float) of flow field interpolation divided by length (float) of inner time step
        do IT = 1, int(INSTP/DTI)

            write(*, "(A)", advance='no') "|"
            IINT = IINT + 1 ! increment total step count
            LAG_TIME = float(IINT)*DTI + float(ISLAG) ! decimal simulation time for output timestamp
            domain%time = LAG_TIME
            domain%daytime = LAG_TIME / 24.0_sp
            domain%clocktime = LAG_TIME - 24.0_sp*floor(domain%daytime)

            if ((domain%clocktime > 6.0_sp) .and. (domain%clocktime < 18.0_sp)) then
                domain%globalIrradiance = 0.5_SP*irradSurf*(1.0_SP + cos(2.0_SP*pi2*domain%daytime))
            else
                domain%globalIrradiance = zero
            end if

            call domain%diffuse() ! vertical diffusion of dissolved toxin
            call particles%movement()
            call particles%random()
            if ( mod(IINT, int(DTOUT/DTI) ) == 0) then
                particles%XPT(:) = particles%XP(:) ! change back to initial coordinate system for output
                particles%YPT(:) = particles%YP(:)
                call particles%writePosition(iocp) ! output position, same for all particle types
                call particles%writeState(iocs)
                write(iotox,"(1f20.3,51F20.6)") domain%time, domain%verticaltox(1:KB)
            end if
        end do
    end do

    write(*,*); write(*,*); write(*,*) "Simulation finished."; write(*,*)
    call random%stats()

end program