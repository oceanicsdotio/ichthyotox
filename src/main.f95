program main

  use variables
  use random, only : random_number_generator
  use simulation, only : domain, topology, getInteger, getString, find_key
  use cyanobacteria, only : CyanobacteriaAgent
  use fish, only : FishAgent

  implicit none

  type :: Container
    type(CyanobacteriaAgent) :: cyanobacteria
    type(FishAgent) :: fish
  end type;
  type(Container) :: agent


  logical, parameter :: &
    & includeAlgae = .true., &
    & includeFish =  .true.

  real(sp), dimension(0:N, KB) :: UNC, UNC1, UNC2, VNC, VNC1, VNC2, WNC, WNC1, WNC2 ! velocity fields start and end of hour
  real(sp), dimension(0:M, KB) :: KHNC, KHNC1, KHNC2, SALNC, SALNC1, SALNC2, TEMPNC, TEMPNC1, TEMPNC2, RHONC, RHONC1, RHONC2 ! diffusion and physical fields start and end of hour
  real(sp), dimension(0:M) :: ELNC, ELNC1, ELNC2 ! free surface height field read, start and end of hour
  real(sp) :: TMP1, TMP2, LAG_TIME

  integer :: NH, I1, I2, IT, HOUR, IINT, ii, index, ionode=100, ioelem=101, exp_type, NCT, ISCAN
  character(len = 100) :: input_file, input, meshfile, foldername, exp_letter, state_format
  logical :: fileExists

  character(len = 120) :: filename

  call get_command_argument(1, input) ! Import casename from command line
  if (len_trim(input) == 0) then
    print *, 'Provide casename as argument 1.'; stop
  end if

  call get_command_argument(2, foldername) ! Import casename from command line
  if (len_trim(foldername) == 0) then
    print *, 'Provide simulation ID (###) as argument 2'; stop
  end if

  folderprefix = adjustl(foldername)
  casename = adjustl(input)

  inquire(file="./"//trim(folderprefix)//"/"//trim(casename)//"_run.dat", exist=fileExists)
  if (.not. fileExists) then
    write(*, *) 'Run file ', filename, ' does not exist'; stop
  end if

  call get_command_argument(3, exp_letter) ! Import casename from command line
  if (len_trim(exp_letter) == 0) then
    print *, 'Provide experiment type on command line'; stop
  else if (len_trim(exp_letter) > 1) then
    print *, 'Unrecognized experiment type'; stop
  end if
  if (exp_letter == 'A') then
    exp_type = 1
  else if (exp_letter == 'B') then
    exp_type = 2
  else if (exp_letter == 'C') then
    exp_type = 3
  else if (exp_letter == 'D') then
    exp_type = 4
  else
    print *, 'Unrecognized experiment type'; stop
  end if

  write(*, *);
  write(*, "(A)", advance='no') 'Importing simulation parameters... '


  ! Read in variables and set values
  filename = "./"//trim(folderprefix)//"/"//trim(CASENAME)//"_run.dat"

  call getInteger(filename, "HOURLAG", HOURLAG) ! day
  call getInteger(filename, "TDRIFT", TDRIFT) ! Total time to move drifters (TDRIFT)

  ! External time step (DTI)
  ISCAN = find_key(trim(filename), "DTI", FSCAL = DTI)
  if (ISCAN /= 0) then
    write (*, *) 'ERROR READING DTI: ', ISCAN
    stop
  end if

  ! Input time step of flow fields (instp)
  ISCAN = find_key(trim(filename),"INSTP", FSCAL = INSTP)
  if (ISCAN /= 0) then
    write(*, *) 'ERROR READING INSTP: ', ISCAN
    stop
  end if

  ! External time step (DTOUT)
  ISCAN = find_key(trim(filename), "DTOUT", FSCAL = DTOUT)
  if (ISCAN /= 0) then
    write (*, *) 'ERROR READING DTOUT: ', ISCAN
    stop
  end if


  ! Horizontal diffusion coefficient (DHOR)
  ISCAN = find_key(trim(filename), "DHOR", FSCAL = DHOR)
  if (ISCAN /= 0) then
    write(*, *) 'ERROR READING DHOR: ', ISCAN
    stop
  end if

  ! Random walk time step
  ISCAN = find_key(trim(filename), "DTRW", FSCAL = DTRW)
  if (ISCAN /= 0) then
    write(*, *) 'ERROR READING DTRW: ', ISCAN
    stop
  end if

  ! Determine number of elements and nodes in the model
  write(*, "(A)", advance='no') "Reading mesh files... "
  meshfile = "./"//trim(folderprefix)//"/"//"mesh_elem.dat"
  inquire(file=trim(meshfile), exist=fileExists)
  if (.not. fileExists) then
    write(*, *) 'Mesh file ', meshfile, ' does not exist'; stop
  end if

  allocate(domain)

  open(ioelem, file = "./"//trim(folderprefix)//"/mesh_elem.dat")
  open(ionode, file = "./"//trim(folderprefix)//"/mesh_node.dat")

  read(ioelem, *) N, KB ! get number of nodes, elements, and sigma layers
  read(ionode, *) M


  domain%nnodes = M; domain%nelements = N; domain%nlayers = KB
  KBM1 = KB - 1
  KBM2 = KB - 2

  write(*, "(A)", advance='no') "Allocating mesh based variables... "

  real(sp), dimension(0:domain%nelements) :: length1, length2, length3, ss, average_thickness
  NCT = N*3

  ! Grid Metrics
  allocate(XC(0:N))            ;XC   = zero   ! X-COORD AT FACE CENTER
  allocate(YC(0:N))            ;YC   = zero   ! Y-COORD AT FACE CENTER
  allocate(VX(0:M))            ;VX   = zero   ! X-COORD AT GRID POINT
  allocate(VY(0:M))            ;VY   = zero   ! Y-COORD AT GRID POINT

  ! Node, Boundary Condition, and Control Volume
  allocate(NV(0:N,4))           ;NV       = 0  ! NODE NUMBERING FOR ELEMENTS
  allocate(NBE(0:N,3))          ;NBE      = 0  ! INDICES OF ELEMENT NEIGHBORS
  allocate(NTVE(0:M))           ;NTVE     = 0
  allocate(ISONB(0:M))          ;ISONB    = 0  ! NODE MARKER = 0,1,2
  allocate(ISBCE(0:N))          ;ISBCE    = 0

  ! 1-d arrays for the sigma coordinate
  allocate(Z(KB))               ; Z      = zero    ! SIGMA COORDINATE VALUE
  allocate(ZZ(KB))              ; ZZ     = zero    ! INTRA LEVEL SIGMA VALUE
  allocate(DZ(KB))              ; DZ     = zero    ! DELTA-SIGMA VALUE

  ! 2-d flow variable arrays at nodes
  allocate(H(0:M))       ;H    = zero       ! BATHYMETRIC DEPTH
  allocate(D(0:M))       ;D    = zero       ! DEPTH
  allocate(EL(0:M))      ;EL   = zero       ! SURFACE ELEVATION
  allocate(ET(0:M))      ;ET  = zero       ! SURFACE ELEVATION PREVIOUS TIMESTEP

  ! internal mode arrays-(element based)
  allocate(U(0:N, KB))       ;U     = zero   ! X-VELOCITY
  allocate(V(0:N, KB))       ;V     = zero   ! Y-VELOCITY
  allocate(WW(0:N, KB))      ;WW    = zero   ! Z-VELOCITY
  allocate(UT(0:N, KB))      ;UT    = zero   ! X-VELOCITY FROM PREVIOUS TIMESTEP
  allocate(VT(0:N, KB))      ;VT    = zero   ! Y-VELOCITY FROM PREVIOUS TIMESTEP
  allocate(WWT(0:N, KB))     ;WWT   = zero   ! Z-VELOCITY FROM PREVIOUS TIMESTEP
  allocate(KH(0:N, KB))     ;KH    = zero   ! TURBULENT QUANTITY

  ! 3d variable arrays-(node based)
  allocate(T1(0:M, KB))       ;T1     = zero  ! TEMPERATURE AT NODES
  allocate(S1(0:M, KB))       ;S1     = zero  ! SALINITY AT NODES
  allocate(R1(0:M, KB))       ;R1   = zero  ! DENSITY AT NODES
  allocate(TT1(0:M, KB))      ;TT1    = zero  ! TEMPERATURE FROM PREVIOUS TIME
  allocate(ST1(0:M, KB))      ;ST1    = zero  ! SALINITY FROM PREVIOUS TIME
  allocate(RT1(0:M, KB))      ;RT1 = zero

  ! Shape coefficient arrays and control volume metrics
  allocate(A1U(0:N, 4))         ;A1U   = zero
  allocate(A2U(0:N, 4))         ;A2U   = zero
  allocate(AWX(0:N, 3))         ;AWX   = zero
  allocate(AWY(0:N, 3))         ;AWY   = zero
  allocate(AW0(0:N, 3))         ;AW0   = zero

  write(*, *) "Finished"

  write(*, "(A)", advance='no') "Getting element data... "
  element_loop: do ii = 1, N
    read(ioelem, *) index, NV(ii, 1), NV(ii, 3), NV(ii, 2) ! get node indices
  end do element_loop
  NV(:, 4) = NV(:, 1) ! duplicate node for computation
  write(*, *) "Finished"

  write(*, "(A)", advance='no') "Getting node data... "
  node_loop: do ii = 1, M
    read(ionode, *) index, VX(ii), VY(ii), H(ii) ! get node position and depth
  end do node_loop
  write(*, *) "Finished"

  close(ioelem)
  close(ionode)

  write(*, *)
  write(*, *) '    Nodes        :', M
  write(*, *) '    Elements     :', N
  write(*, *) '    Sigma layers :', KB
  write(*, *)

  call domain%init(exp_type) ! allocate and initialize global environmental variables, and additional mesh-based variables
  
  write(*, "(A)", advance='no') "Computing mesh topology... "
  call topology
  write(*, *) "Finished"

  state_format="(1F12.6, 9000(I6,3F12.6))"

  write(*, "(A)", advance='no') "Allocating cyanobacteria and fish structures... "
  allocate(cyano)
  allocate(fish_agent)
  allocate(random); call random%init()
  write(*, *) "Finished"

  write(*, "(A)", advance='no') "Loading physical field data... "
  input_file = "./" //trim(folderprefix)//"/"// trim(casename) // "_phys.dat" ! NetCDF file for simulation
  open(unit=iophys, file=input_file, status='old', position='rewind')
  call domain%read(UNC, VNC, WNC, KHNC, ELNC, SALNC, TEMPNC, RHONC)
  write(*, *) "Finished"

  HOUR = HOURLAG ! start time (int) is 0
  ISLAG = HOURLAG ! tracking begin iteration is 0
  IELAG = ISLAG + TDRIFT - 1 ! tracking end iteration

  write(*, "(A)", advance='no') "Initiatializing particle structures... "
  call agent%cyanobacteria%init(exp_type) ! initialize and allocate type specific structures, also reads position and state variables from file
  call agent%fish%init() !
  write(*, *) "Finished"
  write(*, "(A)", advance='no') "Allocating common variables... "
  call agent%cyanobacteria%lag_alloc() ! allocate common variables other than position
  call agent%fish%lag_alloc()

  write(*, *) "Finished"

  write(*, "(A)", advance='no') "Finding host elements... "
  call agent%cyanobacteria%find_host_element(agent%cyanobacteria%XP, agent%cyanobacteria%YP) ! Determine element containing each particlewrite(*, *) "Finished"

  write(*, "(A)", advance='no') "Interpolating physical fields... "
  call agent%cyanobacteria%INTERP_ELH(agent%cyanobacteria%XP, agent%cyanobacteria%YP, H, ELNC, 1) ! interpolate elevation and bathymetry
  call agent%cyanobacteria%INTERP_FIELDS(agent%cyanobacteria%XP, agent%cyanobacteria%YP, agent%cyanobacteria%ZP, SALNC, TEMPNC, RHONC, 0) ! interpolate salinity and temperature
  write(*, *) "Finished"

  write(*, "(A)", advance='no') "Adjusting vertical domain... "
  agent%cyanobacteria%ZP(:) = -abs(agent%cyanobacteria%ZP(:)) ! make depth negative
  agent%cyanobacteria%LAYER(:) = agent%cyanobacteria%zlocate(agent%cyanobacteria%ZP(:)) ! valid when sigma layers are equal thickness

  write(*, "(A)", advance='no') "Writing position and state variables to file... "
  open(unit=iocp, file="./"//trim(folderprefix)//"/"//"cyanobacteria_position.dat", status='replace') ! create new position file
  open(unit=iocs, file="./"//trim(folderprefix)//"/"//"cyanobacteria_state.dat", status='replace') ! create new state variable file
  call agent%cyanobacteria%writeState(iocs)
  call agent%cyanobacteria%writePosition(iocp) ! write particle positions to output file
  write(*, *) "Finished"

  write(*, "(A)", advance='no') "Finding host elements... "
  call agent%fish%find_host_element(agent%fish%XP, agent%fish%YP) ! Determine element containing each particle
  write(*, *) "Finished"

  write(*, "(A)", advance='no') "Interpolating physical fields... "
  call agent%fish%INTERP_ELH(agent%fish%XP, agent%fish%YP, H, ELNC, 1) ! interpolate elevation and bathymetry
  call agent%fish%INTERP_FIELDS(agent%fish%XP, agent%fish%YP, agent%fish%ZP, SALNC, TEMPNC, RHONC, 0) ! interpolate salinity and temperature
  write(*, *) "Finished"

  write(*, "(A)", advance='no') "Adjusting vertical domain... "
  agent%fish%ZP(:) = -1.0_sp*abs(agent%fish%ZP(:)) ! make depth negative
  agent%fish%LAYER(:) = agent%fish%zlocate(agent%fish%ZP(:)) ! valid when sigma layers are equal thickness

  write(*, "(A)", advance='no') "Writing position and state variables to file... "
  open(unit=iofp, file="./"//trim(folderprefix)//"/"//"fish_position.dat", status='replace') ! create new position file
  open(unit=iofs, file="./"//trim(folderprefix)//"/"//"fish_state.dat", status='replace') ! create new state variable file
  call agent%fish%writeState(iofs)
  call agent%fish%writePosition(iofp) ! write particle positions to output file
  write(*, *) "Finished"

  write(*, *) ! Print particle statistics
  write(*, *) '    Tracking Info'
  write(*, *) '        Start iteration :', ISLAG
  write(*, *) '        Final iteration :', IELAG
  write(*, *)
  call agent%cyanobacteria%stats
  call agent%fish%stats

  write(*, *) "Starting simulation loop... "

  HOUR = HOURLAG ! First reading of velocity fields from netcdf file
  input_file = "./" //trim(folderprefix)//"/"// trim(casename) // "_phys.dat" ! NetCDF file for simulation
  call domain%read(UNC1, VNC1, WNC1, KHNC1, ELNC1, SALNC1, TEMPNC1, RHONC1) ! read physical fields


  call domain%geo(VX, VY, NV, H)


  ! subroutine simulation_geometry(self, vertx, verty, node_indices, bathymetry)
  ! calculates area of any triangular mesh or subregion,

  integer :: ii

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

  open(unit=iotox, file="./"//trim(folderprefix)//"/"//"dissolved_toxin.dat", status='replace')

  UT(:,:) = UNC1(:,:) ! assign current u velocity field
  VT(:,:) = VNC1(:,:) ! assign current v velocity field
  WWT(:,:) = WNC1(:,:) ! assign current w velocity field
  ET(:) = ELNC1(:) ! assign current free surface height field
  ST1(:,:) = SALNC1(:,:) !---fish change 10, assign current salinity field
  TT1(:,:) = TEMPNC1(:,:) !--Keeney change, copy temperature field to working array
  RT1(:,:) = RHONC1(:,:) ! copy density field to working array

  IINT = 0
  do NH = ISLAG, IELAG ! timestep units are hours, but not necessarily whole numbers
    write(*,*)
    write(*, "(I4,A,I4,A)", advance='no') NH-ISLAG+1, ' / ', IELAG-ISLAG+1, ' steps: '
    call domain%read(UNC2, VNC2, WNC2, KHNC2, ELNC2, SALNC2, TEMPNC2, RHONC2)

    HOUR = HOUR + 1
    I1 = 1
    I2 = int(INSTP/DTI) ! length (float) of flow field interpolation divided by length (float) of inner time step
    do IT = I1, I2

      write(*,"(A)",advance='no') "|"
      IINT = IINT + 1 ! increment total step count
      LAG_TIME = float(IINT)*DTI + float(ISLAG) ! decimal simulation time for output timestamp
      domain%time = LAG_TIME
      domain%daytime = LAG_TIME / 24.0_sp
      domain%clocktime = LAG_TIME - 24.0_sp*floor(domain%daytime)

      if ((domain%clocktime > 6.0_sp) .and. (domain%clocktime < 18.0_sp)) then
        domain%globalIrradiance = 0.5_SP*irradSurf*(1.0_SP + cos(2.0_SP*pi2*domain%daytime))
      else
        domain%globalIrradiance = ZERO
      end if

      if (includeAlgae) then
        call domain%diffuse() ! vertical diffusion of dissolved toxin
        call agent%cyanobacteria%movement() ! uses runge-kutta integration
        call agent%cyanobacteria%random() ! random walk
        if ( mod(IINT, int(DTOUT/DTI) ) == 0) then
          call agent%cyanobacteria%writePosition(iocp) ! output position, same for all particle types
          call agent%cyanobacteria%writeState(iocs)
          write(iotox,"(1f20.3,51F20.6)") domain%time, domain%verticaltox(1:KB)
        end if
      end if

      if (includeFish) then
        call agent%fish%movement() ! uses runge-kutta integration
        if ( mod(IINT, int(DTOUT/DTI) ) == 0) then
          call agent%fish%writePosition(iofp) ! output position, same for all particle types
          call agent%fish%writeState(iofs) ! output state variables
        end if
      end if

    end do
  end do

  write(*,*); write(*,*); write(*,*) "Simulation finished."; write(*,*)
  call random%stats()

end program