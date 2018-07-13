program PARTICLE_TRAJ
  
  use ALL_VARS, only : casename, folderprefix
  use MOD_SIM, only : domain
  use LIMS, only : N, M, KB
  implicit none
  character(len = 100) :: input, filename, meshfile, foldername, exp_letter
  logical :: fexist
  integer :: exp_type

  call get_command_argument(1, input) ! Import casename from command line
  if (len_trim(input) .eq. 0) then
     print *, 'Please provide casename on command line'; print *, 'Stopping...'; stop
  end if
  
  call get_command_argument(2, foldername) ! Import casename from command line
  if (len_trim(foldername) .eq. 0) then
     print *, 'Please provide simulation ID on command line (###)'; print *, 'Stopping...'; stop
  end if
  folderprefix = adjustl(foldername)
  casename = adjustl(input)
  inquire(file="./"//trim(folderprefix)//"/"//trim(casename)//"_run.dat", exist = fexist)
  if (.not. fexist) then
     write(*, *) 'Run file ', filename, ' does not exist'; write(*, *) 'Stopping...'; stop
  end if
  
  call get_command_argument(3, exp_letter) ! Import casename from command line
  if (len_trim(exp_letter) .eq. 0) then
     print *, 'Please provide experiment type on command line'; print *, 'Stopping...'; stop
  else if (len_trim(exp_letter) .gt. 1) then
     print *, 'Unrecognized experiment type'; print *, 'Stopping...'; stop
  end if
  if (exp_letter .eq. 'A') then
    exp_type = 1
  else if (exp_letter .eq. 'B') then
    exp_type = 2
  else if (exp_letter .eq. 'C') then
    exp_type = 3 
  else if (exp_letter .eq. 'D') then
    exp_type = 4
  else
    print *, 'Unrecognized experiment type'; print *, 'Stopping...'; stop
  end if

  write(*, *); write(*, "(A)", advance='no') 'Importing simulation parameters... '
  call DATA_RUN ! Read parameters controlling model run
  write(*, *) 'Finished'
  
  ! Determine number of elements and nodes in the model
  meshfile = "./"//trim(folderprefix)//"/"//"mesh_elem.dat"
  inquire(file=trim(meshfile), exist=fexist)
  if (.not. fexist) then
     write(*, *) 'Mesh file ', meshfile, ' does not exist'; write(*, *) 'Stopping...'; stop
  end if
  
  allocate(domain)
  call read_mesh(domain) ! set up grid metrics, allocate mesh variables, and read in mesh data
  write(*, *)
  write(*, *) '    Nodes        :', M
  write(*, *) '    Elements     :', N
  write(*, *) '    Sigma layers :', KB
  write(*, *)
  
  call domain%init(exp_type) ! allocate and initialize global environmental variables, and additional mesh-based variables
  write(*, "(A)", advance='no') "Computing mesh topology... "
  call TRIANGLE_GRID_EDGE
  write(*, *) "Finished"
  call SET_LAG(exp_type) ! set up lagrangian particles
  write(*, *) "Starting simulation loop... "
  call LAG_UPDATE ! time step evolution

end program PARTICLE_TRAJ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SET_LAG(exp_type)
  ! Read in lagrangian control parameters and initial lagrangian positions
  use ALL_VARS
  use MOD_TOX, only : cyano
  use MOD_FISH, only : fish
  use MOD_SIM, only : domain
  use MOD_RAND, only : random
  use parameters, only : iocp, iocs, iophys, iofp, iofs
  implicit none

  integer, intent(in) :: exp_type
  character(LEN = 100) :: input_file, state_format
  integer , allocatable, dimension(:) :: INWATER
  integer :: HOUR
  real(sp), dimension(0:N, KB) :: UNC, VNC, WNC ! velocity fields, as read from NetCDF
  real(sp), dimension(0:M, KB) :: KHNC, SALNC, TEMPNC, RHONC ! KH field, as read from NetCDF
  real(sp), dimension(0:M) :: ELNC ! Free surface height field, as read from NetCDF
  
  state_format="(1F12.6, 9000(I6,3F12.6))"
  
  write(*, "(A)", advance='no') "Allocating cyanobacteria and fish structures... "
  allocate(cyano)
  allocate(fish)
  allocate(random); call random%init()
  write(*, *) "Finished"

  write(*, "(A)", advance='no') "Loading physical field data... "
  input_file = "./" //trim(folderprefix)//"/"// trim(casename) // "_phys.dat" ! NetCDF file for simulation
  open(unit=iophys, file=input_file, status='old', position='rewind')
  call domain%load(UNC, VNC, WNC, KHNC, ELNC, SALNC, TEMPNC, RHONC) ! use steady state
  write(*, *) "Finished"
  
  HOUR = HOURLAG ! start time (int) is 0
  ISLAG = HOURLAG ! tracking begin iteration is 0
  IELAG = ISLAG + TDRIFT - 1 ! tracking end iteration
  
  write(*, "(A)", advance='no') "Initiatializing particle structures... "
  call cyano%init(exp_type) ! initialize and allocate type specific structures, also reads position and state variables from file
  call fish%init() !
  write(*, *) "Finished"
  write(*, "(A)", advance='no') "Allocating common variables... "
  call cyano%lag_alloc() ! allocate common variables other than position and itag
  call fish%lag_alloc()
  
  
  allocate( INWATER(cyano%ndrft) ); INWATER(:) = 1 ! allocate inwater flag for only single species
  write(*, *) "Finished"

  cyano%XP(:) = cyano%XPT(:) - VXMIN ! Shift x to model coordinate system
  cyano%YP(:) = cyano%YPT(:) - VYMIN ! Shift y to model coordinate system

  write(*, "(A)", advance='no') "Finding host elements... "
  call cyano%FHE_ROBUST(cyano%XP, cyano%YP, INWATER) ! Determine element containing each particle
  where (cyano%FOUND .eq. 0) cyano%INDOMAIN(:) = 0 ! if not found, particle is not in domain and will not be tracked
  write(*, *) "Finished"
  
  write(*, "(A)", advance='no') "Interpolating physical fields... "
  call cyano%INTERP_ELH(cyano%XP, cyano%YP, H, ELNC, 1) ! interpolate elevation and bathymetry
  call cyano%INTERP_FIELDS(cyano%XP, cyano%YP, cyano%ZP, SALNC, TEMPNC, RHONC, 0) ! interpolate salinity and temperature
  write(*, *) "Finished"

  write(*, "(A)", advance='no') "Adjusting vertical domain... "
  cyano%ZPT(:) = -1.0_sp*abs(cyano%ZPT(:)) ! make depth negative
  cyano%ZP(:) = cyano%sigma(cyano%ZPT(:)) ! convert to sigma coordinate
  cyano%LAYER(:) = cyano%zlocate(cyano%ZP(:)) ! valid when sigma layers are equal thickness

  write(*, "(A)", advance='no') "Writing position and state variables to file... "
  open(unit=iocp, file="./"//trim(folderprefix)//"/"//"cyanobacteria_position.dat", status='replace') ! create new position file
  open(unit=iocs, file="./"//trim(folderprefix)//"/"//"cyanobacteria_state.dat", status='replace') ! create new state variable file
  call cyano%writeState(iocs)
  call cyano%writePosition(iocp) ! write particle positions to output file
  write(*, *) "Finished"

  deallocate(INWATER) ! needs to be resized for new particle structure
  
  allocate( INWATER(fish%ndrft) ); INWATER(:) = 1 ! allocate inwater flag for only single species

  fish%XP(:) = fish%XPT(:) - VXMIN ! Shift x to model coordinate system
  fish%YP(:) = fish%YPT(:) - VYMIN ! Shift y to model coordinate system

  write(*, "(A)", advance='no') "Finding host elements... "
  call fish%FHE_ROBUST(fish%XP, fish%YP, INWATER) ! Determine element containing each particle
  where (fish%FOUND .eq. 0) fish%INDOMAIN(:) = 0 ! if not found, particle is not in domain and will not be tracked
  write(*, *) "Finished"
  
  write(*, "(A)", advance='no') "Interpolating physical fields... "
  call fish%INTERP_ELH(fish%XP, fish%YP, H, ELNC, 1) ! interpolate elevation and bathymetry
  call fish%INTERP_FIELDS(fish%XP, fish%YP, fish%ZP, SALNC, TEMPNC, RHONC, 0) ! interpolate salinity and temperature
  write(*, *) "Finished"

  write(*, "(A)", advance='no') "Adjusting vertical domain... "
  fish%ZPT(:) = -1.0_sp*abs(fish%ZPT(:)) ! make depth negative
  fish%ZP(:) = fish%sigma(fish%ZPT(:)) ! convert to sigma coordinate
  fish%LAYER(:) = fish%zlocate(fish%ZP(:)) ! valid when sigma layers are equal thickness

  write(*, "(A)", advance='no') "Writing position and state variables to file... "
  open(unit=iofp, file="./"//trim(folderprefix)//"/"//"fish_position.dat", status='replace') ! create new position file
  open(unit=iofs, file="./"//trim(folderprefix)//"/"//"fish_state.dat", status='replace') ! create new state variable file
  call fish%writeState(iofs)
  call fish%writePosition(iofp) ! write particle positions to output file
  write(*, *) "Finished"

  deallocate(INWATER) ! needs to be resized for new particle structure
  
  write(*, *) ! Print particle statistics
  write(*, *) '    Tracking Info'
  write(*, *) '        Start iteration :', ISLAG
  write(*, *) '        Final iteration :', IELAG
  write(*, *)
  call cyano%stats
  call fish%stats

end subroutine SET_LAG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LAG_UPDATE
  ! Update particle positions, calculate scalar fields and particle velocities
  use ALL_VARS
  use MOD_TOX, only : cyano
  use MOD_FISH, only : fish
  use MOD_SIM, only : domain
  use MOD_RAND, only : random
  use parameters
  
  implicit none
  
  real(sp), dimension(0:N, KB) :: UNC1, UNC2, VNC1, VNC2, WNC1, WNC2 ! velocity fields start and end of hour
  real(sp), dimension(0:M, KB) :: KHNC1, KHNC2, SALNC1, SALNC2, TEMPNC1, TEMPNC2, RHONC1, RHONC2 ! diffusion and physical fields start and end of hour
  real(sp), dimension(0:M) :: ELNC1, ELNC2 ! free surface height field start and end of hour
  real(sp) :: TMP1, TMP2, LAG_TIME
  integer :: NH, I1, I2, IT, HOUR, IINT
  character(len = 100) :: input_file
  

  HOUR = HOURLAG ! First reading of velocity fields from netcdf file
  input_file = "./" //trim(folderprefix)//"/"// trim(casename) // "_phys.dat" ! NetCDF file for simulation
  call domain%load(UNC1, VNC1, WNC1, KHNC1, ELNC1, SALNC1, TEMPNC1, RHONC1) ! read physical fields
  call domain%geo(VX, VY, NV, H)    
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
  simulation_loop: do NH = ISLAG, IELAG ! timestep units are hours, but not necessarily whole numbers
    write(*,*); write(*, "(I4,A,I4,A)", advance='no') NH-ISLAG+1, ' / ', IELAG-ISLAG+1, ' steps: '
    call domain%load(UNC2, VNC2, WNC2, KHNC2, ELNC2, SALNC2, TEMPNC2, RHONC2) ! use steady state
    
    HOUR = HOUR + 1
    I1 = 1
    I2 = int(INSTP/DTI) ! length (float) of flow field interpolation divided by length (float) of inner time step
    interpolation_loop: do IT = I1, I2
    
      write(*,"(A)",advance='no') "|"
      IINT = IINT + 1 ! increment total step count
      LAG_TIME = float(IINT)*DTI + float(ISLAG) ! decimal simulation time for output timestamp
      domain%time = LAG_TIME
      domain%daytime = LAG_TIME / 24.0_sp
      domain%clocktime = LAG_TIME - 24.0_sp*floor(domain%daytime)
      
      if ((domain%clocktime .gt. 6.0_sp) .and. (domain%clocktime .lt. 18.0_sp)) then
        domain%globalIrradiance = 0.5_SP*irradSurf*(1.0_SP + cos(2.0_SP*pi2*domain%daytime))
      else
        domain%globalIrradiance = ZERO
      end if
      
      if (includeAlgae) then
        call domain%vdiff() ! vertical diffusion of dissolved toxin
        call cyano%movement() ! uses runge-kutta integration
        call cyano%random() ! random walk
        if ( mod(IINT, int(DTOUT/DTI) ) .eq. 0) then
          cyano%XPT(:) = cyano%XP(:) + VXMIN ! change back to initial coordinate system for output
          cyano%YPT(:) = cyano%YP(:) + VYMIN
          call cyano%writePosition(iocp) ! output position, same for all particle types
          call cyano%writeState(iocs)
          write(iotox,"(1f20.3,51F20.6)") domain%time, domain%verticaltox(1:KB)
        end if
      end if
      
      if (includeFish) then
        call fish%movement() ! uses runge-kutta integration
        if ( mod(IINT, int(DTOUT/DTI) ) .eq. 0) then
          fish%XPT(:) = fish%XP(:) + VXMIN ! change back to initial coordinate system for output
          fish%YPT(:) = fish%YP(:) + VYMIN
          call fish%writePosition(iofp) ! output position, same for all particle types
          call fish%writeState(iofs) ! output state variables
        end if
      end if

    end do interpolation_loop
  end do simulation_loop

  write(*,*); write(*,*); write(*,*) "Simulation finished."; write(*,*)
  call random%stats()
  
end subroutine LAG_UPDATE