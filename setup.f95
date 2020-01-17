program setup
  ! program writes initial position and state variables to simulation file for reading by offlag
  use ALL_VARS, only : zero, folderprefix, sp
  use MOD_TOX, only : colonyBaseRadius, densityMin, densitymax, cellFrac, vesicleFrac, vesicleDensity, cellDensityCoefficient
  use ALL_VARS, only : iorun, iophys
  use MOD_RAND, only : random

  implicit none
  
  logical :: makeGrid = .true.
  logical :: writeFields = .true. 
  integer :: fid=666, nnodes, nelements, node_index, element_index, start_pattern, alternate_pattern, x_div, y_div, ii, jj, nlayers
  integer :: ioelem=500, ionode=501, fishid
  character(len=50) :: kb_format, foldername
  real(sp) :: div_inc, x_range, y_range, time, end_time, layerdepth
  real(sp), dimension(:,:), allocatable :: node_pos
  real(sp), dimension(:), allocatable :: diffusivity, temperature, rho_pure, density
  real(sp), dimension(:), allocatable :: cyano_id, cyano_x, cyano_y, cyano_z, cyano_carb, cyano_pro, cyano_mc, cyano_radius
  integer, dimension(:,:), allocatable :: vertices
  real(sp) :: dt = 1.0, day, clock_time, biomass, carbon, carbohydrate, protein, microcystin, ratio, radius, reqCellDensity, initial_position
  real(sp) :: bottom_depth=-5.0_SP, ini_depth
  real(sp) :: salinity=3.0_SP

  !real(sp) :: control_temp=20.0_SP, temp_slope=0.0_SP  ! Experiment A
  real(sp) :: control_temp=20.0_SP, temp_slope=0.00694_SP  ! Experiment B
  !real(sp) :: control_temp=25.0_SP, temp_slope=0.00694_SP  ! Experiment C
  !real(sp) :: control_temp=30.0_SP, temp_slope=0.0_SP  ! Experiment D
  
  integer :: ncolony, nfish, step, x_nodes, y_nodes, ndays = 30
  
  call get_command_argument(1, foldername) ! Import casename from command line
  if (len_trim(foldername) == 0) then
     print *, 'Please provide simulation ID on command line (###)'; print *, 'Stopping...'; stop
  end if
  folderprefix = adjustl(foldername)

  allocate(random);
  call random%init()

  if (makeGrid) then
    write(*, *); write(*, *) "For the simulation we are going to build a simple triangular mesh representing a rectangular reservoir... "
    !write(*, *); write(*, '(A)', advance='no') "    What is the desired X range (meters)? "; read(*, *) x_range
    x_range = 500.0_SP
    y_range = 500.0_SP
    div_inc = 10.0_SP
    nlayers = 26
    !write(*, *); write(*, '(A)', advance='no') "    What is the desired Y range (meters)? "; read(*, *) y_range
    !write(*, *); write(*, '(A)', advance='no') "    About what size should are the elements (meters)? "; read(*, *) div_inc
    !write(*, *); write(*, '(A)', advance='no') "    How many sigma layers are there? "; read(*, *) nlayers
  
    x_div = floor( x_range / div_inc )
    y_div = floor( y_range / div_inc )
    x_nodes = x_div + 1
    y_nodes = y_div + 1
    
    nnodes = (x_div + 1) * (y_div + 1)
    nelements = 2*x_div*y_div
    
    node_index = 1
    element_index = 1
    start_pattern = 1
    
    allocate(node_pos(nnodes,3))
    allocate(vertices(nelements,3))
    
    write(*, *); write(*, "(A)", advance='no') "Building mesh and topology... "
    inc_row: do ii = 1, y_nodes
      alternate_pattern = start_pattern
      inc_column: do jj = 1, x_nodes
        node_pos(node_index, 1) = float(jj-1)*div_inc
        node_pos(node_index, 2) = float(ii-1)*div_inc
        node_pos(node_index, 3) = bottom_depth
        
        if ( (jj < x_nodes) .and. (ii < y_nodes) ) then
          ! bottom row of triangles
          vertices(element_index, 1) = node_index
          vertices(element_index, 2) = node_index + x_nodes ! note index
          vertices(element_index, 3) = node_index + 1
          if (alternate_pattern < 0) vertices(element_index, 2) = vertices(element_index, 2) + 1
          element_index = element_index + 1
          
          ! upper row of triangles
          vertices(element_index, 1) = node_index + x_nodes
          vertices(element_index, 2) = node_index + x_nodes + 1
          vertices(element_index, 3) = node_index + 1
          if (alternate_pattern < 0) vertices(element_index, 3) = vertices(element_index, 3) - 1
          element_index = element_index + 1
    
          alternate_pattern = -alternate_pattern ! switch triangle pattern as you move x=0 to x=x_range
        end if
        
        node_index = node_index + 1 
        
      end do inc_column
      start_pattern = -start_pattern ! switch triangle pattern between rows as grid is built
      
    end do inc_row
    write(*, "(A)") "Finished"
    write(*, *); write(*, "(A)", advance='no') "Writing mesh data to file... "
    
    open(unit=ioelem, file="../"//trim(folderprefix)//"/mesh_elem.dat", status='replace')
    write(*,*) "CHECK"
    write(ioelem, *) nelements, nlayers
    element_loop: do ii = 1, nelements
       write(ioelem, *) ii, vertices(ii, 1), vertices(ii, 2), vertices(ii, 3) ! write node indices
    end do element_loop

    open(unit=ionode, file="../"//trim(folderprefix)//"/mesh_node.dat", status='replace')
    write(ionode, *) nnodes
    node_loop: do ii = 1, nnodes
       write(ionode, *) ii, node_pos(ii, 1), node_pos(ii, 2), node_pos(ii, 3) ! write node position
    end do node_loop
    write(*, "(A)") "Finished "
    
    write(*, *)
    write(*, *) "    Nodes:          ", nnodes
    write(*, *) "    Elements:       ", nelements
  end if

  if (writeFields) then
  
    allocate(diffusivity(nlayers))
    allocate(temperature(nlayers))
    allocate(rho_pure(nlayers))
    allocate(density(nlayers))
      
    write(kb_format, "(A1,I6,A7)") "(", 3*nlayers, "F20.10)"
    write(*, *); write(*, "(A)", advance='no') "Writing physical fields to file... "
    open(unit=iophys, file="../"//trim(folderprefix)//"/ichthyotox_phys.dat", status='replace')
    
    write_loop: do step = 0, 24*ndays
      time = float(step)*dt
      day = time / 24.0_sp
      clock_time = time - floor(day)*24.0_sp
      
      diffusivity(:) = 60.0_SP*60.0_SP*10.0_SP**(-4.0_SP)*0.1
      temperature(:) = control_temp + temp_slope*time

      ! density of pure water from Millero and Poisson, convert to temperature
      density(:) = 999.842594_SP + 0.06793952_SP*temperature(:) - 0.009095290_SP*temperature(:)**(2.0_SP) + 0.0001001685_SP*temperature(:)**(3.0_SP) - (1.120083D-6)*temperature(:)**(4.0_SP) + (6.536332D-9)*temperature(:)**(5.0_SP)
      density(:) = density(:) + salinity*( 0.824493_SP - 0.0040899_SP*temperature(:) + 0.000076438_SP*temperature(:)**(2.0_SP) - (8.2467D-7)*temperature(:)**(3.0_SP) + (5.3875D-9)*temperature(:)**(4.0_SP) )
      density(:) = density(:) + salinity**(1.5_SP)*( -0.00572466_SP + 0.00010227_SP*temperature(:) - (1.6546D-6)*temperature(:)**(2.0_SP) )
      density(:) = density(:) + 0.00048314_SP*salinity**(2.0_SP)
      
      write(iophys, "(1F10.3)", advance='no') time
      write(iophys, kb_format, advance='no') temperature(:), density(:), diffusivity(:)
      write(iophys, *) ! new line
      
    end do write_loop
    close(iophys)
    write(*, "(A)") "Finished"
  end if

  ! write ichthyotox_run.dat
  open(unit=iorun, file="../"//trim(folderprefix)//"/ichthyotox_run.dat", status='replace')
  write(iorun,"(A)") "INFOFILE = screen"
  write(iorun,"(A,F4.2)") "DTI = ", 0.02_SP ! inner interpolation time step, float"
  write(iorun,"(A,F4.2)") "INSTP = ", 1.0_SP ! time step of physical field data, float
  write(iorun,"(A,F4.2)") "DTOUT = ", 0.1_SP ! output time step, >dti
  write(iorun,"(A,F4.2)") "DHOR = ", 0.10_SP ! horizontal diffusion coefficient
  write(iorun,"(A,F4.2)") "DTRW = ", 0.02_SP ! random walk time step
  write(iorun,"(A,I4)") "TDRIFT = ", 24*ndays  ! number of time steps to iterate, int
  write(iorun,"(A,I4)") "YEARLAG = ", 2016
  write(iorun,"(A,I2)") "MONTHLAG = ", 4
  write(iorun,"(A,I2)") "DAYLAG = ", 1
  write(iorun,"(A,I2)") "HOURLAG = ", 0
  write(iorun,"(A,I1)") "IRW = ", 0 ! random walk type
  write(iorun,"(A)") "P_SIGMA = F"
  write(iorun,"(A)") "OUT_SIGMA = F"
  write(iorun,"(A)") "F_DEPTH = F" 
  write(iorun,"(A)") "GEOAREA = box"
  write(iorun,"(A)") "INPDIR=/"
  write(iorun,"(A)") "LAGINI=/"
  write(iorun,"(A)") "OUTDIR=/"

  ! Generate initial values for lagrangian particles and write these to file 
  write(*, *); write(*, *) "We will now generate the initial positions for lagrangian particles..."; 
  !write(*, *); write(*, '(A)', advance = 'no') "    How many cyanobacteria particles should be used? "; read(*, *) ncolony
  ncolony = 100
  write(*, *); write(*, '(A)', advance = 'no') "    What is the total initial protein biomass in grams? "; read(*, *) biomass
  write(*, *); write(*, '(A)', advance = 'no') "    What is the total initial microcystin load in grams? "; read(*, *) microcystin
  cyanobacteria: if (ncolony > 0) then
  
    ! write initial position
    open(unit = fid, file = "../"//trim(folderprefix)//"/cyanobacteria_ini.dat", status = 'replace')
    write(*, *); write(*, *) '        Writing initial positions to ', '../cyanobacteria_ini.dat'
    write(fid, "(I6)") ncolony ! read number of particles
    do ii = 1, ncolony
       initial_position = -1.0*abs(random%uniform()*5.0)
       write(fid, "(I4,3F20.3)") ii, abs(random%uniform()*x_div*div_inc), abs(random%uniform()*y_div*div_inc), initial_position ! write identifier and position for each particle
    end do
    
    ! write state variables
    open(unit = fid, file = "../"//trim(folderprefix)//"/cyanobacteria_var.dat", status = 'replace')
    write(*, *); write(*, *) '        Writing state variables to ', '../cyanobacteria_var.dat'
    write(fid, "(I6)") ncolony
    do ii = 1, ncolony
      radius = colonyBaseRadius + 5.0_SP*10.0_SP**(-6.0_SP)*random%get()
      carbon = biomass/float(ncolony)
      
      reqCellDensity = ((density(1) - (1.0_SP - cellFrac)*(density(1) + 0.7_SP))/cellFrac - vesicleFrac*vesicleDensity)/(1.0_SP - vesicleFrac)
      ratio = log(1.0_SP - (reqCellDensity - densityMin)/(densityMax - densityMin))/(-cellDensityCoefficient)
      protein = carbon
      carbohydrate = protein*ratio
      
      !microcystin = 0.0_sp
      write(fid, "(4F20.6)") radius, carbohydrate, protein, microcystin/float(ncolony)
    end do
    close(fid)
    
  else
    write(*, *); write(*, *) "No cyanobacteria in this simulation..."; write(*, *)
  end if cyanobacteria

  !write(*, *); write(*, *); write(*, '(A)', advance = 'no') "    How many fish particles should be used? "; read(*, *) nfish
  nfish = 200
  if (nfish > 0) then
    ! write initial position
    open(unit = fid, file = "./"//trim(folderprefix)//"/fish_ini.dat", status = 'replace')
    write(*, *); write(*, *) '        Writing initial positions to: ', './fish_ini.dat'
    write(fid, "(I6)") nfish ! read number of particles
    layerdepth = bottom_depth/5.0_sp

    do ii = 1, nfish
      if (ii <= nfish/2) then ! write surface particles
      
        write(fid, "(I4,3F20.6)") ii, abs(random%uniform()*x_range), abs(random%uniform()*y_range), -0.1_SP
        !write(fid, "(I4,3F20.6)") ii, 0.001_SP, abs(random%uniform()*y_range), -0.1_SP ! min growth case
        !write(fid, "(I4,3F20.6)") ii, 250.0_SP, abs(random%uniform()*y_range), -0.1_SP ! max growth case
        
      else ! write demersal particles
      
        write(fid, "(I4,3F20.6)") ii, abs(random%uniform()*x_range), abs(random%uniform()*y_range), -4.9_SP
        !write(fid, "(I4,3F20.6)") ii, 0.001_SP, abs(random%uniform()*y_range), -4.9_SP ! min growth case
        !write(fid, "(I4,3F20.6)") ii, 250.0_SP, abs(random%uniform()*y_range), -4.9_SP ! max growth case
        
      end if
    end do

  else
    write(*, *); write(*, *) "No fish in this simulation..."
  end if

  write(*, *); write(*, *) "Setup finished..."
  call random%stats(); write(*,*)
  deallocate(random)

end program setup
