program grid
    ! This program generates a simple triangular grid for the simulation
    use variables, only : zero, folderprefix, sp
    use random, only : random_number_generator
    use cyanobacteria, only : colonyBaseRadius, densityMin, densityMax, &
        & cellFrac, vesicleFrac, vesicleDensity, cellDensityCoefficient

    implicit none

    integer :: ncolony, nfish
   
    real(sp) :: &
            & biomass, carbon, carbohydrate, protein, microcystin, ratio, radius, &
            & reqCellDensity, &
            & x, y, z

    integer :: &
        & nodes, elements, node_index=1, element_index=1, start_pattern=1, &
        & alternate_pattern, x_div, y_div, ii, jj, layers, fid=500, fid2=501, &
        & x_nodes, y_nodes, step, ndays, steps
    integer, dimension(:, :), allocatable :: vertices
    character(len=50) :: foldername, kb_format
    real(sp) :: bottom_depth, div_inc, x_range, y_range, dt, day, &
        & clock_time, salinity, time, control_temp, temp_slope
    real(sp), dimension(:), allocatable :: diffusivity, temperature, density

    call get_command_argument(1, foldername) ! Import case name from command line
    if (len_trim(foldername) == 0) then
      print *, 'Please provide a simulation ID as command line argument'
      print *, 'Stopping...'; stop
    end if
    folderprefix = "data/"//adjustl(foldername)

    print *, "Build a right-triangulated rectangular domain..."
    print *
    write(*, '(A)', advance='no') "    X range (meters): "; read(*, *) x_range
    write(*, '(A)', advance='no') "    Y range (meters): "; read(*, *) y_range
    write(*, '(A)', advance='no') "    Uniform depth (meters): "; read(*, *) bottom_depth
    write(*, '(A)', advance='no') "    Size of elements (meters): "; read(*, *) div_inc
    write(*, '(A)', advance='no') "    Sigma layers (#): "; read(*, *) layers
    write(*, '(A)', advance='no') "    Length of simulation (days): "; read(*, *) ndays
    write(*, '(A)', advance='no') "    Initial temperature (C): "; read(*, *) control_temp
    write(*, '(A)', advance='no') "    Temperature slope (C/hour): "; read(*, *) temp_slope
    write(*, '(A)', advance='no') "    Uniform salinity (ppt): "; read(*, *) salinity
    write(*, '(A)', advance='no') "    Time step for forcing (hours): "; read(*, *) dt
    write(*, '(A)', advance='no') "    Cyanobacteria particles (#): "; read(*, *) ncolony
    write(*, '(A)', advance='no') "    Total protein biomass (grams): "; read(*, *) biomass
    write(*, '(A)', advance='no') "    Total microcystin load (grams): "; read(*, *) microcystin
    write(*, '(A)', advance='no') "    Fish particles (#): "; read(*, *) nfish
  
    allocate( diffusivity(layers), temperature(layers), density(layers) )
    write(kb_format, "(A1,I6,A7)") "(", 3*layers, "F20.10)"
    x_div = floor( x_range / div_inc )
    y_div = floor( y_range / div_inc )
    x_nodes = x_div + 1
    y_nodes = y_div + 1
    nodes = (x_div + 1) * (y_div + 1)
    elements = 2 * x_div * y_div
    bottom_depth = -abs(bottom_depth) ! force negative depth

    print *
    write(*, '(A,I5)') "    Nodes:    ", nodes
    write(*, '(A,I5)') "    Elements: ", elements
    print *
    
    print *, "Saving mesh data... "
    allocate(vertices(elements,3))
    open(unit=fid, file=trim(folderprefix)//"/nodes.txt", status='replace')
    open(unit=fid2, file=trim(folderprefix)//"/elements.txt", status='replace')
    write(fid, *) nodes
    write(fid2, *) elements, layers
    do ii = 1, y_nodes
        alternate_pattern = start_pattern
        do jj = 1, x_nodes
            write(fid, *) ii, float(jj-1) * div_inc, float(ii-1) * div_inc, bottom_depth
            if ( (jj < x_nodes) .and. (ii < y_nodes) ) then
                ! bottom row of triangles
                vertices(element_index, 1) = node_index
                vertices(element_index, 2) = node_index + x_nodes ! note index
                vertices(element_index, 3) = node_index + 1
                if (alternate_pattern < 0) vertices(element_index, 2) = vertices(element_index, 2) + 1
                write(fid2, *) element_index, vertices(element_index, 1:3)
                element_index = element_index + 1
              
                ! upper row of triangles
                vertices(element_index, 1) = node_index + x_nodes
                vertices(element_index, 2) = node_index + x_nodes + 1
                vertices(element_index, 3) = node_index + 1
                if (alternate_pattern < 0) vertices(element_index, 3) = vertices(element_index, 3) - 1
                write(fid2, *) element_index, vertices(element_index, 1:3)
                element_index = element_index + 1
                alternate_pattern = -alternate_pattern ! switch triangle pattern as you move x=0 to x=x_range
            end if
            node_index = node_index + 1 
        end do
        start_pattern = -start_pattern ! switch triangle pattern between rows as grid is built
    end do
    close(fid)
    close(fid2)

    print *, "Saving forcing data..."
    open(unit=fid, file=trim(folderprefix)//"/forcing.txt", status='replace')
    ! Uniform profile
    diffusivity(:) = 60.0_SP * 60.0_SP * 10.0_SP**(-5.0_SP)
    steps = 24*ndays
    do step = 0, steps
        time = float(step)*dt
        day = time / 24.0_sp
        clock_time = time - floor(day)*24.0_sp
        ! Uniform profile
        temperature(:) = control_temp + temp_slope*time

        ! Density of pure water from Millero and Poisson
        density(:) = 999.842594_SP + 0.06793952_SP*temperature(:) - 0.009095290_SP*temperature(:)**(2.0_SP) +  &
        & 0.0001001685_SP*temperature(:)**(3.0_SP) - (1.120083D-6)*temperature(:)**(4.0_SP) + &
        & (6.536332D-9)*temperature(:)**(5.0_SP)
        density(:) = density(:) + salinity*( 0.824493_SP - 0.0040899_SP*temperature(:) + &
        & 0.000076438_SP*temperature(:)**(2.0_SP) - (8.2467D-7)*temperature(:)**(3.0_SP) + &
        & (5.3875D-9)*temperature(:)**(4.0_SP) )
        density(:) = density(:) + salinity**(1.5_SP)*( -0.00572466_SP + 0.00010227_SP*temperature(:) - &
        & (1.6546D-6)*temperature(:)**(2.0_SP) )
        density(:) = density(:) + 0.00048314_SP*salinity**(2.0_SP)

        write(fid, "(1F10.3)", advance='no') time
        write(fid, kb_format, advance='no') temperature(:), density(:), diffusivity(:)
        write(fid, *) ! new line
    end do
    close(fid)

    allocate(random_number_generator);
    call random_number_generator%init()

    ! write state variables
    open(unit = fid, file = trim(folderprefix)//"/cyanobacteria.txt", status = 'replace')
    print *, 'Saving cyanobacteria state...'
    write(fid, "(I6)") ncolony
    do ii = 1, ncolony
        x = abs(random_number_generator%uniform()*x_div*div_inc)
        y = abs(random_number_generator%uniform()*y_div*div_inc)
        z = -1.0*abs(random_number_generator%uniform()*5.0)
        radius = colonyBaseRadius + 5.0_SP*10.0_SP**(-6.0_SP)*random_number_generator%get()
        carbon = biomass/float(ncolony)
        reqCellDensity = ((density(1) - (1.0_SP - cellFrac)*(density(1) + 0.7_SP)) / &
            & cellFrac - vesicleFrac*vesicleDensity)/(1.0_SP - vesicleFrac)
        ratio = log(1.0_SP - (reqCellDensity - densityMin)/(densityMax - densityMin)) / &
            & (-cellDensityCoefficient)
        protein = carbon
        carbohydrate = protein*ratio

        write(fid, "(7F20.6)") x, y, z, radius, carbohydrate, protein, microcystin/float(ncolony)
    end do
    close(fid)

    open(unit = fid, file = trim(folderprefix)//"/fish.text", status = 'replace')
    print *, 'Saving fish data...'
    write(fid, "(I6)") nfish
    do ii = 1, nfish
        x = abs(random_number_generator%uniform()*x_range)
        y = abs(random_number_generator%uniform()*y_range)
        if (ii <= nfish/2) then ! surface
            z = -0.1_SP
        else ! demersal
            z = -4.9_SP
        end if
        write(fid, "(I4,3F20.6)") ii, x, y, z
    end do
    close(fid)

    call random_number_generator%stats();
    deallocate(random_number_generator)

    print *, "Saving parameters..."
    open(unit=fid, file=trim(folderprefix)//"/parameters.txt", status='replace')
    write(fid,"(A)") "INFOFILE = screen"
    write(fid,"(A,F4.2)") "DTI = ", 0.02_SP ! inner interpolation time step, float"
    write(fid,"(A,F4.2)") "INSTP = ", dt ! time step of physical field data, float
    write(fid,"(A,F4.2)") "DTOUT = ", 0.1_SP ! output time step, >dti
    write(fid,"(A,F4.2)") "DHOR = ", 0.10_SP ! horizontal diffusion coefficient
    write(fid,"(A,F4.2)") "DTRW = ", 0.02_SP ! random walk time step
    write(fid,"(A,I4)") "TDRIFT = ", steps  ! number of time steps to iterate, int
    write(fid,"(A,I4)") "YEARLAG = ", 2016
    write(fid,"(A,I2)") "MONTHLAG = ", 4
    write(fid,"(A,I2)") "DAYLAG = ", 1
    write(fid,"(A,I2)") "HOURLAG = ", 0
    write(fid,"(A,I1)") "IRW = ", 0 ! random walk type
    write(fid,"(A)") "P_SIGMA = F"
    write(fid,"(A)") "OUT_SIGMA = F"
    write(fid,"(A)") "F_DEPTH = F" 
    write(fid,"(A)") "GEOAREA = box"
    write(fid,"(A)") "INPDIR=/"
    write(fid,"(A)") "LAGINI=/"
    write(fid,"(A)") "OUTDIR=/"
    close(fid)

end program
