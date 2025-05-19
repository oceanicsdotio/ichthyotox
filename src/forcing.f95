program forcing
    ! This program generates a simple triangular grid for the simulation
    use variables, only : zero, folderprefix, sp
    use random, only : random_number_generator
    use cyanobacteria, only : waterDensity, writeInitialState

    implicit none
    character(len=50) :: foldername, sigma_data_format
    real(sp), dimension(:), allocatable :: diffusivity, temperature, density
    real(sp) :: &
            & biomass, microcystin, toxin_production, toxin_excretion, &
            & bottom_depth, div_inc, x_range, y_range, dt, day, &
            & clock_time, salinity, time, control_temp, temp_slope
    integer :: &
        & nodes, elements, node_index=1, element_index=1, start_pattern=1, &
        & alternate_pattern, x_div, y_div, ii, jj, layers, fid=500, fid2=501, &
        & x_nodes, y_nodes, step, days, steps, ncolony
    integer, dimension(3) :: vertices

    call get_command_argument(1, foldername) ! Import case name from command line
    if (len_trim(foldername) == 0) then
      print *, 'Please provide a simulation ID as command line argument'
      print *, 'Stopping...'; stop
    end if
    folderprefix = "data/"//adjustl(foldername)

    print *, "Build a right-triangulated rectangular domain..."
    write(*, '(A)', advance='no') "    X range (meters): "; read(*, *) x_range
    write(*, '(A)', advance='no') "    Y range (meters): "; read(*, *) y_range
    write(*, '(A)', advance='no') "    Uniform depth (meters): "; read(*, *) bottom_depth
    write(*, '(A)', advance='no') "    Size of elements (meters): "; read(*, *) div_inc
    write(*, '(A)', advance='no') "    Sigma layers (#): "; read(*, *) layers
    write(*, '(A)', advance='no') "    Length of simulation (days): "; read(*, *) days
    write(*, '(A)', advance='no') "    Initial temperature (C): "; read(*, *) control_temp
    write(*, '(A)', advance='no') "    Temperature slope (C/hour): "; read(*, *) temp_slope
    write(*, '(A)', advance='no') "    Uniform salinity (ppt): "; read(*, *) salinity
    write(*, '(A)', advance='no') "    Time step for forcing (hours): "; read(*, *) dt
    write(*, '(A)', advance='no') "    Cyanobacteria particles (#): "; read(*, *) ncolony
    write(*, '(A)', advance='no') "    Total protein biomass (grams): "; read(*, *) biomass
    write(*, '(A)', advance='no') "    Total microcystin load (grams): "; read(*, *) microcystin
    write(*, '(A)', advance='no') "    Microcystin production rate: "; read(*, *) toxin_production
    write(*, '(A)', advance='no') "    Microcystin excretion rate: "; read(*, *) toxin_excretion
  
    allocate(diffusivity(layers))
    allocate(temperature, density, mold=diffusivity)
    
    x_div = floor( x_range / div_inc )
    y_div = floor( y_range / div_inc )
    x_nodes = x_div + 1
    y_nodes = y_div + 1
    nodes = (x_div + 1) * (y_div + 1)
    elements = 2 * x_div * y_div
    bottom_depth = -abs(bottom_depth) ! force negative depth
    
    print *, "Saving mesh data... "
    write(*, '(A,I5)') "    Nodes:    ", nodes
    write(*, '(A,I5)') "    Elements: ", elements
    open(unit=fid, file=trim(folderprefix)//"/nodes.txt", status='replace')
    open(unit=fid2, file=trim(folderprefix)//"/elements.txt", status='replace')
    write(fid, *) nodes
    write(fid2, *) elements, layers
    do ii = 1, y_nodes
        alternate_pattern = start_pattern
        do jj = 1, x_nodes
            write(fid, *) float(jj-1) * div_inc, float(ii-1) * div_inc, bottom_depth
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

    print *, "Saving forcing data..."
    write(sigma_data_format, "(A1,I6,A7)") "(", 3*layers, "F20.10)"
    open(unit=fid, file=trim(folderprefix)//"/forcing.txt", status='replace')
    ! Uniform profile
    diffusivity(:) = 60.0_SP * 60.0_SP * 10.0_SP**(-5.0_SP)
    steps = 24*days
    do step = 0, steps
        time = float(step)*dt
        day = time / 24.0_sp
        clock_time = time - floor(day)*24.0_sp
        temperature(:) = control_temp + temp_slope*time
        density(:) = waterDensity(temperature, salinity)
        write(fid, "(1F10.3)", advance='no') time
        write(fid, sigma_data_format, advance='no') temperature(:), density(:), diffusivity(:)
        write(fid, *) ! new line
    end do
    close(fid)

    allocate(random_number_generator);
    call random_number_generator%init()
    call writeInitialState(&
        & trim(folderprefix)//"/cyanobacteria.txt", &
        & ncolony, &
        & biomass, &
        & microcystin, &
        & density(1), &
        & random_number_generator)
    call random_number_generator%stats();

    print *, "Saving parameters..."
    open(unit=fid, file=trim(folderprefix)//"/parameters.txt", status='replace')
    write(fid,"(A)") "INFOFILE = screen"
    write(fid,"(A,F4.2)") "DTI = ", 0.02_SP ! inner interpolation time step, float"
    write(fid,"(A,F4.2)") "INSTP = ", dt ! time step of physical field data, float
    write(fid,"(A,F4.2)") "DTOUT = ", 0.1_SP ! output time step, >dti
    write(fid,"(A,F4.2)") "DHOR = ", 0.10_SP ! horizontal diffusion coefficient
    write(fid,"(A,F4.2)") "DTRW = ", 0.02_SP ! random walk time step
    write(fid,"(A,F8.4)") "TOXP = ", 0.0_SP ! toxin production rate
    write(fid,"(A,F8.4)") "TOXE = ", 0.0_SP ! toxin excretion rate
    write(fid,"(A,I4)") "TDRIFT = ", steps  ! number of time steps to iterate, int
    write(fid,"(A,I4)") "YEARLAG = ", 2016
    write(fid,"(A,I2)") "MONTHLAG = ", 4
    write(fid,"(A,I2)") "DAYLAG = ", 1
    write(fid,"(A,I2)") "HOURLAG = ", 0
    write(fid,"(A,I1)") "EXPTYPE = ", 1
    write(fid,"(A)") "GEOAREA = box"
    write(fid,"(A)") "INPDIR=/"
    write(fid,"(A)") "LAGINI=/"
    write(fid,"(A)") "OUTDIR=/"
    close(fid)

end program
