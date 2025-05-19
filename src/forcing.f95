program forcing
    ! This program generates a simple triangular grid for the simulation
    use variables, only : zero, folderprefix, sp
    use random, only : random_number_generator
    use simulation, only : write_mesh_data
    use cyanobacteria, only : waterDensity, writeInitialState

    implicit none
    character(len=50) :: foldername, sigma_data_format
    real(sp), dimension(:), allocatable :: diffusivity, temperature, density
    real(sp) :: &
        & biomass, microcystin, toxin_production, toxin_excretion, &
        & bottom_depth, div_inc, x_range, y_range, dt, day, &
        & clock_time, salinity, time, control_temp, temp_slope
    integer :: &
        & layers, fid=500, step, days, steps, ncolony

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

    call write_mesh_data(folderprefix, bottom_depth, div_inc, x_range, y_range, layers)
    
    allocate(diffusivity(layers))
    allocate(temperature, density, mold=diffusivity)
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
    write(fid,"(A,F4.2)") "DTI = ", 0.02_SP ! inner interpolation time step, float
    write(fid,"(A,F4.2)") "INSTP = ", dt ! time step of physical field data, float
    write(fid,"(A,F4.2)") "DTOUT = ", 0.1_SP ! output time step, >dti
    write(fid,"(A,F4.2)") "DHOR = ", 0.10_SP ! horizontal diffusion coefficient
    write(fid,"(A,F4.2)") "DTRW = ", 0.02_SP ! random walk time step
    write(fid,"(A,F8.4)") "TOXP = ", 0.0_SP ! toxin production rate
    write(fid,"(A,F8.4)") "TOXE = ", 0.0_SP ! toxin excretion rate
    write(fid,"(A,I4)") "TDRIFT = ", steps  ! number of time steps to iterate, int
    write(fid,"(A,I2)") "HOURLAG = ", 0
    write(fid,"(A,I1)") "EXPTYPE = ", 1
    close(fid)

end program
