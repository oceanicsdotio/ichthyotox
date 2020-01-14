FC = gfortran
CFLAGS = -c -std=f2003 -Wextra -Wall -pedantic -ffree-line-length-none -fbounds-check

setup: ichthyotox_setup.o ichthyotox offlag.o mod_var.o mod_lag.o mod_fish_watkins.o
		$(FC) -o setup mod_prec.o parameters.o mod_var.o mod_sim.o ichthyotox_setup.o

ichthyotox_setup.o: mod_prec.o parameters.o mod_var.o mod_lag.o mod_fish_watkins.o mod_sim.o offlag.o ichthyotox ichthyotox_setup.f95
		$(FC) $(CFLAGS) ichthyotox_setup.f95

ichthyotox: offlag.o ncdio.o mod_prec.o mod_var.o parameters.o mod_inp.o util.o mod_sim.o mod_lag.o mod_fish_watkins.o mod_tox.o triangle_grid_edge.o alloc_vars.o data_run.o
		$(FC) -o ichthyotox mod_prec.o ncdio.o parameters.o mod_var.o data_run.o mod_inp.o util.o mod_sim.o mod_lag.o mod_fish_watkins.o mod_tox.o triangle_grid_edge.o alloc_vars.o offlag.o

offlag.o: mod_prec.o mod_var.o data_run.o mod_sim.o mod_lag.o mod_fish_watkins.o mod_tox.o triangle_grid_edge.o alloc_vars.o offlag.f95
		$(FC) $(CFLAGS) offlag.f95 

# multiple dependencies

mod_fish_watkins.o: mod_fish_watkins.f95 mod_var.o mod_tox.o mod_prec.o parameters.o mod_sim.o mod_lag.o
		$(FC) $(CFLAGS) mod_fish_watkins.f95

mod_tox.o: mod_tox.f95 mod_var.o mod_prec.o parameters.o mod_sim.o mod_lag.o
		$(FC) $(CFLAGS) mod_tox.f95
 
triangle_grid_edge.o: triangle_grid_edge.f90 mod_var.o
		$(FC) $(CFLAGS) triangle_grid_edge.f90

ncdio.o: ncdio.f90 mod_var.o alloc_vars.o
		$(FC) $(CFLAGS) ncdio.f90

alloc_vars.o: alloc_vars.f90 mod_var.o
		$(FC) $(CFLAGS) alloc_vars.f90

mod_lag.o: mod_lag.f95 util.f90 mod_prec.o mod_sim.o util.o mod_var.o
		$(FC) $(CFLAGS) mod_lag.f95

mod_sim.o: mod_sim.f95 mod_prec.o parameters.o
		$(FC) $(CFLAGS) mod_sim.f95

data_run.o: data_run.f90 mod_inp.o mod_var.o
		$(FC) $(CFLAGS) data_run.f90

mod_var.o: mod_var.f95 util.f90 mod_prec.o util.o
		$(FC) $(CFLAGS) mod_var.f95

# depends only on precision
parameters.o: parameters.f95 mod_prec.o
		$(FC) $(CFLAGS) parameters.f95

mod_inp.o: mod_inp.f90 mod_prec.o
		$(FC) $(CFLAGS) mod_inp.f90

util.o: util.f90 mod_prec.o
		$(FC) $(CFLAGS) util.f90

# no dependencies
mod_prec.o: mod_prec.f90
		$(FC) $(CFLAGS) mod_prec.f90
		
clean: 
	$(RM) *.o *~ *.mod