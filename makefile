FC = gfortran
CFLAGS = -c -std=f2003 -Wextra -Wall -pedantic -ffree-line-length-none -fbounds-check

setup: setup.o ichthyotox main.o mod_var.o mod_lag.o behavior.o
		$(FC) -o setup parameters.o mod_var.o simulation.o setup.o

setup.o: parameters.o mod_var.o mod_lag.o behavior.o simulation.o main.o ichthyotox setup.f95
		$(FC) $(CFLAGS) setup.f95

ichthyotox: main.o mod_var.o parameters.o mod_inp.o utilities.o simulation.o mod_lag.o behavior.o cyanobacteria.o triangle_grid_edge.o alloc_vars.o
		$(FC) -o ichthyotox parameters.o mod_var.o mod_inp.o utilities.o simulation.o mod_lag.o behavior.o cyanobacteria.o triangle_grid_edge.o alloc_vars.o main.o

main.o: parameters.o mod_var.o mod_inp.o simulation.o mod_lag.o behavior.o cyanobacteria.o triangle_grid_edge.o alloc_vars.o main.f95
		$(FC) $(CFLAGS) main.f95

# multiple dependencies

behavior.o: behavior.f95 mod_var.o cyanobacteria.o parameters.o simulation.o mod_lag.o
		$(FC) $(CFLAGS) behavior.f95

cyanobacteria.o: cyanobacteria.f95 mod_var.o parameters.o simulation.o mod_lag.o
		$(FC) $(CFLAGS) cyanobacteria.f95
 
triangle_grid_edge.o: triangle_grid_edge.f90 mod_var.o
		$(FC) $(CFLAGS) triangle_grid_edge.f90

alloc_vars.o: alloc_vars.f90 mod_var.o
		$(FC) $(CFLAGS) alloc_vars.f90

mod_lag.o: mod_lag.f95 utilities.f90 parameters.o simulation.o utilities.o mod_var.o
		$(FC) $(CFLAGS) mod_lag.f95

simulation.o: simulation.f95 parameters.o
		$(FC) $(CFLAGS) simulation.f95

mod_inp.o: mod_inp.f90 parameters.o mod_var.o
		$(FC) $(CFLAGS) mod_inp.f90

mod_var.o: mod_var.f95 utilities.f90 parameters.o utilities.o
		$(FC) $(CFLAGS) mod_var.f95

# depends only on parameters
utilities.o: utilities.f90 parameters.o
		$(FC) $(CFLAGS) utilities.f90

parameters.o: parameters.f95
		$(FC) $(CFLAGS) parameters.f95

clean: 
	$(RM) *.o *~ *.mod