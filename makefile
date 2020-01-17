FC = gfortran
CFLAGS = -c -std=f2003 -Wextra -Wall -pedantic -ffree-line-length-none -fbounds-check

bin.setup: setup.o ichthyotox main.o variables.o lagrangian.o behavior.o
		$(FC) -o setup variables.o simulation.o setup.o

setup.o: variables.o lagrangian.o behavior.o simulation.o main.o ichthyotox setup.f95
		$(FC) $(CFLAGS) setup.f95

bin/ichthyotox: main.o variables.o simulation.o lagrangian.o behavior.o cyanobacteria.o
		$(FC) -o ichthyotox variables.o simulation.o lagrangian.o behavior.o cyanobacteria.o main.o

main.o: variables.o simulation.o lagrangian.o behavior.o cyanobacteria.o main.f95
		$(FC) $(CFLAGS) main.f95

behavior.o: behavior.f95 variables.o cyanobacteria.o simulation.o lagrangian.o
		$(FC) $(CFLAGS) behavior.f95

cyanobacteria.o: cyanobacteria.f95 variables.o simulation.o lagrangian.o
		$(FC) $(CFLAGS) cyanobacteria.f95

lagrangian.o: lagrangian.f95 simulation.o variables.o
		$(FC) $(CFLAGS) lagrangian.f95

simulation.o: simulation.f95 variables.o
		$(FC) $(CFLAGS) simulation.f95

variables.o: variables.f95
		$(FC) $(CFLAGS) variables.f95

clean: 
	$(RM) *.o *~ *.mod