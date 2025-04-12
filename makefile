FC = gfortran
CFLAGS = -c -std=f2003 -Wextra -Wall -pedantic -ffree-line-length-none -fbounds-check
LDFLAGS = -std=f2003 -Wextra -Wall -pedantic -ffree-line-length-none -fbounds-check
# Object code that contains all global variables
# and is used by all other object code
src/%.o src/%.mod: src/%.f95
	$(FC) $(CFLAGS) -Isrc/ -o $@ -J src $<

setup: bin/variables.o bin/random.o bin/lagrangian.o bin/cyanobacteria.o bin/setup.o
	$(FC) -o setup variables.o simulation.o setup.o

.: src/variables.o src/random.o src/simulation.o src/lagrangian.o src/cyanobacteria.o
	@true

# setup.o: variables.o lagrangian.o behavior.o simulation.o main.o ichthyotox setup.f95
# 	$(FC) $(CFLAGS) setup.f95

# ichthyotox: main.o variables.o simulation.o lagrangian.o behavior.o cyanobacteria.o
# 	$(FC) -o ichthyotox variables.o simulation.o lagrangian.o behavior.o cyanobacteria.o main.o

# main.o: variables.o simulation.o lagrangian.o behavior.o cyanobacteria.o main.f95
# 	$(FC) $(CFLAGS) main.f95

# behavior.o: behavior.f95 variables.o cyanobacteria.o simulation.o lagrangian.o
# 	$(FC) $(CFLAGS) behavior.f95

# cyanobacteria.o: cyanobacteria.f95 variables.o simulation.o lagrangian.o
# 	$(FC) $(CFLAGS) cyanobacteria.f95

# lagrangian.o: lagrangian.f95 simulation.o variables.o
# 	$(FC) $(CFLAGS) lagrangian.f95

# simulation.o: simulation.f95 variables.o
# 	$(FC) $(CFLAGS) simulation.f95

clean: 
	rm src/*.mod
	rm src/*.o