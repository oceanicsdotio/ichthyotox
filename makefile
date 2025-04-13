FC = gfortran
CFLAGS = -c -std=f2008 -Wextra -Wall -pedantic -ffree-line-length-none -fbounds-check
LDFLAGS = -std=f2008 -Wextra -Wall -pedantic -ffree-line-length-none -fbounds-check
OBJ = src/variables.o src/random.o src/simulation.o src/lagrangian.o src/cyanobacteria.o src/io.o

src/%.o src/%.mod: src/%.f95
	$(FC) $(CFLAGS) -Isrc/ -o $@ -J src $<

bin/forcing: src/forcing.f95 $(OBJ) 
	$(FC) -o $@ $^
	chmod +x $@

bin/bloom: src/bloom.f95 $(OBJ) 
	$(FC) -o $@ $^
	chmod +x $@

new: bin/forcing
	mkdir -p data/test
	bin/forcing test

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

clean: 
	rm src/*.mod
	rm src/*.o
	rm bin/*
	rm -rf data/test
	
.PHONY: clean