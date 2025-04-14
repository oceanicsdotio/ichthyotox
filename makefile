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

test.cpython-313-darwin.so: src/variables.o src/random.o
	pixi run f2py -c src/variables.f95 src/random.f95 src/test.f95 -m test

new: bin/forcing
	mkdir -p data/test
	bin/forcing test

clean: 
	@ rm -f src/*.mod
	@ rm -f src/*.o
	@ rm -f bin/*
	@ rm -rf data/test
	@ rm -f test.cpython-313-darwin.so
	
.PHONY: clean