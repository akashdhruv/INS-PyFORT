FF=f2py
FFLAGS = --fcompiler=gnu95 -m


MPIFLAGS = --f90flags=-lmpich --fcompiler=mpif90 -m

all: INS.so POISSON.so HEAT.so

POISSON.so: POISSON.F90
	$(FF) $(FFLAGS) POISSON -c POISSON.F90

INS.so: INS.F90
	$(FF) $(FFLAGS) INS -c INS.F90

HEAT.so: HEAT.F90
	$(FF) $(FFLAGS) HEAT -c HEAT.F90

cleanall:
	rm -f *.so *.dat

clean:
	rm -f *.dat

run: Solver.py
	python Solver.py
