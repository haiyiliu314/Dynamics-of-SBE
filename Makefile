SHELL   = /bin/bash
CMP     = ifort # gfortran 
LIBS    = -mkl #-llapack -lblas #${LINK_LAPACK} #-mkl #-llapack -lblas
FRAMEWORK = # -framework Accelerate
FLAGS   = -O3 #-check all #-no-wrap-margin #-fbounds-check
DEBUG   = #-g -O0 #-ftrapuv #-check all -traceback #-O3
OBJrun  = no_matrices.o constants.o 

run.x: no_matrices.o
	${CMP} -o run.x ${OBJrun} ${LIBS} ${FRAMEWORK} ${OMP} ${DEBUG} ${FLAGS}

no_matrices.o: no_matrices.f90 constants.o 
	${CMP} -c ${DEBUG}  $<

constants.o: constants.f90 
	${CMP} -c ${DEBUG}  $<

clean:
	rm -f *.mod *.o *~ *.exe *.x *.dat

core: run.x
	ulimit -c unlimited; ./run.x	
