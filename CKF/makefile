FC = ifort
FFLAGS = -g -heap-arrays
FOPT = -O3
FLIBS = -qmkl#-llapack -lblas
EXEC = CKF

SRC = main.f90
OBJ = main.o


all: $(OBJ)
	$(FC) $^ -o $(EXEC) $(FFLAGS) $(FLIBS)
$(OBJ): $(SRC)
	$(FC) $(FFLAGS) $(FOPT) $(FLIBS) $^ -c
