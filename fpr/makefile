FC = ifort
FFLAGS = -g
FOPT = -O3
FLIBS = -llapack -lblas
EXEC = fpr

SRC = main.f90
OBJ = main.o


all: $(OBJ)
	$(FC) $^ -o $(EXEC) $(FFLAGS) $(FLIBS)
$(OBJ): $(SRC)
	$(FC) $(FFLAGS) $(FOPT) $^ -c
