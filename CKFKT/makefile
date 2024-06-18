FC = ifort
FFLAGS = -g -heap-arrays
FOPT = -O3
FLIBS = -qmkl#-llapack -lblas
EXEC = CKFKT

SRC = main.f90
OBJ = main.o


run: bld
	@ ./$(EXEC)
bld: $(OBJ)
	@ $(FC) $^ -o $(EXEC) $(FFLAGS) $(FLIBS)
$(OBJ): $(SRC)
	@ $(FC) $(FFLAGS) $(FOPT) $(FLIBS) $^ -c
plt:
	@ gnuplot graphics.gnp
cln:
	@ rm -f $(OBJ) $(EXEC) *.png profiles.dat time.dat vqest.dat vR.dat vTest_c.dat vTest.dat kirch.dat *.mod
