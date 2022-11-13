FC = ifort
FFLAGS = -g #-Wall -ffree-line-length-none -fmax-stack-var-size=32768
FOPT = -O3

all: SSKF1.o SSKF1
SSKF1.o: SSKF1.f90
	$(FC) $(FFLAGS) $(FOPT) SSKF1.f90 -c
SSKF1: SSKF1.o
	$(FC) SSKF1.o -o SSKF1 $(FCCP) $(FFLAGS) $(FOPT)
clean:
	rm -f *.o *.mod SSKF1
cleandat:
	rm -f test.dat residual.dat qest* perfil* tempo.dat
cleanall: clean cleandat
