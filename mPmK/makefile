FC = ifort
FFLAGS = #-Wall -ffree-line-length-none -fmax-stack-var-size=32768
FLIB =  -lslicot -llapack -lblas # -lccplib
FOPT = -O3

all: mPmK1.o mPmK1
mPmK1.o: mPmK1.f90
	$(FC) $(FFLAGS) $(FOPT) mPmK1.f90 -c
mPmK1: mPmK1.o
	$(FC) $(FFLAGS) $(FOPT) mPmK1.o -o mPmK1 $(FLIB)
clean:
	rm -f *.o
	rm -f mPmK1
