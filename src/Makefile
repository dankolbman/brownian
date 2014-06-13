# The makefile for brownian
# Dan Kolbman 2014

# Change if you use a different compiler
FC = gfortran

# Optimization flag
FLAGS = -O2

# The directory to place the binary
BINDIR = ../bin/

NAME = brownian

# Required object files
OBJECTS = grcalc.o grinit.o scale.o formlattice.o movepart.o readcon.o forcecalc.o brownian.o

all: brownian

brownian: $(OBJECTS)
	mkdir -p $(BINDIR)
	$(FC) $(OBJECTS) -o $(BINDIR)$(NAME)

%.o: %.f param
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o *.mod *.MOD
