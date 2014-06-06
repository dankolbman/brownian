Brownian motion simulator for co-culture systems by Julian Butcher

## Install 

First, fetch the source files:

    git clone git://github.com/dankolbman/brownian.git
    cd brownian

The default fortran compiler is set to gfortran. If you are using a different
compiler, change `FC=` in the Makefile. Compiler flags can be set by editing
`FCFLAGS=`. Only optimization flags are used by default.

After configuring the Makefile, build the program:

    cd src/
    make
    cd ../

The binary is now in `bin/`
Run it with `bin/brownian` from the project directory.

### Dependencies

The python scripts in 'scripts/' require python 3 and matplotlib for plotting.

## Scripts

Some python scripts are provided to provide some basic plotting functionality
and post-processing.

#### circplot.py
Useage: `python circplot.py sysparam.dat pos1.dat pos2.dat ...`
Plot a circular boundary simulation with the given position files.

#### grplot.py
Useage: `python grplot.py sysparam.dat gr1.dat gr2.dat ...`
Plot the radial distribution functions.

## Output

Many files will be generated in the directory where the executable was run from.
Here are some of the data files generated:
- sysparam.dat - The system parameters used for the simulation
- config(1,2)(n).dat - The initial positions for species 1 or 2 for the nth run.
- gr(11,12,22)(n).dat - The radial distribution function for species 1, 2, or
 1 and 2 for the nth run.
- fpos(1,2)(n) - The final positions for species 1 or 2 for the nth run.
- fgr(1,2)(n) - The final radial distribution for species 1 or 2 for the nth run
- msd(1,2) - The mean squared displacements for species 1 or 2
- msdave(1,2) - The average mean squared displacements over all iterations
 for species 1 or 2
