CXX=g++

##for Ubuntu
EXTRA_LIB_FLAGS = ##-static -zmuldefs

BOOST_INCLUDE_FLAG = -I /usr/include

##Lapack for Armadillo to do Cholesky decomposition
##GLPK for linear programming procedure
LIB_FLAGS = -larmadillo -llapack -lblas -lglpk -lgfortran $(EXTRA_LIB_FLAGS)

OPT = -O2
## As the Armadillo library uses recursive templates,
## compilation times depend on the level of optimisation:
##
## -O0: quick compilation, but the resulting program will be slow
## -O1: good trade-off between compilation time and execution speed
## -O2: produces programs which have almost all possible speedups,
##      but compilation takes longer

#EXTRA_OPT = -fwhole-program
## Uncomment the above line if you're compiling 
## all source files into one program in a single hit.

#DEBUG = -DARMA_EXTRA_DEBUG
## Uncomment the above line to enable low-level
## debugging.  Lots of debugging information will
## be printed when a compiled program is run.
## Please enable this option when reporting bugs.

FINAL = -DARMA_NO_DEBUG
## Uncomment the above line to disable Armadillo's checks.
## DANGEROUS!  Not recommended unless your code has been
## thoroughly tested.

#
#
#

CXXFLAGS = $(BOOST_INCLUDE_FLAG) $(DEBUG) $(FINAL) $(OPT) $(EXTRA_OPT)

all: main

main: main.cpp vol.cpp vol.h
	$(CXX) $(CXXFLAGS) -o PolyVest main.cpp vol.cpp $(LIB_FLAGS)

clean:
	rm -f PolyVest
