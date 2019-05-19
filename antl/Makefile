#
# This is a non-recursive make file, implemented along the lines of:
# Miller, P. (1997). Recursive Make Considered Harmful.
# 	http://aegis.sourceforge.net/auug97.pdf
# 
# This file was adapted from the example given by:
# Crosby, S. (2001): A build system using "Recursive Make Considered harmful."
#	http://lists.gnu.org/archive/html/autoconf/2001-08/msg00011.html
#
# Run all make commands from the directory containing this makefile.
# 
# To make everything, run:			make all
# To make programs, run:			make appl
# To make library libantl.a, run:	make lib (or just make, with no args)
# To make an individual app, run:  	make appl/my_path_to/my_application
# To make an individual .o, run:	make src/my_path_to/my_obj_module.o
# 
# "make clean" deletes all compiled files.
# "make veryclean" deletes all compiled files and all dependency (.d) files.
#
# Author: John Lees-Miller (2005-06-06).
# Version: $Header: /scrinium/ANTL/ANTL/Makefile,v 1.16 2013/07/27 06:49:20 jacobs Exp $
# 

##
## Project Options
##

# include user-defined flags and paths
include CONFIG

# MPI flags -- cannot compile statically
MFLAGS := $(CFLAGS)
MFLAGS += -DMPICH_IGNORE_CXX_SEEK -DPARALLEL #MPI Stuff
#CFLAGS += -static


#
# Declare all subdirectories w/ files needing attention from make.
# These paths are relative to the location of this Makefile.
# 
DIRS := appl/Tests/Exponentiation
DIRS += appl/Tests/Arithmetic
DIRS += appl/Tests/Quadratic
DIRS += src
DIRS += src/Arithmetic
DIRS += src/Quadratic
DIRS += src/Quadratic/Cube
DIRS += src/Quadratic/Multiply
DIRS += src/Quadratic/Reduce
DIRS += src/Quadratic/Square
DIRS += src/XGCD

##
## Non-Recursive Makefile Implementation
## (You shouldn't need to change anything below this line unless making major
## changes to the structure of the project.)
##

# Make sure that 'lib' is the default rule (first in file).
lib:

# Flush suffixes to prevent make from using implicit suffix rules.
.SUFFIXES:

#
# Configure source groups (group by libs and compile flags).
#

# ANTL_ Source Group: Source files to be compiled into the libantl.a library.
ANTL_SRC := # Empty; added in submakefiles.
ANTL_INC  := -I./include $(GMP_INC) $(NTL_INC) $(BLAS_INC) $(IML_INC) $(GIVARO_INC) $(LINBOX_INC)
ANTL_CFLAGS := $(ANTL_INC) $(CFLAGS) 

# APPL_ Source Group: Source files to be compiled as applications, and linked
# against libantl.a
APPL_SRC := #Empty; added in submakefiles.
APPL_INC := -I./include $(GMP_INC) $(NTL_INC) $(BLAS_INC) $(IML_INC) $(GIVARO_INC) $(LINBOX_INC)
APPL_CFLAGS := $(APPL_INC) $(CFLAGS)
#APPL_LFLAGS := -L./lib/$(ARCH) -lantl $(LINBOX_LIB) $(GIVARO_LIB) $(IML_LIB) $(BLAS_LIB) $(NTL_LIB) $(GMP_LIB) -lm
#APPL_LFLAGS := -L./lib/$(ARCH) -lantl $(LINBOX_LIB) $(GIVARO_LIB) $(IML_LIB) $(BLAS_LIB) $(NTL_LIB) -Wl,--rpath -Wl,/scratch/lib/gmp/lib $(GMP_LIB) -lm
APPL_LFLAGS := -L./lib/$(ARCH) -lantl $(LINBOX_LIB) $(GIVARO_LIB) $(IML_LIB) $(BLAS_LIB) $(NTL_LIB) $(GMP_LIB) -lm 
APPL_MFLAGS := $(APPL_INC) $(MFLAGS)
APPL_LFLAGS += $(LFLAGS)
APPL_MPI := #Empty; added in submakefiles.

#
# Include submakefiles in subdirectories (Makefile.dir files).
# They will add to the CLEANFILES and *_SRC variables, depending on what
# flags are needed by the compiler.
#
CLEANFILES := # Empty; added in submakefiles.

include $(patsubst %,%/Makefile.dir,$(DIRS))

#
# Find object files for each source group, and aggregate them.
#

ANTL_OBJ := $(patsubst %.cpp,%.o,$(filter %.cpp,$(ANTL_SRC)))
APPL_OBJ := $(patsubst %.cpp,%.o,$(filter %.cpp,$(APPL_SRC)))
MPI_OBJ := $(patsubst %.cpp,%.o,$(filter %.cpp,$(APPL_MPI)))

SRC := $(ANTL_SRC) $(APPL_SRC) $(APPL_MPI)
OBJ := $(ANTL_OBJ) $(APPL_OBJ) $(MPI_OBJ)
BIN := $(patsubst %.cpp, %, $(APPL_SRC)) 
MPI := $(patsubst %.cpp, %, $(APPL_MPI))

$(info $(OBJ))
$(info $(BIN))

#
# Define rule for building dependencies, then include .d(ependency) files.
#

# Flags to pass for dependency-finding. All and only include files that are
# _not_ part of 3rd party libraries should be in this path.
DEPFLAGS := -I./include -DMPICH_IGNORE_CXX_SEEK $(IML_INC) 

# This beast is the actual rule that finds the dependencies.
$(patsubst %.cpp,%.d,$(filter %.cpp,$(SRC))) : %.d: %.cpp
	@echo "***** Doing dependencies for $<"
	@set -e; b=`basename $* .cpp` ; d=`dirname $*` ; \
			$(CXX) -MM $(DEPFLAGS) $< \
		| sed "s%\\($$b\\)\\.o[ :]*%$${d}/\\1.o $${d}/\\1.d : %g" > $@; \
		[ -s $@ ] || rm -f $@

# Dependencies are computed for all object files.
DEPS := $(patsubst %.o, %.d, $(OBJ))

# Note: using -include to suppress warnings due to missing .d files.
-include $(DEPS)

$(MPI): CXX := $(MPICXX)
$(MPI): ANTL_CFLAGS := $(ANTL_INC) $(MFLAGS)

#
# Source group build rules.
#

$(ANTL_OBJ): %.o : %.cpp %.d
	$(CXX) $(ANTL_CFLAGS) -o $@ -c $<

#
# Library build rule.
#

LIB_NAME := libantl.a
LIB_PATH := ./lib/$(ARCH)
LIB_FULL_NAME := $(LIB_PATH)/$(LIB_NAME)

$(LIB_FULL_NAME) : $(ANTL_OBJ)
	mkdir -p $(LIB_PATH)
	ar cru $@ $^

#
# Application build rule
#

$(BIN) : % : %.cpp $(LIB_FULL_NAME)
	$(CXX) $(APPL_CFLAGS) $< $(APPL_LFLAGS) -o $@ 

#
# MPI files build rule
#

$(MPI) : % : %.cpp $(LIB_FULL_NAME)
	$(MPICXX) $(APPL_MFLAGS) $< $(APPL_LFLAGS) -o $@


#
# Aggregate rules.
#

all: appl lib
appl: $(BIN)
lib: $(LIB_FULL_NAME)
mpiappl: $(MPI)
obj: $(OBJ)
dep: $(DEPS)

cleandep:
	rm -f $(DEPS)

cleanobj:
	rm -f $(OBJ)

cleanlib:
	rm -f $(LIB_FULL_NAME)

cleanappl:
	rm -f $(BIN) $(MPI)

# Clean objects, binaries, libraries and any submakefile-specific files.
clean: cleanobj cleanappl cleanlib
	rm -f $(CLEANFILES)

# Cleans and destroys dependencies and all libraries.
veryclean: clean cleandep
	rm -rf ./lib/*
