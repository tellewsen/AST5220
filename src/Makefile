# -*- Makefile -*-

###UIO PART###
FC        = ifort     # F90 compiler at uio
OPTIM     = -O3#-g -C# Optimization flags; set to -g -C for useful debug information

# It should hopefully not be necessary to edit anything below
FITSDIR = /mn/stornext/u3/hke/local/lib      # Directory containing libcfitsio.a
LAPACK  = -L/mn/stornext/u3/hke/local/lib -llapack -lblas
HEALPIX = -L/mn/stornext/u3/hke/local/lib -lhealpix
HEALINC = -I/mn/stornext/u3/hke/local/include
OUTPUT  = cmbspec
FFLAGS  = $(HEALPIX) $(LAPACK)
###END UIO PART###


###LAPTOP PART###
#THIS PART IS NOT FUNCTIONAL AT THIS POINT

#FC        = gfortran     # F90 compiler at own laptop
#OPTIM     = -O3       # Optimization flags; set to -g -C for useful debug information

#FITSDIR = ~/usr/lib      # Directory containing libcfitsio.a
#LAPACK  = ~/usr/lib -llapack -lblas 
#HEALPIX = ~/usr/lib -lhealpix
#HEALINC = ~/usr/include
#OUTPUT  = cmbspec
#FFLAGS  = $(HEALPIX) $(LAPACK)
###END LAPTOP PART###






# List of source files to be compiled
OBJS    = math_tools.o spline_1D_mod.o rk_mod.o bs_mod.o ode_solver.o \
	  params.o time_mod.o rec_mod.o cmbspec.o evolution_mod.o \
          spline_2D_mod.o cl_mod.o sphbess_mod.o

# Linking stage
cmbspec: $(OBJS)
	$(FC) $(FFLAGS) $(OPTIM) -o $(OUTPUT) $(OBJS) $(LAPACK)

# Dependencies
cmbspec.o       : time_mod.o rec_mod.o evolution_mod.o cl_mod.o
time_mod.o      : params.o spline_1D_mod.o
evolution_mod.o : spline_2D_mod.o
cl_mod.o        : sphbess_mod.o
# Compilation of source files
%.o : %.f90
	$(FC) $(OPTIM) $(HEALINC) -c $<

# Clean-up command (write "make clean")
.PHONY: clean
clean:
	rm *.mod *.o *~ cmbspec
