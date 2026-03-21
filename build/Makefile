makdep = ./makdep_usrbinperl
precision_mode = r8
CONDA_PREFIX = $(HOME)/miniconda3
usrCPP_INCS = -I$(CONDA_PREFIX)/include
nchdflib = -L$(CONDA_PREFIX)/lib -lnetcdff -lnetcdf -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lsz -lz
usrCPP_DEFS_AER = -DDST -DCAER -DSCYC -DWETDEP -DDRYDEP -DASSIM
r8FLAG =
usrCPP_DEFS_PRC = -DRPREC=12 -DREALSIZE=8
usrCPP_DEFS = $(usrCPP_DEFS_PRC)
ifdef usrCPP_DEFS_AER
  usrCPP_DEFS := $(usrCPP_DEFS) $(usrCPP_DEFS_AER)
endif
ifdef usrCPP_DEFS_TRAC
  usrCPP_DEFS := $(usrCPP_DEFS) $(usrCPP_DEFS_TRAC)
endif
usrFFLAGS  = -O3 -ffixed-line-length-132 -fallow-invalid-boz -fallow-argument-mismatch
specFFLAGS = $(usrFFLAGS)
specFFLAGS2 = $(usrFFLAGS)
usrLDFLAGS = -g

###############################################################################
# User shouldn't need to modify code below this point.
###############################################################################
null :=
space := $(null) $(null)
comma := ,

dirs := $(shell cat ./Filepath)   # Load dependency search path.

VPATH := $(dirs)
VPATH := $(subst $(space),:,$(VPATH))   # Change spaces to colons.

cpp_dirs := ./ $(dirs)    # make cpp start search in ./ when includes use < >
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

srclist_F90 := $(shell grep '\.F90' ./Srcfiles)
srclist_F   := $(shell grep '\.F$$' ./Srcfiles)
srclist_c   := $(shell grep '\.c' ./Srcfiles)

OBJ := \
       $(foreach file,$(srclist_F90),OBJ/$(file:.F90=.o)) \
       $(foreach file,$(srclist_F),OBJ/$(file:.F=.o)) \
       $(foreach file,$(srclist_c),OBJ/$(file:.c=.o))
###############################################################################

CPP_DEFS = -DLINUX
ifdef usrCPP_DEFS
  CPP_DEFS := $(CPP_DEFS) $(usrCPP_DEFS)
endif	
CPP_INCS := $(cpp_path) $(usrCPP_INCS)

FC = /Users/Dfillmor/bin/gfortran
ifeq ($(useF77),TRUE)
  FC = $(fc77)
endif
FFLAGS = $(usrFFLAGS)

CC = /Users/Dfillmor/bin/gcc
CFLAGS = -DLINUX
objflag = -c -o $@

LDFLAGS = $(mpFLAGS) $(nchdflib) $(usrLDFLAGS)

###############################################################################

all : OBJdir match

clean	: 
	rm -rf OBJ
	rm -f *.mod Depends match

OBJdir :
	- if test -d OBJ; then echo OBJ exists; else mkdir OBJ; fi

match : $(OBJ)
	$(FC) -o $@ $(OBJ) $(LDFLAGS) 

# fft can't do bounds checking
OBJ/fft.o : fft.F
	$(FC) $(CPP_INCS) $(CPP_DEFS) $(specFFLAGS) $(r8FLAG)	$(objflag) $<

OBJ/advect.o : advect.F
	$(FC) $(CPP_INCS) $(CPP_DEFS) $(specFFLAGS2) $(r8FLAG)	$(objflag) $<

OBJ/advfix.o : advfix.F
	$(FC) $(CPP_INCS) $(CPP_DEFS) $(specFFLAGS2) $(r8FLAG) $(objflag) $<

OBJ/%.o : %.F
	$(FC) -c $(CPP_INCS) $(CPP_DEFS) $(FFLAGS) $(r8FLAG) $(objflag) $<

OBJ/%.o : %.c
	$(CC) $(CPP_INCS) $(CPP_DEFS) $(CFLAGS) -c -o $@ $<

OBJ/%.o : %.F90
	$(FC) -c $(CPP_INCS) $(CPP_DEFS) $(FFLAGS) $(r8FLAG) $(objflag) $<

Depends : Filepath Srcfiles
	rm -f Depends
	${makdep} -t OBJ Filepath Srcfiles > Depends

###############################################################################
# Source code dependencies:

include Depends
