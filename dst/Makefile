# $Header$

# Purpose: Makefile for aerosol module aer

# Known problems: 
# Solaris f90 does not support -r8 -i4 so netCDF output is garbage
# Solution: None, really

# Use separate directory for aer objects because they are often 
# compiled with -r8 -i4 and thus incompatible with, e.g., map
ifndef PVM_ARCH
 PVM_ARCH := $(shell $(MY_BIN_DIR)/sh/pvmgetarch)
endif
MY_OBJ_DIR := ./obj/$(PVM_ARCH)
MY_BLD_DIR := .
MY_DAT_DIR := .
MY_DOC_DIR := .
MY_DPN_DIR := $(MY_OBJ_DIR)
MY_INC_DIR := .
MY_SRC_DIR := .
ifndef MY_BIN_DIR
 MY_BIN_DIR := ../bin
endif
ifndef MY_LIB_DIR
 MY_LIB_DIR := ../lib
endif
ifndef MY_OBJ_DIR
 MY_OBJ_DIR := ../obj
endif
ifndef OPTS
 OPTS := O
endif
ifndef NETCDF_INC
 NETCDF_INC := /usr/local/include # Directory containing netcdf.h
endif
ifndef NETCDF_LIB
 NETCDF_LIB := /usr/local/lib # Directory containing libnetcdf.a
endif

# NB: Do NOT add comment lines, e.g., # This is a comma, to character definitions
null :=
space := $(null) $(null)
comma := ,
newline := \n
# We replace `/' by `cszzsc', call perl, then replace `cszzsc' by `/'
# Unique character(s) to substitute for $(slash) before passing to perl regex
slash_rx := cszzsc
# Unique character(s) to replace by $(slash_rx) before passing to perl regex
slash := /
MY_OBJ_DIR_RX := $(subst $(slash),$(slash_rx),$(MY_OBJ_DIR))
MY_DPN_DIR_RX := $(subst $(slash),$(slash_rx),$(MY_DPN_DIR))

# Directories to search for source files
MDL_PTH := $(MY_SRC_DIR)
# Find all C++, C, and Fortran files in a given directory
FIND_FNC = $(wildcard $(dir)/*.cc $(dir)/*.c $(dir)/*.F)
# Assemble source files from all directories
SRC_LST = $(foreach dir, $(MDL_PTH),$(FIND_FNC))
# Source file names with directories removed
MDL_SRC := $(notdir $(SRC_LST))
# Dependency list for executable
MDL_OBJ := $(addprefix $(MY_OBJ_DIR)/,$(addsuffix .o, $(basename $(MDL_SRC)))) 
# Dependency (make) file for each object file
MDL_DPN := $(addprefix $(MY_DPN_DIR)/,$(addsuffix .d, $(basename $(MDL_SRC)))) 
# VPATH helps make find dependencies (which are not pathname qualified) in *.d file
VPATH := $(subst $(space),:,$(MDL_PTH))
# Prepend -I to use for compiler argument
CPP_PTH := $(foreach dir,$(MDL_PTH),-I$(dir))

# Variables having to do with binary executables created by module
MDL_BIN_TRG := aer # Binary targets
MDL_BIN_SYM_LNK := # Symbolic links
MDL_BIN_STB := $(MDL_BIN_TRG) $(MDL_BIN_SYM_LNK) # All C++ files in MY_BIN_DIR
MDL_BIN := $(addprefix $(MY_BIN_DIR)/,$(MDL_BIN_STB)) # dst_cln removes these files

# Variables having to do with header files created by module
MDL_INC_TRG := # Include targets
MDL_INC_SYM_LNK := # Symbolic links
MDL_INC_STB := $(MDL_INC_TRG) $(MDL_INC_SYM_LNK) # All header files in MY_INC_DIR
MDL_INC := $(addprefix $(MY_INC_DIR)/,$(MDL_INC_STB)) # dst_cln removes these files

# Redefine default C and C++ pattern rules
$(MY_OBJ_DIR)/%.o : %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $(MY_OBJ_DIR)/$(notdir $@)
$(MY_OBJ_DIR)/%.o : %.cc
	$(C++) $(CPPFLAGS) $(C++FLAGS) -c $< -o $(MY_OBJ_DIR)/$(notdir $@)

# Default Fortran pattern rules: CRAY and RS6K must override these rules
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<

# Automatically generate a dependency file for each source file.
# $* is the stem, e.g., f
# $@ is the filename of the target, e.g., f.d
# Linux gcc may return an extra `.F' on Fortran names, e.g., `hello.F.o: hello.F'
# (.F)? gets rid of this extra `.F'
$(MY_DPN_DIR)/%.d : %.F
# Following command makes, e.g., f.d begin "f.o f.d : f.F ..."
# Since f.o is not preceded by $(MY_OBJ_DIR), objects are not recompiled when sources are touched.
#	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*\.F\.o/$*.o $@/g;' > $@
# Following command makes, e.g., f.d begin "/home/zender/obj/LINUX/f.o f.d : f.F ..."
# This works fairly well, but is a hack.
# First pattern substitutes MY_OBJ_DIR_RX, which has placeholders for slashes.
# Second pattern substitutes slashes for the placeholders.
	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*(\.F)?\.o/$(MY_OBJ_DIR_RX)\/$*.o $(MY_DPN_DIR_RX)\/$(notdir $@)/g;s/$(slash_rx)/\$(slash)/g' > $@
# Following command makes, e.g., f.d begin "$(MY_OBJ_DIR)/f.o f.d : f.F ..."
# This would be the ideal command but I can't get the dollar sign to show up.
#	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*\.F\.o/\$(dollar)MY_OBJ_DIR\/$*.o $@/g;' > $@
$(MY_DPN_DIR)/%.d : %.f
	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*(\.F)?\.o/$(MY_OBJ_DIR_RX)\/$*.o $(MY_DPN_DIR_RX)\/$(notdir $@)/g;s/$(slash_rx)/\$(slash)/g' > $@

$(MY_DPN_DIR)/%.d : %.c
#	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*\.o/$*.o $@/g;' > $@
	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*\.o/$(MY_OBJ_DIR_RX)\/$*.o $(MY_DPN_DIR_RX)\/$(notdir $@)/g;s/$(slash_rx)/\$(slash)/g' > $@

$(MY_DPN_DIR)/%.d : %.cc
#	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*\.o/$*.o $@/g;' > $@
	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*\.o/$(MY_OBJ_DIR_RX)\/$*.o $(MY_DPN_DIR_RX)\/$(notdir $@)/g;s/$(slash_rx)/\$(slash)/g' > $@

# The first LDFLAGS is for typical C programs with netCDF, math, and networking
# The second LDFLAGS enables C/Fortran linking

ifeq ($(PVM_ARCH),AIX)
C++ := xlC
CC := xlc_r
CPP := /lib/cpp -P
CPPFLAGS := -DDST -DBXM -DCCM -D$(PVM_ARCH) $(CPP_PTH) -I$(NETCDF_INC)
FC := xlf90_r
FFLAGS := -qarch=pwr3 -qrealsize=8 -qdpc=e -qsmp=noauto -qspillsize=2500
LD := ld
LDFLAGS := -qsmp=noauto -bmaxdata:0x80000000 -L$(MY_LIB_DIR) -lcsz_f77 -L$(NETCDF_LIB) -lncaru -lnetcdf
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS :=  -O2
 FFLAGS += -O3 -qstrict
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS += -g
endif
$(MY_DPN_DIR)/%.d : %.F
	$(CPP) $(CPPFLAGS) $< > /dev/null;perl -p -e 's/$*(\.F)?\.o/$(MY_OBJ_DIR_RX)\/$*.o $(MY_DPN_DIR_RX)\/$(notdir $@)/g;s/$(slash_rx)/\$(slash)/g' $*.u > $@;/bin/rm $*.u
$(MY_DPN_DIR)/%.d : %.F90
	$(CPP) $(CPPFLAGS) $< > /dev/null;perl -p -e 's/$*(\.F)?\.o/$(MY_OBJ_DIR_RX)\/$*.o $(MY_DPN_DIR_RX)\/$(notdir $@)/g;s/$(slash_rx)/\$(slash)/g' $*.u > $@;/bin/rm $*.u
$(MY_OBJ_DIR)/%.o : %.F90
	$(FC) $(CPPFLAGS) -c $(FFLAGS) -qsuffix=f=f90:cpp=F90 -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) $(CPPFLAGS) -c $(FFLAGS) -qsuffix=f=f:cpp=F -qfixed=132 -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif
# endif AIX

ifeq ($(PVM_ARCH),CRAY)
C++ := g++
CC := cc
CPP := cpp
CPPFLAGS := -DDST -DBXM -DCCM -D$(PVM_ARCH) $(CPP_PTH) -I$(NETCDF_INC)
FC := f90
LD := ld
LDFLAGS := -L$(NETCDF_LIB) -lnetcdf
LDFLAGS += -L/lib -lf
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O
 FFLAGS = -N 132
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS = -g -N 132
endif
ifeq ($(OPTS),X)
 CFLAGS := -g -N 132
 FFLAGS = -g -N 132 -e i
endif
# 1997/10/21 Added -P to suppress #line # directives on Fortran files
# 1999/05/24 Added -C for cpp to retain comments (else cpp thinks Fortran string concatenator // is a C++ comment and removes it)
$(MY_OBJ_DIR)/%.o : %.F
	$(CPP) -C -P $(CPPFLAGS) $< > $(patsubst %.F,%.f,$(notdir $<))
	$(FC) -c $(FFLAGS) $(patsubst %.F,%.f,$(notdir $<)) 
#	$(FC) -e p -c $(FFLAGS) $(notdir $<) 
#	$(FC) -e p -Wp"$(CPPFLAGS)" -c $(FFLAGS) $(notdir $<) 
	-mv -f $(notdir $@) $(MY_OBJ_DIR)
	rm -f $(patsubst %.F,%.f,$(notdir $<)) 
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) $<
	mv -f $(notdir $@) $(MY_OBJ_DIR)
endif
# endif CRAY

ifeq ($(PVM_ARCH),LINUX)
C++ := g++
CC := gcc
CPPFLAGS := -DDST -DBXM -D$(PVM_ARCH) $(CPP_PTH) -I$(NETCDF_INC)
FC := pgf90
LD := ld
LDFLAGS := -mp -L$(NETCDF_LIB) -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq ($(CC),gcc)
ifeq ($(OPTS),O)
 CFLAGS := -O -Wall
endif
ifeq ($(OPTS),D)
 CFLAGS := -g -Wall
endif
ifeq ($(OPTS),R)
 CFLAGS := -Wall
endif
ifeq ($(OPTS),X)
 CFLAGS := -g -O -Wall
endif
 C++FLAGS := $(CFLAGS)
endif
ifeq ($(FC),pgf90)
# NB: -Mrecursive ensures that memory is allocated on the stack
ifeq ($(OPTS),O)
 FFLAGS := -fast -mp -Mr8 -Mi4 -byteswapio -Mextend -Mnosecond_underscore -Mrecursive -Mdalign
endif
ifeq ($(OPTS),D)
 FFLAGS := -g -mp -Mr8 -Mi4 -byteswapio -Mextend -Mnosecond_underscore -Mrecursive -Mdalign
endif
ifeq ($(OPTS),R)
 FFLAGS := -mp -Mr8 -Mi4 -byteswapio -Mextend -Mnosecond_underscore -Mrecursive -Mdalign
endif
ifeq ($(OPTS),X)
 FFLAGS := -g -Mbounds -mp -Mr8 -Mi4 -byteswapio -Mextend -Mnosecond_underscore -Mrecursive -Mdalign
endif
endif
ifeq ($(FC),g77)
ifeq ($(OPTS),O)
 FFLAGS := -O -ffixed-line-length-132 -fno-second-underscore
endif
ifeq ($(OPTS),D)
 FFLAGS := -g -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge
endif
ifeq ($(OPTS),R)
 FFLAGS := -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge
endif
ifeq ($(OPTS),X)
 FFLAGS := -g -O -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge
endif
endif
endif
# endif LINUX

ifeq ($(PVM_ARCH),RS6K)
C++ := g++
CC := gcc -ansi
CPP := /lib/cpp -P
CPPFLAGS := -D$(PVM_ARCH) $(CPP_PTH) -I$(NETCDF_INC)
FC := xlf
LD := ld
LDFLAGS := -L$(MY_LIB_DIR) -lcsz_f77 -L$(NETCDF_LIB) -lncaru -lnetcdf
LDFLAGS += -lxlf90 -lxlf
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 CPP := $(CPP) $(CPPFLAGS)
 PREPROCESS.F := $(CPP) $(CPPFLAGS)
 FFLAGS := -O -NS2000 -qfixed=132
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 CPP := $(CPP) $(CPPFLAGS)
 PREPROCESS.F := $(CPP) $(CPPFLAGS)
 FFLAGS := -g -NS2000 -qfixed=132
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(CPP) $(CPPFLAGS) $< $(MY_OBJ_DIR)/$(basename $<).f 
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $(MY_OBJ_DIR)/$(basename $<).f
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif
# endif RS6K

ifeq ($(PVM_ARCH),SGI5)
C++ := g++
CC := gcc -ansi
CPPFLAGS := -D$(PVM_ARCH) $(CPP_PTH) -I$(NETCDF_INC)
FC := f77
LD := ld
LDFLAGS := -L$(MY_LIB_DIR) -lcsz_f77 -L$(NETCDF_LIB) -lnetcdf
LDFLAGS += -lF77 -lI77 -lU77 -lftn
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 FFLAGS := -O2 -e -Nl99 -extend_source
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS := -g -e -Nl99 -extend_source
endif
endif
# endif SGI5

ifeq ($(PVM_ARCH),SGI64)
C++ := c++
CC := cc
CPP := cpp
CPPFLAGS := -D$(PVM_ARCH) $(CPP_PTH) -I$(NETCDF_INC)
FC := f77
LD := ld
LDFLAGS := -64 -mips4 -L$(MY_LIB_DIR) -lcsz_f77 -L$(NETCDF_LIB) -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq ($(OPTS),O)
 CFLAGS := -64 -mips4 -O2
 FFLAGS := -64 -mips4 -O2 -extend_source
endif
ifeq ($(OPTS),R)
 CFLAGS := -64 -mips4
 FFLAGS := -64 -mips4 -extend_source
endif
ifeq ($(OPTS),D)
 CFLAGS := -64 -mips4 -g
 FFLAGS := -64 -mips4 -g -extend_source
endif
ifeq ($(OPTS),X)
 CFLAGS := -64 -mips4 -g -trapuv
 FFLAGS := -64 -mips4 -g -extend_source -check_bounds -trapuv
endif
endif
# endif SGI64

ifeq ($(PVM_ARCH),SGIMP64)
C++ := CC
CC := cc
CPP := cpp
# 19990525: Added -DSGI to emulate CCM_TRG_ARCH flag
CPPFLAGS := -DDST -DBXM -DCCM -D$(PVM_ARCH) $(CPP_PTH) -I$(NETCDF_INC)
FC := f90 -cpp
#FC := f77
LD := ld
LDFLAGS := -64 -mips4 -L$(NETCDF_LIB) -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq ($(OPTS),O)
 CFLAGS := -64 -mips4 -O2
 FFLAGS := -64 -mips4 -O2 -extend_source
endif
ifeq ($(OPTS),R)
 CFLAGS := -64 -mips4
 FFLAGS := -64 -mips4 -extend_source
endif
ifeq ($(OPTS),D)
 CFLAGS := -64 -mips4 -g
 FFLAGS := -64 -mips4 -g -extend_source
endif
ifeq ($(OPTS),X)
 CFLAGS := -64 -mips4 -g -trapuv
 FFLAGS := -64 -mips4 -g -extend_source -check_bounds -trapuv
endif
 FFLAGS += -r8 -i4
endif
# endif SGIMP64

ifeq ($(PVM_ARCH),SUN4)
C++ := g++
CC := acc
CPP := cpp
CPPFLAGS := -D$(PVM_ARCH) $(CPP_PTH) -I$(NETCDF_INC)
FC := f77
LD := ld
LDFLAGS := -L$(MY_LIB_DIR) -lcsz_f77 -L$(NETCDF_LIB) -lnetcdf
LDFLAGS += -cg92 -L/opt/SUNWspro/SC3.0/lib/cg92 -lF77 -lM77 -lresolv
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 FFLAGS := -fast -e -Nl99
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS := -g -e -Nl99
endif
 FFLAGS += -r8 -i4
endif
# endif SUN4

ifeq ($(PVM_ARCH),SUN4SOL2)
C++ := g++
CC := gcc -ansi
CPP := cpp
CPPFLAGS := -D$(PVM_ARCH) $(CPP_PTH) -I$(NETCDF_INC)
FC := f90
LD := ld
LDFLAGS := -L$(MY_LIB_DIR) -lcsz_f77 -L$(NETCDF_LIB) -lnetcdf -lsunmath -lsocket -lnsl
LDFLAGS += -lF77 -lM77 -lresolv
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 FFLAGS := -O -e
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS := -g -e
endif
ifeq ($(OPTS),X)
 CFLAGS := -g 
 FFLAGS := -g -e
# 1998/06/01 -C (range-checking) is not supported by Sun f90
ifeq ($(FC),f77)
 FFLAGS += -C
endif
endif
# 1999/09/07 -r8 -i4 are not supported by Sun f90
ifeq ($(FC),f77)
 FFLAGS += -r8 -i4
endif
endif
# endif SUN4SOL2

ifeq ($(PVM_ARCH),SUNMP)
C++ := g++
#CC := gcc -ansi -pedantic
CC := gcc 
CPP := cpp
CPPFLAGS := -DDST -DBXM -DCCM -D$(PVM_ARCH) $(CPP_PTH) -I$(NETCDF_INC)
# 1998/09/10 f90 is default compiler for aer routines
FC := f90
LD := ld
LDFLAGS := -L$(NETCDF_LIB) -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq ($(OPTS),O)
 CFLAGS := -O2
 FFLAGS := -O -e
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS := -g -e
endif
ifeq ($(OPTS),X)
 CFLAGS := -g 
 FFLAGS := -g -e
# 1998/06/01 -C (range-checking) is not supported by Sun f90
ifeq ($(FC),f77)
 FFLAGS += -C
endif
endif
# 1999/09/07 -r8 -i4 are not supported by Sun f90
ifeq ($(FC),f77)
 FFLAGS += -r8 -i4
endif
endif
# endif SUNMP

ifeq ($(PVM_ARCH),WIN32)
C_ONLY := Y
C++ := g++
CC := gcc -ansi
# NB: nameser.h needs -Di386, but gcc is sending -Di586 (on pentiums)
CPPFLAGS := -D$(PVM_ARCH) -Di386 $(CPP_PTH) -I$(NETCDF_INC)
FC := g77
LD := ld
LDFLAGS := -L$(MY_LIB_DIR) -lcsz_f77 -L$(NETCDF_LIB) -lnetcdf
LDFLAGS += 
LEX := flex
LINT := lint
YACC := bison
ifeq ($(OPTS),O)
 CFLAGS := -O
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
endif
endif
# endif WIN32

# Define any remaining variables

# Default targets
all: dir $(MDL_BIN_TRG)
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln debug 
# Delete default suffixes---this should increase speed
.SUFFIXES: 
# Define suffixes which matter
.SUFFIXES: .cc .c .o .F .d
# Delete targets which were not successfully made
.DELETE_ON_ERROR:
# Target directories which may not exist
dir: bin_dir obj_dir 
bin_dir:
	-mkdir -p $(MY_BIN_DIR)
obj_dir:
	-mkdir -p $(MY_OBJ_DIR)

# Targets in bin
aer : $(MDL_OBJ)
	$(FC) -o $(MY_BIN_DIR)/$@ $(MDL_OBJ) $(LDFLAGS)
	chmod 755 $(MY_BIN_DIR)/$@
bin : $(MDL_BIN_TRG)
	rm -f $(MDL_BIN)
binclean: bin_cln
bin_cln: 
	rm -f $(MDL_BIN)
# Targets in bld
buildclean: bld_cln
bld_cln :
	cd $(MY_BLD_DIR); rm -f TAGS
# Targets in dat
data:	dat
dat_cln :
	-cd $(MY_DAT_DIR); rm -f aer.nc dst_mss_bdg.nc
# Targets in dpn
depend: dpn
dpn : $(MDL_DPN)
dpn_cln :
	rm -f $(MDL_DPN)
# Targets in inc
include: inc
inc : 
	cp -f -p $(MDL_INC_TRG) $(MY_INC_DIR)
inc_cln :
	rm -f $(MDL_INC)
# Targets in lib
library: lib
lib :
	-cd $(HOME)/f; gmake lib
libclean: lib_cln
lib_cln :
	-cd $(MY_LIB_DIR);
# Targets in obj
object: obj
obj : $(MDL_OBJ)
objclean: obj_cln
obj_cln :
	rm -f $(MDL_OBJ)

# Housekeeping
clean : cln
cln : dat_cln dpn_cln obj_cln
debug: dbg
dbg:
	printf "MDL_SRC = $(MDL_SRC)\n"
	printf "MDL_OBJ = $(MDL_OBJ)\n"
	printf "MDL_DPN = $(MDL_DPN)\n"
	printf "VPATH = $(VPATH)\n"
distclean : dst_cln
dst_cln : cln bld_cln
install : inc 
tags : 
	etags *.com *.F *.h Makefile
uninstall : dst_cln bin_cln inc_cln lib_cln

# Create dependency files only if they will not be immediately deleted
INCLUDE_DPN := TRUE
GOALS_WHICH_DELETE_DEPENDENCY_FILES := cln clean dir distclean dst_cln dpn_cln tags uninstall
ifeq (,$(findstring $(MAKECMDGOALS),$(GOALS_WHICH_DELETE_DEPENDENCY_FILES)))
 INCLUDE_DPN := TRUE
else
 INCLUDE_DPN := FALSE
endif
ifeq ($(INCLUDE_DPN),TRUE)
-include $(MDL_DPN)
endif
