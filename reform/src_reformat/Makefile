PROG = ../reformat_ed2_match.$(CPU).exe


OBJS = sage_clim.o destrat.o sageaot.o readnetcdf_aots.o readnetcdf_vertical.o
OBJS_F = 

TKLIBS = -L$(PGSLIB) -lPGSTK
CLIBS = $(CERESLIB)/data_products.a $(CERESLIB)/cereslib.a
MOD_FLAG = -I. -I$(CERESLIB)/mod
LIBS = 

INC_FLAG = -I$(NETCDFINC) -I$(PGSINC) -I$(HDFINC) -I$(HDF5INC) -I../include -I./include
VPATH = .:$(NETCDFINC).:$(PGSINC):$(HDFINC):$(HDF5INC):../include:./include

HDF_FLAGS = $(NETCDFDIR) $(HDFEOSDIR) $(HDFEOSLIB) $(HDFDIR) $(HDFLIBS) 
HDF5_FLAGS = $(HDFEOS5DIR) $(HDFEOS5LIB) $(HDF5DIR) $(HDF5LIBS) 
DAAC_FLAGS = $(ADD_LFLAGS) $(ADD_LIBS) 

#CC = 
#CFLAGS = 
#F90 = 
#F90COMP = 
#FCOMP = 
#F90LOAD = 

all: $(PROG)

$(PROG): $(OBJS) $(OBJS_F) $(LIBS) $(CLIBS) 
	$(F90) $(F90LOAD) -o $@ $(OBJS) $(OBJS_F) $(NETCDFDIR) $(NETCDFLIBS) $(LIBS) $(CLIBS) \
	$(TKLIBS) $(HDF_FLAGS) $(HDF5_FLAGS) $(DAAC_FLAGS) 

clean:
	-\rm -f $(PROG) $(OBJS) $(OBJS_F) 
	-\rm -f *.mod 

.SUFFIXES: 

.SUFFIXES: .f90 .mod .f .c .o

.f90.o:
	$(F90) $(F90COMP) $(MOD_FLAG) $(INC_FLAG) $< 

.f90.mod:
	$(F90) $(F90COMP) $(MOD_FLAG) $(INC_FLAG) $< 

.f.o:
	$(F90) $(FCOMP) $(MOD_FLAG) $(INC_FLAG) $< 

.c.o:
	$(CC) $(CFLAGS) $(MOD_FLAG) $(INC_FLAG) $< 

readnetcdf_aots.o: readnetcdf_aots.o readnetcdf_vertical.o netcdf.inc
readnetcdf_vertical.o: netcdf.inc
