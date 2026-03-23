src_ccm2nc = ccm2nc.c \
	ccmhdr.c ccm_util.c cr-ieee.c crpack.c fortio.c \
	htdatlen.c iscosb.c nc_util.c parse.c strops.c

obj_ccm2nc = $(src_ccm2nc:.c=.o)

src_hsum = hsum.c \
	ccmhdr.c ccm_util.c cr-ieee.c crpack.c fortio.c \
	htdatlen.c iscosb.c parse.c strops.c

obj_hsum = $(src_hsum:.c=.o)

CC = gcc
NETCDF_DIR = $(CONDA_PREFIX)
cpp_path = -I./ -I$(NETCDF_DIR)/include
CPPFLAGS = -DIEEE_DP $(cpp_path)
CFLAGS = -O
USRLIBS = -I$(NETCDF_DIR)/include -L$(NETCDF_DIR)/lib -lnetcdf

all : ccm2nc hsum

ccm2nc : $(obj_ccm2nc)
	$(CC) -o $@ $(obj_ccm2nc) $(USRLIBS)

hsum : $(obj_hsum)
	$(CC) -o $@ $(obj_hsum) -lm

clean :
	rm *.o ccm2nc hsum

