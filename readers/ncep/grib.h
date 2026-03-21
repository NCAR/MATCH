/* $Id$ */
#ifndef GRIB_WAS_INCLUDED
#define GRIB_WAS_INCLUDED

extern struct Dset_t *GRIBopen_dataset( char *fname, char *options);
extern void GRIBclose_dataset( Dset_t *def);

extern int GRIBget_data( Dset_t *def, ARR_t *arr);
extern double GRIBget_value( Dset_t *def, int time, int species, int xidx, int yidx, 
		int zidx);

extern void GRIBfile_info( FILE *fp, Dset_t *def);

#endif
