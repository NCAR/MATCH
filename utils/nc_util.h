/*$Id: nc_util.h,v 1.4 1998/11/28 22:11:06 eaton Exp $*/

#ifndef NC_UTIL_H
#define NC_UTIL_H

#include "ccmhdr.h"

typedef struct var_id { /* variable ids */

  int day;     /* day of time sample */
  int daysec;  /* seconds relative to day */
  int date;    /* date of time sample */
  int datesec; /* seconds relative to date */
  int time;    /* day + daysec/86400 */
  int hyai;    /* hybrid coef A at interfaces */
  int hybi;    /* hybrid coef B at interfaces */
  int hyam;    /* hybrid coef A at midpoints */
  int hybm;    /* hybrid coef B at midpoints */
  int mlev;    /* 1000*(A+B) at midpoints */
  int lon;     /* Distinct longitudes */
  int lat;     /* Gaussian latitudes */
  int gw;      /* Gaussian weights */
  int nlons;   /* number of longitudes at each latitude (for reduced grid) */

  int nfld;    /* number of requested fields from the CCM file */
  int *fld;    /* ids for all requested fields from the CCM file */
  int *ccmfld; /* index into the Field array in the CCM_HDR structure */
               /* for the requested fields                            */

} var_id;


void
init_nc( CCM_HDR *hdr, int ncid, char* reqvar[], char* cmndline,
	 var_id *vid, int dim_ord, int coards, int doubleOut );
void
write_ncrec( FFILE ff, CCM_HDR *h, int ncid, var_id *vid,
	     size_t its, int dim_ord, int doubleOut );

char*
checkName( char* name );

int
time_units( int date, int sec, char* units );

void
switch_lev_lat( Hfloat* ccmdat, Hfloat* ncdat, int nlon, int nlat,
		int nlev );
void
handleError( int status );

void
appendHistAtt( int ncIn, int ncOut, char* cmndline );

void
ncPutVarHfloat( int ncid, int vid, const Hfloat* data );

void
ncPutVar1Hfloat( int ncid, int vid, const size_t index[], const Hfloat* data );

void
ncPutVaraHfloat( int ncid, int vid, const size_t start[], const size_t count[], const Hfloat* data );

void
ncPutAttHfloat( int ncid, int vid, const char* name, nc_type xtype, size_t len, const Hfloat* data );

#endif
