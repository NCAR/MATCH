/*$Id: ccm_util.h,v 1.1.1.1 1998/02/17 23:44:00 eaton Exp $*/

#ifndef CCM_UTIL_H
#define CCM_UTIL_H

#include "ccmhdr.h"
#include "fortio.h"

extern void print_ccm_header( CCM_HDR *hdr );

extern int print_ccm_stats( FFILE ff, CCM_HDR *hdr );

extern int accum_ccm_stats( Hfloat *buf, int nfldlev, int nlon, int nlonw,
			  double spval, double gw,
                          int *nfldval, double *sum, double *min, double *max );

extern int get_slice( char *fldname, FFILE ff, CCM_HDR *h, void *data,
                      Hfloat *slice );

extern Hfloat cv_spval( char *xrep );

extern int spval_exists( FFILE ff, Hfloat* data, int nvals, Hfloat* spval );

#endif
