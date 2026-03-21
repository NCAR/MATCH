/*$Id: ccm_util.c,v 1.3 1998/04/07 22:32:48 eaton Exp $*/

/* Routines for manipulating ccm history tape format files. */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>

#include "ccm_util.h"
#include "fortio.h"
#include "crpack.h"

#define NINTS 32    /* no. of ints read from header to get size info */
#define SPVAL 1.e36 /* special value used by ccm processor */

/*===========================================================================*/


void print_ccm_header( CCM_HDR *hdr )
{
  /* Print ccm history tape header */

  int i;

  /* If not first time sample, print only day and date info */

  if ( hdr->ts_count > 1 ) {

    printf(
    "\nTime sample no. %d: ndcur = %d, nscur = %d, ncdate = %d, ncsec = %d\n",
     hdr->mfilh, hdr->ndcur, hdr->nscur, hdr->ncdate, hdr->ncsec );

  } else {  /* This branch extends to the end of the procedure */

  /* Case and Title */

  if ( hdr->mftyp%10 == 2 ) {
    printf("\n    *** HEADER FOR CCM1 HISTORY TAPE ***\n\n");
  } else if ( hdr->mftyp%10 == 3 ) {
    printf("\n    *** HEADER FOR CCM2 HISTORY TAPE ***\n\n");
  }
  printf(" CASE:  %.8s\n", hdr->mcase );
  printf("TITLE:  %.80s\n", hdr->mcstit );

  /* Integer header except for field list info */

  printf(
" LENHDI  MFTYP   MFILH   MFILTH  NRBD  MAXSIZ   NDAVU   NSPHER   NLON   NLONW\n" );
  printf( "%5d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
  	  hdr->lenhdi,  hdr->mftyp,  hdr->mfilh,  hdr->mfilth,  hdr->nrbd,
	  hdr->maxsiz,  hdr->ndavu,  hdr->nspher,  hdr->nlon,  hdr->nlonw  );

  printf(
" NOREC    NLEV    NTRM    NTRN   NTRK   NFLDH   NSTEPH  NSTPRH  NITSLF  NDBASE\n");
  printf( "%5d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
	  hdr->norec,  hdr->nlev,  hdr->ntrm,  hdr->ntrn,  hdr->ntrk,
          hdr->nfldh,  hdr->nsteph,  hdr->nstprh,  hdr->nitslf,
	  hdr->ndbase );

  printf(
" NSBASE  NDCUR   NSCUR NBDATE   NBSEC  NCDATE   NCSEC     MDT    MHISF  MFSTRT\n");
  printf( "%5d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
	  hdr->nsbase,  hdr->ndcur,  hdr->nscur,  hdr->nbdate,
	  hdr->nbsec,  hdr->ncdate,  hdr->ncsec,  hdr->mdt,
	  hdr->mhisf,  hdr->mfstrt );

  printf(" LENHDC LENHDR\n");
  printf( "%5d%8d\n\n",  hdr->lenhdc,  hdr->lenhdr );

  /* Character header except for field list info */

  printf("             MSS PATH NAME                             "
	 "DATE     TIME    SEQ NO.\n" );

  if ( hdr->mftyp%10 == 2 ) {

    printf("CURRENT     %-40s %-8s %-8s %-8s\n",
         hdr->lnhstc, hdr->ldhstc, hdr->lthstc, hdr->lshstc );

    printf("PREVIOUS    %-40s %-8s %-8s %-8s\n",
         hdr->lnhstp, hdr->ldhstp, hdr->lthstp, hdr->lshstp );

    printf("FIRST       %-40s %-8s %-8s %-8s\n",
         hdr->lnhstf, hdr->ldhstf, hdr->lthstf, hdr->lshstf );

    printf("INITIALIZED %-40s %-8s %-8s %-8s\n",
         hdr->lnhsti, hdr->ldhsti, hdr->lthsti, hdr->lshsti );

    printf("ANALYZED    %-40s %-8s %-8s %-8s\n",
         hdr->lnhsta, hdr->ldhsta, hdr->lthsta, hdr->lshsta );

  } else if ( hdr->mftyp%10 == 3 ) {

    printf("CURRENT      %-40s %-8s %-8s %-8s\n",
         hdr->lnhstc, hdr->ldhstc, hdr->lthstc, hdr->lshstc );

    printf("FIRST        %-40s %-8s %-8s %-8s\n",
         hdr->lnhstf, hdr->ldhstf, hdr->lthstf, hdr->lshstf );

    printf("INITIAL      %-40s %-8s %-8s %-8s\n",
         hdr->lnhsti, hdr->ldhsti, hdr->lthsti, hdr->lshsti );

    printf("TI BOUNDARY  %-40s %-8s %-8s %-8s\n",
         hdr->lnhstt, hdr->ldhstt, hdr->lthstt, hdr->lshstt );

    printf("SST BOUNDARY %-40s %-8s %-8s %-8s\n",
         hdr->lnhsts, hdr->ldhsts, hdr->lthsts, hdr->lshsts );

    printf("OZ BOUNDARY  %-40s %-8s %-8s %-8s\n",
         hdr->lnhsto, hdr->ldhsto, hdr->lthsto, hdr->lshsto );

  }


  /* Field list */

  printf("\nFIELD LIST\n");
  printf("FLD NO.   NAME       LEV.  TIM.      OFF.  PACK.  UNITS\n");
  for ( i=0; i < hdr->nfldh; ++i )
    printf("   %3d    %-8.8s    %2d    %2d  %8d      %1d    %-8.8s\n",  i+1,
           (hdr->fld+i)->name, (hdr->fld+i)->lev, (hdr->fld+i)->tim,
	   (hdr->fld+i)->off, (hdr->fld+i)->pdens, (hdr->fld+i)->units );


  /* Vertical levels */

  if ( hdr->mphya  ==  NULL ) {                         /* Sigma */
    printf("\nSigma Levels\n");
    for ( i=0; i<2*hdr->nlev+1; ++i )
      printf(" %2d   %8.6f\n", i+1, *(hdr->mphyab+i) );
  } else {                                              /* Hybrid */
    printf("\nHybrid Levels\n");
    printf("  N    A(P0)       B(Ps)       Hlevel=A+B\n");

    for ( i=0; i<2*hdr->nlev+1; ++i )
      printf(" %2d   %8.6f    %8.6f     %8.6fH\n", i+1, *(hdr->mphya+i),
             *(hdr->mphyb+i), *(hdr->mphyab+i) );
  }

  /* Latitudes */

  printf("\nLatitudes\n");

  for ( i=0; i<hdr->norec; ++i ) {
    printf(" %9.5f", *(hdr->mplat+i) );
    if ( (i+1)%8  == 0 ) printf("\n");
  }
  if ( hdr->norec%8 != 0 ) printf("\n");

  /* Gaussian Weights */

  printf("\nGaussian Weights\n");

  for ( i=0; i<hdr->norec; ++i ) {
    printf(" %11.5e", *(hdr->mpwts+i) );
    if ( (i+1)%6  == 0 ) printf("\n");
  }
  if ( hdr->norec%6 != 0 ) printf("\n");

  printf("\n    *** END OF HEADER DATA ***\n");

  }

}                                                        /* print_ccm_header */
/*===========================================================================*/


int print_ccm_stats( FFILE ff, CCM_HDR *h )
{
  /* Read the data records in a ccm history tape for a single time sample
     and print summary of minimum, maximum, and mean values for each field.
     It is assumed that the file is positioned at the beginning of the
     first data record.

     Input args:

     ff    Input stream.
     h     CCM history tape header.

     Return values:

      0  Successful return.
     -1  malloc: Error allocating memory for statistics.
     -2  malloc: Error allocating memory for latitude slice.
     -3  malloc: Error allocating memory for data record.
     -4  read_rec: Error reading data record.
     -5  error returned by get_slice.

  */

  double *stat;       /* Array to hold statistics */
  int *nfldval;       /* number of valid field values for each field */
  Hfloat *fldslice;   /* Array to hold a latitude slice of data. */
  void *drec;         /* data record */
  Field *fldm;        /* pointer to current field structure */
  int nrecbytes;      /* number of bytes in a data record */
  int i, m, nfldlev;
  Hfloat fp_lat;      /* latitude index of a data record */
  int lat;            /* latitude index */
  Hfloat spval;       /* 1.e36 converted from external to internal rep */
  int *oorval;        /* flag out-of-range values for each field */
  int rval;           /* function return value */
  int err = 0;        /* error flag */
  double gw;          /* Gauss weight */

  /* Convert the special value 1.e36 from its external representation, i.e.,
     that used in the data file, to its internal one.  We are assuming that
     if the data contains special values, then it must be unpacked. */

  if ( ff.cosblk == 1 )
    spval = cv_spval( "cray" );
  else if ( ff.fsz == 4 )
    spval = cv_spval( "ieee_sp" );
  else if ( ff.fsz == 8 )
    spval = cv_spval( "ieee_dp" );

  /* Allocate memory for statistics */
  
  if ( (stat = (double *) malloc( sizeof(double)*h->nfldh*3 )) == NULL )
    return -1;

  /* Allocate memory to store number of field values for each field */

  if ( (nfldval = (int *) malloc( sizeof(int)*h->nfldh )) == NULL )
    return -1;

  /* Allocate memory to flag out-of-range values for each field */

  if ( (oorval = (int *) malloc( sizeof(int)*h->nfldh )) == NULL )
    return -1;

  /* Allocate memory for a slice of data (converted to local representaion) */

  fldslice = (Hfloat *) malloc( sizeof(Hfloat)*h->nlonw*(h->nlev+1) );
  if ( fldslice == NULL ) return -2;

  /* Allocate memory for a data record (no conversion from external rep) */

  if ( (h->fld)->pdens == 1 ) {
    nrecbytes = ff.fsz*h->maxsiz;
  } else {
    nrecbytes = 8*h->maxsiz;
  }
  if ( (drec = malloc( nrecbytes )) == NULL ) return -3;

  /* Initialize variables to hold data about each field */

  for ( m = 0; m < h->nfldh; ++m ) {
    *(nfldval+m) = 0;                     /* number of values */
    *(oorval+m) = 0;                      /* out-of-range values flag */
    *(stat+3*m) = 0.0;                    /* sum */
    *(stat+3*m+1) = DBL_MAX;              /* min */
    *(stat+3*m+2) = -DBL_MAX;             /* max */
  }

  /* Loop over latitudes */

  for ( i = 0; i < h->norec; ++i ) {

    if ( (rval=read_rec( ff, drec, nrecbytes )) != nrecbytes ) {
      printf("print_ccm_stats: read_rec returned %d at data record %d\n",
              rval, i+1 );
      printf("                 %d bytes were requested\n", nrecbytes );
      return -4;
    }

    /* Get the latitude index from first word in data record */

    if ( need_fcv( ff ) ) {
      if ( cv_floats( ff, 1, drec, &fp_lat ) != 0 ) return -7;
      lat = fp_lat + .1;
    } else {
      lat = ( *( (Hfloat *)drec ) + .1 );
    }
    gw = *(h->mpwts+lat-1);

    /* Loop over fields */

    for ( m = 0; m < h->nfldh; ++m ) {

      fldm = h->fld + m;

      rval = get_slice( fldm->name, ff, h, drec, fldslice );
      if ( rval == 1 ) ++*(oorval+m);
      else if ( rval < 0 ) return -5;

      /* Determine number of levels in field */

      switch ( fldm->lev ) {
      case 0:
	nfldlev = 1; break;
      case 1:
	nfldlev = h->nlev; break;
      case 2:
	nfldlev = h->nlev; break;
      }

      /* accumulate sum of values and update min and max values */

      accum_ccm_stats( fldslice, nfldlev, h->nlon, h->nlonw, spval, gw,
		       nfldval+m, stat+3*m, stat+3*m+1, stat+3*m+2 );

    } /* end loop over fields */


  } /* end loop over latitudes */


  /* Compute mean field values */

  for ( m = 0; m < h->nfldh; ++m ) {
    *(stat+3*m) /=  2.;              /* mean (2 = sum of gw) */
  }

  /* Print field statistics */

  printf("Field         area mean               min                 max\n");

  for ( m = 0; m < h->nfldh; ++m ) {
    fldm = h->fld + m;
    printf("%-8.8s %19.12e %19.12e %19.12e\n", fldm->name,
	          *(stat+3*m), *(stat+3*m+1), *(stat+3*m+2) );
  }

  /* Print diagnostics */

  for ( m = 0; m < h->nfldh; ++m ) err += *(oorval+m);
  
  if ( err > 0 ) {
    printf(
    "\n\nThe following fields contained values for which errors were encountered"
    "\nin the routines that perform IEEE-CRAY conversions.  These errors have"
    "\nbeen resolved as follows.  When converting from CRAY to IEEE, underflow"
    "\nvalues are set to zero and overflow values are set to the largest possible"
    "\nIEEE value with the correct sign.  When converting from IEEE to CRAY, if"
    "\nIEEE values are found that represent NaN's or +-infinity.  These values"
    "\nare set to the largest possible CRAY representation with the sign bit"
    "\npreserved.  Denormalized IEEE numbers are set to zero.\n");

    for ( m = 0; m < h->nfldh; ++m ) {
      fldm = h->fld + m;
      if ( *(oorval+m) > 0 ) printf("%s\n", fldm->name );
    }
  }

  /* Free memory */

  free( stat );
  free( fldslice );
  free( drec );

  return 0;
}                                                         /* print_ccm_stats */
/*===========================================================================*/


#include <math.h>

int accum_ccm_stats( Hfloat *buf, int nfldlev, int nlon, int nlonw, 
		     double spval, double gw,
		     int *nfldval, double *sum, double *min, double *max )
{
  /* Update the values of sum, min, and max.  Sum is incremented by the sum
     of the values in buf; min and max are updated if a smaller or larger value
     is found among n values in buf.

     It is assumed that buf points to a latitude slice in ccm history tape
     format.  nlonw values are written to each latitude circle, but only
     the first nlon values correspond to model data.
  */

  int i, k;
  double fval, absf;
  double latsum;
  int npts;

  latsum = 0.;
  npts = 0;
  for ( k = 0; k < nfldlev; ++k ) {

    for ( i = 0; i < nlon; ++i ) {

      fval = *(buf+k*nlonw+i);
      if ( fval != spval ) {
	++(*nfldval);
	npts++;
/*	absf = hypot( fval, 0.0 ); */
	latsum += fval;
	if ( fval < *min ) *min = fval;
	if ( fval > *max ) *max = fval;
      }
    }
  }
  if ( npts != 0 ) latsum *= gw/npts;    /* for area weighting */
  *sum += latsum;

  return 0;
}                                                         /* accum_ccm_stats */
/*===========================================================================*/


int get_slice( char *fldname, FFILE ff, CCM_HDR *h, void *data, Hfloat *slice )
{
  /* Extract the data for a specified field from a CCM history tape data
     record.  The data record may contain foreign or packed data.

     Input arguments:

     fldname  Name of field.
     ff       History tape input stream.
     h        History tape header structure.
     data     History tape data record.

     Output arguments:

     slice    Requested field data.

     Return values:

      1  out-of-range values encountered by CRAY-IEEE conversion routine
         (underflow set to zero and overflow set to largest ieee values)
      0  for success.
     -1  error returned by unpack.
     -2  error returned by cv_floats.
     -3  requested field not found.
  */

  Field *fp;    /* pointer to a Field structure in the header */
  char *dp;     /* pointer to start of requested slice in the input data */
  Hfloat *in;   /* pointer to start of requested slice in the input data */
  Hint nfloats; /* number of floats in requested data slice */
  Hfloat *out;  /* pointer into slice */
  int npack;    /* number of 8-byte words in packed latitude circle of data */
  int fldlevs;  /* number of levels in field */
  int n, i, k, stat;
  int oorval = 0; /* cray values found that are out-of-range for ieee */

  for ( n = 0; n < h->nfldh; ++n ) {

    fp = h->fld + n;

    if ( strcmp( fp->name, fldname ) == 0 ) {  /* field found */

      if ( fp->lev == 0 ) {
	nfloats = h->nlon;
	fldlevs = 1;
      } else {
	nfloats = h->nlev * h->nlon;
	fldlevs = h->nlev;
      }

      if ( fp->pdens > 1 ) {  /* Data is packed. */

	/* When data is packed, the offset refers to 8-byte cray words */

	dp = (char *)data + fp->off*8;
	out = slice;
	npack = 3 + (h->nlon-1) / fp->pdens;
	for ( k = 0; k < fldlevs; ++k ) {
#ifdef CRAY_UNPKA
          UNPKA( out, &(h->nlon), dp, &(fp->pdens) );
#elif defined CRAY_UNPKAF
          UNPKAF( out, &(h->nlon), dp, &(fp->pdens) );
#else
	  stat = unpack( out, h->nlon, dp, fp->pdens );
	  if ( stat == 1 ) ++oorval;
	  else if ( stat < 0 ) {
	    printf("get_slice: unpack returned %d\n",stat);
	    return -1;
	  }
#endif
	  dp  += npack * 8;
	  out += h->nlon;
	}

      } else if ( need_fcv( ff ) ) {  /* Foreign data */

	/* Offset in terms of input data floats */

	dp = (char *)data + fp->off*ff.fsz;
	if ( cv_floats( ff, nfloats, dp, slice ) != 0 ) ++oorval;

      } else {  /* no data conversions required */

	dp = (char *)data + fp->off*ff.fsz;
	in = (Hfloat *)dp;
	for ( i = 0; i < nfloats; ++i ) *slice++ = *in++;
      }

      if ( oorval > 0 ) return 1;
      return 0;
    }
  }    /* end loop over fields */

  /* field not found */

  return -3;
}                                                               /* get_slice */
/*===========================================================================*/


Hfloat cv_spval( char *xrep )
{
  /* Convert the special value 1.e36 used by the ccm processor to the local
     representation.

     Input arguments:
  
     xrep   external representation ( "cray", "ieee_sp", or "ieee_dp" )

     Return value:
  
      1.e36 in local float representation
  */

  static unsigned char ieee_sp[] = /* 1.e36 in single precision ieee */
    { (unsigned char)'\x7B', (unsigned char)'\x40', (unsigned char)'\x97', (unsigned char)'\xCE' };

  static unsigned char ieee_dp[] = /* 1.e36 in double precision ieee */
    { (unsigned char)'\x47', (unsigned char)'\x68', (unsigned char)'\x12', (unsigned char)'\xF9',
      (unsigned char)'\xCF', (unsigned char)'\x79', (unsigned char)'\x20', (unsigned char)'\xE3' };

  static unsigned char crayfp[] = /* 1.e36 in cray floating point */
    { (unsigned char)'\x40', (unsigned char)'\x78', (unsigned char)'\xC0', (unsigned char)'\x97',
      (unsigned char)'\xCE', (unsigned char)'\x7B', (unsigned char)'\xC9', (unsigned char)'\x07' };

  Hfloat val;
  float fval;
  double dval;

#ifdef CRAY_FLOAT
  if ( strcmp( xrep, "cray" ) == 0 ) {            /* No conversion */
    memcpy( &val, crayfp, 8 );
    return val;
  } else if ( strcmp( xrep, "ieee_sp" ) == 0 ) {  /* Convert SP IEEE to CRAY */
    if ( sptocf( ieee_sp, &val, 1 ) != 0 ) {
      printf("cv_spval: error return from sptocf\n");
      exit( -1 );
    } else return val;
  } else if ( strcmp( xrep, "ieee_dp" ) == 0 ) {  /* Convert DP IEEE to CRAY */
    if ( dptocf( ieee_dp, &val, 1 ) != 0 ) {
      printf("cv_spval: error return from dptocf\n");
      exit( -1 );
    } else return val;
  }
#elif defined ( IEEE_DP )
  if ( strcmp( xrep, "cray" ) == 0 ) {            /* Convert Cray to DP IEEE */
    if ( ctodpf( crayfp, &val, 1 ) != 0 ) {
      printf("cv_spval: error return from ctodpf\n");
      exit( -1 );
    } else return val;
  } else if ( strcmp( xrep, "ieee_sp" ) == 0 ) {/* Convert SP IEEE to DP IEEE */
    memcpy( &fval, ieee_sp, 4 );
    return (Hfloat)fval;
  } else if ( strcmp( xrep, "ieee_dp" ) == 0 ) {  /* No conversion */
    memcpy( &val, ieee_dp, 8 );
    return val;
  }
#elif defined ( IEEE_SP )
  if ( strcmp( xrep, "cray" ) == 0 ) {            /* Convert Cray to SP IEEE */
    if ( ctospf( crayfp, &val, 1 ) != 0 ) {
      printf("cv_spval: error return from ctospf\n");
      exit( -1 );
    } else return val;
  } else if ( strcmp( xrep, "ieee_sp" ) == 0 ) { /* No conversion */
    memcpy( &val, ieee_sp, 4 );
    return val;
  } else if ( strcmp( xrep, "ieee_dp" ) == 0 ) {/* Convert DP IEEE to SP IEEE */
    memcpy( &dval, ieee_dp, 8 );
    return (Hfloat)dval;
  }
#endif

  printf("cv_spval: unrecognized argument - %s\n",xrep);
  exit( -1 );
  return 0.;
}                                                                /* cv_spval */
/*===========================================================================*/

int spval_exists( FFILE ff, Hfloat* data, int nvals, Hfloat* spval )
{
  /* Check for the special value 1.e36 in the input array.

     Input arguments:
  
     ff     history tape input stream
     data   data that may contain a special value
     nvals  number of data values

     Output argument:

     spval  internal representation of the special value

     Return value:
  
     0 - special value not found
     1 - special value found
  */

  int i;

  /* Convert the special value 1.e36 from its external representation, i.e.,
     that used in the data file, to its internal one.  We are assuming that
     if the data contains special values, then it must be unpacked. */

  if ( ff.cosblk == 1 )
    *spval = cv_spval( "cray" );
  else if ( ff.fsz == 4 )
    *spval = cv_spval( "ieee_sp" );
  else if ( ff.fsz == 8 )
    *spval = cv_spval( "ieee_dp" );

  for ( i=0; i<nvals; ++i ) {
    if ( *(data+i) == *spval ) return 1;
  }

  return 0;
}                                                            /* spval_exists */
/*===========================================================================*/
