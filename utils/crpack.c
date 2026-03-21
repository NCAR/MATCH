/*$Id: crpack.c,v 1.2 1998/04/07 22:32:49 eaton Exp $*/

/* C version of the CRAY unpack routines */

#include <stdio.h>
#include <string.h>

#include "version.h"
#include "crpack.h"
#include "cr-ieee.h"


int unpack( Hfloat *arn, int nrn, void *apn, int npk )
{
  /* Unpack data that was packed using the CRAY packaf, packbf, or packcf
     subroutines.

     Input arguments:

     nrn  number of unpacked floats requested
     apn  packed data
     npk  packing density

     Output arguments:

     arn  unpacked data

     Return values:

      1  out-of-range values encountered by CRAY-IEEE conversion routine
         (underflow set to zero and overflow set to largest ieee values)
      0  success
     -2  illegal packing density specified
  */

  Hfloat amn;      /* min value of real array */
  Hfloat scale;    /* scale factor to convert packed integers to floats */
  float *crayfp;   /* used to interpret values in apn as CRAY floats */
  Hfloat ieee[2];  /* store values in apn as IEEE floats */
  char *strm;      /* interpret apn as byte stream pointer */
  unsigned char (*p4)[4]; /* interpret 4-bytes as an int (npk=2) */
  unsigned char (*p2)[2]; /* interpret 2-bytes as an int (npk=4) */
  int i;
  int noorv = 0;  /* number of out-of-range cray values found by conversion */
                  /* routines.                                              */

#ifndef CRAY_FLOAT
  union {
    unsigned int pkint;  /* packed integer */
    unsigned char c4[4]; /* 4 bytes */
  } u;
#else
  union {
    unsigned int pkint;  /* packed integer */
    unsigned char c8[8]; /* 8 bytes */
  } u;
#endif

  /* Make sure size of int matches number of chars specified in union */

#ifndef CRAY_FLOAT
  if ( sizeof( int ) != 4 ) {
    printf("unpack: sizeof( int ) must be 4 bytes\n");
    exit( -1 );
  }
#endif

  u.pkint = 0;  /* Zero all bits in the union */

  /* The 1st 2 values in the packed data are CRAY floats.  Do conversions
     if necessary */

#ifdef CRAY_FLOAT
  crayfp = (Hfloat *)apn;
  amn = *crayfp++;
  scale = *crayfp;
#elif defined ( IEEE_DP )
  if ( ctodpf( apn, ieee, 2 ) != 0 ) ++noorv;
  amn = ieee[0];
  scale = ieee[1];
#elif defined ( IEEE_SP )
  if ( ctospf( apn, ieee, 2 ) != 0 ) ++noorv;
  amn = ieee[0];
  scale = ieee[1];
#endif

  strm = (char *)apn + 16;  /* advance pointer to packed integers */

  if ( npk == 2 ) {

    p4 = (unsigned char (*)[4])strm;   /* init pointer to 4-byte ints */
    for ( i = 0; i < nrn; ++i ) {

      /* copy 4 bytes into union */

#ifndef CRAY_FLOAT
      memcpy( u.c4, p4++, 4 );
#else
      memcpy( &(u.c8[4]), p4++, 4 );
#endif
      arn[i] = amn + u.pkint/scale;
    }

  } else if ( npk == 4 ) {

    p2 = (unsigned char (*)[2])strm;   /* init pointer to 2-byte ints */
    for ( i = 0; i < nrn; ++i ) {

      /* copy 2 bytes into union */

#ifndef CRAY_FLOAT
      memcpy( &(u.c4[2]), p2++, 2 );
#else
      memcpy( &(u.c8[6]), p2++, 2 );
#endif
      arn[i] = amn + u.pkint/scale;

    }

  } else
    return -2;

  if ( noorv == 0 ) return 0;
  return 1;
}                                                                  /* unpack */
/*===========================================================================*/
