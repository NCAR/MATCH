/*$Id: cr-ieee.c,v 1.2 1998/03/30 19:59:42 eaton Exp $*/

#include "cr-ieee.h"

static const unsigned char signmask = (unsigned char)'\200';  /* Mask to get 1st of 8 bits */

/* Cray constants */

static const int cexpbias = 040000;      /* exponent bias */
static const unsigned char crbig[] =    /* largest positive float */
 { (unsigned char)'\137', (unsigned char)'\376', (unsigned char)'\377', (unsigned char)'\377',
   (unsigned char)'\377', (unsigned char)'\377', (unsigned char)'\377', (unsigned char)'\377' };
static const unsigned char crbigneg[] = /* largest negative float */
 { (unsigned char)'\337', (unsigned char)'\376', (unsigned char)'\377', (unsigned char)'\377',
   (unsigned char)'\377', (unsigned char)'\377', (unsigned char)'\377', (unsigned char)'\377' };

/* Single precision (SP) IEEE constants */

static const int spexpbias = 127;       /* exponent bias */
static const int spbexpmx = 255;       /* max value of biased exponent */
static const unsigned char spbig[] =    /* largest positive float */
 { (unsigned char)'\177', (unsigned char)'\177', (unsigned char)'\377', (unsigned char)'\377' };
static const unsigned char spbigneg[] = /* largest negative float */
 { (unsigned char)'\377', (unsigned char)'\177', (unsigned char)'\377', (unsigned char)'\377' };

/* Double precision (DP) IEEE constants */

static const int dpexpbias = 1023;       /* exponent bias */
static const int dpbexpmx = 2047;       /* max value of biased exponent */
static const unsigned char dpbig[] =    /* largest positive float */
 { (unsigned char)'\177', (unsigned char)'\357', (unsigned char)'\377', (unsigned char)'\377',
   (unsigned char)'\377', (unsigned char)'\377', (unsigned char)'\377', (unsigned char)'\377' };
static const unsigned char dpbigneg[] = /* largest negative float */
 { (unsigned char)'\377', (unsigned char)'\357', (unsigned char)'\377', (unsigned char)'\377',
   (unsigned char)'\377', (unsigned char)'\377', (unsigned char)'\377', (unsigned char)'\377' };


int
ctospi( void *in, void *out, int n )
{
  /* Convert 8-byte CRAY integers to 4-byte IEEE integers.

     Input arguments:

     in    Data stream of CRAY integers.
     n     Number of integers to convert.

     Output arguments:

     out  Data stream of IEEE integers.

     Return value:

      0 if successful,
     -1 if CRAY integers are found that are out-of-range for 4-byte
        IEEE representation.
 */

  int i, j;
  unsigned char sign;         /* sign bit */
  unsigned char *pin;         /* current position in input stream */
  unsigned char *pout;        /* current position in output stream */

  pin = (unsigned char *)in;
  pout = (unsigned char *)out;

  for ( i = 0; i < n; i++ ) {

    sign = *pin & signmask;	/*  sign bit  */

    /* Check for Cray integer that is out-of-range in 4-byte IEEE. */

    if ( sign == 0 ) {    /*  int >= 0; 1st 33 bits must be 0.  */

      if ( *pin != 0 || *(pin+1) != 0 || *(pin+2) != 0 || *(pin+3) != 0 ||
	   (*(pin+4) & signmask) != 0 ) return -1;

    } else {              /*  int < 0; 1st 33 bits must be 1.  */

      if ( *pin != 0377 || *(pin+1) != 0377 || *(pin+2) != 0377 ||
	   *(pin+3) != 0377 || (*(pin+4) & signmask) < 0200 ) return -1;
    }

    /* Skip 4 bytes in the Cray integer, then assign 4 low order bytes to
       output stream.  The range check just completed insures that the
       sign bit is correct.
    */

    pin += 4;
    for ( j = 0; j < 4; ++j ) *pout++ = *pin++;
  }
  return 0;
}
/*===========================================================================*/

void
sptoci( void *in, void *out, int n )
{
  /* Convert 4-byte IEEE integers to 8-byte CRAY integers.  

     Input arguments:

     in   Data stream of IEEE integers.
     n    Number of integers to convert.

     Output arguments:

     out  Data stream of Cray integers.
 */

  int i, j;
  unsigned char sign;    /* sign bit */
  unsigned char *pin;    /* current position in input stream */
  unsigned char *pout;   /* current position in output stream */

  pin = (unsigned char *)in;
  pout = (unsigned char *)out;

  for ( i = 0; i < n; i++ ) {

    sign = *pin & signmask;	/*  sign bit  */

    if ( sign != 0x00 ) sign = 0xff;

    for ( j = 0; j < 4; ++j ) *pout++ = sign;
    for ( j = 0; j < 4; ++j ) *pout++ = *pin++;
  }
}
/*===========================================================================*/

int
ctospf( void *in, void *out, int n )
{
  /* Convert 8-byte CRAY floats to 4-byte IEEE floats.

     Cray 64 bit float representation:
 	bits 0-47:	mantissa
 	bits 48-62:	exponent
 	bit 63:		sign

     IEEE 32 bit float representation:
 	bits 0-22:	mantissa
 	bits 23-30:	exponent
 	bit 31:		sign
 
     Input arguments:

     in    Data stream of CRAY floats.
     n     Number of floats to convert.

     Output arguments:

     out  Data stream of IEEE floats.

     Return value:

      0 if successful,
     -1 if CRAY floats are found that are out of range for 4-byte
        IEEE representation.  In this case, the overflow values are set to
	the largest possible IEEE reps, and the underflow values are set to
	zero.
 */

  int i, j;
  unsigned char sign;   /* sign bit for mantissa */
  int exp;              /* biased exponent */
  unsigned char spexp;  /* biased exponent for SP IEEE */
  unsigned char *pin;   /* current position in input stream */
  unsigned char *pout;  /* current position in output stream */
  int outofrange;       /* count out of range Cray floats */

  outofrange = 0;
  pin = (unsigned char *)in;
  pout = (unsigned char *)out;

  for ( i = 0; i < n; i++ ) {

    sign = *pin & signmask;

    exp = (*(pin+1)) | ((*pin & 0x7f) << 8);

    if ( exp == 0 ) {

      /* The zero value does not have a biased exponent, so a zero exponent
	 implies a zero value. */

      for ( j = 0; j < 4; ++j ) *pout++ = 0;

    } else {

      /* Convert biased Cray exponent to biased IEEE exponent.  Subtract 1
	 from the cray exponent because we must shift the cray mantissa 1
	 bit to the left (multiply it by 2) to obtain the IEEE normal form
	 (assuming that we are starting from a Cray normal form). */

      exp = exp - cexpbias -1 + spexpbias;

      if ( exp >= spbexpmx ) {

	/* Overflow.  Assign largest possible float with correct sign. */

	if ( sign == 0 )
	  for ( j = 0; j < 4; ++j ) *pout++ = spbig[j];
	else
	  for ( j = 0; j < 4; ++j ) *pout++ = spbigneg[j];

	++outofrange;

      } else if ( exp < 0 ) {

	/* Underflow.  Assign zero. */

	for ( j = 0; j < 4; ++j ) *pout++ = 0;

	++outofrange;

      } else {

	/* In-range or denormalized. */

	spexp = exp;
	*pout++ = sign | (spexp >> 1);
	*pout++ = ( (spexp & 1) << 7 ) | ( *(pin+2) & ~signmask );
	*pout++ = *(pin+3);
	*pout++ = *(pin+4);
      }
    }
    pin += 8;
  }
  if ( outofrange > 0 ) return -1;
  return 0;
}
/*===========================================================================*/

int
sptocf( void *in, void *out, int n )
{
  /* Convert SP IEEE (4-byte) floats to CRAY (8-byte) floats.  Ignore the
     representations of denormalized numbers in the IEEE format and just
     set these to zero.  (This is because these numbers should be converted
     to a normalized Cray representation, and the trouble of doing this
     doesn't seem worth the effort.)

     IEEE 32 bit float representation:
 	bits 0-22:	mantissa
 	bits 23-30:	exponent
 	bit 31:		sign
 
     Cray 64 bit float representation:
 	bits 0-47:	mantissa
 	bits 48-62:	exponent
 	bit 63:		sign

     Input arguments:

     in    Data stream of IEEE floats.
     n     Number of floats to convert.

     Output arguments:

     out  Data stream of CRAY floats.

     Return value:

      0 if successful,
     -1 if IEEE values are found that represent NaN's or +-infinity.  These
        values are set to the largest possible CRAY representation with
	the sign bit preserved.
 */

  int i, j;
  unsigned char sign;   /* sign bit for mantissa */
  unsigned char spexp;  /* biased exponent for SP IEEE */
  int exp;              /* cray biased exponent */
  unsigned char *pin;   /* current position in input stream */
  unsigned char *pout;  /* current position in output stream */
  int outofrange;       /* count IEEE NaN's and +-infinities */

  outofrange = 0;
  pin = (unsigned char *)in;
  pout = (unsigned char *)out;

  for ( i = 0; i < n; i++ ) {

    sign = *pin & signmask;

    spexp = *(pin+1)>>7  |  *pin<<1;

    if ( spexp == 0 ) {

      /* Set zero and denormalized IEEE numbers to zero (unsigned). */

      for ( j = 0; j < 8; ++j ) *pout++ = 0;
      pin += 4;
      continue;
    }

    if ( spexp == spbexpmx ) {

      /* Convert NaN or infinity to largest possible float with correct sign.*/

      if ( sign == 0 )
	for ( j = 0; j < 8; ++j ) *pout++ = crbig[j];
      else
	for ( j = 0; j < 8; ++j ) *pout++ = crbigneg[j];

      ++outofrange;
      pin += 4;
      continue;
    }

    /* Convert biased IEEE exponent to biased CRAY exponent.  The extra
       1 comes from the differenct between Cray and IEEE normal forms. */

    exp = spexp + cexpbias + 1 - spexpbias;

    /* Convert in-range number to Cray representation. */

    *pout++ = sign | exp>>8;
    *pout++ = exp;
    *pout++ = 0x80 | ( 0x7f & *(pin+1) );
    *pout++ = *(pin+2);
    *pout++ = *(pin+3);
    *pout++ = 0;
    *pout++ = 0;
    *pout++ = 0;

    pin += 4;
  }

  if ( outofrange > 0 ) return -1;
  return 0;
}
/*===========================================================================*/

int
ctodpf( void *in, void *out, int n )
{
  /* Convert 8-byte CRAY floats to 8-byte (DP) IEEE floats.

     Cray 64 bit float representation:
 	bits 0-47:	mantissa
 	bits 48-62:	exponent
 	bit 63:		sign

     IEEE 64 bit float representation:
 	bits 0-51:	mantissa
 	bits 52-62:	exponent
 	bit 63:		sign
 
     Input arguments:

     in    Data stream of CRAY floats.
     n     Number of floats to convert.

     Output arguments:

     out  Data stream of IEEE floats.

     Return value:

      0 if successful,
     -1 if CRAY floats are found that are out of range for 8-byte
        IEEE representation.  In this case, the overflow values are set to
	the largest possible IEEE reps, and the underflow values are set to
	zero.
 */

  int i, j;
  unsigned char sign;   /* sign bit for mantissa */
  int exp;              /* biased exponent */
  unsigned char *pin;   /* current position in input stream */
  unsigned char *pout;  /* current position in output stream */
  int outofrange;       /* count out of range Cray floats */

  outofrange = 0;
  pin = (unsigned char *)in;
  pout = (unsigned char *)out;

  for ( i = 0; i < n; i++ ) {

    sign = *pin & signmask;

    exp = (*(pin+1)) | ((*pin & 0x7f) << 8);

    if ( exp == 0 ) {

      /* The zero value does not have a biased exponent, so a zero exponent
	 implies a zero value. */

      for ( j = 0; j < 8; ++j ) *pout++ = 0;

    } else {

      /* Convert biased Cray exponent to biased IEEE exponent.  Subtract 1
	 from the cray exponent because we must shift the cray mantissa 1
	 bit to the left (multiply it by 2) to obtain the IEEE normal form
	 (assuming that we are starting from a Cray normal form). */

      exp = exp - cexpbias -1 + dpexpbias;

      if ( exp >= dpbexpmx ) {

	/* Overflow.  Assign largest possible float with correct sign. */

	if ( sign == 0 )
	  for ( j = 0; j < 8; ++j ) *pout++ = dpbig[j];
	else
	  for ( j = 0; j < 8; ++j ) *pout++ = dpbigneg[j];

	++outofrange;

      } else if ( exp < 0 ) {

	/* Underflow.  Assign zero. */

	for ( j = 0; j < 8; ++j ) *pout++ = 0;

	++outofrange;

      } else {

	/* In-range or denormalized. */

	*pout++ = sign        | exp>>4;
	*pout++ = exp<<4      | ( *(pin+2) & 0x7f )>>3;
	*pout++ = *(pin+2)<<5 | *(pin+3)>>3;
	*pout++ = *(pin+3)<<5 | *(pin+4)>>3;
	*pout++ = *(pin+4)<<5 | *(pin+5)>>3;
	*pout++ = *(pin+5)<<5 | *(pin+6)>>3;
	*pout++ = *(pin+6)<<5 | *(pin+7)>>3;
	*pout++ = *(pin+7)<<5;
      }
    }
    pin += 8;
  }
  if ( outofrange > 0 ) return -1;
  return 0;
}
/*===========================================================================*/

int
dptocf( void *in, void *out, int n )
{
  /* Convert DP IEEE (8-byte) floats to CRAY (8-byte) floats.  Ignore the
     representations of denormalized numbers in the IEEE format and just
     set these to zero.  (This is because these numbers should be converted
     to a normalized Cray representation, and the trouble of doing this
     doesn't seem worth the effort.)

     IEEE 64 bit float representation:
 	bits 0-51:	mantissa
 	bits 52-62:	exponent
 	bit 63:		sign
 
     Cray 64 bit float representation:
 	bits 0-47:	mantissa
 	bits 48-62:	exponent
 	bit 63:		sign

     Input arguments:

     in    Data stream of IEEE floats.
     n     Number of floats to convert.

     Output arguments:

     out  Data stream of CRAY floats.

     Return value:

      0 if successful,
     -1 if IEEE values are found that represent NaN's or +-infinity.  These
        values are set to the largest possible CRAY representation with
	the sign bit preserved.
 */

  int i, j;
  unsigned char sign;   /* sign bit for mantissa */
  int dpexp;            /* biased exponent for DP IEEE */
  int crexp;            /* cray biased exponent */
  unsigned char *pin;   /* current position in input stream */
  unsigned char *pout;  /* current position in output stream */
  int outofrange;       /* count IEEE NaN's and +-infinities */

  outofrange = 0;
  pin = (unsigned char *)in;
  pout = (unsigned char *)out;

  for ( i = 0; i < n; i++ ) {

    sign = *pin & signmask;

    dpexp = (*pin & 0x7f)<<4 | *(pin+1)>>4;

    if ( dpexp == 0 ) {

      /* Set zero and denormalized IEEE numbers to zero (unsigned). */

      for ( j = 0; j < 8; ++j ) *pout++ = 0;
      pin += 8;
      continue;
    }

    if ( dpexp == dpbexpmx ) {

      /* Convert NaN or infinity to largest possible float with correct sign.*/

      if ( sign == 0 )
	for ( j = 0; j < 8; ++j ) *pout++ = crbig[j];
      else
	for ( j = 0; j < 8; ++j ) *pout++ = crbigneg[j];

      ++outofrange;
      pin += 8;
      continue;
    }

    /* Convert biased IEEE exponent to biased CRAY exponent.  The extra
       1 comes from the differenct between Cray and IEEE normal forms. */

    crexp = dpexp + cexpbias + 1 - dpexpbias;

    /* Convert in-range number to Cray representation. */

    *pout++ = sign | crexp>>8;
    *pout++ = crexp;
    *pout++ = 0x80 | ( 0xf & *(pin+1) )<<3 | *(pin+2)>>5;
    *pout++ = *(pin+2)<<3 | *(pin+3)>>5;
    *pout++ = *(pin+3)<<3 | *(pin+4)>>5;
    *pout++ = *(pin+4)<<3 | *(pin+5)>>5;
    *pout++ = *(pin+5)<<3 | *(pin+6)>>5;
    *pout++ = *(pin+6)<<3 | *(pin+7)>>5;

    pin += 8;
  }

  if ( outofrange > 0 ) return -1;
  return 0;
}
/*===========================================================================*/
