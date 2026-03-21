/*$Id: htdatlen.c,v 1.3 1998/04/07 22:32:51 eaton Exp $*/

#include <stdio.h>

#include "version.h"
#include "htdatlen.h"
#include "cr-ieee.h"

#define NINTS 32       /* no. of ints read from header to get size info */

int htdatlen( const char *fpath, int *isz, int *fsz )
{
/* Determine the lengths of integers and floats in a binary CCM history file.
   Assume the file contains f77/unix type records.
   If the integer header record contains 4-byte integers the data representation
   is assumed to be ieee.  If the integer header record contains 8-byte integers
   then the data format is assumed to be Cray.

   Currently we assume the integers are 4 bytes.  If this is not the case
   an error is flagged.

   Input argument:

   fpath   file pathname of a CCM history tape.

   Output arguments:

   isz     Size in bytes of integers in the history tape.
   fsz     Size in bytes of floats in the history tape.

   Return values:

    0  successful return
   -1  error encountered trying to open the file
   -2  error reading initial control word
   -3  error reading first 32 integers from header
   -4  inconsistency between length of integer header and info in RCW
       i.e., it looks like the integers are neither 4 nor 8 bytes.
   -5  error skipping unread part of integer record
   -6  error reading RCW at end of integer record
   -7  initial and final RCWs don't match for integer record
   -8  error skipping character record
   -9  error reading RCW at beginning of float record

  */

  FILE *fp;
  const int ncw=1;          /* no. of control words */
  int nbrec;                /* number of bytes in record (initial RCW) */
  int nbrecf;               /* number of bytes in record (final RCW) */
  int in[NINTS*2];          /* storage for NINTS integers (possibly 8-byte */
  int out[NINTS*2];         /* storage for integer conversion if necessary*/
  int *ih;                  /* part of integer header */
  int lenhdi;               /* number of ints in header */
  int lenhdc;               /* number of char[8]s in header */
  int lenhdr;               /* number of floats in header */

  /* Open the file. */
  if ( (fp=fopen( fpath, "r" )) == NULL ) return -1;

  /* Read 4-byte record control word */
  if ( fread( in, 4, 1, fp ) != 1 ) return -2;

  /* Convert ieee integer to Cray if necessary. */
#ifdef CRAY_FLOAT
  sptoci( in, out, 1 );
  nbrec = *out;
#else
  nbrec = *in;
#endif

  /* Read first NINTS integers - assuming they are 4 bytes */
  if ( fread( in, 4, NINTS, fp ) != NINTS ) return -3;
#ifdef CRAY_FLOAT
  sptoci( in, out, NINTS );
  ih = out;
#else
  ih = in;
#endif

  lenhdi = ih[0];
  lenhdc = ih[30];
  lenhdr = ih[31];

  /* Check to see if integer header length from header record agrees with
     record control word size.  If it does then the data is ieee.  If not,
     try reading the data as 8-byte integers.  If the check works then we
     assume the data is cray binary.
  */
  if ( ih[0]*4 != nbrec ) {
    /* Position file after first 4-byte RCW */
    if ( fseek( fp, 4L, SEEK_SET ) != 0 ) return -5;
    /* Read first NINTS integers - assuming they are 8 bytes */
    if ( fread( in, 8L, NINTS, fp ) != NINTS ) return -3;
#ifdef CRAY_FLOAT
    ih = in;
#else
    ctospi( in, out, NINTS );
    ih = out;
#endif
    /* Does integer header length correspond to record length? */
    if ( ih[0]*8 == nbrec ) {
      /* assume if integers are 8-bytes, that records contain cray binary */
      *isz = 8;
      *fsz = 8;
      fclose( fp );
      return 0;
    } else {
      return -4;
    }
  }

  *isz = 4;

  /* skip through rest of integer header */
  if ( fseek( fp, 4*(lenhdi-NINTS), SEEK_CUR ) != 0 ) return -5;

  /* Read RCW at end of integer record */
  if ( fread( in, 4, 1, fp ) != 1 ) return -6;
#ifdef CRAY_FLOAT
  sptoci( in, out, 1 );
  nbrecf = *out;
#else
  nbrecf = *in;
#endif

  /* make sure initial and final RCWs match */
  if ( nbrec != nbrecf ) return -7;

  /* skip over character record + 2 RCWs */
  if ( fseek( fp, 8*(lenhdc+1), SEEK_CUR ) != 0 ) return -8;

  /* Read RCW at begining of float record */
  if ( fread( in, 4, 1, fp ) != 1 ) return -9;
#ifdef CRAY_FLOAT
  sptoci( in, out, 1 );
  nbrec = *out;
#else
  nbrec = *in;
#endif

  *fsz = nbrec / lenhdr;

  /* return error if float size not 4 or 8 bytes */
  if ( *fsz != 4  &&  *fsz != 8 ) return -10;

  fclose( fp );
  return 0;
}            /* htdatlen */
