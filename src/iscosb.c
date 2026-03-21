/* $Id$ */
/* Try to determine if a Fortran binary sequential file is COS blocked. */
/* This version to be called from FORTRAN code */

#include <stdio.h>
#include "fca2cs.h"

typedef struct {
  unsigned int type  :  4;
  unsigned int junk  : 27;
  unsigned int bnhi  :  1;
  unsigned int bnlo  : 23;
  unsigned int fwi   :  9;
} BCW;
/*
  The BCW bitfield is used to extract the type of control word and the
  block number from a block control word.  In an 8 byte Cray control word
  the block number would be read as a single 24 bit piece of data.  In
  BCW the block number is split into two pieces, bnhi and bnlo, so that
  the block number field does not span a word boundary on a machine with
  4-byte words.

  See the Cray manpage BLOCKED(4F) for more information.
*/

#ifdef CRAY
int ISCOSB( _fcd ffpath )
#endif
#if ( defined AIX ) || ( defined HP )
int iscosb( const char *ffpath, int nc1 )
#endif
#if ( defined SUN ) || ( defined SGI ) || ( defined DEC ) || ( defined LINUX )
#endif
int iscosb_( const char *ffpath, int nc1 )
{
  /* Try to identify a COS blocked binary file by looking at the data where
     first two block control words would be and making sure that the 
     information contained in them is consistent, i.e., that their types
     are block control words, and that the block numbers are 0 and 1.

     Input arguement:

     fpath   file pathname

     Return values:

     1   file is COS blocked
     0   file NOT COS blocked
     -1  error encountered trying to open the file
     -2  error reading initial control word
     -3  error reading first block
     -4  error reading second control word

  */

  FILE *fp;
  BCW cw;
  char bytes[8][511];
  char* fpath;

  /* Convert Fortran character array to C string */
#ifdef CRAY
  fpath = fca2cs( ffpath );
#else
  fpath = fca2cs( ffpath, nc1 );
#endif

  /* Open the file. */

  if ( (fp=fopen( fpath, "r" )) == NULL ) {
    printf( "iscosb: can\'t open \"%s\"\n", fpath );
    return -1;
  }

  /* Read and interpret initial control word. */

  if ( fread( &cw, sizeof( BCW ), 1, fp ) != 1 ) return -2;

  if ( cw.type != 0 || cw.bnhi != 0 || cw.bnlo != 0 ) return 0;

  /* Skip over the rest of the first block. */

  if ( fread( bytes, 8, 511, fp ) != 511 ) return -3;

  /* Read and interpret second control word. */

  if ( fread( &cw, sizeof( BCW ), 1, fp ) != 1 ) return -4;

  fclose( fp );

  if ( cw.type != 0 || cw.bnhi != 0 || cw.bnlo != 1 ) return 0; 

  return 1;  /* file is COS blocked */
}
