/*$Id: iscosb.c,v 1.1.1.1 1998/02/17 23:44:02 eaton Exp $*/

/* Try to determine if a Fortran binary sequential file is COS blocked. */

#include <stdio.h>
#include "iscosb.h"

int iscosb( const char *fpath )
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

  /* Open the file. */

  if ( (fp=fopen( fpath, "r" )) == NULL ) return -1;

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
