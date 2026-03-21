/* $Id$ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fca2cs.h"

#ifdef CRAY
char* fca2cs( _fcd fca )
#else
char* fca2cs( const char* fca, int nc )
#endif
{
  /* Convert a Fortran character array to a C string */

  char *cs, *end;
  int fcalen;
#ifdef CRAY
  int nc;

  nc = _fcdlen( fca );     /* Lengths of FORTRAN character arrays */
#endif

  if ( nc < 1 ) return 0;

  /* Start by copying fca into local storage */

  if ( (cs=malloc( nc+1 )) == NULL ) {
    printf( "fca2cs: malloc failed\n" );
    return 0;
  }
#ifdef CRAY
  strncpy( cs, _fcdtocp( fca ), nc );
#else
  strncpy( cs, fca, nc );
#endif

  /* Start at end of input FORTRAN character array, and determine position
     of last significant (non-blank or non-null) character.
  */

  end = cs + nc - 1;
  for ( fcalen = nc; fcalen > 0; --fcalen, --end ) {
    if ( *end != ' ' && *end != '\0' ) break;
    if ( fcalen == 1 ) return 0;
  }

  /* Add null byte after last significant character */

  *(cs+fcalen) = '\0';

  return cs;
}
