/* $Id$ */
#include <stdlib.h>
#include <string.h>

#ifdef CRAY
#include <fortran.h>
int FGETENV( _fcd name, _fcd value )
#endif
#if ( defined AIX ) || ( defined HP )
int fgetenv( const char *name, char *value, int nc1, int nc2 )
#endif
#if ( defined SUN ) || ( defined SGI ) || ( defined DEC ) || ( defined LINUX )
#endif
int fgetenv_( const char *name, char *value, int nc1, int nc2 )
{
  /* FORTRAN interface for the standard C library function getenv.

     Input argument:

     name   Name of environment variable.

     Output argument:

     value  Environment string associated with name (if found).

     FORTRAN parameters:  (Non-Cray versions)
     nc1 and nc2 are the lengths of the FORTRAN character arrays that are
     passed as the first two arguments.  This data is provide by the 
     FORTRAN compiler, and should not appear as actual arguments when
     fgetenv is called from the users FORTRAN code.

     On return:
      0 => environment variable string with 0 < length <= nc2 was found
     -1 => input string of length zero or no storage allocated for return value
     -2 => error return from malloc
     -3 => environment variable not found
     -4 => environment variable string of zero length was found
     -5 => environment variable string of length > nc2 was found
  */

  char *cname, *end, *cvalue;
  int namlen, vlen;
#ifdef CRAY
  int nc1, nc2;

  nc1 = _fcdlen( name );     /* Lengths of FORTRAN character arrays */
  nc2 = _fcdlen( value );
#endif

  if ( nc1 < 1 ) return -1;
  if ( nc2 < 1 ) return -1;

  /* Convert name to a C string */
  /* Start by copying name into local storage */

  if ( (cname=malloc( nc1+1 )) == NULL ) return -2;
#ifdef CRAY
  strncpy( cname, _fcdtocp( name ), nc1 );
#else
  strncpy( cname, name, nc1 );
#endif

  /* Start at end of input FORTRAN string, and determine position of last
     significant (non-blank or non-null) character.
  */

  end = cname + nc1 - 1;
  for ( namlen = nc1; namlen > 0; --namlen, --end ) {
    if ( *end != ' ' && *end != '\0' ) break;
    if ( namlen == 1 ) return -1;
  }

  /* Add null byte after last significant character */

  *(cname+namlen) = '\0';

  /* Get the value of the environment variable */

  if ( (cvalue=getenv( cname )) == NULL ) return -3;

  vlen = strlen( cvalue );
  
  if ( vlen == 0 ) return -4;

#ifdef CRAY
  strncpy( _fcdtocp( value ), cvalue, nc2 );
#else
  strncpy( value, cvalue, nc2 );
#endif

  if ( vlen > nc2 )
    return -5;
  else
    return 0;
}
