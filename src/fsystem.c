/* $Id$ */
#include <stdlib.h>
#include <string.h>

#ifdef CRAY
#include <fortran.h>
int FSYSTEM( _fcd cmnd )
#endif
#if ( defined AIX ) || ( defined HP )
int fsystem( const char *cmnd, int nc1 )
#endif
#if ( defined SUN ) || ( defined SGI ) || ( defined DEC) || ( defined LINUX )
#endif
int fsystem_( const char *cmnd, int nc1 )
{
  /* FORTRAN interface for the standard C library function system.

     Input argument:

     cmnd   command to be passed to the operating system for execution

     FORTRAN parameters:  (Non-Cray versions)
     nc1 is the length of the FORTRAN character array that is
     passed as the first argument.  This data is provide by the 
     FORTRAN compiler, and should not appear as an actual argument when
     fsystem is called from the users FORTRAN code.

     On return:
     fsystem takes the return value of the system function
  */

  char *ccmnd, *end;
  int len;
#ifdef CRAY
  int nc1;

  nc1 = _fcdlen( cmnd );     /* Length of FORTRAN character array */
#endif

  /* Convert cmnd to a C string */
  /* Start by copying cmnd into local storage */

  ccmnd = malloc( nc1+1 );
#ifdef CRAY
  strncpy( ccmnd, _fcdtocp( cmnd ), nc1 );
#else
  strncpy( ccmnd, cmnd, nc1 );
#endif

  /* Start at end of input FORTRAN string, and determine position of last
     significant (non-blank or non-null) character.
  */

  end = ccmnd + nc1 - 1;
  for ( len = nc1; len > 0; --len, --end ) {
    if ( *end != ' ' && *end != '\0' ) break;
    if ( len == 1 ) {
      len = 0;
      break;
    }
  }

  /* Add null byte after last significant character */

  *(ccmnd+len) = '\0';

  /* Pass command to system */

  return system( ccmnd );
}
