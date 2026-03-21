/* $Id$ */
#ifndef FCA2CS_H
#define FCA2CS_H

#ifdef CRAY
#include <fortran.h>
extern char* fca2cs( _fcd name );
#else
extern char* fca2cs( const char* name, int nc );
#endif

#endif
