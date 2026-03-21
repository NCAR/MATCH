/*$Id: strops.h,v 1.1.1.1 1998/02/17 23:44:02 eaton Exp $*/

/* String operations not in the standard C library. */

#ifndef STROPS_H
#define STROPS_H

extern char *strcdel( char *s, const char *ct, int c );
extern char *strshl( char *s, int n );
extern char *cmd_ln_sng(int argc,char **argv);

#endif
