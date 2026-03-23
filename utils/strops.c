/*$Id: strops.c,v 1.1.1.1 1998/02/17 23:44:02 eaton Exp $*/

/* String operations not in the standard C library. */

#include <stdlib.h>
#include <string.h>

#include "strops.h"


char *strcdel( char *s, const char *ct, int c )
{
  /* Remove all occurances of the character c from the string ct and place
     the resulting string in s.  s is also the return value.
  */

  char *ssave; /* save pointer to beginning of s */

  ssave = s;

  while ( *ct != '\0' ) {

    if ( *ct != c )
      *s++ = *ct++;
    else
      ++ct;
  }
  *s = '\0';

  return ssave;

}   /* strcdel */
/*===========================================================================*/


char *strshl( char *s, int n )
{
  /* Shift the string s, n characters to the left */

  char *ssave; /* save pointer to beginning of s */
  char *t;     /* part of string being shifted */

  if ( (int)strlen( s ) <= n ) {    /* return null string */
    *s = '\0';
    return s;
  }

  ssave = s;
  t = s + n;

  while ( *t != '\0' ) *s++ = *t++;

  *s = '\0';
  return ssave;

}   /* strshl */
/*===========================================================================*/


char *
cmd_ln_sng(int argc,char **argv)
/* 
   Concatenate the command line args into a " " delimited string.
   From Charlie Zender, 28 September 1995.

   int argc: input argument count
   char **argv: input argument list
   char *cmd_ln_sng(): output command line
*/ 
{
  char *cmd_ln;
  
  int cmd_ln_sz=0;
  int idx;
 
  for(idx=0;idx<argc;idx++){
    cmd_ln_sz+=(int)strlen(argv[idx])+1;
  } /* end loop over args */ 
  cmd_ln=(char *)malloc(cmd_ln_sz*sizeof(char));
  if(argc <= 0){
    cmd_ln[0]='\0';
  }else{
    (void)strcpy(cmd_ln,argv[0]);
    for(idx=1;idx<argc;idx++){
      (void)strcat(cmd_ln," ");
      (void)strcat(cmd_ln,argv[idx]);
    } /* end loop over args */ 
  } /* end else */ 
 
  return cmd_ln;
} /* end cmd_ln_sng() */ 
