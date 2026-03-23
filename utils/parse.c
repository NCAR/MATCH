/*$Id: parse.c,v 1.1.1.1 1998/02/17 23:43:58 eaton Exp $*/

/* some parsing functions */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parse.h"

char** str2str_array( char* s, const char* delim )
{
  /* Parse delimited string.  Return null terminated array of 
     substrings.

     Input arguments:
     s             delimited string
     delim         characters to use as delimiters

     Return value:
     Array of strings.  Null if input is null or delimiter string is null.
  */

  char* s_beg = s;          /* save pointers to beginning of input strings */
  const char* delim_beg = delim;

  char** str_list;
  char** str_list_beg;
  char* tok;
  int ndchar;
  int ndelim = 0;

  if ( s == NULL || delim == NULL ) {     /* check for null input */
    return (char**)0;
  }

  for (; *s != NULL; s++, delim=delim_beg ) {    /* count delimiters */
    while ( *delim ) {
      if ( *s == *delim++ ) { ndelim++; break; }
    }
  }

  /* allocate array of pointers */

  str_list_beg = (char**)malloc( (ndelim+2)*sizeof(char*) );

  /* create list from delimited substrings */

  str_list = str_list_beg;
  s = s_beg;
  for ( tok=strtok( s, delim ); tok != NULL; tok=strtok( (char*)0, delim ) ) {
    *str_list++ = tok;
  }
  *str_list = (char*)0;

  return str_list_beg;
}

/*===========================================================================*/
int test_str2str_array()
{
  char str[] = "a, ,,,,b#c, ";
  char** str_list;
  str_list = str2str_array( str,", #" );
  while ( *str_list ) printf("%s\n",*str_list++);
  return 0;
}
