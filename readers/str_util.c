/* $Id$ */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "str_util.h"

#define MAX_TOKEN 40

/*
 ************** Utilities
 */
int STRblnk(char *s)
   /* Remove leading and trailing blanks from the string s              */
   /* return the length of the blanked out string */
   {
      int         first,
                  last,
                  i;
      int         lengthi;

      /* find first non blank */
      for (i=0; s[i] && isspace( (int) s[i]); i++);
      first = i;

      /* find last non blank */
      for (i= (int) strlen(s)-1; i>0 && isspace( (int) s[i]); i--);
      last = i;

      lengthi = last - first + 1;
      if (lengthi < 0)
         lengthi = 0;

      memmove( s, &s[first], lengthi);
      s[lengthi] = 0;

      return( lengthi);
   }

static char *Smax_copy(char *s1, const char *s2, int max_froms2, int maxs1)
   {
      int froms2;

      froms2 = MIN( max_froms2, maxs1-1);
      if (froms2 > 0)
	 {
         strncpy(s1, s2, (size_t) froms2);
         s1[froms2] = 0;
	 }

      return(s1);
   }

#include <ctype.h>
int Spos(const char *s, const char* c)
   {
      register int i, j, l;

      l = (int) strlen(c);

      for (i=0; s[i] != 0; i++)
         {
         for (j=0;j<l;j++)
            if ((s[i] == c[j]) || ((c[j] == ' ') && isspace( s[i])))
               {
               return i;
               }
         }

      return -1;
   }

char *STRremove( char *s, char c)
   {
      int   di = 0;
      int   si;

      for (si=0; s[si] != 0; si++)
         {
         if (s[si] != c)
            s[di++] = s[si];
         }
      s[di] = 0;

      return (s);
   }
   

char  *STRtokn( char **str_ptr, char *token, int max_toksiz, char *delim)
/* char  **str_ptr;     current point in the string to parse  */
/* char  *token;        put the found token here             */
/* int   max_toksiz;    maximum size of token (include zero byte) */
/* char  *delim;        delimiter for tokens                 */

/* Returns pointer to the original string where the token was found,
 *  and NULL if no token was found.
 *  str_ptr is updated past the token and delim(s), or set to NULL
 *  when theres nothing more to scan.
 */
{
    int   	delim_index;   /* index in string of the delimiter    */
    char 	*retptr;        /* return value */

    if (NULL == *str_ptr)
	return (NULL);

    /* skip leading delimiter characters */
    while (0 == (delim_index = Spos( *str_ptr, delim)))
    {
	(*str_ptr)++;
    }

    if (delim_index == -1 )
    {
	/* didnt find the delimiter */
	if ((*str_ptr != NULL) && (*str_ptr[0] != 0))
	{
	    /* copy nbytes or as much as caller has allocated for
	     * (which ever is less) to string, i.e. to EOS
	     */
	    Smax_copy( token, *str_ptr, 
			(int) strlen( *str_ptr), max_toksiz);
	    
	    /* scanning pointer set to null - we're done */
	    retptr = (*str_ptr);
	    *str_ptr = NULL;
	    return (retptr);
	}
	else
	{
	    *str_ptr = NULL;
	    return NULL;
	}
    }

    /* copy this token out */
    Smax_copy( token, *str_ptr, delim_index, max_toksiz);
    
    /* increment scanning pointer past the token and delimiter */
    retptr = (*str_ptr);
    (*str_ptr) += delim_index+1;

    return (retptr);
}
