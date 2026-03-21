/* $Id$ */
#ifndef STR_UTIL_WAS_INCLUDED
#define STR_UTIL_WAS_INCLUDED

/* #include <malloc.h> */

#ifndef MAX
#define MAX(a,b)        (((a) > (b)) ?  (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b)        (((a) < (b)) ?  (a) : (b))
#endif
#define RAD_PER_DEG 0.017453293 /* radians per degree */

#define MALLOC(s) (s *) calloc(1,sizeof(s))
#define CALLOC(n,s) (s *) calloc(n,sizeof(s))
#define REALLOC(p,n,s) (s *) realloc(p,(n)*sizeof(s))
#define FREE(p) if(p)free(p)

extern int STRblnk(char *s);
extern char *STRremove( char *s, char c);
extern char  *STRtokn( char **str_ptr, char *token, int max_toksiz, char *delim);

#endif
