/*$Id: fortio.h,v 1.1.1.1 1998/02/17 23:43:59 eaton Exp $*/

/* Interface for functions to handle i/o on binary files written using
   FORTRAN i/o. */

#ifndef FORTIO_H
#define FORTIO_H

#include <stdio.h>

typedef union STRM {
  FILE *fp;           /* FILE pointer for standard C i/o */
  int fd;             /* file descriptor for Flexible File i/o */
} STRM;

typedef struct FFILE {
  STRM strm;          /* i/o stream */
  int cosblk;         /* 1 => file is COS blocked, 0 => file NOT COS blocked */
  int isz;            /* length in bytes of integers in the file */
  int fsz;            /* length in bytes of floats in the file */
} FFILE;


typedef struct {
  unsigned int type  :  4;
  unsigned int junk1 : 28;
  unsigned int junk2 : 23;
  unsigned int fwi   :  9;
} CWORD;

/*
  The CWORD bitfield is used to extract information from the control words
  in COS blocked files.  It extracts the type of control word and the number of
  words to the next control word (fwi) from either a block control word
  or a record control word.  The "junk" in the middle is different for the
  two types of control word.  It is split so that the field does not span
  a word boundary on a machine with 4-byte words.

  type = 0           => block control word
         10 (octal)  => EOR
	 16 (octal)  => EOF
	 17 (octal)  => EOD

  See the Cray manpage BLOCKED(4F) for more information.
*/


/* "User Interface" level functions */

extern int attach_ccm( const char *locpath, const char *mspath, int bfsiz,
                       FFILE *ff );
extern int close_ccm( FFILE ff );
extern int read_rec( FFILE ff, void *buf, int nbytes );
extern int rdcv_irec( FFILE ff, int ilen, void *buf, void *irec );
extern int rdcv_frec( FFILE ff, int flen, void *buf, void *frec );
extern int cv_floats( FFILE ff, int flen, void *in, void *out );
extern int skip_rec( FFILE ff );
extern int prev_rec( FFILE ff );
extern int rewind_file( FFILE ff );
extern int need_icv( FFILE ff );
extern int need_fcv( FFILE ff );

/* Functions that use Cray's flexible file i/o */

extern int ffread_rec( int fd, void *buf, int nbytes );
extern int ffskip_rec( int fd );
extern int ffprev_rec( int fd );

/* Functions that use standard C library i/o on files with f77/UNIX
   record structure */

extern int fread_f77rec( FILE *fp, void *buf, int nbytes );
extern int fskip_f77rec( FILE *fp );
extern int fprev_f77rec( FILE *fp );


/* Functions that use standard C library i/o on COS blocked files */

extern int fread_cosrec( FILE *fp, void *buf, int nbytes );
extern int fskip_cosrec( FILE *fp );

#endif
