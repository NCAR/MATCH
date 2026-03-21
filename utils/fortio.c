/*$Id: fortio.c,v 1.4 1999/01/21 18:11:28 eaton Exp $*/

/* Functions to handle i/o on binary files written using FORTRAN i/o. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* includes for the stat system call */
#include <sys/types.h>
#include <sys/stat.h>

#include "fortio.h"
#include "strops.h"
#include "iscosb.h"
#include "htdatlen.h"
#include "cr-ieee.h"
#include "version.h"

#if defined CRAY_FFIO       /* includes for Cray Flexible File i/o */
#include <fcntl.h>
#include <ffio.h>
#endif


/*
Function definitions:
*/


int attach_ccm( const char *locpath, const char *mspath, int bfsiz, FFILE *ff )
{
  /* 
     Open the CCM1 or CCM2 file locpath for reading.  If locpath is not on the 
     local file system, then it is acquired from mspath on the MSS.
     locpath is checked to see if it is a COS blocked file before it is
     opened.  Also determine if the real header record contains 4 or 8 byte
     floats.

     If locpath exists on the local file system, then mspath need not be
     specified.

     Input arguments:

     locpath   Pathname of the file on the local file system.
               If locpath is specified but does not exist, then the MSS
	       file mspath is copied to locpath.
     mspath    Pathname of the file on the MSS.  If locpath exists then
               mapath is not used.
     bfsiz     Size to use for the i/o buffer. (disabled)

     Output arguments:

     ff        Stream for opened file.

     Return values:

     0 for successful return, -1 otherwise.

    */

  char str[120];     /* string for system commands */
  char lpath[96];    /* local pathname */
  char ldir[80];     /* local directory */
  char lfn[80];      /* local filename */
  int ncshl;         /* number of characters to shift left */
  char *last_slash;  /* pointer to last slash in char string */
  int ncldir;        /* number of characters in local directory name */
  int st;            /* status from some function calls */
  struct stat stbuf;

#if defined CRAY_FFIO
  struct ffsw stffio;  /* status of flexible i/o request */
  char layer[16];      /* flexible file i/o layer specification */
#endif


  /* Set local directory and filename */

  strcpy( lpath, locpath );               /* Local pathname */

  last_slash = strrchr( lpath, '/' );

  if ( last_slash == NULL ) {
    strcpy( ldir, "./" );
    strcpy( lfn, lpath );
  } else {
    ncldir = ( last_slash - lpath ) + 1;
    strncpy( ldir, lpath, ncldir );
    strcpy( lfn, last_slash+1 );
  }

  /* Check for file on local system. */

  if ( (ff->strm.fp=fopen( lpath, "r" )) == NULL ) {

    /* File not found... make local directory if it doesn't exist. */
    if ( stat( ldir, &stbuf) != 0 ) {
      strcpy( str, "mkdir -m 777 -p " );
      strcat( str, ldir );
      fprintf(stderr,"attach_ccm: Issue system command -- %s\n", str );
      system( str );
    }

    /* Acquire file from MSS */
    strcpy( str, "msread " );
    strcat( str, lpath );
    strcat( str, " " );
    strcat( str, mspath );
    fprintf(stderr,"attach_ccm: Issue system command -- %s\n", str );
    if ( system( str ) != 0 ) return -1;

  } else {

    fprintf(stderr,"attach_ccm: using local file %s\n", lpath );
    fclose( ff->strm.fp );

  }

  /* Set COS blocking flag */

  if ( (ff->cosblk=iscosb( lpath )) < 0 ) {
    fprintf(stderr,"attach_ccm - iscosb: error return\n");
    return -1;
  }

  /* Set int and float data sizes */

  if ( ff->cosblk == 0 ) {
    if ( (st=htdatlen( lpath, &ff->isz, &ff->fsz )) < 0 ) {
      fprintf(stderr,"attach_ccm - htdatlen: status = %d\n", st );
      return -1;
    }
  } else { /* if file is COS blocked, assume cray internal representations */
    ff->isz = 8;
    ff->fsz = 8;
  }

  /* Open file */

#if defined CRAY_FFIO

  if ( ff->cosblk == 1 ) {
    strcpy( layer, "cos" );
  } else {
    strcpy( layer, "f77" );
  }

  if ( (ff->strm.fd = ffopens( lpath, O_RDONLY, 0, 0, &stffio, layer ))
                                                                    == -1 ) {
    fprintf(stderr,"attach_ccm - ffopens: error opening %s\n", lpath );
    return -1;
  }

#else

  if ( (ff->strm.fp=fopen( lpath, "rb" )) == NULL ) {
    fprintf(stderr,"attach_ccm - fopen: error opening %s\n", lpath );
    return -1;
  }

#endif

  return 0;
}                                                              /* attach_ccm */
/*===========================================================================*/

int close_ccm( FFILE ff )
{
  /* 
     Close the CCM1 or CCM2 file.

     Input arguments:

     ff        Stream for opened file.

     Return values:

     0 for successful return, -1 otherwise.

    */

#if defined CRAY_FFIO
  return ffclose( ff.strm.fd );
#else
  return fclose( ff.strm.fp );
#endif
}                                                               /* close_ccm */
/*===========================================================================*/


int attach_fort( const char *locpath, const char *mspath, int bfsiz,
                 FFILE *ff )
{
  /* 
     Open a binary Fortran file locpath for reading.  If locpath is not on the 
     local file system, then it is acquired from mspath on the MSS.
     locpath is checked to see if it is a COS blocked file before it is
     opened.

     The difference between this routine and attach_ccm is that for a
     general binary Fortran file it is not possible to distinquish between
     single and double precision ieee data.  So we make the decision that
     float data will be assumed to be the same precision as that being
     used to internally represent the data as set in the Hfloat typedef.

     If locpath exists on the local file system, then mspath need not be
     specified.

     Input arguments:

     locpath   Pathname of the file on the local file system.
               If locpath is specified but does not exist, then the MSS
	       file mspath is copied to locpath.
     mspath    Pathname of the file on the MSS.  If locpath exists then
               mapath is not used.
     bfsiz     Size to use for the i/o buffer. (used for Cray version only)

     Output arguments:

     ff        Stream for opened file.

     Return values:

     0 for successful return, -1 otherwise.

    */

  char str[120];     /* string for system commands */
  char lpath[96];    /* local pathname */
  char ldir[80];     /* local directory */
  char lfn[80];      /* local filename */
  int ncshl;         /* number of characters to shift left */
  char *last_slash;  /* pointer to last slash in char string */
  int ncldir;        /* number of characters in local directory name */
  char cbfsiz[8];    /* bfsiz as a character string */
  int st;            /* status from some function calls */

#if defined CRAY_FFIO
  struct ffsw stat;  /* status of flexible i/o request */
  char layer[16];   /* flexible file i/o layer specification */
#endif


  /* Set local directory and filename */

  strcpy( lpath, locpath );               /* Local pathname */

  last_slash = strrchr( lpath, '/' );

  if ( last_slash == NULL ) {
    strcpy( ldir, "./" );
    strcpy( lfn, lpath );
  } else {
    ncldir = ( last_slash - lpath ) + 1;
    strncpy( ldir, lpath, ncldir );
    strcpy( lfn, last_slash+1 );
  }

  /* Check for file on local system. */

  if ( (ff->strm.fp=fopen( lpath, "r" )) == NULL ) {

    /* Not found... Make sure local directory exits. */

    strcpy( str, "mkdir -m 777 -p " );
    strcat( str, ldir );
    printf("attach_fort: Issue system command -- %s\n", str );
    system( str ); /* Don't check return status because an error is returned */
                   /* if the directory already exits. */

    /* Acquire file from MSS */
    
    strcpy( str, "msread " );
    strcat( str, lpath );
    strcat( str, " " );
    strcat( str, mspath );
    printf("attach_fort: Issue system command -- %s\n", str );
    if ( system( str ) != 0 ) return -1;

  } else {

    printf("attach_fort: using local file %s\n", lpath );
    fclose( ff->strm.fp );

  }

  /* Set COS blocking flag */

  if ( (ff->cosblk=iscosb( lpath )) < 0 ) {
    printf("attach_fort - iscosb: error return");
    return -1;
  }

  /* Set int and float data sizes */

  if ( ff->cosblk == 0 ) {  /* ieee data */
    ff->isz = 4;
    ff->fsz = sizeof( Hfloat );
  } else { /* if file is COS blocked, assume cray internal representations */
    ff->isz = 8;
    ff->fsz = 8;
  }

  /* Open file */

#if defined CRAY_FFIO

  sprintf( cbfsiz, "%d", bfsiz );

  if ( ff->cosblk == 1 ) {
    strcpy( layer, "cos:" );
    strcat( layer, cbfsiz );
  } else {
    strcpy( layer, "f77" );
  }

  if ( (ff->strm.fd = ffopens( lpath, O_RDONLY, 0, 0, &stat, layer ))
                                                                    == -1 ) {
    printf("attach_fort - ffopens: error opening %s\n", lpath );
    return -1;
  }

#else

  if ( (ff->strm.fp=fopen( lpath, "rb" )) == NULL ) {
    printf("attach_fort - fopen: error opening %s\n", lpath );
    return -1;
  }

#endif

  return 0;
}                                                             /* attach_fort */
/*===========================================================================*/


int read_rec( FFILE ff, void *buf, int nbytes )
{
  /* Read specified number of bytes from the next Fortran record.
  
     Input args:
  
     ff      Input stream descriptor and flag for COS blocking.
     nbytes  Number of bytes to read.

     Output arguments:

     buf     Data read from record.
  
     Return value:
  
     number of bytes actually read, or a negative value indicating the
     error returned from the function called to read the data.
  */

#if defined CRAY_FFIO

  return ffread_rec( ff.strm.fd, buf, nbytes );

#else

  if ( ff.cosblk == 0 ) {
    return fread_f77rec( ff.strm.fp, buf, nbytes );
  } else {
    return fread_cosrec( ff.strm.fp, buf, nbytes );
  }

#endif

}                                                                /* read_rec */
/*===========================================================================*/


#if defined CRAY_FFIO
int ffread_rec( int fd, void *buf, int nbytes )
{
  /* Read specified number of bytes from the next Fortran record.
     Use ffread from CRAY's flexible file i/o library.  The file can be
     COS blocked, or can contain records in f77/UNIX format.  This
     distinction is made as a "layer" specification to the ffopens
     function.
  
     Input args:
  
     fd      Input stream.
     nbytes  Number of bytes to read.

     Output arguments:

     buf     Data read from record.
  
     Return value:
  
     number of bytes actually read, or -1 for error in ffread.  ffread will
     return 0 bytes read when it hits EOF rather than return a -1 error.
  */

  return ffread( fd, buf, nbytes );

}                                                              /* ffread_rec */
#endif
/*===========================================================================*/


int fread_f77rec( FILE *fp, void *buf, int nbytes )
{
  /* Read specified number of bytes from the next Fortran record.
     Use fread from the standard C library, and assume the records are
     in f77/UNIX non-VAX format, i.e., each record has a leading and
     trailing integer that contains the number of bytes in the
     record.
  
     Input args:
  
     fp      Input stream.
     nbytes  Number of bytes to read.

     Output arguments:

     buf     Data read from record.
  
     Return value:
  
     number of bytes actually read, or
     -1 for error in fread reading initial integer
     -2 for error in fread reading data
     -3 for error in fseed skipping unread part of record
     -4 for error in fread reading final integer
     -5 if initial and final integers don't match
  */

  int nbytes1, nbytes2;

  /* Read integer that tells how long the FORTRAN record is. */

  if ( fread( &nbytes1, sizeof( int ), 1, fp ) != 1 )  return -1;

  if ( nbytes > nbytes1 ) nbytes = nbytes1;

  /* Read "data" part of record. */

  if ( fread( buf, 1, (size_t)nbytes, fp ) != nbytes )  return -2;

  /* skip unread part of record */

  if ( nbytes1 > nbytes ) {
    if ( fseek( fp, nbytes1-nbytes, SEEK_CUR ) != 0 ) return -3;
  }

  /* check that last integer in record matches first */

  if ( fread( &nbytes2, sizeof( int ), 1, fp ) != 1 ) return -4;

  if ( nbytes1 != nbytes2 ) return -5;

  return nbytes;
}                                                            /* fread_f77rec */
/*===========================================================================*/


int fread_cosrec( FILE *fp, void *buf, int nbytes )
{
  /* Read specified number of bytes of data from the next Fortran record.
     Use fread from the standard C library, and assume the records are
     in COS blocked format.  On return the file will be positioned at the
     beginning of the next Fortran record.
  
     Input args:
  
     fp      Input stream.
     nbytes  Number of bytes of data to read (don't include number of bytes
             for control words).

     Output arguments:

     buf     Data read from record.
  
     Return value:
  
     number of bytes actually read, or
     -1 for error in fread reading initial control word
     -2 for error in fread reading data
     -3 for error in fseed skipping unread part of record
     -4 for error in fread reading control word
  */

  CWORD cw;        /* control word structure */
  int nbytestogo;  /* number of bytes that still need to be read and written */
  char *bufp;      /* keep track of where we are in the buffer being written */

  nbytestogo = nbytes;
  bufp = (char *)buf;

  /* Read initial control word. */

  if ( fread( &cw, sizeof( CWORD ), 1, fp ) != 1 ) return( -1 );

  do {

    if ( nbytestogo >= (int)cw.fwi*8 ) {

      /* read data to next control word */

      if ( fread( bufp, 1, cw.fwi*8, fp ) != cw.fwi*8 )
	return( -2 );

      nbytestogo -= cw.fwi*8;
      bufp += cw.fwi*8;
      
    } else {

      /* read amount of data requested then skip to next control word */

      if ( fread( bufp, 1, nbytestogo, fp ) != nbytestogo )
	return( -2 );

      if ( fseek( fp, cw.fwi*8-nbytestogo, SEEK_CUR ) != 0 ) return( -3 );
      nbytestogo = 0;
    }

    /* Read next control word. */

    if ( fread( &cw, sizeof( CWORD ), 1, fp ) != 1 ) return( -4 );

  }
  while ( cw.type == 0  &&  nbytestogo > 0 );

  /* We drop out of loop when we come to a control word that is not a 
     block control word (BCW) or nbytestogo == 0.
  */

  /* Reposition file in front of previous control word. */
  fseek( fp, -sizeof( CWORD ), SEEK_CUR );

  if ( cw.type == 0 ) {
    /* A block control word was encountered, so nbytestogo == 0.
       Skip to the next record control word. */
    fskip_cosrec( fp );
  } else {
    /* A record control word has been encountered.  Either we are done, or
       less than nbytes of data was found. */
    if ( nbytestogo > 0 )
      return nbytes-nbytestogo;
  }

  return nbytes;
}                                                            /* fread_cosrec */
/*===========================================================================*/


int skip_rec( FFILE ff )
{
  /* Skip the next Fortran record.
  
     Input args:
  
     ff      Input stream.
  
     Return value:
  
      0 for successful return or negative for error return.
  */

  int stat;

#if defined CRAY_FFIO

  if ( (stat=ffskip_rec( ff.strm.fd )) < 0 )
    fprintf(stderr,"skip_rec: ffskip_rec returned %d\n", stat );

#else

  if ( ff.cosblk == 0 ) {
    if ( (stat=fskip_f77rec( ff.strm.fp )) < 0 )
      fprintf(stderr,"skip_rec: fskip_f77rec returned %d\n", stat );
  } else {
    if ( (stat=fskip_cosrec( ff.strm.fp )) < 0 )
      fprintf(stderr,"skip_rec: fskip_cosrec returned %d\n", stat );
  }

#endif

  return stat;
}                                                                /* skip_rec */
/*===========================================================================*/

#if defined CRAY_FFIO
int ffskip_rec( int fd )
{
  /* Skip next FORTRAN record.
     Use ffread from CRAY's flexible file i/o library.  The file can be
     COS blocked, or can contain records in f77/UNIX format.  This
     distinction is made as a "layer" specification to the ffopens
     function.

     Input args:
  
     fd      Input stream.
  
     Return value:
  
     0  Successful return
     -1 for error in ffread
  */

  char buf[8];

  return ffread( fd, buf, 8 );

}                                                              /* ffskip_rec */
#endif
/*===========================================================================*/


int fskip_f77rec( FILE *fp )
{
  /* Skip the next Fortran record.
     Use fread and fseek from the standard C library, and assume the
     records are in f77/UNIX non-VAX format, i.e., each record has a
     leading and trailing 4-byte integer that contains the number of
     bytes in the record.
  
     Input args:
  
     fp      Input stream.
  
     Return value:
  
      0 for successful return
     -1 for error in fread reading initial integer
     -3 for error in fseed skipping data part of record
     -4 for error in fread reading final integer
     -5 if initial and final integers don't match
  */

  int nbytes1, nbytes2;

  /* Read integer that tells how long the FORTRAN record is. */

  if ( fread( &nbytes1, sizeof( int ), 1, fp ) != 1 )  return -1;

  /* Skip data part of record */

  if ( fseek( fp, nbytes1, SEEK_CUR ) != 0 ) return -3;

  /* check that last integer in record matches first */

  if ( fread( &nbytes2, sizeof( int ), 1, fp ) != 1 ) return -4;

  if ( nbytes1 != nbytes2 )  return -5;

  return 0;
}                                                            /* fskip_f77rec */
/*===========================================================================*/


int fskip_cosrec( FILE *fp )
{
  /* Reposition file in front of next record control word encountered.
     Use fread and fseek from the standard C library and assume the
     file is COS blocked.
     Note that on entry the file is not necessarily positioned at the
     end of a record.  It can be at any control word.

     Input args:
  
     fp      Input stream.
  
     Return value:
  
      0  Successful return
     -1  error in fread reading initial control word
     -2  error returned from fseek
     -3  error in fread reading control word
  */

  CWORD cw;

  /* Read and interpret initial control word. */

  if ( fread( &cw, sizeof( CWORD ), 1, fp ) != 1 )  return -1;

  do {

    /* skip to next control word */

    if ( fseek( fp, cw.fwi*8, SEEK_CUR ) != 0 )  return -2;

    /* read next control word */

    if ( fread( &cw, sizeof( CWORD ), 1, fp ) != 1 )  return -3;

  }
  while ( cw.type == 0 );

  /* A non-block control word has been found.  Position file in front of it. */

  fseek( fp, -sizeof( CWORD ), SEEK_CUR );

  return 0;
}                                                            /* fskip_cosrec */
/*===========================================================================*/


int prev_rec( FFILE ff )
{
  /* Reposition file to beginning of the previous FORTRAN record.

     Input args:
  
     ff      Input stream.
  
     Return value:
  
     0 for successful return or negative value for error return.
  */

  int stat;

#if defined CRAY_FFIO

  if ( ff.cosblk == 1 ) {
    if ( (stat=ffprev_rec( ff.strm.fd )) < 0 ) {
      fprintf(stderr,"prev_rec: ffprev_rec returned %d\n", stat );
    }
  } else {
    printf("prev_rec: backspacing f77/UNIX data not yet implemented\n");
    return -1;
  }

#else

  if ( ff.cosblk == 0 ) {
    if ( (stat=fprev_f77rec( ff.strm.fp )) < 0 ) {
      fprintf(stderr,"prev_rec: fprev_f77rec returned %d\n", stat );
    }
  } else {
    printf("prev_rec: backspacing COS blocked data not yet implemented\n");
    return -1;
  }

#endif

  return stat;
}                                                                /* prev_rec */
/*===========================================================================*/


#if defined CRAY_FFIO
int ffprev_rec( int fd )
{
  /* Reposition file to beginning of the previous FORTRAN record.
     Use ffbksp from CRAY's flexible file i/o library.  The file can be
     COS blocked, or can contain records in f77/UNIX format.  This
     distinction is made as a "layer" specification to the ffopens
     function.

     Input args:
  
     fd      Input stream.
  
     Return value:
  
     0  Successful return
     -1 for error in ffbksp
  */

  return ffbksp( fd );
}                                                              /* ffprev_rec */
#endif
/*===========================================================================*/


int fprev_f77rec( FILE *fp )
{
  /* Reposition file to beginning of the previous FORTRAN record.
     Use fread and fseek from the standard C library, and assume the
     records are in f77/UNIX non-VAX format, i.e., each record has a
     leading and trailing 4-byte integer that contains the number of
     bytes in the record.

     Input args:
  
     fp      Input stream.
  
     Return value:
  
     0  Successful return
     -1 error in fseek backing up past control word
     -2 error in fread reading control word
     -3 error in fseek backing up past record
  */

  int nbytes;

  /* Back up and read integer that gives length of previous Fortran record. */

  if ( fseek( fp, -sizeof( int ), SEEK_CUR ) != 0 ) return -1;
  if ( fread( &nbytes, sizeof( int ), 1, fp ) != 1 )  return -2;

  /* Back up past the data + two control words */

  if ( fseek( fp, -(nbytes+2*sizeof( int )), SEEK_CUR ) != 0 )
    return -3;

  return 0;
}                                                            /* fprev_f77rec */
/*===========================================================================*/


int rewind_file( FFILE ff )
{
  /* Reposition file to beginning.

     Input args:
  
     ff      Input stream.
  
     Return value:
  
     0  Successful return
     -1 for error
  */

  int stat;

#if defined CRAY_FFIO

  if ( (stat=ffseek( ff.strm.fd, 0, 0 )) < 0 ) {
    fprintf(stderr,"rewind_file: ffrewind returned %d\n", stat );
  }

#else

  if ( (stat=fseek( ff.strm.fp, 0, 0 )) < 0 ) {
    fprintf(stderr,"rewind_file: frewind returned %d\n", stat );
  }

#endif

  return stat;
}                                                             /* rewind_file */
/*===========================================================================*/


int rdcv_irec( FFILE ff, int ilen, void *buf, void *irec )
{
  /* Read and convert data in an integer record.  Partial records may be
     requested.

     Input arguments:
  
     ff      Input stream.
     ilen    Length of integer record (no. of int).
     buf     Storage used to read data from file.

     Output arguments:

     irec    Integer record.
  
     Return values:
  
      number of integers read and converted if successful.
     -1  Error returned by read_rec reading data.
     -2  No. of bytes read corresponds to a non-integral number of integers.
     -3  Error converting Cray to IEEE integers.
  */

  int ihbytes;  /* No. of bytes in integer record before conversion */
  int nbread;   /* No. of bytes read by read_rec */

  ihbytes = ilen*ff.isz;

  if ( (nbread=read_rec( ff, buf, ihbytes )) < 0 ) return -1;

  if ( nbread != ihbytes  &&  nbread%ff.isz != 0 ) return -2;

  ilen = nbread / ff.isz;

#ifdef CRAY_FLOAT
  sptoci( buf, irec, ilen );               /* convert ieee to cray */
#else
  if ( ctospi( buf, irec, ilen ) != 0 )    /* convert cray to ieee */
    return -3;
#endif

  return ilen;
}                                                               /* rdcv_irec */
/*===========================================================================*/


int rdcv_frec( FFILE ff, int flen, void *buf, void *frec )
{
  /* Read and convert data in an float record.  Partial records may be
     requested.

     Input arguments:
  
     ff      Input stream.
     flen    Length of float record (no. of floats).
     buf     Storage used to read data from file.

     Output arguments:

     frec    Float record.
  
     Return values:
  
     Number of floats read and converted if successful, or
     -1  Error returned by read_rec reading data.
     -2  No. of bytes read corresponds to a non-integral number of floats.
     -3  Error converting between Cray and IEEE floats.
  */

  int fhbytes;  /* No. of bytes requested in float record */
  int nbread;   /* No. of bytes read by read_rec */

  fhbytes = flen*ff.fsz;

  if ( (nbread=read_rec( ff, buf, fhbytes )) < 0 ) return -1;

  if ( nbread != fhbytes  &&  nbread%ff.fsz != 0 ) return -2;

  flen = nbread / ff.fsz;

  if ( cv_floats( ff, flen, buf, frec ) != 0 ) return -3;

  return flen;
}                                                               /* rdcv_frec */
/*===========================================================================*/


int cv_floats( FFILE ff, int flen, void *in, void *out )
{
  /* Convert input float data to machine representation of specified
     VERSION.

     Input arguments:
  
     ff      Input stream structure.
     flen    No. of floats to convert.
     in      Input float data.

     Output arguments:

     out     Converted floats.
  
     Return values:
  
      0 if successful, or
     -1  Error converting between Cray and IEEE floats.
  */

  int i;
  float *psp;   /* pointer into SP IEEE float stream */
  double *pdp;  /* pointer into DP IEEE float stream */

#ifdef CRAY_FLOAT
  if ( ff.fsz == 4 ) {
    /* Convert SP IEEE to CRAY */
    if ( sptocf( in, out, flen ) != 0 ) return -1;
  } else {
    /* Convert DP IEEE to CRAY */
    if ( dptocf( in, out, flen ) != 0 ) return -1;
  }
#elif defined(IEEE_DP)
  if ( ff.isz == 8 ) {
    /* Convert Cray to DP IEEE */
    if ( ctodpf( in, out, flen ) != 0 ) return -1;
  } else {
    /* Convert SP IEEE to DP IEEE (use implicit type conversion) */
    psp = (float *)in;
    pdp = (double *)out;
    for ( i = 0; i < flen; ++i ) *pdp++ = *psp++;
  }
#elif defined(IEEE_SP)
  if ( ff.isz == 8 ) {
    /* Convert Cray to SP IEEE */
    if ( ctospf( in, out, flen ) != 0 ) return -1;
  } else {
    /* Convert DP IEEE to SP IEEE (use implicit type conversion) */
    psp = (float *)out;
    pdp = (double *)in;
    for ( i = 0; i < flen; ++i ) *psp++ = *pdp++;
  }
#endif

  return 0;
}                                                               /* cv_floats */
/*===========================================================================*/


int need_icv( FFILE ff )
{
  /* Need integer conversion?

     Input argument:

     ff    Input stream.

     Return values:

     0 => no conversion necessary
     1 => conversion necessary
  */

#ifdef CRAY_FLOAT
  if ( ff.isz == 4 ) return 1;
#elif defined(INT64)
  if ( ff.isz == 4 ) return 1;
#else
  if ( ff.isz == 8 ) return 1;
#endif

  return 0;
}                                                                /* need_icv */
/*===========================================================================*/


int need_fcv( FFILE ff )
{
  /* Need float conversion?

     Input argument:

     ff    Input stream.

     Return values:

     0 => no conversion necessary
     1 => conversion necessary
  */

#ifdef CRAY_FLOAT
  if ( ff.isz == 4 ) return 1;
#else
  if ( ff.isz == 8 ) return 1;
  if ( ff.fsz != sizeof( Hfloat ) ) return 1;
#endif

  return 0;
}                                                                /* need_fcv */
/*===========================================================================*/

