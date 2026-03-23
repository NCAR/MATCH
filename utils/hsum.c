/*$Id: hsum.c,v 1.1.1.1 1998/02/17 23:43:58 eaton Exp $*/

/*
NAME
     hsum - Print summary of a CCM1 or CCM2 file (history tape).

SYNOPSIS
     hsum [-h] [-l locpath] ccm_file

DESCRIPTION
     Print summary of information contained in history tape header records.
     Then print the mean absolute, minimum, and maximum values of each
     field for each time sample in the tape.

     The specified file can be either a local file, or a file on the
     MSS.  In the latter case, the full MSS pathname must be specified.
     If the file is an MSS file, that file will be copyed to
     locpath which may be specified with the -l option.  If locpath is
     not specified, then the default is to use the MSS pathname as the
     local one.  The directory part of the MSS pathname is created as
     a subdirectory of the current directory, then the MSS file is copied
     to that directory.

OPTIONS
     -h
          Print summary of header information only.  (The default is to
	  print a summary of header and field values.)

     -l locpath
          Pathname in local file system used to store the history tapes.
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ccmhdr.h"
#include "ccm_util.h"
#include "fortio.h"

extern void printusage( void );
extern void print_ccm_type( FFILE ff );

int main( int argc, char **argv )
{

  char *curarg;        /* current command line argument */
  FFILE ff;            /* stream for ccm file */
  CCM_HDR hdr;         /* generic ccm header */
  int i, n, status;
  int hdronly = 0;     /* Set to 1 to have only header info output */
  int set_locpath = 0; /* flag that user has set locpath */
  char *locpath;       /* local pathname for the ccm file */
  char *mspath;        /* MSS pathname for the ccm file */

  /* Process command line arguments */

  if ( argc == 1 ) {
    printusage();
    exit( -1 );
  }

  while ( --argc > 0 ) {

    curarg = *++argv;
    if ( *curarg == '-' ) {
      switch ( curarg[1] ) {
      case 'h':
	hdronly = 1;
	continue;
      case 'l':
        set_locpath = 1;
	locpath = *++argv;
	--argc;
	continue;
      default:
	fprintf(stderr,"Unrecognized argument: %s\n", curarg );
        printusage();
	exit( -1 );
      }
    }

    /* Arguments not starting with '-' are assumed to be file.
       If specified filename is not an MSS filename then assume it is a
       file on the local system.  */

    if ( curarg[0] == '/' && curarg[1] >= 'A' && curarg[1] <= 'Z' ) {
      mspath = curarg;
      if ( ! set_locpath ) locpath = mspath+1;
    } else {
      locpath = curarg;
    }

    if ( attach_ccm( locpath, mspath, 5000, &ff ) == -1 ) {
      printf("main: Can't open %s\n", curarg );
      exit( -1 );
    }

    print_ccm_type( ff );

    /* Initialize header structure. */

    if ( (status = init_hdr( ff, &hdr )) != 0 ) {
      printf("main: error %d returned from init_hdr\n", status );
      exit( -1 );
    }

    /* Print full header first time only */

    print_ccm_header( &hdr );

    /* Loop over time samples in file */

    do {

      /* Print time info for current data records */

      printf("\nDate = %6d, seconds = %5d, day = %6d, seconds = %5d\n",
	     hdr.ncdate, hdr.ncsec, hdr.ndcur, hdr.nscur );

      /* Print data stats */

      if ( hdronly == 0 ) {
	if ( (status=print_ccm_stats( ff, &hdr )) != 0 ) {
          printf("main: error %d returned from print_ccm_stats\n", status );
          exit(-1);
        }
      } else {
	for ( i = 1 ; i <= hdr.norec ; ++i ) skip_rec( ff );
      }

    } while ( next_hdr( ff, &hdr ) == 0 );

  }  /* End of command line argument processing */

  return 0;
}                                                          /* main */
/*===========================================================================*/

void printusage( void )
{
  fprintf(stderr,"\nUsage: hsum [options] ccm_file \n");
  fprintf(stderr,"\tOptions are: \n");
  fprintf(stderr,"\t\t-h\t\tPrint summary of header information only. \n");
  fprintf(stderr,"\t\t-l locpath\tLocal path to use for file on MSS. \n");
}
/*===========================================================================*/

void print_ccm_type( FFILE ff )
{
  /*
     Print characteristics of ccm file.
  */

  if ( ff.cosblk ) {
    printf("This CCM file has COS blocked structure.\n");
    printf("The real headers contain 8-byte cray floats.\n");
  } else {
    printf("This CCM file has f77/UNIX Fortran record blocking.\n");
    if ( ff.isz == 4 ) {
      if ( ff.fsz == 4 )
	printf("The real headers contain 4-byte ieee floats.\n");
      else
	printf("The real headers contain 8-byte ieee floats.\n");
    } else {
	printf("The real headers contain 8-byte cray floats.\n");
    }
  }
}
/*===========================================================================*/
