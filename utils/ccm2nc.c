/*$Id: ccm2nc.c,v 1.4 1998/08/19 19:57:01 eaton Exp $*/

/*
  ccm2nc - convert a CCM history format file to netCDF format

  Author - Brian Eaton
           National Center for Atmospheric Research
	   Climate and Global Dynamics Division
	   Climate Modeling Section
	   P.O. Box 3000
	   Boulder, CO 80307-3000
	   eaton@ucar.edu
*/

#include <netcdf.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* includes for the stat system call */
#include <sys/types.h>
#include <sys/stat.h>

#include "fortio.h"
#include "ccmhdr.h"
#include "nc_util.h"
#include "parse.h"
#include "strops.h"

void printusage( void );

int main( int argc, char **argv )
{

  char *curarg;          /* current command line argument */
  char* cmndline;        /* command line as single string */
  char** reqvar = NULL;  /* requested variable names */
  char *ccm_fname;       /* file name for CCM input file */
  char *nc_fname;        /* file name for netCDF output file */
  FFILE ff;              /* stream for history tape */
  CCM_HDR hdr;           /* CCM header info */
  char *msspath = NULL;  /* NCAR MSS pathname for the history tape file */
  char *locpath = NULL;  /* local pathname for the history tape file */
  int ncid;              /* i.d. for netCDF file */
  var_id vid;            /* variable ids */
  int its = 0;           /* time sample index */
  int status;            /* function return status */
  int first_ccm = 1;     /* flag for first CCM file processed */
  int coards = 0;        /* 0 => CSM conventions, 1 => COARDS conventions */
  int dim_ord = 0;       /* 0 => lon,lat,lev,tim; 1 => lon,lev,lat,tim */
  int doubleOut = 0;     /* 0 => variables output as NC_FLOAT, 1 => NC_DOUBLE */
  struct stat stbuf;
  int nFileArgs = argc-1;  /* number of arguments that are files */
  /*--------------------------------------------------------------------------*/

  /* ccm2nc assumes that floating point data is 8 bytes, i.e., float on crays and
     double on ieee machines.  Check that this is true. */
  if ( sizeof( Hfloat ) != 8 ) {
    fprintf( stderr, "main: ccm2nc must be compiled so that the Hfloat type is 8 bytes\n" );
    exit( EXIT_FAILURE );
  }

  if ( argc < 3 ) {
    printusage();
    exit( -1 );
  }

  cmndline = cmd_ln_sng( argc, argv );

  while ( --argc > 1 ) {  /* don't process last arg which we assume is */
                          /* the output netCDF filename */
    curarg = *++argv;
    if ( *curarg == '-' ) {
      --nFileArgs;
      switch ( curarg[1] ) {
      case 'C':
	coards = 1;
	continue;
      case 'v':
	--argc;
	reqvar = str2str_array( *++argv, "," );
	if ( reqvar == NULL ) {
	  fprintf( stderr, "main: error parsing argument: %s\n", curarg );
	  printusage();
	  exit( -1 );
	}
	continue;
      case 'x':
	if ( curarg[2] == 'z' ) {
	  dim_ord = 1;
	} else {
	  fprintf( stderr, "main: error parsing argument: %s\n", curarg );
	  printusage();
	  exit( -1 );
	}
	continue;
      case 'r':
	if ( curarg[2] == '8' ) {
	  doubleOut = 1;
	} else {
	  fprintf( stderr, "main: error parsing argument: %s\n", curarg );
	  printusage();
	  exit( -1 );
	}
	continue;
      default:
	fprintf( stderr, "main: error parsing argument: %s\n", curarg );
        printusage();
	exit( -1 );
      }
    }

    /* If COARDS conventions requested then override -xz flag */
    if ( coards == 1 ) dim_ord = 0;

    /* assume the remaining arguments are the required filenames */
    ccm_fname = curarg;

    /* If specified ccm filename is an MSS filename then will attempt to read
       it from the mss if it is not found on the local disk with a relative
       pathname that is the same as the MSS pathname without the initial slash.
    */
    if ( ccm_fname[0] == '/' && ccm_fname[1] >= 'A' && ccm_fname[1] <= 'Z' ) {
      msspath = ccm_fname;
      locpath = ccm_fname+1;
    } else {
      msspath = NULL;
      locpath = ccm_fname;
    }
    if ( attach_ccm( locpath, msspath, 0, &ff ) == -1 ) {
      fprintf( stderr, "main: Can't open %s\n", locpath );
      exit( -1 );
    }

    if ( first_ccm ) {  /* initialization done for first CCM file only */
      first_ccm = 0;

      /* Initialize CCM header structure. */
      if ( (status = init_hdr( ff, &hdr )) != 0 ) {
	fprintf( stderr, "main: init_hdr returned %d\n", status );
	exit( -1 );
      }

      /* Create netCDF file.  If the specified output file exists, and it
	 is a netCDF file, then clobber it with the new output file.  If
	 the output file exists, but it's not a netCDF file, assume the
	 user forgot to specify the output file and exit.
      */
      nc_fname = argv[argc-1];
      if ( stat( nc_fname, &stbuf) == 0 ) {  /* file exists */
	if ( (status = nc_open( nc_fname, NC_NOWRITE, &ncid )) != NC_NOERR ) { /* it's not a netCDF file */
	  fprintf( stderr, "main: output file %s exists\n", nc_fname );
	  exit( -1 );
	}
      }
      if ( (status = nc_create( nc_fname, NC_CLOBBER, &ncid )) != NC_NOERR ) {
	fprintf( stderr, "%s\n", nc_strerror(status) );
	fprintf( stderr, "main: nc_create failed on file %s\n", nc_fname );
	exit( -1 );
      }

      /* Initialize netCDF file */
      init_nc( &hdr, ncid, reqvar, cmndline, &vid, dim_ord, coards, doubleOut );

    } else {  /* Read the first 3 header records of new CCM file */

      if ( (status=next_hdr( ff, &hdr )) != 0 ) {
	fprintf( stderr, "main: next_hdr returned %d\n", status );
	exit( -1 );
      }
      ++its;
    }

    do {  /* Loop over time samples in file */

      /* Put record data into file */
      write_ncrec( ff, &hdr, ncid, &vid, its, dim_ord, doubleOut );

      /* Read next CCM header records if they exist */
      if ( next_hdr( ff, &hdr ) != 0 ) break;

      ++its;

    } while ( 1 );

    if ( close_ccm( ff ) == -1 ) {
      fprintf( stderr, "main: Can't close %s\n", locpath );
      exit( -1 );
    }
  }     /* process next CCM file */

  if ( nFileArgs < 2 ) {
    printusage();
    exit( -1 );
  }

  /* Close netCDF file */
  handleError( nc_close( ncid ) );

  return 0;
}
/*===========================================================================*/

void printusage( void )
{
  fprintf(stderr,"\nUsage: ccm2nc [options] ccm_file[s] netcdf_file \n");
  fprintf(stderr,"Options: \n");
  fprintf(stderr,"   -C \n\toutput netCDF file conforms to COARDS conventions \n");
  fprintf(stderr,"   -r8 \n\toutput variables as NC_DOUBLE \n");
  fprintf(stderr,"   -v var1[,var2...]\n\tvariable names (case sensitive)\n");
  fprintf(stderr,"   -xz \n\tuse (lon,lev,lat) order for 3D fields \n");
}
/*===========================================================================*/
