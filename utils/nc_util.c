/*$Id: nc_util.c,v 1.5 1998/11/28 22:11:05 eaton Exp $*/

/* Utilities for creating netCDF files from CCM files */

#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <netcdf.h>

#include "nc_util.h"
#include "ccmhdr.h"
#include "ccm_util.h"

/*===========================================================================*/

void
init_nc( CCM_HDR *hdr, int ncid, char* reqvar[], char* cmndline,
	 var_id *vid, int dim_ord, int coards, int doubleOut )
{
  /*
    Input arguments:
    hdr         CCM header structure.
    ncid        netCDF file id.
    reqvar      requested fields from the CCM file.
    cmndline    command line that invoked this code
    dim_ord     dimension ordering flag
    coards      1 => follow COARDS conventions
    doubleOut   0 => NC_FLOAT output, 1 => NC_DOUBLE output

    Output arguments:
    vid         structure containing IDs for variables with a record dimension
  */    

  /* IDs for variables that are set in this routine from info in the CCM header records */
  int lon_id;    /* number of distinct longitudes */
  int mlev_id;   /* number of levels for fields at layer midpoints */
  int ilev_id;   /* number of levels for fields at layer interfaces */
  int lev_id;    /* temp var */
  int lat_id;    /* number of Gaussian latitudes */
  int ntime_id;  /* number of time samples */
  int str_id;    /* length of character info strings */
  int p0_id;
  int ntrm_id, ntrn_id, ntrk_id;
  int ndbase_id, nsbase_id, nbdate_id, nbsec_id, mdt_id, mhisf_id;
  int curr_id, first_id, init_id;
  int tibnd_id, sst_id, ozone_id;
  int ilev_vid;

  size_t k;
  int count;
  int i, j, m;
  int rc;
  int oldFillMode;  /* old fill mode */
  size_t itmp;
  Field *fp;   /* pointer to field info in CCM header */

  int ndim;      /* number of dimensions for a field */
  int dim_id[4]; /* array to hold dimension ids */

  char lev_att[40]; /* string for level attribute values */
  char tim_att[40]; /* string for time attribute values */

  Hfloat mlev_val, ilev_val;  /* hybrid coordinate values in the form 1000*(A+B) */
  Hfloat lon_val, dlon;
  Hfloat* lat_val;

  const Hfloat zero = 0.;

  nc_type FPTypeOut = NC_FLOAT;

  char* sourceAtt = "Data converted from CCM History Tape Format";

  char* pressAtt = "\n"
    "Pressure at a grid point (lon(i),lat(j),lev(k)) is computed    \n"
    "using the formula:                                             \n"
    "           p(i,j,k) = A(k)*PO + B(k)*PS(i,j)                   \n"
    "where A, B, PO, and PS are contained in the variables whose    \n"
    "names are given by the attributes of the vertical coordinate   \n"
    "variable A_var, B_var, P0_var, and PS_var respectively.        \n";
  /*--------------------------------------------------------------------------------------*/

  if ( doubleOut == 1 ) FPTypeOut = NC_DOUBLE;

  /* Write Global Attributes */
  {
    char* convAtt;
    if ( coards == 0 ) {
      convAtt = "NCAR-CSM";
    } else {
      convAtt = "COARDS";
    }
    handleError( nc_put_att_text( ncid, NC_GLOBAL, "Conventions", strlen(convAtt)+1, convAtt ) );
  }
  handleError( nc_put_att_text( ncid, NC_GLOBAL, "source", strlen(sourceAtt)+1, sourceAtt ) );
  handleError( nc_put_att_text( ncid, NC_GLOBAL, "case", strlen(hdr->mcase)+1, hdr->mcase ) );
  handleError( nc_put_att_text( ncid, NC_GLOBAL, "title", strlen(hdr->mcstit)+1, hdr->mcstit ) );
  handleError( nc_put_att_text( ncid, NC_GLOBAL, "hybrid_sigma_pressure", strlen(pressAtt)+1, pressAtt ) );
  appendHistAtt( ncid, ncid, cmndline );

  /* Define dimensions. */
  handleError( nc_def_dim( ncid, "lon", hdr->nlon, &lon_id ) );
  handleError( nc_def_dim( ncid, "lev", hdr->nlev, &mlev_id ) );
  handleError( nc_def_dim( ncid, "ilev", hdr->nlev+1, &ilev_id ) );
  handleError( nc_def_dim( ncid, "lat", hdr->norec, &lat_id ) );
  handleError( nc_def_dim( ncid, "nchar", 81, &str_id ) );
  handleError( nc_def_dim( ncid, "time", NC_UNLIMITED, &ntime_id ) );

  /* non-record variables not in the CCM header */
  {
    char* name = "reference pressure";
    char* units = "Pa";
    handleError( nc_def_var( ncid, "P0", FPTypeOut, 0, 0, &p0_id ) );
    handleError( nc_put_att_text( ncid, p0_id, "long_name", strlen(name)+1, name ) );
    handleError( nc_put_att_text( ncid, p0_id, "units", strlen(units)+1, units ) );
  }

  /* non-record variables from the 1st set of header records */
  {
    char* name = "spectral truncation parameter M";
    handleError( nc_def_var( ncid, "ntrm", NC_INT, 0, 0, &ntrm_id ) );
    handleError( nc_put_att_text( ncid, ntrm_id, "long_name", strlen(name)+1, name ) );
  }
  {
    char* name = "spectral truncation parameter N";
    handleError( nc_def_var( ncid, "ntrn", NC_INT, 0, 0, &ntrn_id ) );
    handleError( nc_put_att_text( ncid, ntrn_id, "long_name", strlen(name)+1, name ) );
  }
  {
    char* name = "spectral truncation parameter K";
    handleError( nc_def_var( ncid, "ntrk", NC_INT, 0, 0, &ntrk_id ) );
    handleError( nc_put_att_text( ncid, ntrk_id, "long_name", strlen(name)+1, name ) );
  }
  {
    char* name = "base day for this case";
    handleError( nc_def_var( ncid, "ndbase", NC_INT, 0, 0, &ndbase_id ) );
    handleError( nc_put_att_text( ncid, ndbase_id, "long_name", strlen(name)+1, name ) );
  }
  {
    char* name = "seconds to complete base day";
    char* units = "s";
    handleError( nc_def_var( ncid, "nsbase", NC_INT, 0, 0, &nsbase_id ) );
    handleError( nc_put_att_text( ncid, nsbase_id, "long_name", strlen(name)+1, name ) );
    handleError( nc_put_att_text( ncid, nsbase_id, "units", strlen(units)+1, units ) );
  }
  {
    char* name = "base date as 8 digit integer (YYYYMMDD)";
    handleError( nc_def_var( ncid, "nbdate", NC_INT, 0, 0, &nbdate_id ) );
    handleError( nc_put_att_text( ncid, nbdate_id, "long_name", strlen(name)+1, name ) );
  }
  {
    char* name = "seconds to complete base date";
    char* units = "s";
    handleError( nc_def_var( ncid, "nbsec", NC_INT, 0, 0, &nbsec_id ) );
    handleError( nc_put_att_text( ncid, nbsec_id, "long_name", strlen(name)+1, name ) );
    handleError( nc_put_att_text( ncid, nbsec_id, "units", strlen(units)+1, units ) );
  }
  {
    char* name = "model timestep";
    char* units = "s";
    handleError( nc_def_var( ncid, "mdt", NC_INT, 0, 0, &mdt_id ) );
    handleError( nc_put_att_text( ncid, mdt_id, "long_name", strlen(name)+1, name ) );
    handleError( nc_put_att_text( ncid, mdt_id, "units", strlen(units)+1, units ) );
  }
  {
    char* name = "frequency of model writes (timesteps)";
    handleError( nc_def_var( ncid, "mhisf", NC_INT, 0, 0, &mhisf_id ) );
    handleError( nc_put_att_text( ncid, mhisf_id, "long_name", strlen(name)+1, name ) );
  }
  {
    char* name = "MSS pathname of this file";
    handleError( nc_def_var( ncid, "current_mss", NC_CHAR, 1, &str_id, &curr_id ) );
    handleError( nc_put_att_text( ncid, curr_id, "long_name", strlen(name)+1, name ) );
    handleError( nc_put_att_text( ncid, curr_id, "date", strlen(hdr->ldhstc)+1, hdr->ldhstc ) );
    handleError( nc_put_att_text( ncid, curr_id, "time", strlen(hdr->lthstc)+1, hdr->lthstc ) );
    handleError( nc_put_att_text( ncid, curr_id, "seq", strlen(hdr->lshstc)+1, hdr->lshstc ) );
  }
  {
    char* name = "MSS pathname of first file for this case";
    handleError( nc_def_var( ncid, "first_mss", NC_CHAR, 1, &str_id, &first_id ) );
    handleError( nc_put_att_text( ncid, first_id, "long_name", strlen(name)+1, name ) );
    handleError( nc_put_att_text( ncid, first_id, "date", strlen(hdr->ldhstf)+1, hdr->ldhstf ) );
    handleError( nc_put_att_text( ncid, first_id, "time", strlen(hdr->lthstf)+1, hdr->lthstf ) );
    handleError( nc_put_att_text( ncid, first_id, "seq", strlen(hdr->lshstf)+1, hdr->lshstf ) );
  }
  {
    char* name = "MSS pathname of initial data";
    handleError( nc_def_var( ncid, "init_mss", NC_CHAR, 1, &str_id, &init_id ) );
    handleError( nc_put_att_text( ncid, init_id, "long_name", strlen(name)+1, name ) );
    handleError( nc_put_att_text( ncid, init_id, "date", strlen(hdr->ldhsti)+1, hdr->ldhsti ) );
    handleError( nc_put_att_text( ncid, init_id, "time", strlen(hdr->lthsti)+1, hdr->lthsti ) );
    handleError( nc_put_att_text( ncid, init_id, "seq", strlen(hdr->lshsti)+1, hdr->lshsti ) );
  }
  {
    char* name = "MSS pathname of time-invariant boundary data";
    handleError( nc_def_var( ncid, "tibnd_mss", NC_CHAR, 1, &str_id, &tibnd_id ) );
    handleError( nc_put_att_text( ncid, tibnd_id, "long_name", strlen(name)+1, name ) );
    handleError( nc_put_att_text( ncid, tibnd_id, "date", strlen(hdr->ldhstt)+1, hdr->ldhstt ) );
    handleError( nc_put_att_text( ncid, tibnd_id, "time", strlen(hdr->lthstt)+1, hdr->lthstt ) );
    handleError( nc_put_att_text( ncid, tibnd_id, "seq", strlen(hdr->lshstt)+1, hdr->lshstt ) );
  }
  {
    char* name = "MSS pathname of SST boundary data";
    handleError( nc_def_var( ncid, "sst_mss", NC_CHAR, 1, &str_id, &sst_id ) );
    handleError( nc_put_att_text( ncid, sst_id, "long_name", strlen(name)+1, name ) );
    handleError( nc_put_att_text( ncid, sst_id, "date", strlen(hdr->ldhsts)+1, hdr->ldhsts ) );
    handleError( nc_put_att_text( ncid, sst_id, "time", strlen(hdr->lthsts)+1, hdr->lthsts ) );
    handleError( nc_put_att_text( ncid, sst_id, "seq", strlen(hdr->lshsts)+1, hdr->lshsts ) );
  }
  {
    char* name = "MSS pathname of ozone boundary data";
    handleError( nc_def_var( ncid, "ozone_mss", NC_CHAR, 1, &str_id, &ozone_id ) );
    handleError( nc_put_att_text( ncid, ozone_id, "long_name", strlen(name)+1, name ) );
    handleError( nc_put_att_text( ncid, ozone_id, "date", strlen(hdr->ldhsto)+1, hdr->ldhsto ) );
    handleError( nc_put_att_text( ncid, ozone_id, "time", strlen(hdr->lthsto)+1, hdr->lthsto ) );
    handleError( nc_put_att_text( ncid, ozone_id, "seq", strlen(hdr->lshsto)+1, hdr->lshsto ) );
  }
  {
    char* name = "hybrid A coefficient at layer interfaces";
    handleError( nc_def_var( ncid, "hyai", FPTypeOut, 1, &ilev_id, &vid->hyai ) );
    handleError( nc_put_att_text( ncid, vid->hyai, "long_name", strlen(name)+1, name ) );
  }
  {
    char* name = "hybrid B coefficient at layer interfaces";
    handleError( nc_def_var( ncid, "hybi", FPTypeOut, 1, &ilev_id, &vid->hybi ) );
    handleError( nc_put_att_text( ncid, vid->hybi, "long_name", strlen(name)+1, name ) );
  }
  {
    char* name = "hybrid A coefficient at layer midpoints";
    handleError( nc_def_var( ncid, "hyam", FPTypeOut, 1, &mlev_id, &vid->hyam ) );
    handleError( nc_put_att_text( ncid, vid->hyam, "long_name", strlen(name)+1, name ) );
  }
  {
    char* name = "hybrid B coefficient at layer midpoints";
    handleError( nc_def_var( ncid, "hybm", FPTypeOut, 1, &mlev_id, &vid->hybm ) );
    handleError( nc_put_att_text( ncid, vid->hybm, "long_name", strlen(name)+1, name ) );
  }
  {
    char* units;
    char* name = "hybrid level at layer midpoints (1000*(A+B))";
    char* at1 = "hyam";
    char* at2 = "hybm";
    char* at3 = "P0";
    char* at4 = "PS";
    char* at5 = "down";
    char* at6 = "ilev";
    handleError( nc_def_var( ncid, "lev", FPTypeOut, 1, &mlev_id, &vid->mlev ) );
    handleError( nc_put_att_text( ncid, vid->mlev, "long_name", strlen(name)+1, name ) );
    if ( coards == 0 ) {
      units = "hybrid_sigma_pressure";
    } else {
      units = "level";
    }
    handleError( nc_put_att_text( ncid, vid->mlev, "units", strlen(units)+1, units ) );
    handleError( nc_put_att_text( ncid, vid->mlev, "positive", strlen(at5)+1, at5 ) );
    handleError( nc_put_att_text( ncid, vid->mlev, "A_var", strlen(at1)+1, at1 ) );
    handleError( nc_put_att_text( ncid, vid->mlev, "B_var", strlen(at2)+1, at2 ) );
    handleError( nc_put_att_text( ncid, vid->mlev, "P0_var", strlen(at3)+1, at3 ) );
    handleError( nc_put_att_text( ncid, vid->mlev, "PS_var", strlen(at4)+1, at4 ) );
    handleError( nc_put_att_text( ncid, vid->mlev, "bounds", strlen(at6)+1, at6 ) );
  }
  {
    char* units;
    char* name = "hybrid level at layer interfaces (1000*(A+B))";
    char* at1 = "hyai";
    char* at2 = "hybi";
    char* at3 = "P0";
    char* at4 = "PS";
    char* at5 = "down";
    handleError( nc_def_var( ncid, "ilev", FPTypeOut, 1, &ilev_id, &ilev_vid ) );
    handleError( nc_put_att_text( ncid, ilev_vid, "long_name", strlen(name)+1, name ) );
    if ( coards == 0 ) {
      units = "hybrid_sigma_pressure";
    } else {
      units = "level";
    }
    handleError( nc_put_att_text( ncid, ilev_vid, "units", strlen(units)+1, units ) );
    handleError( nc_put_att_text( ncid, ilev_vid, "positive", strlen(at5)+1, at5 ) );
    handleError( nc_put_att_text( ncid, ilev_vid, "A_var", strlen(at1)+1, at1 ) );
    handleError( nc_put_att_text( ncid, ilev_vid, "B_var", strlen(at2)+1, at2 ) );
    handleError( nc_put_att_text( ncid, ilev_vid, "P0_var", strlen(at3)+1, at3 ) );
    handleError( nc_put_att_text( ncid, ilev_vid, "PS_var", strlen(at4)+1, at4 ) );
  }
  {
    char* name = "latitude";
    char* units = "degrees_north";
    handleError( nc_def_var( ncid, "lat", FPTypeOut, 1, &lat_id, &vid->lat ) );
    handleError( nc_put_att_text( ncid, vid->lat, "long_name", strlen(name)+1, name ) );
    handleError( nc_put_att_text( ncid, vid->lat, "units", strlen(units)+1, units ) );
  }
  {
    char* name = "gauss weights";
    handleError( nc_def_var( ncid, "gw", FPTypeOut, 1, &lat_id, &vid->gw ) );
    handleError( nc_put_att_text( ncid, vid->gw, "long_name", strlen(name)+1, name ) );
  }
  {
    char* name = "longitude";
    char* units = "degrees_east";
    handleError( nc_def_var( ncid, "lon", FPTypeOut, 1, &lon_id, &vid->lon ) );
    handleError( nc_put_att_text( ncid, vid->lon, "long_name", strlen(name)+1, name ) );
    handleError( nc_put_att_text( ncid, vid->lon, "units", strlen(units)+1, units ) );
  }
  {
    char* name = "number of longitudes at each latitude";
    handleError( nc_def_var( ncid, "nlons", NC_INT, 1, &lat_id, &vid->nlons ) );
    handleError( nc_put_att_text( ncid, vid->nlons, "long_name", strlen(name)+1, name ) );
  }

  /* record variables */
  {
    char* name = "current date as 8 digit integer (YYYYMMDD)";
    handleError( nc_def_var( ncid, "date", NC_INT, 1, &ntime_id, &vid->date ) );
    handleError( nc_put_att_text( ncid, vid->date, "long_name", strlen(name)+1, name ) );
  }
  {
    char* name = "seconds to complete current date";
    char* units = "s";
    handleError( nc_def_var( ncid, "datesec", NC_INT, 1, &ntime_id, &vid->datesec ) );
    handleError( nc_put_att_text( ncid, vid->datesec, "long_name", strlen(name)+1, name ) );
    handleError( nc_put_att_text( ncid, vid->datesec, "units", strlen(units)+1, units ) );
  }
  {
    char* name = "time";
    char units[31];
    time_units( hdr->nbdate, hdr->nbsec, units );
    handleError( nc_def_var( ncid, "time", NC_DOUBLE, 1, &ntime_id, &vid->time ) );
    handleError( nc_put_att_text( ncid, vid->time, "long_name", strlen(name)+1, name ) );
    handleError( nc_put_att_text( ncid, vid->time, "units", strlen(units)+1, units ) );
  }

  /* variables for the model fields */
  fp = hdr->fld;

  /* number of requested fields */
  if ( reqvar == NULL ) {
    vid->nfld = hdr->nfldh;
    /* create reqvar list containing all fields from CCM header */
    reqvar = (char**)malloc( sizeof(char*) * hdr->nfldh );
    for ( i=0; i < hdr->nfldh; ++i ) reqvar[i] = (fp+i)->name;
  } else {
    for ( vid->nfld=0, i=0; reqvar[i] != NULL; vid->nfld++, i++ );
  }

  /* allocate int array for field variable ids */
  vid->fld = (int *)malloc( vid->nfld * sizeof( int ) );
  vid->ccmfld = (int *)malloc( vid->nfld * sizeof( int ) );

  /* loop over requested fields */
  for ( i=0; i < vid->nfld; ++i ) {

    /* Find requested field in CCM header */
    for ( m=0; m < hdr->nfldh; ++m ) {

      if ( strcmp( reqvar[i], (fp+m)->name ) != 0  ) {
	if ( m == hdr->nfldh-1 ) {
	  fprintf( stderr, "init_nc: requested field %s not found\n", reqvar[i] );
	  exit( EXIT_FAILURE );
	} else {
	  continue;
	}
      } else {

        /* save index for requested field's data in the list of fields
           in the CCM history file */
	vid->ccmfld[i] = m;

	/* Set dimension id array */
	dim_id[0] = ntime_id;
	if ( (fp+m)->lev == 0 ) {
	  ndim = 3;
	  dim_id[1] = lat_id;
	  dim_id[2] = lon_id;
	} else if ( (fp+m)->lev == 1 ) {
	  ndim = 4;
	  lev_id = ilev_id;
	} else if ( (fp+m)->lev == 2 ) {
	  ndim = 4;
	  lev_id = mlev_id;
	}
	if ( ndim == 4 ) {
	  if ( dim_ord == 0 ) {     /* (lon,lat,lev,tim) */
	    dim_id[1] = lev_id;
	    dim_id[2] = lat_id;
	    dim_id[3] = lon_id;
	  } else {                  /* (lon,lev,lat,tim) */
	    dim_id[1] = lat_id;
	    dim_id[2] = lev_id;
	    dim_id[3] = lon_id;
	  }
	}

	/* define the field variable */
	if ( strcmp( (fp+m)->name, checkName( (fp+m)->name ) ) != 0 ) {
	  printf("CCM variable %s is not a legal netCDF name.  Changed to %s\n",
		 (fp+m)->name, checkName( (fp+m)->name ) );
	}
	handleError( nc_def_var( ncid, checkName( (fp+m)->name ), FPTypeOut, ndim, dim_id, &vid->fld[i] ) );
	handleError( nc_put_att_text( ncid, vid->fld[i], "long_name", strlen( (fp+m)->long_name )+1, (fp+m)->long_name ) );
	if ( (fp+m)->udunit != NULL ) {
	  itmp = strlen( (fp+m)->udunit );
	  handleError( nc_put_att_text( ncid, vid->fld[i], "units", itmp+1, (fp+m)->udunit ) );
	}

	/* time representation */
	if ( (fp+m)->tim == 1 ) {
	  strcpy( tim_att, "average" );
	} else if ( (fp+m)->tim == 2 ) {
	  strcpy( tim_att, "minimum" );
	} else if ( (fp+m)->tim == 3 ) {
	  strcpy( tim_att, "maximum" );
	}
	if ( (fp+m)->tim != 0 ) {
	  handleError( nc_put_att_text( ncid, vid->fld[i], "time_op", strlen( tim_att )+1, tim_att ) );
	}

	break;  /* requested field found... break out of loop over CCM fields */
      }

    }  /* end loop over CCM fields */

  }  /* end loop over requested fields */

  /* Turn off prefilling of values (for efficiency) */
  handleError( nc_set_fill( ncid, NC_NOFILL, &oldFillMode ) );

  /* End define mode */
  handleError( nc_enddef( ncid ) );

  { /* reference pressure (Pa) for hybrid coordinate system */
    Hfloat val = 1.e5;
    ncPutVarHfloat( ncid, p0_id, &val );
  }

  /* Set non-record variables from 1st set of headers. */
  handleError( nc_put_var_int( ncid, ntrm_id, &hdr->ntrm ) );
  handleError( nc_put_var_int( ncid, ntrn_id, &hdr->ntrn ) );
  handleError( nc_put_var_int( ncid, ntrk_id, &hdr->ntrk ) );
  handleError( nc_put_var_int( ncid, ndbase_id, &hdr->ndbase ) );
  handleError( nc_put_var_int( ncid, nsbase_id, &hdr->nsbase ) );
  handleError( nc_put_var_int( ncid, nbdate_id, &hdr->nbdate ) );
  handleError( nc_put_var_int( ncid, nbsec_id, &hdr->nbsec ) );
  handleError( nc_put_var_int( ncid, mdt_id, &hdr->mdt ) );
  handleError( nc_put_var_int( ncid, mhisf_id, &hdr->mhisf ) );
  handleError( nc_put_var_text( ncid, curr_id, hdr->lnhstc ) );
  handleError( nc_put_var_text( ncid, first_id, hdr->lnhstf ) );
  handleError( nc_put_var_text( ncid, init_id, hdr->lnhsti ) );
  handleError( nc_put_var_text( ncid, tibnd_id, hdr->lnhstt ) );
  handleError( nc_put_var_text( ncid, sst_id, hdr->lnhsts ) );
  handleError( nc_put_var_text( ncid, ozone_id, hdr->lnhsto ) );

  if ( hdr->mphya != NULL ) {
    for ( k=0; k < hdr->nlev; ++k ) {
      ncPutVar1Hfloat( ncid, vid->hyam, &k, hdr->mphya+2*k+1 );
      ncPutVar1Hfloat( ncid, vid->hybm, &k, hdr->mphyb+2*k+1 );
      mlev_val = 1000.*( *(hdr->mphya+2*k+1) + *(hdr->mphyb+2*k+1) );
      ncPutVar1Hfloat( ncid, vid->mlev, &k, &mlev_val );
    }
    for ( k=0; k <= hdr->nlev; ++k ) {
      ncPutVar1Hfloat( ncid, vid->hyai, &k, hdr->mphya+2*k );
      ncPutVar1Hfloat( ncid, vid->hybi, &k, hdr->mphyb+2*k );
      ilev_val = 1000.*( *(hdr->mphya+2*k) + *(hdr->mphyb+2*k) );
      ncPutVar1Hfloat( ncid, ilev_vid, &k, &ilev_val );
    }
  } else {
    for ( k=0; k < hdr->nlev; ++k ) {
      ncPutVar1Hfloat( ncid, vid->hyam, &k, &zero );
      ncPutVar1Hfloat( ncid, vid->hybm, &k, hdr->mphyb+2*k+1 );
      mlev_val = 1000.*( *(hdr->mphyb+2*k+1) );
      ncPutVar1Hfloat( ncid, vid->mlev, &k, &mlev_val );
    }
    for ( k=0; k <= hdr->nlev; ++k ) {
      ncPutVar1Hfloat( ncid, vid->hyai, &k, &zero );
      ncPutVar1Hfloat( ncid, vid->hybi, &k, hdr->mphyb+2*k );
      ilev_val = 1000.*( *(hdr->mphyb+2*k) );
      ncPutVar1Hfloat( ncid, ilev_vid, &k, &ilev_val );
    }
  }

  /* Add latitude variable */
  count = hdr->norec;
  /* Make sure the latitudes are ordered S -> N */
  if ( hdr->mplat[0] > hdr->mplat[1] ) {
    lat_val = (Hfloat*)malloc( count * sizeof( Hfloat ) );
    for ( j=0; j<count; ++j ) lat_val[j] = hdr->mplat[count-j-1];
    ncPutVarHfloat( ncid, vid->lat, lat_val );
    free( lat_val );
  } else {
    ncPutVarHfloat( ncid, vid->lat, hdr->mplat );
  }
  ncPutVarHfloat( ncid, vid->gw, hdr->mpwts );

  /* Add longitude variable */
  dlon = 360./hdr->nlon;
  for ( k=0; k < hdr->nlon; k++ ) {
    lon_val = k * dlon;
    ncPutVar1Hfloat( ncid, vid->lon, &k, &lon_val );
  }

}
/*===========================================================================*/

void
write_ncrec( FFILE ff, CCM_HDR* h, int ncid, var_id* vid, size_t its,
	     int dim_ord, int doubleOut )
{
  /* Read the data records in a ccm history tape for a single time sample
     and put the data for each field into a netCDF file.
     It is assumed that the CCM file is positioned at the beginning of the
     first data record.

     Input arguments:
     ff         Input stream.
     h          CCM history file header.
     ncid       CDF file id.
     vid        variable id's
     its        index for current time sample
     dim_ord    dimension ordering flag
     doubleOut  0 => NC_FLOAT output, 1 => NC_DOUBLE output
  */

  int nfld;                     /* number of fields */
  int nlon, nlat, nlev, nlevi;  /* field dimensions */
  int len2D, len3D, len3Di;     /* field lengths */

  Hfloat** datap;     /* pointers to data for each requested field */
  Hfloat* tmp3D;      /* temp data storage to rearrange 3D field at layer midpoints */
  Hfloat* tmp3Di;     /* temp data storage to rearrange 3D field at layer interfaces */
  Hfloat* tmpp;
  Hfloat* fldslice;   /* pointer to latitude slice of data for a field. */
  Hfloat drec12[2];   /* latitude index and number of longitudes in a data record */
  void* drec;         /* data record */
  int nrecbytes;      /* number of bytes in a data record */
  Field* fldm;        /* pointer to current CCM field structure */
  int fldlen;         /* field length */

  int i, j, stat;
  int lat;            /* latitude index */
  int* nlons;         /* number of longitudes in latitude slice */
  int rval;           /* function return value */
  Hfloat time;        /* day no. + fraction */
  Hfloat spval;       /* internal rep of missing value used in CCM history files */
  size_t start[] = { 0, 0, 0, 0 };
  size_t count[] = { 1, 1, 1, 1 };

  nc_type FPTypeOut = NC_FLOAT;

  char* mallocErrorMsg = "  malloc unable to allocate requested memory.  Must be able to allocate\n"
                         "  enough memory (using 8-byte words) for a global copy of each requested\n"
                         "  field for a single time sample.\n";
  /*----------------------------------------------------------------------------------------------*/

  if ( doubleOut == 1 ) FPTypeOut = NC_DOUBLE;

  nfld = vid->nfld;
  nlon = h->nlon;
  nlat = h->norec;
  nlev = h->nlev;
  nlevi = nlev + 1;
  len2D = nlon*nlat;
  len3D = nlon*nlat*nlev;
  len3Di = nlon*nlat*nlevi;

  /* Allocate memory for array of pointers to the data for all requested fields. */
  datap = (Hfloat**)malloc( sizeof(Hfloat*)*nfld );

  if ( dim_ord == 0 ) {
    /* Allocate memory for 3D fields for rearranging data */
    if ( (tmp3D=(Hfloat*)malloc( sizeof(Hfloat)*len3D )) == NULL ) {
      fprintf( stderr, "write_ncrec: FATAL ERROR:\n %s", mallocErrorMsg );
      exit( EXIT_FAILURE );
    }
    if ( (tmp3Di=(Hfloat*)malloc( sizeof(Hfloat)*len3Di )) == NULL ) {
      fprintf( stderr, "write_ncrec: FATAL ERROR:\n %s", mallocErrorMsg );
      exit( EXIT_FAILURE );
    }
  }

  /* Allocate memory for global data for each requested field */
  for ( i=0; i<nfld; ++i ) {
    fldm = h->fld + vid->ccmfld[i];
    switch ( fldm->lev ) {
    case 0:
      if ( (datap[i]=(Hfloat*)malloc( sizeof(Hfloat)*len2D )) == NULL ) {
	fprintf( stderr, "write_ncrec: FATAL ERROR:\n %s", mallocErrorMsg );
	exit( EXIT_FAILURE );
      }
      break;
    case 1:
      if ( (datap[i]=(Hfloat*)malloc( sizeof(Hfloat)*len3Di )) == NULL ) {
	fprintf( stderr, "write_ncrec: FATAL ERROR:\n %s", mallocErrorMsg );
	exit( EXIT_FAILURE );
      }
      /* init this data to zero since the history tape is missing the data for
	 one level */
      for ( j=0; j<len3Di; ++j ) *(datap[i]+j) = 0.;
      break;
    case 2:
      if ( (datap[i]=(Hfloat*)malloc( sizeof(Hfloat)*len3D )) == NULL ) {
	fprintf( stderr, "write_ncrec: FATAL ERROR:\n %s", mallocErrorMsg );
	exit( EXIT_FAILURE );
      }
      break;
    }
  }

  /* Allocate memory for a CCM data record */
  if ( (h->fld)->pdens == 1 ) {
    nrecbytes = ff.fsz*h->maxsiz;
  } else {
    nrecbytes = 8*h->maxsiz;
  }
  if ( (drec=malloc( nrecbytes )) == NULL ) {
    fprintf( stderr, "write_ncrec: FATAL ERROR:\n %s", mallocErrorMsg );
    exit( EXIT_FAILURE );
  }

  /* Allocate memory for nlons */
  if ( (nlons=(int*)malloc( sizeof(int)*nlat )) == NULL ) {
    fprintf( stderr, "write_ncrec: FATAL ERROR:\n %s", mallocErrorMsg );
    exit( EXIT_FAILURE );
  }

  /* Loop over latitudes */
  for ( j=0; j<nlat; ++j ) {

    if ( (rval=read_rec( ff, drec, nrecbytes )) != nrecbytes ) {
      fprintf( stderr, "write_ncrec: read_rec read %d bytes from data record %d\n", rval, j+1 );
      fprintf( stderr, "             %d bytes were requested\n", nrecbytes );
      exit( EXIT_FAILURE );
    }

    if ( h->mftyp%10 == 3 ) { /* Get latitude index from first word in data record */
      if ( need_fcv( ff ) ) {
        if ( cv_floats( ff, 2, drec, &drec12 ) != 0 ) {
	  fprintf( stderr, "write_ncrec: cv_floats encountered error converting between\n"
		           "             CRAY and IEEE floating point representations\n" );
	  exit( EXIT_FAILURE );
	}
        lat = drec12[0] + .1;
        nlons[lat-1] = drec12[1] + .1;
      } else {
        lat = *( (Hfloat*)drec ) + .1;
        nlons[lat-1] = *( ((Hfloat*)drec+1) ) + .1;
      }
    } else {                  /* ping-pong order, N -> S */
      if ( j%2 == 0 ) {
        lat = nlat - j/2;
	nlons[lat-1] = nlon;
      } else {
        lat = j/2 + 1;
	nlons[lat-1] = nlon;
      }
    }

    /* Loop over fields */
    for ( i=0; i<nfld; ++i ) {
      fldm = h->fld + vid->ccmfld[i];
      switch ( fldm->lev ) {
      case 0:
	fldslice = datap[i] + (lat-1)*nlon;
	break;
      case 1:
	/* assume that the data missing from the history tape is the top interface */
	fldslice = datap[i] + (lat-1)*nlon*nlevi + nlon;
	break;
      case 2:
	fldslice = datap[i] + (lat-1)*nlon*nlev;
	break;
      }
      /* Copy a latitude slice of data for the current field into the global array */
      if ( (stat=get_slice( fldm->name, ff, h, drec, fldslice )) < 0 ) {
	fprintf( stderr, "write_ncrec: get_slice returned %d\n", stat );
	exit( EXIT_FAILURE );
      }
    }
  }        /* end loop over latitudes */

  /* Write nlons array for 1st time sample only */
  if ( its == 0 ) {
    handleError( nc_put_var_int( ncid, vid->nlons, nlons ) );
  }

  /* Check for special values in any field if the history tape data is
     unpacked, and if a special value has not already been found in that
     field.  If a special value is found, then add the special value 
     attribute for that field. */
  for ( i=0; i<nfld; ++i ) {
    fldm = h->fld + vid->ccmfld[i];
    if ( fldm->pdens == 1 ) {
      if ( fldm->spval_flg == 0 ) {
	switch ( fldm->lev ) {
	case 0:
	  fldlen = len2D;
	  break;
	case 1:
	  fldlen = len3Di;
	  break;
	case 2:
	  fldlen = len3D;
	  break;
	}
	if ( spval_exists( ff, datap[i], fldlen, &spval ) ) {
	  fldm->spval_flg = 1;
	  handleError( nc_redef( ncid ) );
	  ncPutAttHfloat( ncid, vid->fld[i], "missing_value", FPTypeOut, 1, &spval );
	  handleError( nc_enddef( ncid ) );
	}
      }
    }
  }

  /* Reorder 3D output fields if [lev][lat][lon] order is requested */
  if ( dim_ord == 0 ) {
    for ( i=0; i<nfld; ++i ) {
      fldm = h->fld + vid->ccmfld[i];
      switch ( fldm->lev ) {
      case 0:
	break;
      case 1:
	switch_lev_lat( datap[i], tmp3Di, nlon, nlat, nlevi );
	tmpp = datap[i];
	datap[i] = tmp3Di;
	tmp3Di = tmpp;
	break;
      case 2:
	switch_lev_lat( datap[i], tmp3D, nlon, nlat, nlev );
	tmpp = datap[i];
	datap[i] = tmp3D;
	tmp3D = tmpp;
	break;
      }
    }
  }

  /* Put data from CCM header records into netCDF file */
  handleError( nc_put_var1_int( ncid, vid->date, &its, &h->ncdate ) );
  handleError( nc_put_var1_int( ncid, vid->datesec, &its, &h->ncsec ) );
  time = h->ndcur + h->nscur/86400.;
  ncPutVar1Hfloat( ncid, vid->time, &its, &time );

  /* Put field data into netCDF file */
  start[0] = its;
  count[3] = nlon;
  for ( i=0; i<nfld; ++i ) {
    fldm = h->fld + vid->ccmfld[i];
    switch ( fldm->lev ) {
    case 0:
      count[1] = nlat;
      count[2] = nlon;
      break;
    case 1:
      count[1] = nlevi;
      count[2] = nlat;
      if ( dim_ord == 1 ) {
	count[1] = nlat;
	count[2] = nlevi;
      }
      break;
    case 2:
      count[1] = nlev;
      count[2] = nlat;
      if ( dim_ord == 1 ) {
	count[1] = nlat;
	count[2] = nlev;
      }
      break;
    }
    ncPutVaraHfloat( ncid, vid->fld[i], start, count, datap[i] );
  }

  /* Free memory */
  for ( i=0; i<nfld; ++i ) free( datap[i] );
  free( datap );
  free( drec );
  if ( dim_ord == 0 ) {
    free( tmp3D );
    free( tmp3Di );
  }
  free( nlons );
}
/*===========================================================================*/

char* checkName( char* name )
{
  /* Check that name is a legal netCDF variable or attribute name.
     Replace illegal characters with '_', unless it is the first character in
     which case replace it with an 'x'.
  */

  char* newStr;
  char* s;

  /* names must be non-null strings */
  if ( strlen(name) == 0 ) return "x";

  newStr = (char*)malloc( sizeof(char) * (strlen(name)+1) );
  s = newStr;

  /* First character must be alphabetic */
  if ( isalpha( *name ) ) {
    *s++ = *name++;
  } else {
    *s++ = 'x';
    ++name;
  }

  /* Subsequent characters can be alphanumeric, '_', or '-'. */
  while( *name != '\0' ) {
    if ( isalnum( *name ) || strncmp( name, "_", 1 )==0 || strncmp( name, "-", 1 )==0 ) {
      *s++ = *name++;
    } else {
      *s++ = '_';
      ++name;
    }
  }
  *s = '\0';
  return newStr;
}
/*===========================================================================*/

int
time_units( int date, int sec, char* units )
{
  /* construct time units string

     Input args:

     date    yymmdd
     sec     seconds relative to date at 0Z

     Output arg:

     units   string in the form "days since 0000-00-00 00:00:00"

     Return values:

      0  Successful return.
  */

  int year, month, day, hours, minutes, seconds;

  year = date/10000;
  month = (date%10000)/100;
  day = date%100;

  hours = sec/3600;
  minutes = (sec - 3600*hours)/60;
  seconds = sec - 3600*hours - 60*minutes;

  sprintf( units, "days since %04i-%02i-%02i %02i:%02i:%02i", year, month,
                   day, hours, minutes, seconds );

  return 0;
}
/*===========================================================================*/

void
switch_lev_lat( Hfloat* ccmdat, Hfloat* ncdat, int nlon, int nlat,
		int nlev )
{
  /* Reorder the input array from [nlat][nlev][nlon] to [nlev][nlat][nlon]

     Input arguments:
     ccmdat - array ordered [nlat][nlev][nlon]
     nlon   - length of longitude dimension
     nlat   - length of latitude dimension
     nlev   - length of level dimension

     Output arguments:
     ncdat  - array ordered [nlev][nlat][nlon]
  */

  int i, j, k;
  int joff, koff;      /* offsets for the output array */
  int nlonlat;
  
  nlonlat = nlon*nlat;

  /* Loop over the elements of the input array in storage order and copy them
     to the appropriate place in the output array. */
  for ( j = 0; j < nlat; ++j ) {
    joff = j*nlon;
    for ( k = 0; k < nlev; ++k ) {
      koff = k*nlonlat;
      for ( i = 0; i < nlon; ++i ) {
	*( ncdat + i + joff + koff ) = *ccmdat++;
      }
    }
  }
}
/*===========================================================================*/

void
handleError( int status )
{
  if ( status != NC_NOERR ) {
    fprintf(stderr, "%s\n", nc_strerror(status));
    abort();
  }
}
/*===========================================================================*/

void
appendHistAtt( int ncIn, int ncOut, char* cmndline )
{
  /* Create a history entry for the current operator/utility by appending the
     input command line to the current timestamp.  If the input netCDF file has
     a history attribute then append the current history entry to it and write
     that to the output file.  If the input netCDF file has no history attribute 
     then write the current history entry as the history attribute to the output
     netCDF file.  In the case where a new netCDF file is being created then
     setting ncIn=ncOut will write the history entry to the output file.
  */

  char* histLine;  /* history entry for current operator/utility */
  char* ctimeStr;
  char timeStampStr[25];
  time_t clock;
  size_t dataLen;
  nc_type dataType;
  char* histAttIn;
  char* histAttOut;  /* history attribute for input netCDF file with current history
			entry appended */

  /* create history entry for current operator/utility */
  clock = time( (time_t*)NULL );
  ctimeStr = ctime( &clock );
  /* Get rid of the carriage return in the ctimeStr */ 
  (void)strncpy( timeStampStr, ctimeStr, 24 );
  timeStampStr[24] = '\0';
  histLine = (char*)malloc( (30+strlen(cmndline))*sizeof(char) );
  strcpy( histLine, "\n" );
  strcat( histLine, timeStampStr );
  strcat( histLine, "> " );
  strcat( histLine, cmndline );

  /* Check to see if history attribute is present */
  if ( nc_inq_att( ncIn, NC_GLOBAL, "history", &dataType, &dataLen ) == NC_NOERR ) {
    histAttIn = (char*)malloc( (dataLen+1)*sizeof(char) );
    /* Get history att from input file and make sure it's a string */
    handleError( nc_get_att_text( ncIn, NC_GLOBAL, "history", histAttIn ) );
    histAttIn[dataLen] = '\0';
    histAttOut = (char*)malloc( (strlen(histAttIn)+strlen(histLine)+1)*sizeof(char) );
    strcpy( histAttOut, histAttIn );
    strcat( histAttOut, histLine );
  } else {
    histAttOut = histLine;
  }
  handleError( nc_put_att_text( ncOut, NC_GLOBAL, "history", strlen(histAttOut)+1, histAttOut ) );
}
/*===========================================================================*/

void
ncPutVarHfloat( int ncid, int vid, const Hfloat* data )
{
#ifdef CRAY_FLOAT
  handleError( nc_put_var_float( ncid, vid, data ) );
#else
  handleError( nc_put_var_double( ncid, vid, data ) );
#endif
}
/*===========================================================================*/

void
ncPutVar1Hfloat( int ncid, int vid, const size_t index[], const Hfloat* data )
{
#ifdef CRAY_FLOAT
  handleError( nc_put_var1_float( ncid, vid, index, data ) );
#else
  handleError( nc_put_var1_double( ncid, vid, index, data ) );
#endif
}
/*===========================================================================*/

void
ncPutVaraHfloat( int ncid, int vid, const size_t start[], const size_t count[], const Hfloat* data )
{
#ifdef CRAY_FLOAT
  handleError( nc_put_vara_float( ncid, vid, start, count, data ) );
#else
  handleError( nc_put_vara_double( ncid, vid, start, count, data ) );
#endif
}
/*===========================================================================*/

void
ncPutAttHfloat( int ncid, int vid, const char* name, nc_type xtype, size_t len, const Hfloat* data )
{
#ifdef CRAY_FLOAT
  handleError( nc_put_att_float( ncid, vid, name, xtype, len, data ) );
#else
  handleError( nc_put_att_double( ncid, vid, name, xtype, len, data ) );
#endif
}
/*===========================================================================*/
