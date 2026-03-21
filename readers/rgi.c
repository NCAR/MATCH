/* $Id$ */
/* rgi.c   regrid interpolation only

   this is for putting just the interpolation part into
   another program, without dragging vmd code into it

   see vmd/regrid.c for algorithm info

   04/07/97: split from regrid.c
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

#ifdef VMD
#include <toolsa/MEMM.h>
#else
#include "str_util.h"
#endif

#define FRIEND_OF_RGI
#include "rgi.h"

#if (defined  SGI ) || ( defined LINUX )
#define TINTERP tinterp_
#define AREAAVG areaavg_
#elif ( defined AIX ) 
#define TINTERP tinterp
#define AREAAVG areaavg
#endif

#if (REALSIZE == 4)
#define REAL float
#else
#define REAL double
#endif

static char *scat( char *s1, char *s2)
    {
	static char s[200];
	strcpy(s, s1);
	strcat(s,s2);
	return s;
    }

/****************************************************/

static void WriteCoordVar( int ncid, char *name, long size, double *v)
    {
        int dimid[1], varid;
        long start[1], count[1];

	ncredef( ncid);
	
        dimid[0] = ncdimdef( ncid, name, size);
        varid = ncvardef( ncid, name, NC_DOUBLE, 1, dimid);

	ncendef( ncid);

        start[0] = 0;
        count[0] = size;
        ncvarput( ncid, varid, start, count, v);
    }

static void WriteGrid( RGgrid_t *rgg, int ncid, char *prefix)
    {
        int dimid[1], varid;
        long start[1], count[1];

        WriteCoordVar( ncid, scat(prefix, "lon"), rgg->nlon, rgg->lon);
        WriteCoordVar( ncid, scat(prefix, "lone"), rgg->nlon+1, rgg->lone);
        WriteCoordVar( ncid, scat(prefix, "lat"), rgg->nlat, rgg->lat);
        WriteCoordVar( ncid, scat(prefix, "late"), rgg->nlat+1, rgg->late);

        ncredef( ncid);
        dimid[0] = ncdimid( ncid, scat(prefix, "lat"));
        varid = ncvardef( ncid, scat(prefix, "gw"), NC_DOUBLE, 1, dimid);
        ncendef( ncid);

	if (rgg->gw != NULL)
	    {
            start[0] = 0;
            count[0] = rgg->nlat;
            ncvarput( ncid, varid, start, count, rgg->gw);
	    }
    }

int RGIwrite( RG_t *rgn, char *filename)
    {
        int ncid;
        int fromidx_varid, pct_varid, start_varid;

        int dimid[2];
        long start[2], count[2];

        ncid = nccreate( filename, NC_CLOBBER);

	ncendef( ncid);

        WriteGrid( rgn->from, ncid, "f");
        WriteGrid( rgn->to, ncid, "t");

	ncredef( ncid);

        dimid[0] = ncdimdef( ncid, "isect", (long) rgn->n_isects);
        fromidx_varid = ncvardef( ncid, "from_idx", NC_LONG, 1, dimid);
        pct_varid = ncvardef( ncid, "pct", NC_FLOAT, 1, dimid);

        dimid[0] = ncdimid( ncid, "tlat");
        dimid[1] = ncdimid( ncid, "tlon");
        start_varid = ncvardef( ncid, "start_isect", NC_LONG, 2, dimid);

	ncendef( ncid);

        start[0] = 0;
        count[0] = rgn->n_isects;
        ncvarput( ncid, fromidx_varid, start, count, rgn->fromidx);
        ncvarput( ncid, pct_varid, start, count, rgn->pct);

        start[0] = 0;
        count[0] = rgn->to->nlat;
        start[1] = 0;
        count[1] = rgn->to->nlon;
        ncvarput( ncid, start_varid, start, count, rgn->start_isect);

        ncclose( ncid);
        printf("wrote regrid file %s\n", filename);
        return 1;
    }

/****************************************************/

static void ReadCoordVar( int ncid, char *name, int *n, double **v)
    {
        int dimid, varid;
	long start[1], count[1];
	long size;

        dimid = ncdimid( ncid, name);
        ncdiminq( ncid, dimid, NULL, &size);
        varid = ncvarid( ncid, name);

        start[0] = 0;
        count[0] = size;
        *v = CALLOC( size, double);
        ncvarget( ncid, varid, start, count, *v);
	*n = size;
    }

static RGgrid_t *ReadGrid( int ncid, char *prefix)
    {
	int dumm, varid;
        long start[1], count[1];

	RGgrid_t *rgg = MALLOC( RGgrid_t);

	ReadCoordVar( ncid, scat(prefix, "lon"), &rgg->nlon, &rgg->lon);
	ReadCoordVar( ncid, scat(prefix, "lone"), &dumm, &rgg->lone);
	ReadCoordVar( ncid, scat(prefix, "lat"), &rgg->nlat, &rgg->lat);
	ReadCoordVar( ncid, scat(prefix, "late"), &dumm, &rgg->late);

	varid = ncvarid( ncid, scat(prefix, "gw"));
        start[0] = 0;
        count[0] = rgg->nlat;
        rgg->gw = CALLOC( rgg->nlat, double);
        ncvarget( ncid, varid, start, count, rgg->gw);

	return rgg;
    }

RG_t *RGIread( char *filename)
    {
	extern int ncopts;

        int ncid;
        int isect_dimid;
        int fromidx_varid, pct_varid, start_varid;

        long start[2], count[2], size;
        RG_t *rgn;

	ncopts = 0;
        ncid = ncopen( filename, NC_NOWRITE);
	if (ncerr < 0)
	    return NULL;
	ncopts = NC_VERBOSE | NC_FATAL;

        rgn = MALLOC( RG_t);
        rgn->from = ReadGrid( ncid, "f");
        rgn->to = ReadGrid( ncid, "t");

        isect_dimid = ncdimid( ncid, "isect");
        ncdiminq( ncid, isect_dimid, NULL, &size);
        rgn->n_isects = size;

        fromidx_varid = ncvarid( ncid, "from_idx");
        pct_varid = ncvarid( ncid, "pct");
        start_varid = ncvarid( ncid, "start_isect");

        start[0] = 0;
        count[0] = rgn->n_isects;

        rgn->fromidx = CALLOC( rgn->n_isects, int);
        ncvarget( ncid, fromidx_varid, start, count, rgn->fromidx);

        rgn->pct = CALLOC( rgn->n_isects, float);
        ncvarget( ncid, pct_varid, start, count, rgn->pct);

        start[0] = 0;
        count[0] = rgn->to->nlat;
        start[1] = 0;
        count[1] = rgn->to->nlon;
        rgn->start_isect = CALLOC( rgn->to->nlat * rgn->to->nlon, int);
        ncvarget( ncid, start_varid, start, count, rgn->start_isect);

        ncclose( ncid);
        printf("read regrid file %s\n", filename);
        return rgn;
    }

/* this does the interpolation. It assumes the src and dest is in
   canonical order (x,y,z), with standard strides. also assumes
   that the data are both float.  See RGconvertArr() to relax those
   assumptions
 */
int RGInterp( RG_t *rg, REAL *from_data, int from_nx, int from_ny,
	REAL *to_data, int to_nx, int to_ny, int nz)
    {
	int i, j;

	int from_sz = from_nx * from_ny;

	int to_sz = to_nx * to_ny;

	if (from_nx != rg->from->nlon)
	    {
	    printf("RGInterp bad from nlons %d != %d\n", from_nx, rg->from->nlon);
	    return 0;
	    }
	if (from_ny != rg->from->nlat)
	    {
	    printf("RGInterp bad from nlats %d != %d\n", from_ny, rg->from->nlat);
	    return 0;
	    }
	if (to_nx != rg->to->nlon)
	    {
	    printf("RGInterp bad to nlons %d != %d\n", to_nx, rg->to->nlon);
	    return 0;
	    }
	if (to_ny != rg->to->nlat)
	    {
	    printf("RGInterp bad to nlats %d != %d\n", to_ny, rg->to->nlat);
	    return 0;
	    }

        /* loop over the "to" points */
        for (i=0; i< to_nx*to_ny; i++)
            {
            int start = rg->start_isect[i];
            int end = (i == to_nx*to_ny-1) ? rg->n_isects : rg->start_isect[i+1];
            int npts = end - start;
            int idx[6];
            REAL pct[6];
        
            memset( idx, 0, sizeof(idx));
            memset( pct, 0, sizeof(pct));
            for (j=0; j<npts; j++)
                {       
                pct[j] = *(rg->pct + start + j);
                idx[j] = *(rg->fromidx + start + j) + 1;
                }

	    /* this allows us to vectorize */
            AREAAVG( &nz, from_data, &from_sz, to_data, &to_sz, pct, idx);
	    to_data++;
            }

#ifdef OLD_WAY
        /* loop over the "to" points */
        for (i=0; i< to_nx*to_ny; i++)
            {
	    int x, y, z;
	    int from_sx = 1;
	    int from_sy = from_nx;
	    int to_sx = 1;
	    int to_sy = to_nx;
            int start = rg->start_isect[i];
            int end = (i == to_nx*to_ny-1) ? rg->n_isects : rg->start_isect[i+1];

	    /* loop over levels */
	    for (z=0; z<nz; z++)
	      {
	      float val, sum = 0.0;

	      /* weighted sum of "from" points */
	      for (j=start; j<end; j++)
		{
		x = rg->fromidx[j] % from_nx;
		y = rg->fromidx[j] / from_nx;
		
		val = *(from + x * from_sx + y * from_sy + z * from_sz);
                sum += rg->pct[j] * val;
		}
	    
	      x = i % to_nx;
	      y = i / to_nx;

	      *(to + x * to_sx + y * to_sy + z * to_sz) = sum;

	      /* if (0 != malloc_check(-1))
		printf("corrupted at %d %d %d %d\n", i, x, y, z); */
	      } 
	    } 
#endif

	return 1;
    }

void RGIfree( RG_t *rg)
    {
        FREE( rg->fromidx);
        FREE( rg->pct);
        FREE( rg->start_isect);
        FREE(rg);
    }

