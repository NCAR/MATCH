/* $Id$ */
#ifndef RGI_WAS_INCLUDED
#define RGI_WAS_INCLUDED

#ifdef FRIEND_OF_RGI
typedef struct RGgrid_t {
    int nlon;
    double *lon;        /* lon[nlon] */
    double *lone;       /* lone[nlon+1] */

    int nlat;
    double *lat;        /* lat[nlat] */
    double *late;       /* late[nlat+1] */
    double *gw;         /* gw[nlat] */
} RGgrid_t;

typedef struct RG_t {
    RGgrid_t    *from;
    RGgrid_t    *to;

        /* intersections */
    int         malloc_isects;  /* number of malloced intersections */
    int         n_isects;
    int         *fromidx;       /* fromidx[n_isects] */
    float       *pct;           /* pct[n_isects] */

    int         *start_isect;   /* start_isect[to->nlat, to->nlon] : index into intersection */
} RG_t;
#else
typedef void *RG_t;
#endif

#if (REALSIZE == 4)
#define REAL float
#else
#define REAL double
#endif

extern RG_t *RGIread( char *filename);
extern int RGIwrite( RG_t *rg, char *filename);

extern int RGInterp( RG_t *rg, REAL *from, int from_nx, int from_ny,
        REAL *to, int to_nx, int to_ny, int nz);
extern void RGIfree( RG_t *rg);

#endif
