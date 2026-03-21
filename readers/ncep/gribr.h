/* $Id$ */
#ifndef GRIBR_WAS_INCLUDED
#define GRIBR_WAS_INCLUDED

#include <time.h>


#if (REALSIZE == 4)
#define REAL float
#else
#define REAL double
#endif

/********** grib special projections  */
typedef struct GRIBRgds_latlon_t
    {
    unsigned char length[3];
    unsigned char nv;           /* num vert coord params */
    unsigned char pv;           /* location of vert coord or num points */
    unsigned char type;         /* projection type */

    unsigned char num_lon_points[2];
    unsigned char num_lat_points[2];
    unsigned char la1[3];
    unsigned char lo1[3];
    unsigned char res_flags;
    unsigned char la2[3];
    unsigned char lo2[3];
    unsigned char di[2];
    unsigned char dj[2];
    unsigned char scan_mode;
    unsigned char reserved[4];
    } GRIBRgds_latlon_t;

/* binary data section */
typedef struct GRIBRbds_t
{
    unsigned char length[3];
    unsigned char flag;
    unsigned char binary_scale_factor[2];
    unsigned char reference_value[4];
    unsigned char number_of_bits;
    unsigned char *data;
} GRIBRbds_t;

typedef struct GRIBRbms_t
    {
    unsigned char length[3];
    unsigned char unused_bits;
    unsigned char numeric_id[2];
    unsigned char *bit_map;
} GRIBRbms_t;

/* product descriptor section */
typedef struct GRIBRpds_t
{
    unsigned char length[3];
    unsigned char parameter_table_ver;
    unsigned char center_id;
    unsigned char process_id;
    unsigned char grid_id;
    unsigned char gds_bms_flag;
    unsigned char parameter;
    unsigned char level_type;
    unsigned char level_params[2];
    unsigned char year;
    unsigned char month;
    unsigned char day;
    unsigned char hour;
    unsigned char minute;
    unsigned char time_unit;
    unsigned char time_p1;
    unsigned char time_p2;
    unsigned char time_range;
    unsigned char number_in_average[2];
    unsigned char number_missing;
    unsigned char century;
    unsigned char reserved0;
    unsigned char scale_factor[2];
} GRIBRpds_t;

typedef struct GRIBR_t
{
    int total_len;
    long file_offset;

    GRIBRpds_t pds;
    GRIBRgds_latlon_t gds;
    GRIBRbms_t bms;
    GRIBRbds_t bds;

} GRIBR_t;

typedef struct Param_t {
    char *long_name;
    char *unit;
    char *short_name;
} Param_t;

typedef struct Time_t {
    time_t utime;
    long offset;
} Time_t;

typedef struct Level_t {
    int id;
    int param;
    int ntim_alloc;
    int ntimes;
    Time_t *time;
} Level_t;

typedef struct Field_t {
    char *name;
    char *long_name;
    char *units;

    int nlev_alloc;
    int nlevels;
    Level_t *level;
} Field_t;

typedef struct GRIBRhead_t
{
    int nx, ny, nz;

    double *lat;
    double *lat_edge;
    double *gaussw;

    double *lev;
    double *lev_edge;

    int   nflds;
    Field_t *fld;

	/* private */
    Param_t *table2;

    FILE *fp;
} GRIBRhead_t;

extern GRIBRhead_t *GRIBRopen(char *filename, int force, int print);
/* open a grib file, and constuct index of all products if it doesnt
   exist. force = 1: always construct index. print=1: print products
   to stdout. return NULL if fails.
 */ 

extern void GRIBRclose( GRIBRhead_t *head);

extern void GRIBRprint( FILE *fp, GRIBRhead_t *gh);

extern int GRIBRread_data( GRIBRhead_t *gh, int field_idx, int level_idx, 
	int time_idx, REAL *outdata);
/* read the product data defined by field_idx, level_idx, time_idx. 
   level_idx = -1 means read all levels. 
   Put into outdata in natural data order (x, y, then z varies fastest). 
   Caller must ensure that enough space is allocated.

   Note these are indices, ie:
	gh->fld[ field_idx].lev[ lev_idx].tim[ time_idx]
 */


#endif
