/* $Id$ */
/*
   grib.c

  4/11/97 : add level
  4/14/97 : add regridding
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

#include <ctype.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "gribr.h"
#include "mss.h"
#include "str_util.h"
#include "utim.h"

#define FRIEND_OF_RGI
#include "rgi.h"

static int Debug = 0;
static int DebugM = 0;
#define DPRINT if (Debug) fprintf(stderr,
#define DPRINTM if (DebugM) fprintf(stderr,

#if defined( SGI  )
#define TINTERP tinterp_
#define AREAAVG areaavg_
#elif defined( LINUX  )
#define TINTERP tinterp_
#define AREAAVG areaavg_
#elif defined( AIX )
#define TINTERP tinterp
#define AREAAVG areaavg
#endif

#if (REALSIZE == 4)
#define REAL float
#else
#define REAL double
#endif

extern void TINTERP(REAL*, REAL*, REAL*, REAL*, int*);
extern void AREAAVG(int*, REAL*, int*, REAL*, int*, REAL*, int*);

static unsigned long Reader_cl = 0, Regrid_cl = 0, Grib_cl = 0, Interp_cl = 0;
static UTIMdiff_t *Total_ud, *Reader_ud, *Regrid_ud, *Grib_ud, *Interp_ud;

#define DEG2RAD 0.017453293       /* radians per degree */

/**************** local structures ******************/
typedef struct DataVol {
    int	    timeidx;
    time_t  utime;
    REAL  *data;
    REAL  *regrid_data;
    int  *data_guard;	/* debugging */
} DataVol;

typedef struct Cfield_t {
    char    *name;
    int     grib_id;
    int     is_3d;

    int     last_used;
    DataVol vol[2];
} Cfield_t;

typedef struct Ctime {
   time_t   utime;	/* secs since 1970 */

   MSSfile_t *file;	/* grib file */	
   GRIBRhead_t 	*head;  /* grib header */
   int	    gindex;	/* time index into grib header */
} Ctime;

#define MAX_FILES 4
typedef struct Dataset {
   MSS_t	*mss;		/* mss files */
   int		mss_current;	/* current mss file index */
   
   GRIBRhead_t 	*head[MAX_FILES];  	/* list of active grib headers */
   MSSfile_t 	*file[MAX_FILES];  	/* their mss file */
   int		last_timeidx[MAX_FILES];/* their last time index in time array */

   int      	ntimes;     	/* number of time steps */
   int      	ntimes_malloc; 	/* number of time steps malloced */
   Ctime    	*time;      	/* time steps */

   int		nflds;
   Cfield_t	*fld;

   int		nlon,  nlat,  nlev;	/* result size: may not match grib if regridding */

   int          remove;         /* remove old files */
} Dataset;

static Dataset *DS2, *DS3;	/* 2d, 3d datasets */
static int Regrid = 0;
static RG_t *Rgi = NULL;

static int malloc_check( int dumm)
{  return 0; }

#define MAX_TOKEN 40

static char *gmtime_str(time_t *t) 
    {
	static char buff1[100], buff2[100];
	static int idx = 0;
	char *b = idx ? buff1 : buff2;
	idx = (idx +1) % 2;

	strcpy(b, asctime( gmtime( t)));
	b[strlen(b)-1] = 0;
	return b;
    }

/* on Cray, the char length doesnt work, so we just use
   the first blank as terminator (!!!) */
static char *f2c_blank(char *src)
    {
        static char buff1[400], buff2[400];
        static int use = 0;
        char *dst, *buff = (use % 2) ? buff1 : buff2;
        use++;

	dst = buff;
        while (*src && (*src != ' '))
            *dst++ = *src++;
        *dst = 0;
        return buff;
    }

/* convert integer YYMMDD secs -> time_t */
static time_t ConvertTime2( int date, int secs)
    {
        int y,  m,  d,  w;

        w = date;
        DPRINT "ConvertTime2: date is = %d \n",date);

        if (date < 1000000) 
         {
           y = (w / 10000) + 1900;
           w += 19000000;
           DPRINT "ConvertTime2: passed 2 digit year, y = %d\n",y);
         }     
        else                  
         {     
           y = w / 10000;
           DPRINT "ConvertTime2: passed 4 digit year, y = %d\n",y);
         }

        w -= y * 10000;
        m = w / 100;
        d = w % 100;

        DPRINT "y in ConvertTime2 = %d \n",y);
        DPRINT "d in ConvertTime2 = %d \n",d);
        DPRINT "m in ConvertTime2 = %d \n",m);

       return UTIMdate2unix( y,  m,  d,  0,  0,  secs);   
    }



/*----------------------------------------------------------------------*/
#if defined( SGI ) 
void modify_date_datesec_
#elif defined( LINUX  )
void modify_date_datesec_
#elif defined( AIX )
void modify_date_datesec
#elif defined( CRAY )
void MODIFY_DATE_DATESEC
#endif
     (int *date, int *secs, int *offset, int *offset_date,
      int *offset_secs)

/*   *offset is the amount to shift the day by in seconds       
 *   *offset_date is the new offsetted date
 *   *offset_sec is the new offsetted date in seconds
 *
 *   This function handles both a 2 digit year and a 
 *   four digit year (as of 20010817, jdf).
 *   For compatibility with later scripts, the same number 
 *   of digits is returned as were passed to the script. */

   {

     long int  orig_time;
     long int offset_time;
     long int ref_day,temp_day;
     long int offset_year,offset_month,offset_day;
     int year_digits;

     DPRINT "\nmodify_date_datesec: date passed =  %d \n",*date);
         if (*date < 1000000) 
         {
	   year_digits = 2;
         }     
        else                  
         {     
	   year_digits = 4;
         }
     DPRINT "passed %d digit year\n",year_digits);


       /* if sec - offset < 0 then subtract day */
       if (*offset > *secs)
       { 
         /*orig_time is seconds since 1970/1/1 UTC until date specified*/
	 orig_time = ConvertTime2( *date, *secs);

	 DPRINT "orig_time (long) after ConvertTime2: %ld \n",orig_time); 
	 DPRINT "\nmodify_date_datesec: Original time: \n");
	 DPRINT "date, secs: %d, %d \n",*date,*secs);
	 DPRINT "long : %ld \n", orig_time);

	 DPRINT "\nmodify_date_datesec: Offset: \n");
	 DPRINT "offset (in secs): %d  \n",*offset);

	 offset_time = orig_time - *offset;
	 DPRINT   "long modified: %li \n", offset_time);

	 ref_day = UTIMjulian_date(1,1,1970);
	 temp_day = (offset_time / 86400);

	 UTIMcalendar_date((ref_day + temp_day), &offset_day, &offset_month, &offset_year);
	 DPRINT   "\nmodify_date_datesec: Modifed: \n");
	 DPRINT   "day, month, year: " );
	 DPRINT   "%ld %ld %ld \n", offset_day, offset_month, offset_year);
	 DPRINT   "long modified: %li \n", offset_time);
	 
	 *offset_secs = offset_time % 86400;

         /*Now we need to give back a 4 digit year if modify_date_datesec (this function) was given a 4 digit year.  Likewise, return a 2 digit year if passed a 2 digit year
	  */
         if (year_digits == 2) 
         {
         offset_year = offset_year - 1900;
	 *offset_date = 10000*offset_year + 100*offset_month + offset_day;
         }     
        else                  
         {     
	 *offset_date = 10000*offset_year + 100*offset_month + offset_day;
         }
         DPRINT   "Used previous day \n");
	 DPRINT   "Offset_date %d \n",*offset_date);
	 DPRINT   "Offset_secs %d \n",*offset_secs);
       }
     else
       {
	 *offset_secs = *secs - *offset;
	 *offset_date = *date;
	 DPRINT "  didn't use previous day \n");
	 DPRINT "  Offset_date %d \n",*offset_date);
	 DPRINT "  Offset_secs %d \n",*offset_secs);

       }
   }


/*----------------------------------------------------------------------*/



/* check for old files and delete */
static int CheckOld( Dataset *ds, int timidx)
    {
	int i;
	static int count = 0;
	if (++count < 100)
	    return 1;
	count = 0;

        for (i=0; i<MAX_FILES; i++)
            if ((ds->head[i] != NULL) && (ds->last_timeidx[i] < timidx))
                {
		if ( ds->remove ) MSSdelete( ds->mss, ds->file[i]);
		GRIBRclose( ds->head[i]);

		MSSget_next(ds->mss);

                ds->file[i] = NULL;
                ds->head[i] = NULL;
                ds->last_timeidx[i] = 0;
                }

        return 0;
    }

/* opens the file in ds->mss given by index */
static int OpenOneFile( Dataset *ds, int index)
    {
	int i, idx;
	char pathname[500];
	MSSfile_t *file;
	GRIBRhead_t *head;

	if (index < 0)
	    return -1;
	if (index >= ds->mss->nfiles)
	    {
	    fprintf(stderr,"OpenOneFile: file index %d beyond end of file list (%d files)\n",
		index, ds->mss->nfiles);
	    return 0;
	    }
	file = &ds->mss->file[index];

        /* open the file */
	if (0 > MSSget( ds->mss, file->name, 1))
	    return 0;

	sprintf( pathname, "%s/%s", ds->mss->lpath, file->name);
	if (NULL == (head = GRIBRopen(pathname, 0, 0)))
	    return 0;
	DPRINT "    OpenOneFile: file open ok on <%s>\n", pathname);

        /* assume all fields have the same times; find first non-empty one */
	for (idx=0; idx< head->nflds; idx++)
	    if (head->fld[idx].nlevels > 0)
		break;

        {
        Level_t *lev = &head->fld[idx].level[0];

        while (ds->ntimes + lev->ntimes >= ds->ntimes_malloc)
            {
            ds->ntimes_malloc *= 2;
            ds->time = REALLOC( ds->time, ds->ntimes_malloc, Ctime);
            }

        for (i=0; i < lev->ntimes; i++)
            {
	    Time_t *tim = &lev->time[i];
	    Ctime *time = &ds->time[ds->ntimes];
	    ds->ntimes++;

	    time->utime = tim->utime;
	    time->file = file;
	    time->head = head;
	    time->gindex = i;

	    DPRINT "   add time %ld %s\n", (long)tim->utime,  gmtime_str( &time->utime));
	    }
	}

	/* find a spot for it in the ds */
	for (i=0; i<MAX_FILES; i++)
	    if (ds->head[i] == NULL)
		{
		ds->head[i] = head;
		ds->file[i] = file;
		ds->last_timeidx[i] = ds->ntimes-1;
		return 1;
		}

	fprintf(stderr, "ERROR: no room in ds->head array!!!\n");

	return 0;
    }

/* find the file whose starting time is the largest value less than or equal to
   the want time.   */
static int FindFileByTime( Dataset *ds, time_t want_time)
    {
        int idx;

        /* LOOK: wont detect going off the end, always end at n-1 */
        idx = 0;
        while ((idx < ds->mss->nfiles-1) && (ds->mss->file[idx+1].utime <= want_time))
            idx++;

        if (idx >= ds->mss->nfiles)
            {
            fprintf(stderr,"GRIB FindFilebyTime: cant find file for time = %s\n", 
		gmtime_str( &want_time));
            return -1;
            }

	if (Debug)
           fprintf(stderr,"GRIB FindFilebyTime found file %s for time %s\n", 
		ds->mss->file[idx].name, gmtime_str( &want_time));

	return idx;
}


Dataset *DatasetInit( char *lpath, char *rpath, char *list, int start_time, 
	int nlon, int nlat, int nlev, int remove)
    {
	Dataset *ds = MALLOC (Dataset);
	int i, count, err = 0;
	GRIBRhead_t *head;

	ds->nlat = nlat;
	ds->nlon = nlon;
	ds->nlev = nlev;
	ds->remove = remove;
	    	
        ds->ntimes_malloc = 20;
        ds->time = CALLOC( ds->ntimes_malloc, Ctime);

	/* initialize the mss list */
	if (NULL == (ds->mss = MSSinit( rpath, lpath, mss_explicit, list)))
	    return NULL;

	/* open the first file, as defined by start_time */
	ds->mss_current = FindFileByTime( ds, start_time);
	if (0 == OpenOneFile(ds, ds->mss_current))
	    return NULL;

       	/* set up the field array. This assumes all files have the same fields. */
        count = 0;
	head = ds->head[0];
        for (i=0; i<head[0].nflds; i++)
            if (head->fld[i].nlevels > 0)
                count++;
        ds->nflds = count;
        if (ds->nflds > 0)
            ds->fld = CALLOC(ds->nflds, Cfield_t);


        count = 0;
        for (i=0; i<head->nflds; i++)
            if (head->fld[i].nlevels > 0)
                {
                ds->fld[count].name = head->fld[i].name;
                ds->fld[count].grib_id = i;
                ds->fld[count].is_3d = (head->fld[i].nlevels == head->nz);
                ds->fld[count].vol[0].timeidx = -1;
                ds->fld[count].vol[1].timeidx = -1;
                count++;
                }

        /* check that dimensions are compatible */
        err = 0;


	if (Regrid)
	  {
          if (ds->nlon != Rgi->to->nlon)
            {
            printf("GRIB dataset regrid has %d nlon, wanted %d\n", Rgi->to->nlon, ds->nlon);
            err = 1;
            }
          if (ds->nlat != Rgi->to->nlat)
            {
            printf("GRIB dataset regrid has %d nlat, wanted %d\n", Rgi->to->nlat, ds->nlat);
            err = 1;
            }
	  }
	else
	  {
          if (ds->nlon != head->nx)
            {
            printf("GRIB dataset has %d nlon, wanted %d\n", head->nx,ds->nlon);
            err = 1;
            }
          if (ds->nlat != head->ny)
            {
            printf("GRIB dataset has %d nlat, wanted %d\n", head->ny,ds->nlat);
            err = 1;
            }
	  }

        if ((ds->nlev != head->nz) && (ds->nlev != 0))
            {
            printf("GRIB datset has %d nlev, wanted %d\n", head->nz,ds->nlev);
            err = 1;
            }

	return err ? NULL : ds;
    }

#if defined( SGI  )
int grib_open_
#elif defined( LINUX  )
int grib_open_
#elif defined( AIX )
int grib_open
#elif defined( CRAY )
int GRIB_OPEN
#endif
	( char *lpth, char *rpth, int *date, int *datesec, 
	   int *plon, int *plat, int *plev, int *remove,
	   REAL *phi, REAL *latwts,  REAL *hyam, REAL *hybm, REAL *hyai, 
	   REAL *hybi, int *regrid,
	   int lpth_len, int rpth_len )
    {
	int y, z;
	time_t start_time;
	char *lpath, *rpath;
	GRIBRhead_t *head;

	/* timing info */
	Total_ud = MALLOC( UTIMdiff_t);
	Reader_ud = MALLOC( UTIMdiff_t);
	Regrid_ud = MALLOC( UTIMdiff_t);
	Grib_ud = MALLOC( UTIMdiff_t);
	Interp_ud = MALLOC( UTIMdiff_t);
	UTIMdiff(Total_ud);

	lpath = strdup(f2c_blank(lpth));
	rpath = strdup(f2c_blank(rpth));
	Regrid = *regrid;

	fprintf(stderr,"Open NMC GRIB files in %s %s\n", lpath, rpath);

	if (Regrid)
	    Rgi = RGIread("T62n.T42n.avg.nc");

	start_time = ConvertTime2( *date, *datesec);

	fprintf(stderr,"grib_open, grb2d modify %d %d %d \n", start_time,*date,*datesec);


	if (NULL == (DS2 = DatasetInit(lpath, rpath, "grb2d.list", start_time, *plon, 
		*plat, 0, *remove)))
	    return 0;

	if (NULL==(DS3=DatasetInit(lpath, rpath, "grbsanl.list", start_time, *plon, 
		*plat, *plev, *remove)))
	    return 0;

	/* match wants levels from top to bottom.
	   it wants lats ordered south to north

	   in order to make the regridding transparent, when we regrid, we
	   preserve the original data order (north to south, top to bottom),
	   and put it into match order when we do the interpolation  
  	 */

	head = DS3->head[0];
	if (Regrid)
	  {
	  RGgrid_t *to = Rgi->to;
          for (y=0; y<to->nlat; y++)
            {
            phi[y] = DEG2RAD * to->lat[DS3->nlat - 1 - y];
            latwts[y] = to->gw[DS3->nlat - 1 - y];
            }
	  }
	else
	  {
	  for (y=0; y<DS3->nlat; y++)
	    {
	    phi[y] = DEG2RAD * head->lat[DS3->nlat - 1 - y];
	    latwts[y] = head->gaussw[DS3->nlat - 1 - y];
	    fprintf(stderr,"grib_open phi %d %g \n", y, phi[y]);
	    }
	  }

	for (z=0; z<DS3->nlev+1; z++)
	    {
	    hybi[z] = (REAL) head->lev_edge[DS3->nlev - z];
	    hyai[z] = 0.0;
	    fprintf(stderr,"grib_open hybi %d %g \n", z, hybi[z]);
	    }
	for (z=0; z<DS3->nlev; z++)
	    {
	    hybm[z] = (REAL) head->lev[DS3->nlev - 1 - z];
	    fprintf(stderr,"grib_open hybm %d %g \n", z, hybm[z]);
	    hyam[z] = 0.0;
	    }
	    
	return 1;
    }


#if defined( SGI  )
void grib_time_ ( void)
#elif defined( LINUX  )
void grib_time_ ( void)
#elif defined( AIX )
void grib_time ( void)
#elif defined( CRAY )
void GRIB_TIME (void)
#endif
    {
	printf("TIME (millisecs) total %u reader %lu regrid %lu interp %lu gribr %lu\n", 
		UTIMdiff( Total_ud), Reader_cl, Regrid_cl, Interp_cl, Grib_cl);
    }

static 	DataVol *GetData( Dataset *ds, int fldidx, int level, int timidx)
    {
	Cfield_t *fld = &ds->fld[fldidx];
	DataVol *vol;
	Ctime *time = &ds->time[timidx];
	int i, ret, nz, size;
	
	/* see if we already got it LOOK: 2.5d fields !!!! */
	for (i=0; i< 2; i++)
	    if (fld->vol[i].timeidx == timidx)
	        {
		fld->last_used = i;
		return &fld->vol[i];
		}
		
	/* nope - keep track of last used */
	fld->last_used = (fld->last_used + 1) % 2;
	vol = &fld->vol[ fld->last_used];
	vol->timeidx = timidx;
	vol->utime = time->utime;
		
	if (vol->data == NULL)
	    {
	    nz = fld->is_3d ? time->head->nz : 1;
	    size = nz * time->head->nx * time->head->ny;

	    if (Debug)
	        fprintf(stderr,"  malloc %d for fld %s\n", size, fld->name);

	    DPRINTM "calloc data = %d\n", malloc_check(-1));
	    vol->data = CALLOC(size + 1, REAL);
	    vol->data_guard = (int *) vol->data + size; /* barf */
	    *vol->data_guard = 12345678;
	    }
	    
	Grib_ud->init = 0;
        UTIMdiff( Grib_ud);
	if (Debug)
	    printf("calling GRIBR for %s\n", fld->name);
	ret = GRIBRread_data(time->head, fld->grib_id, level, time->gindex, vol->data);
	Grib_cl += UTIMdiff( Grib_ud);

       	if (Regrid)
          {
	  Regrid_ud->init = 0;
	  UTIMdiff( Regrid_ud);

          nz = fld->is_3d ? time->head->nz : 1;
	  if (vol->regrid_data == NULL)
            {
            size = nz * Rgi->to->nlon * Rgi->to->nlat;
	    DPRINTM "calloc regrid data = %d\n", malloc_check(-1));
            vol->regrid_data = CALLOC(size, REAL);
            if (Debug)
               fprintf(stderr,"  malloc %d for fld %s regrid\n", size, fld->name);
            }

            DPRINTM "before interp data = %d\n", malloc_check(-1));
          ret &= RGInterp( Rgi, vol->data, Rgi->from->nlon, Rgi->from->nlat,
		vol->regrid_data, Rgi->to->nlon, Rgi->to->nlat, nz);
            DPRINTM "interp data = %d\n", malloc_check(-1));

	  Regrid_cl += UTIMdiff( Regrid_ud);
          }

	return (ret == 1) ? vol : NULL;
    }
    
#if defined( SGI  )
int grib_getdata_
#elif defined( LINUX )
int grib_getdata_
#elif defined( AIX )
int grib_getdata
#elif defined( CRAY )
int GRIB_GETDATA
#endif
   ( char *fld, int *is_2d, int *level, int *date, int *datesec, 
     REAL *result, int *stridex, int *stridey, int *stridez, 
     int *last_timidx, int fld_len)

   {
	time_t want_time = ConvertTime2( *date,  *datesec);
	int fldidx,  timidx;
	char fldname[300];
	Dataset *ds = (*is_2d == 1) ? DS2 : DS3;
	DataVol *this;

	Reader_ud->init = 0;
        UTIMdiff( Reader_ud);

        strcpy(fldname, f2c_blank(fld));

	DPRINTM "before GetData = %d\n", malloc_check(-1));

        DPRINT "GRIB_GETDATA %d %d %d %s %g %d %d %d %d %d\n", *date, *datesec, *is_2d, fldname, *result,
            *stridex, *stridey, *stridez, *last_timidx, fld_len);
        DPRINT "   want time %ld %s\n", (long)want_time,  gmtime_str( &want_time));

	for (fldidx=0; fldidx<ds->nflds; fldidx++)
	    if (0 == strcasecmp( fldname, ds->fld[fldidx].name))
		goto gotit;
		
	fprintf(stderr,"ERROR grib_getdata cant find fld = <%s> want time = %ld\n", fldname,
		(long)want_time);
	return 0;
	
    gotit:
	/* find the bracketing time indices */
	if (want_time < ds->time[0].utime)
	    {
	    fprintf(stderr,"WARNING grib_getdata: time %d %d == %s before first data time = %s, clamping\n",
		*date, *datesec,  gmtime_str( &want_time), gmtime_str( &ds->time[0].utime));
	    want_time = ds->time[0].utime;
	    }
	 if (want_time >= ds->time[ds->ntimes-1].utime)
	    {
	    /* here we force it to get the next file in the sequence */
	    int index = ds->mss_current + 1;
	    if ((index >= ds->mss->nfiles) || (1 != OpenOneFile( ds, index)))
	    	{
	    	fprintf(stderr,"WARNING grib_getdata: time %d %d == %s after last data time = %s, clamping\n",
			*date, *datesec,  gmtime_str( &want_time), gmtime_str( &ds->time[ds->ntimes-1].utime));
		want_time = ds->time[ds->ntimes-1].utime;
		}
	    else
		{
		ds->mss_current = index;
		}
	    }
	    
	/* find the largest timidx thats less than or equal to want_time */
	timidx = 0;
        while ((timidx < ds->ntimes-1) && (want_time >= ds->time[timidx+1].utime))
	    timidx++;

	/* if its the same as last time, do nothing */
	if (timidx == (*last_timidx))
	    return 2;
		
        DPRINT "  grib_get_data fld %d = %s level %d time %s, use %d = %s\n", fldidx, fldname,
             *level, gmtime_str( &want_time), timidx, gmtime_str( &ds->time[timidx].utime));

        /* Fetch the data */
        if (NULL == (this = GetData( ds, fldidx, *level, timidx)))
            return 0;
	DPRINTM "check malloc after GetData = %d\n", malloc_check(-1));

	/* reorder to what match wants */
	{
	int x, y, z;
	int rxs = *stridex,  rys = *stridey,  rzs = *stridez;
	int nx = ds->nlon, ny = ds->nlat, nz = (rzs == 0) ? 1 : ds->nlev;
	int dxs = 1,  dys = nx,  dzs = (rzs == 0) ? 0 : nx * ny;
	REAL *data = Regrid ? this->regrid_data : this->data;

	/* these indices are the "result" grid */
	for (z=0; z < nz; z++)
	 for (y=0; y < ny; y++)
	  for (x=0; x < nx; x++)
	    {
	    int result_elem = x * rxs + (ny-1-y) * rys + (nz-1-z) * rzs;
	    int data_elem = x * dxs + y * dys + z * dzs;
	  
	    *(result + result_elem) = (REAL) *(data + data_elem);
	    }
	}
	
	DPRINTM "check malloc after grib_getdata = %d\n", malloc_check(-1));
        /* delete old files */
        CheckOld( ds, timidx);

	*last_timidx = timidx;
	Reader_cl += UTIMdiff( Reader_ud);
	return 1;
    }

#if defined( SGI  )
int grib_interp_
#elif defined( LINUX )
int grib_interp_
#elif defined( AIX )
int grib_interp
#elif defined( CRAY )
int GRIB_INTERP
#endif
    ( char *fld, int *is_2d, int *level, int *date, int *datesec, 
      REAL *result, int *stridex, int *stridey, int *stridez, int fld_len)

    {
	time_t want_time = ConvertTime2( *date,  *datesec);
	int fldidx,  timidx;
	char fldname[300];
        Dataset *ds = (*is_2d == 1) ? DS2 : DS3;
	
 	Reader_ud->init = 0;
        UTIMdiff( Reader_ud);

        strcpy(fldname, f2c_blank(fld));


	for (fldidx=0; fldidx<ds->nflds; fldidx++)
	    if (0 == strcasecmp( fldname, ds->fld[fldidx].name))
		goto gotit;
		
	fprintf(stderr,"ERROR grib_interp cant find fld = <%s> (%s %d) want time = %ld\n", 
		fldname, fld, fld_len, (long)want_time);
	return 0;
	
    gotit:
	/* find the bracketing time indices */
	if (want_time < ds->time[0].utime)
	    {
	    /* Clamp to first available time step instead of failing.
	       This allows running without the previous month's data;
	       the first few hours use uninterpolated (held) values. */
	    fprintf(stderr,"WARNING grib_interp: time %d %d == %s before first data time = %s, clamping\n",
		*date, *datesec,  gmtime_str( &want_time), gmtime_str( &ds->time[0].utime));
	    want_time = ds->time[0].utime;
	    }
         if (want_time >= ds->time[ds->ntimes-1].utime)
            {
            int index = ds->mss_current + 1;
            if ((index >= ds->mss->nfiles) || (1 != OpenOneFile( ds, index)))
                {
                /* Clamp to last available time step instead of failing.
                   This allows running without the next month's data;
                   the final few hours use uninterpolated (held) values. */
                fprintf(stderr,"WARNING grib_interp: time %d %d == %s after last data time = %s, clamping\n",
                        *date, *datesec,  gmtime_str( &want_time), gmtime_str( &ds->time[ds->ntimes-1].utime));
                want_time = ds->time[ds->ntimes-1].utime;
		}
	    else
		{
		ds->mss_current = index;
		}
            }
	    
        /* find the largest timidx thats less than or equal to want_time */
	timidx = 0;
        while ((timidx < ds->ntimes-1) && (want_time >= ds->time[timidx+1].utime))
	    timidx++;

        DPRINT "  grib_interp fld %d = %s level %d time %s, use %d = %s\n", fldidx, fldname,
	     *level, gmtime_str( &want_time), timidx, gmtime_str( &ds->time[timidx].utime));
	
	{
        int x, y, z;
        int rxs = *stridex,  rys = *stridey,  rzs = *stridez;
        int nx = ds->nlon, ny = ds->nlat, nz = (rzs == 0) ? 1 : ds->nlev;
        int dxs = 1,  dys = nx,  dzs = (rzs == 0) ? 0 : nx * ny;
	REAL dt, dtm1;
	REAL *adata, *bdata;

	/* Fetch the data; clamp after-index to avoid out-of-bounds
	   when want_time was clamped to the last available time step */
	int timidx_after = (timidx+1 < ds->ntimes) ? timidx+1 : timidx;
	DataVol *before = GetData( ds, fldidx, *level, timidx);
	DataVol *after = GetData( ds, fldidx, *level, timidx_after);
	if ((before == NULL) || (after == NULL))
	    return 0;
	    
        bdata = Regrid ? before->regrid_data : before->data;
        adata = Regrid ? after->regrid_data : after->data;

	/* interpolate (guard divide-by-zero when clamped to single time step) */
	if (after->utime == before->utime)
	    { dt = 1.0; dtm1 = 0.0; }
	else
	    { dt = (REAL)(after->utime-want_time)/(REAL)(after->utime - before->utime);
	      dtm1 = 1.0 - dt; }
	
	DPRINT "    %ld = %g * %ld + %g * %ld\n",
		(long)want_time, dt, (long)before->utime, dtm1, (long)after->utime);
 
        Interp_ud->init = 0;
        UTIMdiff( Interp_ud);

	if (rxs == 1) 
	  {
	  /* can vectorize over x dimension */
          for (z=0; z < nz; z++)
           for (y=0; y < ds->nlat; y++)
             {
             int result_elem = (ny-1-y) * rys + (nz-1-z) * rzs;
             int data_elem = y * dys + z * dzs;

	     TINTERP( &dt, bdata + data_elem, adata + data_elem, result+result_elem,
		 &ds->nlon);

                if ((y==0) && (z==0))
                     DPRINT "grib_interp, check result:  %g  = bdata %g * %g + adata %g * %g\n", 
		       *(result + result_elem),
		       *(bdata + data_elem), dt, *(adata + data_elem), dtm1);
	      }
	  }
	else
	  {
	  /* these indices are the "result" grid */
	  for (z=0; z < nz; z++)
	   for (y=0; y < ds->nlat; y++)
	    for (x=0; x < ds->nlon; x++)
	     {
             int result_elem = x * rxs + (ny-1-y) * rys + (nz-1-z) * rzs;
             int data_elem = x * dxs + y * dys + z * dzs;
	    
	     *(result + result_elem) = (REAL)
	        *(bdata + data_elem) * dt + *(adata + data_elem) * dtm1;

	     if (Debug && (x == 0) && (y==0))
		printf("%d == result %d %g  data %d %g %g\n", z, result_elem, 
			*(result + result_elem),
			data_elem, *(bdata + data_elem), *(adata + data_elem));

		if ((x == 1) && (y == 1))
		DPRINT "%d == result %d %g  data %d %g %g\n", z, result_elem, 
			*(result + result_elem),
			data_elem, *(bdata + data_elem), *(adata + data_elem));
		
	     }
	  }

        Interp_cl += UTIMdiff( Interp_ud);
	}
	
        /* delete old files */
        CheckOld( ds, timidx);

	Reader_cl += UTIMdiff( Reader_ud);
	return 1;
    }
