/* $Id$ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#include "utim.h"

#ifdef VMD
#include <toolsa/MEMM.h>
#include <toolsa/str.h>
#include "arr.h"
extern int DBG_file;
#else
#include "str_util.h"
#define ARRmissing_f -9999.99	/* lame */
static int DBG_file = 0;
#endif

#include "gribr.h"

static int Debug = 0;	/* general, very small */
static int Debug2 = 0;	/* each time plan is read */
static int Debug3 = 0;	/* full - crashes on suns atleast drb */  

static int Table2_count = 100;
static Param_t *Table2 = NULL;

#define LEN2B(b)	((int) ((b[0]<<8)+b[1]))
#define LEN3B(b)	((int) ((b[0]<<16)+(b[1]<<8)+b[2]))
#define PDS_HAS_GDS(pds)        ((pds.gds_bms_flag & 128) != 0)
#define PDS_HAS_BMS(pds)        ((pds.gds_bms_flag & 64) != 0)

#define MAX_TOKEN 100

static double sigma_lev[] = {
           0.9950,    0.9821,    0.9644,    0.9425,    0.9159, 
           0.8838,    0.8458,    0.8014,    0.7508,    0.6943, 
           0.6329,    0.5681,    0.5017,    0.4357,    0.3720, 
           0.3125,    0.2582,    0.2101,    0.1682,    0.1326, 
           0.1028,    0.0782,    0.0580,    0.0418,    0.0288, 
           0.0183,    0.0101,    0.0027};

static double sigma_thick[] = {   
	   0.01,       0.01582,    0.01959,    0.02405,    0.02919, 
           0.03493,    0.04115,    0.04754,    0.05372,    0.05919, 
           0.06347,    0.06606,    0.06669,    0.06526,    0.06197, 
           0.05716,    0.05135,    0.04503,    0.03867,    0.03262, 
           0.02709,    0.02222,    0.01803,    0.01451,    0.0116, 
           0.00923,    0.00729,    0.00657 };

static double  gaussw[]  =  {
     .000831,        .001933,        .003035,        .004134,        .005228,
     .006316,        .007397,        .008471,        .009534,        .010588,
     .011629,        .012658,        .013673,        .014672,        .015656,
     .016622,        .017570,        .018498,        .019406,        .020292,
     .021156,        .021997,        .022813,        .023604,        .024370,
     .025108,        .025818,        .026500,        .027152,        .027775,
     .028367,        .028927,        .029456,        .029952,        .030415,
     .030844,        .031239,        .031600,        .031925,        .032216,
     .032471,        .032690,        .032873,        .033019,        .033129,
     .033203,        .033239,        .033239,        .033203,        .033129,
     .033019,        .032873,        .032690,        .032471,        .032216,
     .031925,        .031600,        .031239,        .030844,        .030415,
     .029952,        .029456,        .028927,        .028367,        .027775,
     .027152,        .026500,        .025818,        .025108,        .024370,
     .023604,        .022813,        .021997,        .021156,        .020292,
     .019406,        .018498,        .017570,        .016622,        .015656,
     .014672,        .013673,        .012658,        .011629,        .010588,
     .009534,        .008471,        .007397,        .006316,        .005228,
     .004134,        .003035,        .001933,        .000831  };

static double lat[] =  {
   88.541950,      86.653167,      84.753230,      82.850772,      80.947359,
   79.043483,      77.139349,      75.235055,      73.330657,      71.426186,
   69.521663,      67.617101,      65.712509,      63.807894,      61.903260,
   59.998611,      58.093950,      56.189279,      54.284599,      52.379912,
   50.475219,      48.570520,      46.665817,      44.761110,      42.856399,
   40.951685,      39.046969,      37.142250,      35.237529,      33.332806,
   31.428081,      29.523355,      27.618628,      25.713900,      23.809170,
   21.904439,      19.999708,      18.094976,      16.190243,      14.285510,
   12.380776,      10.476042,       8.571308,       6.666573,       4.761838,
    2.857103,        .952368,       -.952368,      -2.857103,      -4.761838,
   -6.666573,      -8.571308,     -10.476042,     -12.380776,     -14.285510,
  -16.190243,     -18.094976,     -19.999708,     -21.904439,     -23.809170,
  -25.713900,     -27.618628,     -29.523355,     -31.428081,     -33.332806,
  -35.237529,     -37.142250,     -39.046969,     -40.951685,     -42.856399,
  -44.761110,     -46.665817,     -48.570520,     -50.475219,     -52.379912,
  -54.284599,     -56.189279,     -58.093950,     -59.998611,     -61.903260,
  -63.807894,     -65.712509,     -67.617101,     -69.521663,     -71.426186,
  -73.330657,     -75.235055,     -77.139349,     -79.043483,     -80.947359,
  -82.850772,     -84.753230,     -86.653167,     -88.541950 };

static double lat_edge[] =  {
   90.000000,      87.664205,      85.738977,      83.826518,      81.917702,
   80.010430,      78.103958,      76.197953,      74.292245,      72.386738,
   70.481372,      68.576109,      66.670925,      64.765802,      62.860726,
   60.955688,      59.050682,      57.145701,      55.240742,      53.335800,
   51.430875,      49.525962,      47.621060,      45.716169,      43.811286,
   41.906410,      40.001541,      38.096678,      36.191820,      34.286966,
   32.382117,      30.477271,      28.572429,      26.667589,      24.762752,
   22.857917,      20.953084,      19.048253,      17.143424,      15.238596,
   13.333769,      11.428943,       9.524118,       7.619293,       5.714469,
    3.809646,       1.904823,        .000000,      -1.904823,      -3.809646,
   -5.714469,      -7.619293,      -9.524118,     -11.428943,     -13.333769,
  -15.238596,     -17.143424,     -19.048253,     -20.953084,     -22.857917,
  -24.762752,     -26.667589,     -28.572429,     -30.477271,     -32.382117,
  -34.286966,     -36.191820,     -38.096678,     -40.001541,     -41.906410,
  -43.811286,     -45.716169,     -47.621060,     -49.525962,     -51.430875,
  -53.335800,     -55.240742,     -57.145701,     -59.050682,     -60.955688,
  -62.860726,     -64.765802,     -66.670925,     -68.576109,     -70.481372,
  -72.386738,     -74.292245,     -76.197953,     -78.103958,     -80.010430,
  -81.917702,     -83.826518,     -85.738977,     -87.664205,     -90.000000 };

/************************ grid reading and decoding **************/
/** find next grib header
 *
 * pos = initial position to start looking at  ( = 0 for 1st call)
 *       returns with position of next grib header (units=bytes)
 * len_grib = length of the grib record (bytes)
 * buffer[buf_len] = buffer for reading/writing
 *
 * returns (char *) to start of GRIB header+PDS
 *         NULL if not found
 */

#define NTRY 3
#define LEN_HEADER_PDS (28+8)
#define MSEEK 1024

static unsigned char *Seek_grib(FILE *file, long *pos, long *len_grib,
        unsigned char *buffer) 
{

    int i, len, nread, count = 0;

    while (1) {
        if (fseek(file, *pos, SEEK_SET) == -1) {
            *len_grib = 0;
	    fprintf(stderr, " end of file\n");
            return (unsigned char *) NULL;
        }

        nread = fread(buffer, sizeof (unsigned char), MSEEK, file);

	len = nread  - LEN_HEADER_PDS;

	if (len == 0) {
            *len_grib = 0;
            fprintf(stderr, " end of file\n");
            return (unsigned char *) NULL;
        }
 

        for (i = 0; i < len; i++) {
            if (buffer[i] == 'G' && buffer[i+1] == 'R' && buffer[i+2] == 'I'
                && buffer[i+3] == 'B' && buffer[i+7] == 1) {
                    *pos += i;
                    *len_grib = (buffer[i+4] << 16) + (buffer[i+5] << 8) +
                            buffer[i+6];

		    if (Debug2 && (count+i > 0))
		 	fprintf(stderr, " %ld %ld %ld: that took %d bytes to find next record\n",
			    *pos, *len_grib, (*pos) + (*len_grib), count+i);

                    return (buffer+i);
            }
        }

        *pos += len;
	count += len;
    }
}

static int Read_grib(FILE *file, long pos, long len_grib, unsigned char *buffer) 
{
    int i;

    if (fseek(file, pos, SEEK_SET) == -1) {
            return 0;
    }

    i = fread(buffer, sizeof (unsigned char), len_grib, file);
    return (i == len_grib);
}

static int GetGridParams (GRIBR_t *g, int *nx, int *ny)
{
     if (g->pds.grid_id != 255)
	return -1;
    switch (g->gds.type) {
      case 4:			/* Gaussian Lat/Lon grid */
	*nx = LEN2B(g->gds.num_lon_points);	
	*ny = LEN2B(g->gds.num_lat_points);	
	break;
      case 0:			/* Gaussian Lat/Lon grid */
	*nx = LEN2B(g->gds.num_lon_points);	
	*ny = LEN2B(g->gds.num_lat_points);	
	break;
    }

    return 1;
}

static float CvtDeg3( unsigned char *b)
{
    int sign = b[0] & 0x80;
    int val = LEN3B(b);

    val &= 0x7fffff;

    return (sign ? -1.0 : 1.0) * val * .001;
}

static float CvtS2( unsigned char *b)
{
    int sign = b[0] & 0x80;
    int val = LEN2B(b);

    val &= 0x7fff;

    return (sign ? -1.0 : 1.0) * val;
}

static float CvtDeg2( unsigned char *b)
{
	return CvtS2(b) * .001;
}

static int PrintGridParams (FILE *fp, GRIBR_t *g)
{
    if (g->pds.grid_id != 255)
        return -1;

    switch (g->gds.type) {
      case 4:                   /* Gaussian Lat/Lon grid */
        fprintf(fp,"Nlons = %d\n", LEN2B(g->gds.num_lon_points));
        fprintf(fp,"Nlats = %d\n", LEN2B(g->gds.num_lat_points));
        fprintf(fp,"First lat = %6.3f\n", CvtDeg3(g->gds.la1));
        fprintf(fp,"Last  lat = %6.3f\n", CvtDeg3(g->gds.la2));
        fprintf(fp,"First lon = %6.3f\n", CvtDeg3(g->gds.lo1));
        fprintf(fp,"Last  lon = %6.3f\n", CvtDeg3(g->gds.lo2));
	fprintf(fp,"Res flags = %0x\n", g->gds.res_flags);
        fprintf(fp,"delta lon = %6.3f\n", CvtDeg2(g->gds.di));
        fprintf(fp,"nlat circles = %d\n", LEN2B(g->gds.dj));
	fprintf(fp,"Scan mode = %x\n", g->gds.scan_mode);
        break;
    }

    return 1;
}


static char *Level_desc(int lev_id, int lev_param)
{
    int o11 =  lev_param / 256;
    int o12 =  lev_param % 256;
    static char stemp[100];

        switch (lev_id) {

        case 1: return ("surface");
        case 2: return ("cloud base");
        case 3: return ("cloud top");
        case 4: return ("0C isotherm");
        case 5: return ("cond lev");
        case 6: return ("max wind lev");
        case 7: return ("tropopause");
        case 8: return ("nom. top");
        case 9: return ("sea bottom");
        case 200:
        case 10: return ("atmos col");
        case 12:
        case 212: return ("low cld bot");
        case 13:
        case 213: return ("low cld top");
        case 14:
        case 214: return ("low cld lay");
        case 22:
        case 222: return ("mid cld bot");
        case 23:
        case 223: return ("mid cld top");
        case 24:
        case 224: return ("mid cld lay");
        case 32:
        case 232: return ("high cld bot");
        case 33:
        case 233: return ("high cld top");
        case 34:
        case 234: return ("high cld lay");


        case 100: sprintf(stemp, "%d mb",lev_param);
                return stemp;
        case 101: sprintf(stemp, "%d-%d mb",o11*10,o12*10);
                return stemp;
        case 103: sprintf(stemp, "%d m above MSL",lev_param);
                return stemp;
        case 104: sprintf(stemp, "%d-%d m above msl",o11*100,o12*100);
                return stemp;
        case 105: sprintf(stemp, "%d m above gnd",lev_param);
                return stemp;
        case 106: sprintf(stemp, "%d-%d m above gnd",o11*100,o12*100);
                return stemp;
        case 107: sprintf(stemp, "sigma=%.4f",lev_param/10000.0);
                return stemp;
        case 108: sprintf(stemp, "sigma %.2f-%.2f",o11/100.0,o12/100.0);
                return stemp;
        case 109: sprintf(stemp, "hybrid lev %d",lev_param);
                return stemp;
        case 110: sprintf(stemp, "hybrid %d-%d",o11,o12);
                return stemp;
        case 111: sprintf(stemp, "%d cm down",lev_param);
                return stemp;
        case 112: sprintf(stemp, "%d-%d cm down",o11,o12);
                return stemp;
        case 113: sprintf(stemp, "%dK",lev_param);
                return stemp;
        case 114: sprintf(stemp, "%d-%dK",475-o11,475-o12);
                return stemp;
        case 115: sprintf(stemp, "%d mb above gnd",lev_param);
                return stemp;
        case 116: sprintf(stemp, "%d-%d mb above gnd",o11,o12);
                return stemp;
        case 121: sprintf(stemp, "%d-%d mb",1100-o11,1100-o12);
                return stemp;
        default:
                return "unknown level type";
        }
}


#define BUFF_ALLOC0 40000

static GRIBR_t *Read_Gribr(FILE *input, long pos)
{
    static GRIBR_t gribr;
    static unsigned char *buffer = NULL;
    static int buffer_size = 0;

    unsigned char *mpos;
    long len_grib;
    unsigned long len_pds, len_gds, len_bms, len_bds;


    	if (NULL == buffer) {
            buffer_size = BUFF_ALLOC0;
	    buffer = CALLOC(buffer_size, unsigned char);
	}

	/* find the next grib record, return its position and its length */
	mpos = Seek_grib(input, &pos, &len_grib, buffer);
	if (mpos == NULL) {
	    return NULL;
	}
	gribr.total_len = len_grib;
	gribr.file_offset = pos;

	/* make sure we got enough room for the grib record */
        if (len_grib + mpos - buffer > buffer_size) {
            buffer_size = len_grib + mpos - buffer + 1000;
            buffer = (unsigned char *) realloc( buffer, buffer_size);
            if (buffer == NULL) {
                fprintf(stderr,"ran out of memory\n");
                return NULL;
            }
        }

        /* read entire grib record */
        if (!Read_grib(input, pos, len_grib, buffer))
	    {
	    fprintf(stderr, "end of file\n");
	    return NULL;
	    }

	/* parse grib message */
	mpos = buffer + 8;
        memcpy( &gribr.pds, mpos, sizeof(GRIBRpds_t));
        len_pds = LEN3B(gribr.pds.length);
        mpos += len_pds;

        if (PDS_HAS_GDS( gribr.pds)) {
	    memcpy( &gribr.gds, mpos, sizeof(GRIBRgds_latlon_t));
            len_gds = LEN3B(gribr.gds.length);
  	    mpos += len_gds;
        }
        else {
	    memset( &gribr.gds, 0, sizeof(GRIBRgds_latlon_t));
	    len_gds = 0;
        }

        if (PDS_HAS_BMS( gribr.pds)) {
	    memcpy( &gribr.bms, mpos, 6);
	    gribr.bms.bit_map = mpos + 6;
            len_bms = LEN3B(gribr.bms.length);
	    mpos += len_bms;
        }
        else {
	    memset( &gribr.bms, 0, sizeof(GRIBRbms_t));
	    gribr.bms.bit_map = NULL;
	    len_bms = 0;
        }

        memcpy( &gribr.bds, mpos, 11);
	gribr.bds.data = mpos + 11;
        len_bds = LEN3B(gribr.bds.length);
        mpos += len_bds;

        /* end section - "7777" in ascii */
        if (mpos[0] != 0x37 || mpos[1] != 0x37 ||
            mpos[2] != 0x37 || mpos[3] != 0x37) {
            fprintf(stderr,"  missing 7777 end gribr start = %lx len = %ld\n", 
		pos, mpos - buffer);
	    fprintf(stderr,"   lens = %ld %lu %lu %lu %lu \n", len_grib, 
		len_pds, len_gds, len_bms, len_bds); 
        }


    return &gribr;
}

/************* data unpacking ************************/

/* ibm = sign * 2^-24 * mant * 16^(exp-64) */
static double ibm2flt(unsigned char *ibm) {

        int positive, power;
        unsigned int abspower;
        long int mant;
        double value, exp;

        positive = (ibm[0] & 0x80) == 0;
        mant = (ibm[1] << 16) + (ibm[2] << 8) + ibm[3];
        power = (int) (ibm[0] & 0x7f) - 64;

        abspower = power > 0 ? power : -power;

        /* calc exp */
        exp = 16.0;
        value = 1.0;
        while (abspower) {
                if (abspower & 1) {
                        value *= exp;
                }
                exp = exp * exp;
                abspower >>= 1;
        }

        if (power < 0) value = 1.0 / value;
        value = value * mant / 16777216.0;
        if (positive == 0) value = -value;
        return value;
}

static int Read_data( GRIBR_t *gr, REAL *outdata, int nout)
    {
	unsigned char *indata = gr->bds.data;
	int n_bits = (int) gr->bds.number_of_bits;
	int dec_scale = CvtS2(gr->pds.scale_factor);
	int bin_scale = CvtS2(gr->bds.binary_scale_factor);
	double  ref = ibm2flt(gr->bds.reference_value);
	double idscale, bscale;
	unsigned char *bitmap = gr->bms.bit_map;
	unsigned char map_mask, got;

    	int i, k, wbits;
    	unsigned int pval, hi;
	unsigned char work;

        bscale = pow(2.0, bin_scale);
        idscale = 1.0/pow(10.0, dec_scale);

    /* decode numeric data */
    /* value * 10 ^ dec_scale = ref + pack * 2 ^ bin_scale */
    /* value = (ref + pack * bscale) * idscale  */

      work = *indata;
      wbits = 8;
      map_mask = 0x80;

      for (i = 0; i < nout; i++) {

	if (bitmap) {

            got = (*bitmap & map_mask);
            if ((map_mask >>= 1) == 0) {
                map_mask = 128;
                bitmap++;
            }
            if (got == 0) {
               *outdata++ = ARRmissing_f;
                continue;
            }
        }

        pval = 0;
        k = n_bits;

        while (k) {
            if (k >= wbits) {
                pval <<= wbits;
		hi = work;
		hi >>= (8-wbits);		/* get hi wbits from work */
		pval += hi;
                k -= wbits;

                indata++;
		work = *indata;
		wbits = 8;
            }
            else {
		pval <<= k;
		hi = work;
		hi >>= (8-k);		/* get hi k bits from work */
		pval += hi;

		work <<= k; 		/* finished k bits */
		wbits -= k;		
		k = 0;			/* in any case, this pval is done */

	        if (wbits == 0)		/* advance to new indata char */
		    {
                    indata++;
		    work = *indata;
		    wbits = 8;
		    }
            }
        }
        *outdata++ = (REAL) ((ref + bscale*pval) * idscale);
      }

    if (Debug2) {
      fprintf(stderr, "GRIBR read_data:  dec_scale = %d, bin_scale = %d, ref = %g",
	dec_scale, bin_scale, ref);
      fprintf(stderr, " nbits = %d, consumed %d bytes, expected %d, grib_len %d\n", n_bits,
	(int)(indata - gr->bds.data), (nout * n_bits+7)/8, gr->total_len);
    }

    return 1;
    }

int GRIBRread_data( GRIBRhead_t *gh, int field_idx, int level_idx,
        int time_idx, REAL *outdata)
    {
    	Field_t *fld;
    	Level_t *lev;
    	Time_t *tim;
	GRIBR_t *gr;
	int i, data_size = gh->nx * gh->ny;

	if (field_idx >= gh->nflds)
	    {
	    fprintf(stderr, "GRIBRread_data bad field index %d\n", field_idx);
	    return 0;
	    }
        fld = &gh->fld[field_idx];

	/**** read all levels ******/
	if (level_idx < 0)
	  {
          for (i=0; i< fld->nlevels; i++)
            {
            lev = &fld->level[i];

            if (time_idx >= lev->ntimes)
                {
                fprintf(stderr,"GRIBRread_data fld %d lev %d does not have time %d\n", 
		    field_idx, i, time_idx);
                return 0;
                }

            tim = &lev->time[time_idx];

            if (NULL == (gr = Read_Gribr( gh->fp, tim->offset)))
                return 0;
            if (0 == Read_data( gr, outdata, data_size))
                return 0;

            outdata += data_size;
            }
	  return 1;
	  }

	/**** read one level ******/
	if (level_idx >= fld->nlevels)
	    {
	    fprintf(stderr, "GRIBRread_data bad level index %d\n", level_idx);
	    return 0;
	    }
	lev = &fld->level[level_idx];

	if (time_idx >= lev->ntimes)
	    {	
	    fprintf(stderr,"fld %d lev %d does not have time %d\n", field_idx,
		    level_idx, time_idx);
	    return 0;
	    }

	tim = &lev->time[time_idx];

	if (NULL == (gr = Read_Gribr( gh->fp, tim->offset)))
	    return 0;
	if (0 == Read_data( gr, outdata, data_size))
	    return 0;

	return 1;
    }

/************** Keep track of Products available: fld / lev / time ************/

static Level_t *AddLevel( Field_t *fld, int level, int param)
    {
	Level_t *lev;

	if (fld->nlevels >= fld->nlev_alloc)
	    {
	    fld->level = REALLOC( fld->level, fld->nlev_alloc+10, Level_t);
	    memset( fld->level + fld->nlev_alloc, 0, 10 * sizeof(Level_t));
	    fld->nlev_alloc += 10;
	    }

	lev = fld->level + fld->nlevels;
	lev->id = level;
	lev->param = param;
	fld->nlevels++;
	return lev;
    }

static Time_t *AddTime( Level_t *lev, time_t time)
    {
        Time_t *tim;

        if (lev->ntimes >= lev->ntim_alloc)
            {
            lev->time = REALLOC( lev->time, lev->ntim_alloc+10, Time_t);
	    memset( lev->time + lev->ntim_alloc, 0, 10 * sizeof(Time_t));
            lev->ntim_alloc += 10;
            }

        tim = lev->time + lev->ntimes;
        tim->utime = time;
        lev->ntimes++;
        return tim;
    }

static Level_t *FindLevel( Field_t *fld, int level, int param)
    {
	int i;
	for (i=0; i<fld->nlevels; i++)
	    if ((fld->level[i].id == level) &&
	       (fld->level[i].param == param))
		return fld->level + i;
	return NULL;
    }

static Time_t *FindTime( Level_t *lev, time_t time)
    {
        int i;
        for (i=0; i<lev->ntimes; i++)
            if (lev->time[i].utime == time)
                return lev->time + i;
        return NULL;
    }


static int AddProduct( GRIBRhead_t *gh, int field, int level, int lev_param, 
	time_t time, long pos)
{
    Field_t *fld;
    Level_t *lev;
    Time_t *tim;

	if (field >= gh->nflds)
	    {
	    int add = field + 100 - gh->nflds;
	    gh->fld = REALLOC( gh->fld, gh->nflds + add, Field_t);
	    memset( gh->fld + gh->nflds, 0, add * sizeof(Field_t));
	    gh->nflds += add;
	    }

	fld = &gh->fld[field];
	fld->name = gh->table2[field].short_name;
	fld->long_name = gh->table2[field].long_name;
	fld->units = gh->table2[field].unit;

	lev = FindLevel( fld, level, lev_param);
	if (NULL == lev)
	   lev = AddLevel( fld, level, lev_param);
	lev->param = lev_param;

        if (FindTime( lev, time))
            {
            printf("duplicate grib record skipped at pos %ld: fld %d %s lev %d %d time %d !!!\n",
                pos, field, fld->name, level, lev_param, time);
            return 0;
            }
	tim = AddTime( lev, time);
	tim->offset = pos;
	
    	return 1;	
    }

static void PrintIndex( FILE *fp, GRIBR_t *g, time_t time)
{
    GRIBRpds_t *p = &g->pds;

    fprintf(fp, "%3d %3d %3d %3d %9ld %9ld # %10s %10s\n",
                p->parameter, p->level_type, p->level_params[0], p->level_params[1],
                (long)time,
                g->file_offset,
                Table2[p->parameter].short_name,
                Level_desc( (int) p->level_type, LEN2B(p->level_params)));

}

static int ReadIndex(GRIBRhead_t *gh, char *filename)
{
    FILE *fp = fopen(filename, "r");
    char line[1000], token[MAX_TOKEN];
    int i;
    long time, pos;

    if (NULL == fp)
        return 0;

    /* make sure it finished last time */
    while (NULL != fgets( line, 1000, fp));
    if (NULL == strstr(line,"DONE-DONE-DONE"))
    	{
	fclose(fp);
	return 0;
	}

    fseek( fp, SEEK_SET, 0);
    while (NULL != fgets( line, 1000, fp))
        {
        char *ptr = line;
        int p[4];

        if (NULL != strstr(line,"DONE-DONE-DONE"))
	    break;

        for (i=0; i<4; i++)
            {
            STRtokn( &ptr, token, MAX_TOKEN, " ");
            p[i] =  atoi(token);
            }
        STRtokn( &ptr, token, MAX_TOKEN, " ");
        time =  atol(token);
        STRtokn( &ptr, token, MAX_TOKEN, " ");
        pos =  atol(token);

        AddProduct( gh, p[0], p[1], 256 * p[2] + p[3], time, pos);
        }

    fclose(fp);
    return 1;
}

static void PrintProducts( FILE *fp, GRIBRhead_t *gh)
{
    int f, l;
    int need_times = 0;
    Param_t *table = gh->table2;

    fprintf(fp, "Products Available by Level\n");
    for (f=0; f<gh->nflds; f++)
	{
	Field_t *fld = &gh->fld[f];
	if (fld->nlevels == 0)
	    continue;

	fprintf(fp, "   %3d %s %s %s nlev = %d\n", f, table[f].short_name,
                table[f].unit, table[f].long_name, fld->nlevels);

	if (DBG_file)
	  for (l=0;l<fld->nlevels;l++)
	    {
	    Level_t *lev = &fld->level[l];

	    if (need_times == 0)
		need_times = lev->ntimes;
	    else if (need_times != lev->ntimes)
		fprintf(fp," HEYY! %d != %d", need_times, lev->ntimes);
		
	    fprintf(fp, "      %3d %3d %s\n", 
		lev->id, lev->param,
		Level_desc(lev->id, lev->param));
	    }
	}
}

static int TimeCmp( const void *e1, const void *e2)
    {
        Time_t *t1 = (Time_t *) e1;
        Time_t *t2 = (Time_t *) e2;

        if (t1->utime == t2->utime)
            return 0;
        else if (t1->utime > t2->utime)
            return 1;
        else
            return -1;
    }

static int LevCmp( const void *e1, const void *e2)
    {
	Level_t *lev1 = (Level_t *) e1;
	Level_t *lev2 = (Level_t *) e2;

	if (lev1->id == lev2->id)
	    {
	    if (lev1->param == lev2->param)
		return 0;	
	    else if (lev1->param > lev2->param)
		return -1;	
	    else
		return 1;
	    }
	else
	    {
	    if (lev1->id == lev2->id)
		return 0;	
	    else if (lev1->id > lev2->id)
		return -1;	
	    else
		return 1;
	    }
    }


/* sort levels and times, so theres no dependency on order in the file */
/* might as well do some consistency checks here, also */
static void SortProducts( GRIBRhead_t *gh)
{
    int f, l;
    int need_times = 0;

    for (f=0; f<gh->nflds; f++)
        {
        Field_t *fld = &gh->fld[f];
        if (fld->nlevels == 0)
            continue;

	qsort( fld->level, fld->nlevels, sizeof( Level_t), LevCmp);


        for (l=0;l<fld->nlevels;l++)
            {
            Level_t *lev = &fld->level[l];

	    qsort( lev->time, lev->ntimes, sizeof( Time_t), TimeCmp);

            if (need_times == 0)
		{
                need_times = lev->ntimes;
		if (Debug)
		    fprintf(stderr," GRIB SortProducts: file has %d times\n", need_times);
		}
            else if (need_times != lev->ntimes)
                fprintf(stderr," POSSIBLE PROBLEM: GRIB SortProducts fld %s #times %d != %d", 
			fld->name, need_times, lev->ntimes);
            }
        }
}


static int ReadTable2(char *filename)
{
    FILE *fp = fopen(filename, "r");
    char line[1000], token[MAX_TOKEN];
    int id;

    if (NULL == fp)
	{
	fprintf(stderr,"cant open parameter file %s\n", filename);
	return 0;
	}

    Table2 = CALLOC( Table2_count, Param_t);

    while (NULL != fgets( line, 1000, fp))
	{
	char *ptr = line;

	STRtokn( &ptr, token, MAX_TOKEN, ":");   
	id = atoi(token);

	if (id >= Table2_count)
	    {
	    Table2 = REALLOC( Table2, (Table2_count+100), Param_t);
	    memset(Table2 + Table2_count, 0, 100 * sizeof(Param_t));
	    Table2_count = id + 100;
	    }

	if (NULL == STRtokn( &ptr, token, MAX_TOKEN, ":"))
	    goto next;   
	Table2[id].long_name = strdup(token);
	STRblnk( Table2[id].long_name);

	if (NULL == STRtokn( &ptr, token, MAX_TOKEN, ":"))
	    goto next;   
	Table2[id].unit = strdup(token);
	STRblnk( Table2[id].unit);

	if (NULL == STRtokn( &ptr, token, MAX_TOKEN, ":"))
	    goto next;   
	Table2[id].short_name = strdup(token);
	STRblnk( Table2[id].short_name);

    next:
	if (Debug3)
	    fprintf(stderr,"%d: <%s> <%s> <%s>\n", id, Table2[id].long_name,
		Table2[id].unit, Table2[id].short_name);
	}

    fclose(fp);
    return 1;
}

void GRIBRprint( FILE *fp, GRIBRhead_t *head)
    {
        GRIBR_t *g;

   	/* read in first record to get global info from */
   	if (NULL == (g = Read_Gribr(head->fp, 0)))
            return;
	PrintGridParams(fp, g);

	PrintProducts( fp, head);
    }

void GRIBRclose( GRIBRhead_t *head)
    {
    	int f, l;

    	for (f=0; f<head->nflds; f++)
          {
          Field_t *fld = &head->fld[f];
          for (l=0;l<fld->nlevels;l++)
              {
              Level_t *lev = &fld->level[l];
	      FREE(lev->time);
              }
	  FREE( fld->level);
          }
	FREE(head->fld);

	fclose( head->fp);
	FREE (head->lev_edge);
	FREE (head);
    }
	
GRIBRhead_t *GRIBRopen(char *filename, int force, int print)
{
    FILE *input;
    char stemp[400];
    GRIBR_t *g;
    GRIBRhead_t *ghead = MALLOC( GRIBRhead_t);
    long pos;
    int i;

    if (Table2 == NULL)
#ifdef VMD
        ReadTable2("grib/grib.table2.rean");
#else
        ReadTable2("grib.table2.rean");
#endif

    ghead->table2 = Table2;

    if ((input = fopen(filename,"rb")) == NULL) {
        fprintf(stderr,"GRIBRopen: could not open file: %s\n", filename);
	return NULL;
    }

    /* read in index file if you can */
    sprintf(stemp, "%s.index", filename);
    if (force || (0 == ReadIndex( ghead, stemp)))
	{
	FILE *index = fopen( stemp, "w");
	long pos = 0;

	/* otherwise construct it */
        while (NULL != (g = Read_Gribr(input, pos)))
            {
            time_t time = UTIMdate2unix( 1900+ (int) g->pds.year,  (int) g->pds.month,  
		    (int) g->pds.day, (int) g->pds.hour, (int) g->pds.minute,  0);

            if (Debug3)
		PrintIndex( stderr, g, time);

            PrintIndex( index, g, time);

            AddProduct( ghead, (int) g->pds.parameter,
                (int) g->pds.level_type,
                LEN2B(g->pds.level_params),
		time, 
                g->file_offset);

	    pos = g->file_offset + g->total_len;
            }

	fprintf(index,"DONE-DONE-DONE");
	fclose(index);
	}
    SortProducts( ghead);

   /* read in first record to get global info from */
   pos = 0;
   if (NULL == (g = Read_Gribr(input, pos)))
	return NULL;

    GetGridParams( g, &ghead->nx, &ghead->ny);
    ghead->nz = 28;

    ghead->lat = lat;
    ghead->lat_edge = lat_edge;
    ghead->lev = sigma_lev;
    ghead->gaussw = gaussw;

    ghead->lev_edge = CALLOC( ghead->nz+1, double);
    ghead->lev_edge[0] = ghead->lev[0] + .5 * sigma_thick[0];
    for (i=1; i<ghead->nz; i++)
        ghead->lev_edge[i] = ghead->lev_edge[i-1] - sigma_thick[i-1];
    ghead->lev_edge[ghead->nz] = 0.0;

    ghead->fp = input;
    return ghead;
}

#ifdef TEST
int main(int argc, char *argv[])
{
	char *filename = "/mop3/caron/data/nmc/A04860";
	GRIBRhead_t *gh = NULL;
	long pos = 0;
	int i, force = 0, print = 0;

	for (i=1; i<argc; i++)
	    if (0 == strcmp(argv[i], "-file"))
		filename = argv[++i];
	    else if (0 == strcmp(argv[i], "-force"))
		force = 1; 
	    else if (0 == strcmp(argv[i], "-print"))
		print = 1; 

	ReadTable2("grib.table2.rean");
	if (NULL == (gh = GRIBRopen(filename, force, print)))
	    return -1;

	if (print)
	    GRIBRprintProducts(stderr, gh);

	fclose (gh->fp);
}
#endif
