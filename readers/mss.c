/* $Id$ */
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "mss.h"
#include "str_util.h"
#include "utim.h"

#define MAX_BUFF 200
#define MAX_COMMAND 500
#define MAX_PATHNAME 200
#define MAX_TOKEN 100

#define gmtime_str(t) asctime( gmtime( t))

static int debug = 1;
static int debug2 = 0;

/* convert YYYYMMDDHHMM secs -> time_t */
/* modified to require a 4 digit year 19 Sept 2001 bundy */
static time_t ConvertTime( char *yyyymmddhhmm)
    {
        int y,  m,  d,  h,  min;
	char s[4], *in = yyyymmddhhmm;
	memcpy(s, in, 4);
	y = atoi(s);
	s[2] = s[3] = 0;
	memcpy(s, in+4, 2);
	m = atoi(s);
	memcpy(s, in+6, 2);
	d = atoi(s);
	memcpy(s, in+8, 2);
	h = atoi(s);
	memcpy(s, in+10, 2);
	min = atoi(s);

	return UTIMdate2unix( y,  m,  d, h, min,  0);
    }

/* set up management of mss <--> local files */
MSS_t *MSSinit( char *rpath, char *lpath, MSSlist_t type, char *filelist)
    {
	MSS_t *mss;
	MSSfile_t *file;

	FILE *fp;
	char buff[MAX_BUFF], token[MAX_TOKEN], *ptr;
	int err = 0, count;

	if (debug)
	    fprintf(stderr, "MSSinit %s\n", filelist);

	if (NULL == (fp = fopen(filelist, "r")))
	    {
	    printf("MSSinit: Cant open filelist %s\n", filelist);
	    return NULL;
	    }
	
	count = 0;
	while (NULL != fgets( buff, MAX_BUFF, fp))
	    {
	    if (0 == STRblnk( buff))
		continue;
	    if (buff[0] == '#')
		continue;
	    count++;
	    }

	mss = MALLOC( MSS_t);
	mss->filelist = strdup(filelist);
	mss->rpath = strdup(rpath);
	mss->lpath = strdup(lpath);

	mss->nfiles = count;
	mss->file = CALLOC( mss->nfiles, MSSfile_t);
	file = mss->file;

	fseek(fp, 0, SEEK_SET);
        while (NULL != fgets( buff, MAX_BUFF, fp))
            {
            if (0 == STRblnk( buff))
                continue;
            if (buff[0] == '#')
                continue;

	    ptr = buff;
	    if (type == mss_msls) 
		{
	    	STRtokn( &ptr, token, MAX_TOKEN, " ");
	    	STRtokn( &ptr, token, MAX_TOKEN, " ");
	    	STRtokn( &ptr, token, MAX_TOKEN, " ");
	    	file->size = atol(token);
	    	STRtokn( &ptr, token, MAX_TOKEN, " ");
	    	file->name = strdup(token);
	       	} 
	    else  if (type == mss_explicit)
	 	{
	    	STRtokn( &ptr, token, MAX_TOKEN, " ");
	    	file->name = strdup(token);
	    	STRtokn( &ptr, token, MAX_TOKEN, " ");
	    	file->size = atol(token);
	    	STRtokn( &ptr, token, MAX_TOKEN, " ");
		file->utime = ConvertTime(token);
		}
	    else
		{
		printf("MSSinit: unknown list file type %d\n", type);
		err = 1;
		break;
		}

	    if (debug2)
		fprintf(stderr,"  %s %ld %ld %s", file->name, file->size, (long)file->utime,
		    gmtime_str( &file->utime));

	    file++;
	    }

	fclose( fp);

	return err ? NULL : mss;
    }

int MSSfile2idx( MSS_t *mss, char *filename)
    {
	int idx;

        /* find the file in the list */
        for (idx=0; idx< mss->nfiles; idx++)
            if (0 == strcasecmp(filename, mss->file[idx].name))
                break;

        if (idx >= mss->nfiles)
            {
            printf("MSSget: file %s not in list\n", filename);
            return -1;
            }
	return idx;
    }

static int Check( MSS_t *mss, int idx)
    {
        MSSfile_t *file = &mss->file[idx];
        char pathname[MAX_PATHNAME];
        struct stat fs;
	int ret = 0;

        /* see if it exists, what its' size */
        sprintf( pathname, "%s/%s", mss->lpath, file->name);

        errno = 0;
        stat( pathname, &fs);
	
        if (0 == errno)
            {
	      file->size = 0;
	      fprintf(stderr,"WARNING: disabled file size check");
            if (fs.st_size >= file->size)
                {
		  /*	      fprintf(stderr,"DEBUG2 Check: errno %d fs.st_size %d >= file->size %d \n",
			      errno, fs.st_size , file->size);*/
                file->has_local = 1;
                file->in_progress = 0;
                mss->current = idx;	/* for get_next */

                if (debug)
                    fprintf(stderr," got file %s\n", file->name);
                ret = 1;
                }
	    else 
	      {
		fprintf(stderr,"DRBDBG: Warning, grib file size not as big as grblist claims");
		fprintf(stderr,"            fs.st_size %d < file->size %d \n",  
			fs.st_size , file->size);
	      }
            }
        else if (errno != ENOENT)
            {
            printf("MSSget: file stat on %s failed = %s\n", pathname,
                strerror(errno));
            ret = -1;
            }

        if (debug)
            fprintf(stderr," check %s = %d\n", pathname, ret);


	return ret;
    }

/* check if we already have it. if not, fetch from mass store.
   if wait = 1, sleep till file is ready 
   return -1: ERR, 0: in progress, 1 : ready */
int MSSget( MSS_t *mss, char *filename, int wait)
    {
	MSSfile_t *file;
	char cmd[MAX_COMMAND];
	int idx, ret, fail = 0;

	if (debug)
	    fprintf(stderr,"MSSget on %s", filename);

	/* find the file in the list */
	if (0 > (idx = MSSfile2idx(mss, filename)))
	    return -1;
	file = &mss->file[idx];

	/* see if it exists, what its' size */
	ret = Check( mss, idx);

	if ((ret == 0) && !file->in_progress)
	    {
	    /* send a message to the mss to get the file */
#if defined(HAS_MSREAD)
	    sprintf( cmd, "msread -f BI -R %s/%s %s/%s",
	        mss->lpath, filename, mss->rpath, filename);
#else
	    /*
	    sprintf( cmd, "nrnet msget %s/%s l df=bi r flnm=%s/%s",
	        mss->lpath, filename, mss->rpath, filename);
	    */
	    sprintf( cmd, "msrcp mss:%s/%s %s/%s; cosconvert -b %s/%s",
	         mss->rpath, filename, mss->lpath, filename, mss->lpath, filename);
#endif
	    if (debug)
	        fprintf(stderr," send mss get = <%s>\n", cmd);

	    fail = (0 != system( cmd));
	    file->in_progress = 1;
	    }

	if (wait)
	    {
	    if (fail)
		return -1;
	    while (0 == (ret = Check(mss, idx))) 
		sleep(60);
	    return ret;
	    }

	return 0;
    }


/* this is for sequentially accessing the files.
   when its called, the next file in sequence is fetched from mss.
   call it way before you actually need the file.  */
void MSSget_next( MSS_t *mss)
    {
	if (mss->current >= mss->nfiles-1)
	    {
	    printf("MSSget_next : no more files: %d, %d\n",
		mss->current+1, mss->nfiles);
	    }

	MSSget( mss, mss->file[mss->current+1].name, 0);
    }

/* delete local copy */
void MSSdelete( MSS_t *mss, MSSfile_t *file)
    {
        char cmd[MAX_COMMAND];

        if (file->has_local)
            {
            /* delete the file */

            sprintf( cmd, "rm %s/%s", mss->lpath, file->name);
            system( cmd);
  
            file->has_local = 0;
            printf("MSSclear: deleted file %s\n",file->name);
            }

    }

/* delete local copy */
void MSSclean( MSS_t *mss, int keep)
    {
	int i;

	if (debug)
	    fprintf(stderr,"MSSclean\n");

	for (i=mss->current - keep; i >= 0; i--)
	    {
	    MSSfile_t *file = &mss->file[i];
	    MSSdelete( mss, file);
	    }

    }

