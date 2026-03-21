/* $Id$ */
#ifndef MSS_WAS_INCLUDED
#define MSS_WAS_INCLUDED

/* filelist format:
   skip blanks or lines tarting with '#'

    mss_explicit:
	<filename> <size in bytes> <time=YYMMDDHHMM>
 */

typedef enum {mss_msls, mss_explicit} MSSlist_t;

typedef struct MSSfile_t {
    char 	*name;
    long 	size;
    time_t	utime;	/* unix time */

    int  has_local;	/* is on local file system */
    int  in_progress;	/* waiting for mss to finish */

	/* used only by calling program */
    FILE	*fp;
} MSSfile_t;

typedef struct MSS_t {
    char *filelist;  	/* file with list of ms files in msls format */

    char *rpath;	/* mss path */
    char *lpath;	/* local path */

    int nfiles;
    MSSfile_t *file;

    int current;
} MSS_t;


/* set up management of mss <--> local files */
extern MSS_t *MSSinit( char *rpath, char *lpath, MSSlist_t type, char *filelist);

/* check if we already have it. if not, fetch from mass store.
   return -1: ERR, 0: in progress, 1 : ready */
extern int MSSget( MSS_t *mss, char *filename, int wait);

/* this is for sequentially accessing the files.
   when its called, the next file in sequence is fetched from mss.
   call it way before you actually need the file.  */
extern void MSSget_next( MSS_t *mss);

/* delete local copies more than <keep> before current */
extern void MSSclean( MSS_t *mss, int keep);
extern void MSSdelete( MSS_t *mss, MSSfile_t *file);


extern int MSSfile2idx( MSS_t *mss, char *filename);

#endif
