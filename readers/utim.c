/* $Id$ */
/* UTIM.c **/

#include <time.h>
#include <sys/time.h>
#include <math.h>
#include "utim.h"

#if defined(SUNOS4) || defined(SUNOS5) || defined(IRIX5)
extern int gettimeofday(struct timeval *tp,struct timezone *tzp);
#endif

/*************************************************************************
 *	JULIAN_DATE: Calc the Julian calandar Day Number
 *	As Taken from Computer Language- Dec 1990, pg 58
 */

long UTIMjulian_date(int day, int month, int year)
{
    long 	a,b;
    double	yr_corr;

    /* correct for negative year */
    yr_corr = (year > 0? 0.0: 0.75);
    if(month <=2) {
	year--;
	month += 12;
    }
    b=0;

    /* Cope with Gregorian Calandar reform */
    if(year * 10000.0 + month * 100.0 + day >= 15821015.0) {
	a = year / 100;
	b = 2 - a + a / 4;
    }
    return((365.25 * year - (double)yr_corr) + 
	   (long) (30.6001 * ((long)month +1)) + (long)day + 1720994L + b);
}

/*************************************************************************
 *	CALENDAR_DATE: Calc the calandar Day from the Julian date
 *	As Taken from Computer Language- Dec 1990, pg 58
 */

void UTIMcalendar_date(long jdate,long *day, long *month, long *year)
{
	long	a,b,c,d,e,z,alpha;

	z = jdate +1;

	/* Gregorian reform correction */
	if (z < 2299161) { 
		a = z; 
	} else {
		alpha = (long) ((z - 1867216.25) / 36524.25);
		a = z + 1 + alpha - alpha / 4;
	}

	b = a + 1524;
	c = (long) ((b - 122.1) / 365.25);
	d = (long) (365.25 * c);
	e = (long) ((b - d) / 30.6001);
	*day = (int) b - d - (long) (30.6001 * e);
	*month = (int) (e < 13.5)? e - 1 : e - 13;
	*year = (int) (*month > 2.5)? (c - 4716) : c - 4715;
}


time_t UTIMdate2unix( int year, int month, int day, int hour, int min, 
	int sec)
{
	long	u_day,this_day,days;
	long	u_time;

	/* quick fix for hardcoded 1900 in gribr.c, time_t function,
	   Del Corral Nov 2001*/

	if (year < 1950) year=year+100; 

	u_day = UTIMjulian_date(1,1,1970);
	this_day = UTIMjulian_date(day,month,year);

	days = this_day - u_day;

	u_time = (days * 86400) + (hour * 3600) + (min * 60) + sec;

	return u_time;
}

#ifdef NOT_NOW
/*************************************************************************
 *
 */

    /* days in the month */
static int Dom[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static char *Month_name[] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
	"Aug", "Sep", "Oct", "Nov", "Dec"};

void UTIMunix_to_date(time_t utime, UTIMstruct *date)
{
	long	u_day, day;

	u_day = UTIMjulian_date(1,1,1970);

	day = (utime / 86400);

	UTIMcalendar_date((u_day + day), &date->day, &date->month, &date->year);

	day = (utime % 86400);
	date->hour = day / 3600;
	date->min = (day / 60) - (date->hour * 60);
	date->sec = day % 60;
	date->unix_time = utime;

}
 
int UTIMdayofyear( UTIMstruct *date)
    /* return the day of the year for this date */
    {
	int i, doy;

	doy = date->day;
	for (i=0; i<date->month-1; i++)
	    doy += Dom[i];
	
	/* leap year */
	if (((date->year % 4) == 0) && (date->month > 2))
	   doy++;

	return doy;
    }

void UTIMmonthday( int year, int doy, int *month, int *day)
    /* return the month, day knowing the year and the day of year */
    {
	int m = 1;

	/* leap year adjustment */
	if ((year % 4) == 0) 
	    Dom[1] = 29;

	while (doy > Dom[m-1])
	    {
	    doy -= Dom[m-1];
	    m++;
	    }
	*day = doy;
	*month = m;

	Dom[1] = 28;
    }


/* given some # days since 0/0/0, format nice string */
char *UTIMday2str( double day)
    {
        static char result[100];
        char stemp[20];
        int year, month, daym;

        year = day / 365;
        day -= year*365;
	UTIMmonthday( 0, day, &month, &daym);	

        sprintf(result,"%s/%d/%d",Month_name[month-1],daym+1,year);

        return result;
    }

/* given some # seconds, format nice string */
char *UTIMsec2str( double time)
    {
        static char result[100];
        char stemp[20];
        long year, month, day, hour, min, sec;

        day = (long) floor(time / 86400.0);
        time -= day * 86400.0;

        year = day / 365;
        day -= year*365;
        month = day/30;
        day -= month*30;

        sec = (long) time;
        hour = sec / 3600;
        sec -= hour * 3600;
        min = sec / 60;

        result[0] = 0;
        if (year > 0)
            {
            sprintf(stemp,"%ld years ",year);
            strcat(result,stemp);
            }
        if (month > 0)
            {
            sprintf(stemp,"%ld months ",month);
            strcat(result,stemp);
            }
        if (day > 0)
            {
            sprintf(stemp,"%ld days ",day);
            strcat(result,stemp);
            }
        if (hour > 0)
            {
            sprintf(stemp,"%ld hours ",hour);
            strcat(result,stemp);
            }
        if (min > 0)
            {
            sprintf(stemp,"%ld minutes ",min);
            strcat(result,stemp);
            }

        return result;
    }

double UTIMtoTAI( UTIMstruct *date)
    /* convert to eos TAI format: (real continuous seconds since 12AM UTC 1-1-93)
       you must fill out all fields in *date except the unix_time.
     */
    {
	UTIMstruct tai_base;
	
	memset( &tai_base, 0, sizeof(UTIMstruct));
	tai_base.year = 1993;	
	tai_base.month = 1;
	tai_base.day = 1;
	UTIMdate_to_unix( &tai_base);

	UTIMdate_to_unix( date);
	
	/* LOOK: 6 leap seconds for 1999 !! */
	return (double) (date->unix_time - tai_base.unix_time + 6);
    }

time_t UTIM_TAI_to_unix( double tai)
    /* convert from eos TAI format: (real continuous seconds since 12AM UTC 1-1-93)
       to unix time
     */
    {
        UTIMstruct ut;

	ut.unix_time = (time_t) tai - 6;
	UTIMunix_to_date(  ut.unix_time, &ut);
	ut.year += 23;	/* offset from unix time */	
	return UTIMdate_to_unix( &ut);
    }

#endif

/********************************************************************/

unsigned long UTIMdiff( UTIMdiff_t *tdiff)
   /* returns time differences in millisecs.
      caller initializes tdiff->init to 0 to start the clock.
      subsequent calls returns the time difference since the start
    */      
    {
        struct timeval  *start_tv = &tdiff->start_tv;
        struct timeval  tv;
        long diff;

        if (tdiff->init == 0)
            {
            gettimeofday( start_tv, NULL);
	    tdiff->init = 1;
            return 0;
            }

        if (0 != gettimeofday( &tv, NULL))
            return 0;

        diff = (tv.tv_sec - start_tv->tv_sec) * 1000;
        diff += (tv.tv_usec - start_tv->tv_usec) / 1000;
        return (unsigned long) diff;
    }

unsigned long UTIMtime_diff( int start)
   /* returns time differences in millisecs.
      call with start = 1 to begin, then subsequent calls
      return the (wall clock) time difference since that starting
      point in millisecs.
    */
    {
	static UTIMdiff_t tdiff;

	if (start)
	    tdiff.init = 0;
	return UTIMdiff( &tdiff);
    }
