/* $Id$ */
#ifndef UTIM_WAS_INCLUDED
#define UTIM_WAS_INCLUDED

#include <time.h>
#include <sys/time.h>

typedef struct UTIMdiff_t {
        int     init;
        struct timeval start_tv;
} UTIMdiff_t;

extern time_t UTIMdate2unix( int year, int month, int day, int hour, int min, 
	int sec);

extern long UTIMjulian_date(int day, int month, int year);

extern void UTIMcalendar_date(long jdate,long *day, long *month, long *year);

#endif
