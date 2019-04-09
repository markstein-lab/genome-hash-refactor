//copyright (c) in silico Labs, LLC 2006

/*
 * This is a timer interface for a millisecond resolution
 * timer that works on both windows and unix.
 *
 * use as in:
 *
 *	 TIMEVAL time1, time2;
 *	 getTime(&time1);
 *	 < do something that takes awhile >
 *	 getTime(&time2);
 *	 cout << "That took: " << setprecision(3) <<
 *			 (getDiffMillisecs(&time1, &time2)/1000.0) <<
 *			 " sec." << endl;
 *
 */

#include "timer.h"

void getTime(TIMEVAL *t)
{
#ifdef WIN32
		clock_t ticks = clock();
		t->tv_sec = ticks / CLOCKS_PER_SEC;
		ticks -= t->tv_sec * CLOCKS_PER_SEC;
		t->tv_usec = ticks * (1000000 / CLOCKS_PER_SEC);
#else
		gettimeofday(t, (void *) 0);
#endif
}

int getDiffMillisecs(TIMEVAL *first, TIMEVAL *second)
{
		return (second->tv_sec - first->tv_sec) * 1000 +
				(second->tv_usec - first->tv_usec) / 1000;
}
