#ifndef _TIMER_INCLUDED_
#define _TIMER_INCLUDED_

#ifdef WIN32
#include <time.h>
typedef struct timeval {
		unsigned long tv_sec;
		long tv_usec;
} TIMEVAL;
#else
/* #include <sys/time.h> */
/* typedef struct timeval TIMEVAL; */
#include <time.h>
typedef struct timespec TIMEVAL;
#endif

#if defined(__cplusplus)
extern "C"
{
#endif

void getTime(TIMEVAL *t);
int getDiffMillisecs(TIMEVAL *t1, TIMEVAL *t2);
long getDiffNanosecs(TIMEVAL *t1, TIMEVAL *t2);

#if defined(__cplusplus)
}
#endif

#endif
