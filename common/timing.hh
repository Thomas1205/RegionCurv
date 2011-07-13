/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef TIMING_HH
#define TIMING_HH

#include <ctime>

#ifndef WIN32
#include <sys/time.h>
#else
typedef struct timeval {
  long tv_sec;
  long tv_usec;
} timeval;
#endif

double diff_mseconds(const timeval& end, const timeval& start);

double diff_seconds(const timeval& end, const timeval& start);

double diff_seconds(const std::clock_t& end, const std::clock_t& start);

#endif
