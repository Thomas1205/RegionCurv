/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#include "timing.hh"
#include <algorithm>
#include <iostream>

double diff_mseconds(const timeval& end, const timeval& start) {

//   std::cerr << "end= " << end.tv_usec << " micro + " << end.tv_sec << " sec" << std::endl;
//   std::cerr << "start= " << start.tv_usec << " micro + " << start.tv_sec << " sec" << std::endl;

  double end_mseconds = (end.tv_usec * 0.001)  + (end.tv_sec * 1000);
  double start_mseconds = (start.tv_usec * 0.001) + (start.tv_sec * 1000);

  return end_mseconds - start_mseconds;
}

double diff_seconds(const timeval& end, const timeval& start) {
  return diff_mseconds(end,start) / 1000.0;
}


double diff_seconds(const std::clock_t& end, const std::clock_t& start) {
  return ((double) (end - start)) / ((double) (CLOCKS_PER_SEC));
}
