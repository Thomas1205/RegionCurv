/******* written by Thomas Schoenemann during his Ph.D. thesis at the University of Bonn ********/

#ifndef COLOR_CONVERSION_HH
#define COLOR_CONVERSION_HH

#include "makros.hh"

void yuv2rgba(uchar y, uchar u, uchar v, uchar& r, uchar& g, uchar& b);

void rgb2yuv(uchar r, uchar g, uchar b, uchar& y, uchar& u, uchar& v);

#endif
