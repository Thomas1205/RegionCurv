/******* partiallly written by Thomas Schoenemann during his Ph.D. thesis at the University of Bonn ********/
/******* the other part was written by Thomas Schoenemann as an employee of Lund University, Sweden, 2011 *******/

#include "color_conversion.hh"

void yuv2rgba(uchar y, uchar u, uchar v, uchar& r, uchar& g, uchar& b) {

    double dy = (double) y;
    double du = (double) u;    
    du -= 128;
    du *= 224.0 / 255.0;
    double dv = (double) v;
    dv -= 128;
    dv *= 314.0 / 255.0;
    //    v -= 156;

    // properly implemented ?
    int red = ((uint) (dy + dv * 1.14));
    if (red > 255)
	red = 255;
    if (red < 0)
	red = 0;
    int green = ((int) ( dy - du * 0.394 - dv * 0.581));
    if (green < 0)
	green = 0;
    if (green > 255)
	green = 255;
    int blue = ((int)  (dy + du * 2.028) ) ;
    if (blue > 255)
	blue = 255;
    if (blue < 0)
	blue = 0;

    r = red;
    g = green;
    b = blue;
}

void rgb2yuv(uchar r, uchar g, uchar b, uchar& y, uchar& u, uchar& v) {

  double dr = r;
  double dg = g;
  double db = b;

  double dy =  (0.299*dr + 0.587*dg + 0.114*db);
  if (dy < 0.0)
    dy = 0.0;
  if (dy > 255.0)
    dy = 255.0;

  y = (uchar) dy;
  
  double du = (db - dy) * 0.493 + 128.0;
  du /= 224.0 / 255.0;
  double dv = (dr - dy) * 0.877 + 128.0;
  dv /= 314.0 / 255.0;

  if (du < 0.0)
    du = 0.0;
  if (du > 255.0)
    du = 255.0;

  u = (uchar) du;

  if (dv < 0.0)
    dv = 0.0;
  if (dv > 255.0)
    dv = 255.0;

  v = (uchar) dv;
}
