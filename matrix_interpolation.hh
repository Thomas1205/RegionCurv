#ifndef MATRIX_INTERPOLATION_HH
#define MATRIX_INTERPOLATION_HH

#include "matrix.hh"

template<typename T>
double bilinear_interpolation(const Math2D::Matrix<T>& matrix, double pos_x, double pos_y);


/********************** implementation ***********************/

template<typename T>
double bilinear_interpolation(const Math2D::Matrix<T>& matrix, double pos_x, double pos_y) {

  const int x_bound = matrix.xDim()-1;
  const int y_bound = matrix.yDim()-1;

  assert(x_bound > 0 && y_bound > 0);

  if (pos_x < 0)
    pos_x = 0;
  else if (pos_x > x_bound)
    pos_x = x_bound;
  
  if (pos_y < 0)
    pos_y = 0;
  else if (pos_y > y_bound)
    pos_y = y_bound;

  const int lx = (uint) (pos_x);
  const int ly = (uint) (pos_y);

  const int hx = (lx == x_bound) ? lx : lx + 1;
  const int hy = (ly == y_bound) ? ly : ly + 1;

  const double fac_x = hx - pos_x;
  const double fac_y = hy - pos_y;

  const double neg_fac_x = 1.0 - fac_x;
  const double neg_fac_y = 1.0 - fac_y;

  double val = matrix(lx,ly) * fac_x * fac_y
    + matrix(hx,ly) * neg_fac_x * fac_y
    + matrix(lx,hy) * fac_x * neg_fac_y
    + matrix(hx,hy) * neg_fac_x * neg_fac_y;
  
  return val;
}



#endif
