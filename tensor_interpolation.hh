#ifndef TENSOR_INTERPOLATION_HH
#define TENSOR_INTERPOLATION_HH

//interpolation in the x-y plane, but NOT in z-direction
template<typename T>
double bilinear_interpolation(const Math3D::Tensor<T>& tensor, double pos_x, double pos_y, size_t z);



/********************** implementation ***********************/

template<typename T>
double bilinear_interpolation(const Math3D::Tensor<T>& tensor, double pos_x, double pos_y, size_t z) {

  assert(z < tensor.zDim());

  const int x_bound = tensor.xDim()-1;
  const int y_bound = tensor.yDim()-1;

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

  double val = tensor(lx,ly,z) * fac_x * fac_y
    + tensor(hx,ly,z) * neg_fac_x * fac_y
    + tensor(lx,hy,z) * fac_x * neg_fac_y
    + tensor(hx,hy,z) * neg_fac_x * neg_fac_y;
  
  return val;
}




#endif
