/*** written by Thomas Schoenemann as a private person without employment, October 2009 ***/

#ifndef SAMPLING_HH
#define SAMPLING_HH

#include "matrix.hh"
#include "tensor.hh"

#include "matrix_interpolation.hh"
#include "tensor_interpolation.hh"

template <typename T>
void upsample_matrix(const Math2D::Matrix<T>& source, Math2D::Matrix<T>& target);

template<typename T>
void downsample_matrix(const Math2D::Matrix<T>& source, Math2D::Matrix<T>& target);

template <typename T>
void upsample_tensor(const Math3D::Tensor<T>& source, Math3D::Tensor<T>& target);

template<typename T>
void downsample_tensor(const Math3D::Tensor<T>& source, Math3D::Tensor<T>& target);


/************************************** implementation **************************************/

template <typename T>
void upsample_matrix(const Math2D::Matrix<T>& source, Math2D::Matrix<T>& target) {


  assert(target.xDim() >= source.xDim());
  assert(target.yDim() >= source.yDim());

  const double scale_x = ((double) source.xDim()) / ((double) target.xDim());
  const double scale_y = ((double) source.yDim()) / ((double) target.yDim());

  const uint x_bound = source.xDim() - 1;
  const uint y_bound = source.yDim() - 1;

  for (uint y=0; y < target.yDim(); y++) {

    double pos_y = y*scale_y;
    const uint ly = (uint) (pos_y);
    const uint hy = (ly == y_bound) ? ly : ly + 1;
    const double fac_y = hy - pos_y;
    const double neg_fac_y = 1.0 - fac_y;

    for (uint x=0; x < target.xDim(); x++) {

      double pos_x = x*scale_x;
      const uint lx = (uint) (pos_x);
      const uint hx = (lx == x_bound) ? lx : lx + 1;
      
      const double fac_x = hx - pos_x;
      const double neg_fac_x = 1.0 - fac_x;

      target(x,y) = (source(lx,ly) * fac_x
		     + source(hx,ly) * neg_fac_x) * fac_y
	+ (source(lx,hy) * fac_x
	   + source(hx,hy) * neg_fac_x) * neg_fac_y;
    }
  }
}


template <typename T>
void upsample_tensor(const Math3D::Tensor<T>& source, Math3D::Tensor<T>& target) {

  assert(target.xDim() >= source.xDim());
  assert(target.yDim() >= source.yDim());
  assert(target.zDim() == source.zDim());

  double scale_x = ((double) source.xDim()) / ((double) target.xDim());
  double scale_y = ((double) source.yDim()) / ((double) target.yDim());

  const uint x_bound = source.xDim() - 1;
  const uint y_bound = source.yDim() - 1;

  for (uint y=0; y < target.yDim(); y++) {

    double pos_y = y*scale_y;
    const uint ly = (uint) (pos_y);
    const uint hy = (ly == y_bound) ? ly : ly + 1;
    const double fac_y = hy - pos_y;
    const double neg_fac_y = 1.0 - fac_y;

    for (uint x=0; x < target.xDim(); x++) {

      double pos_x = x*scale_x;
      const uint lx = (uint) (pos_x);
      const uint hx = (lx == x_bound) ? lx : lx + 1;
      
      const double fac_x = hx - pos_x;
      const double neg_fac_x = 1.0 - fac_x;
      
      for (uint z=0; z < target.zDim(); z++) {

	target(x,y,z) = (source(lx,ly,z) * fac_x
			 + source(hx,ly,z) * neg_fac_x) * fac_y
	  + (source(lx,hy,z) * fac_x
	     + source(hx,hy,z) * neg_fac_x) * neg_fac_y;
      }
    }
  }


}


template<typename T>
void downsample_matrix(const Math2D::Matrix<T>& source, Math2D::Matrix<T>& target) {

  const uint txDim = uint(target.xDim());
  const uint tyDim = uint(target.yDim());
  
  const double sxDim = double(source.xDim());
  const double syDim = double(source.yDim());
  
  assert(txDim <= sxDim);
  assert(tyDim <= syDim);
  
  const double scale_x = ((double) sxDim) / ((double) txDim);
  const double scale_y = ((double) syDim) / ((double) tyDim);
  const double pixel_area = scale_x*scale_y; 

  for (uint y=0; y < tyDim; y++) {

    const double ly = y*scale_y;
    const double hy = std::min(ly + scale_y,syDim-1.0);
    const uint uly = (uint) ly;
    const uint uhy = (uint) hy;

    for (uint x=0; x < txDim; x++) {

      const double lx = x*scale_x;
      const double hx = std::min(lx + scale_x,sxDim-1.0);
      const uint ulx = (uint) lx;
      const uint uhx = (uint) hx;

      double sum = 0.0;
      double weight_sum = 0.0;

      for (uint yy= uly; yy <= uhy; yy++) {
	for (uint xx= ulx; xx <= uhx; xx++) {
	  
	  double len_x = std::min(1.0,xx+1-lx);
	  len_x = std::min(len_x,hx-xx);
	  
	  double len_y = std::min(1.0,yy+1-ly);
	  len_y = std::min(len_y,hy-yy);
	  
	  sum += source(xx,yy) * len_x * len_y;
	  weight_sum += len_x*len_y;
	}
      }
      
      if (x+1 < txDim && y+1 < tyDim)
	assert(std::abs(weight_sum-pixel_area) < 1e-4);
      
      target(x,y) = T(sum / weight_sum);
    }
  }

}


template<typename T>
void downsample_tensor(const Math3D::Tensor<T>& source, Math3D::Tensor<T>& target) {

  const uint txDim = target.xDim();
  const uint tyDim = target.yDim();
  
  const double sxDim = source.xDim();
  const double syDim = source.yDim();

  const uint zDim = target.zDim();
  
  assert(txDim <= sxDim);
  assert(tyDim <= syDim);
  assert(zDim == source.zDim());
  
  const double scale_x = ((double) sxDim) / ((double) txDim);
  const double scale_y = ((double) syDim) / ((double) tyDim);
  const double pixel_area = scale_x*scale_y; 

  for (uint y=0; y < tyDim; y++) {

    const double ly = y*scale_y;
    const double hy = std::min(ly + scale_y,syDim-1.0);
    const uint uly = (uint) ly;
    const uint uhy = (uint) hy;

    for (uint x=0; x < txDim; x++) {

      const double lx = x*scale_x;
      const double hx = std::min(lx + scale_x,sxDim-1.0);
      const uint ulx = (uint) lx;
      const uint uhx = (uint) hx;

      double sum = 0.0;
      double weight_sum = 0.0;

      for (uint z=0; z < zDim; z++) {
	for (uint yy= uly; yy <= uhy; yy++) {
	  for (uint xx= ulx; xx <= uhx; xx++) {
	    
	    double len_x = std::min(1.0,xx+1-lx);
	    len_x = std::min(len_x,hx-xx);
	    
	    double len_y = std::min(1.0,yy+1-ly);
	    len_y = std::min(len_y,hy-yy);
	    
	    sum += source(xx,yy,z) * len_x * len_y;
	    weight_sum += len_x*len_y;
	  }
	}

	if (x+1 < txDim && y+1 < tyDim)
	  assert(std::abs(weight_sum-pixel_area) < 1e-4);

	target(x,y,z) = sum / weight_sum;
      }
    }
  }
}



#endif
