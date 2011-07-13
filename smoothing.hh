/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef SMOOTHING_HH
#define SMOOTHING_HH

#include "tensor.hh"
#include "vector.hh"
#include "convolution.hh"
#include "matrix.hh"

template<typename T>
void smooth_hor_binomial(Math3D::Tensor<T>& image);

template<typename T>
void smooth_vert_binomial(Math3D::Tensor<T>& image);

template<typename T>
void smooth_binomial(Math3D::Tensor<T>& image);

template<typename T>
void smooth_isotropic_gauss(Math2D::Matrix<T>& image, double sigma);

template<typename T>
void smooth_isotropic_gauss(Math3D::Tensor<T>& image, double sigma);

template<typename T>
void vol_smooth_isotropic_gauss(Math3D::Tensor<T>& image, double sigma);

/*********** implementation ************/
template<typename T>
void smooth_hor_binomial(Math3D::Tensor<T>& image) {

  uint xDim = image.xDim();
  uint yDim = image.yDim();
  uint zDim = image.zDim();

  for (uint z=0; z < zDim; z++) {

    for (uint y=0; y < yDim; y++) {

      T last = image(0,y,z);
      
      for (uint x=0; x < xDim; x++) {
	
	T cur = image(x,y,z);
	T next = ((x+1) < xDim) ? image(x,y,z) : cur;
	
	image(x,y,z) = 0.25 * (last + 2.0*cur + next);
	
	last = cur;
      }
    }
  }
}


template<typename T>
void smooth_vert_binomial(Math3D::Tensor<T>& image) {

  uint xDim = image.xDim();
  uint yDim = image.yDim();
  uint zDim = image.zDim();

  for (uint z=0; z < zDim; z++) {

    for (uint x=0; x < xDim; x++) {

      T last = image(x,0,z);
      
      for (uint y=0; y < yDim; y++) {
	
	T cur = image(x,y,z);
	T next = ((y+1) < yDim) ? image(x,y+1,z) : cur;
	
	image(x,y,z) = 0.25 * (last + 2.0*cur + next);
	
	last = cur;
      }
    }
  }

}


template<typename T>
void smooth_binomial(Math3D::Tensor<T>& image) {

  smooth_hor_binomial(image);
  smooth_vert_binomial(image);
}


template<typename T>
void smooth_isotropic_gauss(Math3D::Tensor<T>& image, double sigma) {

  uint width = (uint) ceil(sigma * 2.14); 
  uint size = 2*width+1;

  double variance = sigma*sigma;

  T norm = 1.0;
  Math1D::Vector<T> gauss_kernel(size);
  gauss_kernel[width] = 1.0;
  for (uint i=1; i <= width; i++) {
    T val = std::exp((-1.0*i*i)/variance);
    //std::cerr << "val: " << val << std::endl;
    norm += 2.0 * val;
    gauss_kernel[width-i] = val;
    gauss_kernel[width+i] = val;
  }

  for (uint i=0; i < size; i++)
    gauss_kernel[i] /= norm;

  //std::cerr << "kernel: " << gauss_kernel << std::endl; 

  convolute_2D_separable(image,gauss_kernel,width);
}

template<typename T>
void smooth_isotropic_gauss(Math2D::Matrix<T>& image, double sigma) {

  uint width = (uint) ceil(sigma * 2.14); 
  uint size = 2*width+1;

  double variance = sigma*sigma;

  T norm = 1.0;
  Math1D::Vector<T> gauss_kernel(size);
  gauss_kernel[width] = 1.0;
  for (uint i=1; i <= width; i++) {
    T val = std::exp((-1.0*i*i)/variance);
    //std::cerr << "val: " << val << std::endl;
    norm += 2.0 * val;
    gauss_kernel[width-i] = val;
    gauss_kernel[width+i] = val;
  }

  for (uint i=0; i < size; i++)
    gauss_kernel[i] /= norm;

  //std::cerr << "kernel: " << gauss_kernel << std::endl; 

  convolute_2D_separable(image,gauss_kernel,width);
}

template<typename T>
void vol_smooth_isotropic_gauss(Math3D::Tensor<T>& image, double sigma) {

  uint width = (uint) ceil(sigma * 2.14); 
  uint size = 2*width+1;

  double variance = sigma*sigma;

  T norm = 1.0;
  Math1D::Vector<T> gauss_kernel(size);
  gauss_kernel[width] = 1.0;
  for (uint i=1; i <= width; i++) {
    T val = std::exp((-1.0*i*i)/variance);
    //std::cerr << "val: " << val << std::endl;
    norm += 2.0 * val;
    gauss_kernel[width-i] = val;
    gauss_kernel[width+i] = val;
  }

  for (uint i=0; i < size; i++)
    gauss_kernel[i] /= norm;

  //std::cerr << "kernel: " << gauss_kernel << std::endl; 

  convolute_3D_separable(image,gauss_kernel,width);
}

#endif

