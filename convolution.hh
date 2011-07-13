/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef CONVOLUTE_HH
#define CONVOLUTE_HH

#include "vector.hh"
#include "matrix.hh"
#include "tensor.hh"

template <typename T>
void convolute(const Math1D::Vector<T>& signal, const Math1D::Vector<T>& kernel, uint kernel_zero_pos,
	       Math1D::Vector<T>& result);

template <typename T>
void convolute_2D_separable(Math2D::Matrix<T>& signal, const Math1D::Vector<T>& kernel, 
			    uint kernel_zero_pos);

template <typename T>
void convolute_2D_separable(Math3D::Tensor<T>& signal, const Math1D::Vector<T>& kernel, 
			    uint kernel_zero_pos);


template <typename T>
void convolute_3D_separable(Math3D::Tensor<T>& signal, const Math1D::Vector<T>& kernel, 
			    uint kernel_zero_pos);


/************************************** implementation *******************************************/


template <typename T>
void convolute(const Math1D::Vector<T>& signal, const Math1D::Vector<T>& kernel, uint kernel_zero_pos,
	       Math1D::Vector<T>& result) {

  int signal_size = signal.size(); 

  result.resize(signal.size());

  for (int x=0; x < signal_size; x++) {
    
    T sum = 0.0;
    for (int k=0; k < ((int) kernel.size()); k++) {

      int xx = x - k + kernel_zero_pos;
      if (xx < 0) //mirror conditions at the border
	xx = -(xx-1);
      if (xx >= signal_size)
	xx = xx - (xx - signal_size  + 1);
      sum += signal[xx]*kernel[k];
    }
    
    result[x] = sum;
  }
}


template <typename T>
void convolute_2D_separable(Math2D::Matrix<T>& signal, const Math1D::Vector<T>& kernel, 
			    uint kernel_zero_pos) {

  //horizontal convolution

  Math1D::Vector<T> temp1(signal.xDim());
  Math1D::Vector<T> temp2(signal.xDim());
  for (uint y=0; y < signal.yDim(); y++) {

    for (uint x=0; x < signal.xDim(); x++) {
      temp1[x] = signal(x,y);
    }

    convolute(temp1,kernel,kernel_zero_pos,temp2);

    for (uint x=0; x < signal.xDim(); x++) {
      signal(x,y) = temp2[x];
    }
  }

  //vertical convolution
  temp1.resize(signal.yDim());
  temp2.resize(signal.yDim());
  for (uint x=0; x < signal.xDim(); x++) {
    
    for (uint y=0; y < signal.yDim(); y++)
      temp1[y] = signal(x,y);
    
    convolute(temp1,kernel,kernel_zero_pos,temp2);
    
    for (uint y=0; y < signal.yDim(); y++)
      signal(x,y) = temp2[y];
  }
}


template <typename T>
void convolute_2D_separable(Math3D::Tensor<T>& signal, const Math1D::Vector<T>& kernel, 
			    uint kernel_zero_pos) {

  for (uint z=0; z < signal.zDim(); z++) {

    //horizontal convolution
    Math1D::Vector<T> temp1(signal.xDim());
    Math1D::Vector<T> temp2(signal.xDim());
    for (uint y=0; y < signal.yDim(); y++) {
      
      for (uint x=0; x < signal.xDim(); x++) {
	temp1[x] = signal(x,y,z);
      }
      
      convolute(temp1,kernel,kernel_zero_pos,temp2);
      
      for (uint x=0; x < signal.xDim(); x++) {
	signal(x,y,z) = temp2[x];
      }
    }
    
    //vertical convolution
    temp1.resize(signal.yDim());
    temp2.resize(signal.yDim());
    for (uint x=0; x < signal.xDim(); x++) {
      
      for (uint y=0; y < signal.yDim(); y++)
	temp1[y] = signal(x,y,z);
      
      convolute(temp1,kernel,kernel_zero_pos,temp2);
      
      for (uint y=0; y < signal.yDim(); y++)
	signal(x,y,z) = temp2[y];
    }
  }
}
template <typename T>
void convolute_3D_separable(Math3D::Tensor<T>& signal, const Math1D::Vector<T>& kernel, 
			    uint kernel_zero_pos) {
  
  //horizontal convolution
  Math1D::Vector<T> temp1(signal.xDim());
  Math1D::Vector<T> temp2(signal.xDim());
  for (uint z=0; z < signal.zDim(); z++) {
    for (uint y=0; y < signal.yDim(); y++) {
      
      for (uint x=0; x < signal.xDim(); x++) {
	temp1[x] = signal(x,y,z);
      }
      
      convolute(temp1,kernel,kernel_zero_pos,temp2);
      
      for (uint x=0; x < signal.xDim(); x++) {
	signal(x,y,z) = temp2[x];
      }
    }
  }
  
  //vertical convolution
  temp1.resize(signal.yDim());
  temp2.resize(signal.yDim());
  
  for (uint z=0; z < signal.zDim(); z++) {
    
    for (uint x=0; x < signal.xDim(); x++) {
      
      for (uint y=0; y < signal.yDim(); y++) {
	temp1[y] = signal(x,y,z);
      }

      convolute(temp1,kernel,kernel_zero_pos,temp2);
	
      for (uint y=0; y < signal.yDim(); y++)
	signal(x,y,z) = temp2[y];
    }
  }

  //depth convolution
  temp1.resize(signal.zDim());
  temp2.resize(signal.zDim());
  
  for (uint y=0; y < signal.yDim(); y++) {
      
    for (uint x=0; x < signal.xDim(); x++) {

      for (uint z=0; z < signal.zDim(); z++) 
	temp1[z] = signal(x,y,z);

      convolute(temp1,kernel,kernel_zero_pos,temp2);

      for (uint z=0; z < signal.zDim(); z++) 
	signal(x,y,z) = temp1[z];
    }
  }
}



#endif
