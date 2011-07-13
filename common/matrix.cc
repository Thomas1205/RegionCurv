/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#include "matrix.hh"

namespace Math2D {

  template<>   
  double Matrix<uint>::norm_l1() const {
    
    uint result = 0;
    for (size_t i=0; i < Storage2D<uint>::size_; i++) {
      result += Storage2D<uint>::data_[i];
    }
    
    return result;    
  }

  template<>   
  double Matrix<ushort>::norm_l1() const {
    
    uint result = 0;
    for (size_t i=0; i < Storage2D<ushort>::size_; i++) {
      result += Storage2D<ushort>::data_[i];
    }
    
    return result;    
  }

  template<>   
  double Matrix<uchar>::norm_l1() const {
    
    uint result = 0;
    for (size_t i=0; i < Storage2D<uchar>::size_; i++) {
      result += Storage2D<uchar>::data_[i];
    }
    
    return result;    
  }

  template<>
  float Matrix<float>::max() const {
    return Makros::max(Storage2D<float>::data_, Storage2D<float>::size_);
  }

  template<>
  float Matrix<float>::min() const {
    return Makros::min(Storage2D<float>::data_, Storage2D<float>::size_);
  }

  template<>
  void Matrix<float>::operator*=(const float scalar) {
    Makros::mul_array(Storage2D<float>::data_, Storage2D<float>::size_, scalar);
  }

  template<>
  void Matrix<double>::operator*=(const double scalar) {
    Makros::mul_array(Storage2D<double>::data_, Storage2D<double>::size_, scalar);
  }


}
