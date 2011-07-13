/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#include "vector.hh"

namespace Math1D {

  template<>
  double Vector<uint>::norm_l1() const {
    
    uint result = 0;
    for (size_t i=0; i < Storage1D<uint>::size_; i++) {
      result += Storage1D<uint>::data_[i];
    }
    
    return result;    
  }


  template<>
  double Vector<uchar>::norm_l1() const {
    
    uint result = 0;
    for (size_t i=0; i < Storage1D<uchar>::size_; i++) {
      result += Storage1D<uchar>::data_[i];
    }
    
    return result;    
  }


  template<>
  double Vector<ushort>::norm_l1() const {
    
    uint result = 0;
    for (size_t i=0; i < Storage1D<ushort>::size_; i++) {
      result += Storage1D<ushort>::data_[i];
    }
    
    return result;    
  }

  template<>
  float Vector<float>::max() const {
    return Makros::max(Storage1D<float>::data_, Storage1D<float>::size_);
  }

  template<>
  float Vector<float>::min() const {
    return Makros::min(Storage1D<float>::data_, Storage1D<float>::size_);
  }

  template<>
  void Vector<float>::operator*=(const float scalar) {
    Makros::mul_array(Storage1D<float>::data_, Storage1D<float>::size_, scalar);
  }

  template<>
  void Vector<double>::operator*=(const double scalar) {
    Makros::mul_array(Storage1D<double>::data_, Storage1D<double>::size_, scalar);
  }


}
