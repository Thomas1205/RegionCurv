/**** written by Thomas Schoenemann as a private person without employment, October 2009 ***/

#include "makros.hh"
#include "tensor.hh"

namespace Math3D {

  template<>
  float Tensor<float>::max() const {
    return Makros::max(Storage3D<float>::data_, Storage3D<float>::size_);
  }
  
  template<>
  float Tensor<float>::min() const {
    return Makros::min(Storage3D<float>::data_, Storage3D<float>::size_);
  }

  template<>
  void Tensor<float>::operator*=(const float scalar) {
    Makros::mul_array(Storage3D<float>::data_, Storage3D<float>::size_, scalar);
  }

  template<>
  void Tensor<double>::operator*=(const double scalar) {
    Makros::mul_array(Storage3D<double>::data_, Storage3D<double>::size_, scalar);
  }

}

