/*** written by Thomas Schoenemann as an employee of Lund University, August 2010 ***/

#include "storage2D.hh"
#include <cstring>

template<>
void Storage2D<int>::operator=(const Storage2D<int>& toCopy) {

  if (size_ != toCopy.size()) {
    if (data_ != 0)
      delete[] data_;

    size_ = toCopy.size();
    data_ = new int[size_];
  }

  xDim_ = toCopy.xDim();
  yDim_ = toCopy.yDim();
  assert(size_ == xDim_*yDim_);

  memcpy(data_,toCopy.direct_access(),size_*sizeof(int));
}

template<>
void Storage2D<uint>::operator=(const Storage2D<uint>& toCopy) {

  if (size_ != toCopy.size()) {
    if (data_ != 0)
      delete[] data_;

    size_ = toCopy.size();
    data_ = new uint[size_];
  }

  xDim_ = toCopy.xDim();
  yDim_ = toCopy.yDim();
  assert(size_ == xDim_*yDim_);

  memcpy(data_,toCopy.direct_access(),size_*sizeof(uint));
}

template<>
void Storage2D<float>::operator=(const Storage2D<float>& toCopy) {

  if (size_ != toCopy.size()) {
    if (data_ != 0)
      delete[] data_;

    size_ = toCopy.size();
    data_ = new float[size_];
  }

  xDim_ = toCopy.xDim();
  yDim_ = toCopy.yDim();
  assert(size_ == xDim_*yDim_);

  memcpy(data_,toCopy.direct_access(),size_*sizeof(float));
}

template<>
void Storage2D<double>::operator=(const Storage2D<double>& toCopy) {

  if (size_ != toCopy.size()) {
    if (data_ != 0)
      delete[] data_;

    size_ = toCopy.size();
    data_ = new double[size_];
  }

  xDim_ = toCopy.xDim();
  yDim_ = toCopy.yDim();
  assert(size_ == xDim_*yDim_);

  memcpy(data_,toCopy.direct_access(),size_*sizeof(double));
}

template<>
void Storage2D<long double>::operator=(const Storage2D<long double>& toCopy) {

  if (size_ != toCopy.size()) {
    if (data_ != 0)
      delete[] data_;

    size_ = toCopy.size();
    data_ = new long double[size_];
  }

  xDim_ = toCopy.xDim();
  yDim_ = toCopy.yDim();
  assert(size_ == xDim_*yDim_);

  memcpy(data_,toCopy.direct_access(),size_*sizeof(long double));
}
