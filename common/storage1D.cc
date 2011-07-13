/*** written by Thomas Schoenemann as a private person without employment, December 2009 ***/

#include "storage1D.hh"

/** template specializations for the copy operator **/
template<>
void Storage1D<uint>::operator=(const Storage1D<uint>& toCopy) {

  if (size_ != toCopy.size()) {
    
    if (data_ != 0)
      delete[] data_;
    
    size_ = toCopy.size();
    data_ = new uint[size_];
  }

  memcpy(data_,toCopy.direct_access(),size_*sizeof(uint));
}


template<>
void Storage1D<int>::operator=(const Storage1D<int>& toCopy) {

  if (size_ != toCopy.size()) {
    
    if (data_ != 0)
      delete[] data_;
    
    size_ = toCopy.size();
    data_ = new int[size_];
  }

  memcpy(data_,toCopy.direct_access(),size_*sizeof(int));
}

template<>
void Storage1D<float>::operator=(const Storage1D<float>& toCopy) {

  if (size_ != toCopy.size()) {
    
    if (data_ != 0)
      delete[] data_;
    
    size_ = toCopy.size();
    data_ = new float[size_];
  }

  memcpy(data_,toCopy.direct_access(),size_*sizeof(float));
}

template<>
void Storage1D<double>::operator=(const Storage1D<double>& toCopy) {

  if (size_ != toCopy.size()) {
    
    if (data_ != 0)
      delete[] data_;
    
    size_ = toCopy.size();
    data_ = new double[size_];
  }

  memcpy(data_,toCopy.direct_access(),size_*sizeof(double));
}

template<>
void Storage1D<long double>::operator=(const Storage1D<long double>& toCopy) {

  if (size_ != toCopy.size()) {
    
    if (data_ != 0)
      delete[] data_;
    
    size_ = toCopy.size();
    data_ = new long double[size_];
  }

  memcpy(data_,toCopy.direct_access(),size_*sizeof(long double));
}


template<>
Storage1D<int>::Storage1D(const Storage1D<int>& toCopy) {

  size_ = toCopy.size();
  data_ = new int[size_];

  memcpy(data_,toCopy.direct_access(),size_*sizeof(int));
}

template<>
Storage1D<uint>::Storage1D(const Storage1D<uint>& toCopy) {

  size_ = toCopy.size();
  data_ = new uint[size_];

  memcpy(data_,toCopy.direct_access(),size_*sizeof(uint));
}

template<>
Storage1D<float>::Storage1D(const Storage1D<float>& toCopy) {

  size_ = toCopy.size();
  data_ = new float[size_];

  memcpy(data_,toCopy.direct_access(),size_*sizeof(float));
}

template<>
Storage1D<double>::Storage1D(const Storage1D<double>& toCopy) {

  size_ = toCopy.size();
  data_ = new double[size_];

  memcpy(data_,toCopy.direct_access(),size_*sizeof(double));
}
