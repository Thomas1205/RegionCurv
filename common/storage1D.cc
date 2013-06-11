/*** first version written by Thomas Schoenemann as a private person without employment, September 2009 ***/
/*** much refined by Thomas Schoenemann  at Lund University, Sweden, the University of Pisa, Italy, ***
 *** and the University of Düsseldorf, Germany 2010 - 2012 **/
/*** if you desire the checked version, make sure your compiler defines the option SAFE_MODE on the command line ***/

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
void Storage1D<ushort>::operator=(const Storage1D<ushort>& toCopy) {

  if (size_ != toCopy.size()) {
    
    if (data_ != 0)
      delete[] data_;
    
    size_ = toCopy.size();
    data_ = new ushort[size_];
  }

  memcpy(data_,toCopy.direct_access(),size_*sizeof(ushort));
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
Storage1D<ushort>::Storage1D(const Storage1D<ushort>& toCopy) {

  size_ = toCopy.size();
  data_ = new ushort[size_];

  memcpy(data_,toCopy.direct_access(),size_*sizeof(ushort));
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


//maintains the values of existing positions, new ones are undefined 
template<>
void Storage1D<float>::resize(size_t new_size) {

  if (data_ == 0) {
    data_ = new float[new_size];
  }
  else if (size_ != new_size) {
    float* new_data = new float[new_size];

    memcpy(new_data,data_,std::min(size_,new_size)*sizeof(float));
    
    delete[] data_;
    data_ = new_data;
  }

  size_ = new_size;
}

template<>
void Storage1D<double>::resize(size_t new_size) {

  if (data_ == 0) {
    data_ = new double[new_size];
  }
  else if (size_ != new_size) {
    double* new_data = new double[new_size];

    memcpy(new_data,data_,std::min(size_,new_size)*sizeof(double));
    
    delete[] data_;
    data_ = new_data;
  }

  size_ = new_size;
}

template<>
void Storage1D<long double>::resize(size_t new_size) {

  if (data_ == 0) {
    data_ = new long double[new_size];
  }
  else if (size_ != new_size) {
    long double* new_data = new long double[new_size];

    memcpy(new_data,data_,std::min(size_,new_size)*sizeof(long double));
    
    delete[] data_;
    data_ = new_data;
  }

  size_ = new_size;
}

template<>
void Storage1D<int>::resize(size_t new_size) {

  if (data_ == 0) {
    data_ = new int[new_size];
  }
  else if (size_ != new_size) {
    int* new_data = new int[new_size];

    memcpy(new_data,data_,std::min(size_,new_size)*sizeof(int));
    
    delete[] data_;
    data_ = new_data;
  }

  size_ = new_size;
}

template<>
void Storage1D<uint>::resize(size_t new_size) {

  if (data_ == 0) {
    data_ = new uint[new_size];
  }
  else if (size_ != new_size) {
    uint* new_data = new uint[new_size];

    memcpy(new_data,data_,std::min(size_,new_size)*sizeof(uint));
    
    delete[] data_;
    data_ = new_data;
  }

  size_ = new_size;
}


//maintains the values of existing positions, new ones are filled with <code> fill_value </code>
template<>
void Storage1D<float>::resize(size_t new_size, float fill_value) {

  if (data_ == 0) {
    data_ = new float[new_size];
    for (size_t i=0; i < new_size; i++)
      data_[i] = fill_value;
  }
  else if (size_ != new_size) {
    float* new_data = new float[new_size];
    // for (size_t i=size_; i < new_size; i++)
    //   new_data[i] = fill_value;
    if (new_size > size_)
      std::fill_n(new_data+size_,new_size-size_,fill_value);


    memcpy(new_data,data_,std::min(size_,new_size)*sizeof(float));
    
    delete[] data_;
    data_ = new_data;
  }

  size_ = new_size;
}

template<>
void Storage1D<double>::resize(size_t new_size, double fill_value) {

  if (data_ == 0) {
    data_ = new double[new_size];
    for (size_t i=0; i < new_size; i++)
      data_[i] = fill_value;
  }
  else if (size_ != new_size) {
    double* new_data = new double[new_size];
    // for (size_t i=size_; i < new_size; i++)
    //   new_data[i] = fill_value;
    if (new_size > size_)
      std::fill_n(new_data+size_,new_size-size_,fill_value);


    memcpy(new_data,data_,std::min(size_,new_size)*sizeof(double));
    
    delete[] data_;
    data_ = new_data;
  }

  size_ = new_size;
}

template<>
void Storage1D<long double>::resize(size_t new_size, long double fill_value) {

  if (data_ == 0) {
    data_ = new long double[new_size];
    for (size_t i=0; i < new_size; i++)
      data_[i] = fill_value;
  }
  else if (size_ != new_size) {
    long double* new_data = new long double[new_size];
    // for (size_t i=size_; i < new_size; i++)
    //   new_data[i] = fill_value;
    if (new_size > size_)
      std::fill_n(new_data+size_,new_size-size_,fill_value);

    memcpy(new_data,data_,std::min(size_,new_size)*sizeof(long double));
    
    delete[] data_;
    data_ = new_data;
  }

  size_ = new_size;
}

template<>
void Storage1D<int>::resize(size_t new_size, int fill_value) {

  if (data_ == 0) {
    data_ = new int[new_size];
    for (size_t i=0; i < new_size; i++)
      data_[i] = fill_value;
  }
  else if (size_ != new_size) {
    int* new_data = new int[new_size];
    // for (size_t i=size_; i < new_size; i++)
    //   new_data[i] = fill_value;
    if (new_size > size_)
      std::fill_n(new_data+size_,new_size-size_,fill_value);

    memcpy(new_data,data_,std::min(size_,new_size)*sizeof(int));
    
    delete[] data_;
    data_ = new_data;
  }

  size_ = new_size;
}

template<>
void Storage1D<uint>::resize(size_t new_size, uint fill_value) {

  if (data_ == 0) {
    data_ = new uint[new_size];
    for (size_t i=0; i < new_size; i++)
      data_[i] = fill_value;
  }
  else if (size_ != new_size) {
    uint* new_data = new uint[new_size];
    // for (size_t i=size_; i < new_size; i++)
    //   new_data[i] = fill_value;
    if (new_size > size_)
      std::fill_n(new_data+size_,new_size-size_,fill_value);

    memcpy(new_data,data_,std::min(size_,new_size)*sizeof(uint));
    
    delete[] data_;
    data_ = new_data;
  }

  size_ = new_size;
}



template<>
FlexibleStorage1D<uint>::FlexibleStorage1D(const FlexibleStorage1D<uint>& toCopy) {

  size_ = toCopy.size();
  reserved_size_ = toCopy.reserved_size();
  
  data_ = new uint[reserved_size_];

  memcpy(data_,toCopy.direct_access(),size_*sizeof(uint));
}

template<>
void FlexibleStorage1D<uint>::operator=(const FlexibleStorage1D<uint>& toCopy) {
  
    uint new_res = toCopy.reserved_size();
  if (new_res != reserved_size_) {
    reserved_size_ = new_res;

    if (data_ != 0)
      delete[] data_;
    data_ = new uint[reserved_size_];
  }

  size_ = toCopy.size();

  memcpy(data_,toCopy.direct_access(),size_*sizeof(uint));
}
