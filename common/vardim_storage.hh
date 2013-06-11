/**** written by Thomas Schoenemann as a private person without employment, July 2011 ****/

#ifndef VAR_DIM_STORAGE_HH
#define VAR_DIM_STORAGE_HH

#include "vector.hh"

template <typename T>
class VarDimStorage {
public:
  VarDimStorage(const Math1D::Vector<size_t>& dim);

  VarDimStorage(const Math1D::Vector<size_t>& dim, T fill);

  //copy constructor
  VarDimStorage(const VarDimStorage& toCopy);

  ~VarDimStorage();

  size_t dim(uint n) const;
  
  size_t size() const;

  size_t nDims() const;

  const Math1D::Vector<size_t>& dim_vector() const;

  const T& operator()(Math1D::Vector<size_t>& pos) const;

  T& operator()(Math1D::Vector<size_t>& pos);

  T data(uint pos) const;

  void operator=(const VarDimStorage& toCopy);

protected:

  Math1D::Vector<size_t> dim_;
  T* data_;

  uint size_;
};


/*********** implementation *******/


template <typename T>
VarDimStorage<T>::VarDimStorage(const Math1D::Vector<size_t>& dim) : dim_(dim) {

  size_ = (dim.size() == 0) ? 0 : 1;
  
  for (uint k=0; k < dim.size(); k++)
    size_ *= dim[k];

  data_ = new T[size_];
}

template <typename T>
VarDimStorage<T>::VarDimStorage(const Math1D::Vector<size_t>& dim, T fill) : dim_(dim) {

  size_ = (dim.size() == 0) ? 0 : 1;
  
  for (size_t k=0; k < dim.size(); k++)
    size_ *= dim[k];

  data_ = new T[size_];

  for (size_t k=0; k < dim.size(); k++) 
    data_[k] = fill;
}

template <typename T>
VarDimStorage<T>::VarDimStorage(const VarDimStorage& toCopy) {

  size_ = toCopy.size();
  dim_ = toCopy.dim_vector();
  
  data_ = new T[size_];
  for (uint k=0; k < size_; k++) {
    data_[k] = toCopy.data(k);
  }
}

template <typename T>
VarDimStorage<T>::~VarDimStorage() {
  delete[] data_;
}

template <typename T>
size_t VarDimStorage<T>::size() const {
  return size_;
}


template <typename T>
void VarDimStorage<T>::operator=(const VarDimStorage& toCopy) {

  size_ = toCopy.size();
  dim_ = toCopy.dim_vector();
  
  data_ = new T[size_];
  for (uint k=0; k < size_; k++) {
    data_[k] = toCopy.data(k);
  }
}

template <typename T>
size_t VarDimStorage<T>::dim(uint n) const {
  return dim_[n];
}

template <typename T>
size_t VarDimStorage<T>::nDims() const {
  return dim_.size();
}

template <typename T>
T VarDimStorage<T>::data(uint pos) const {

  return data_[pos];
}

template <typename T>
const Math1D::Vector<size_t>& VarDimStorage<T>::dim_vector() const {
  return dim_;
}

template <typename T>
const T& VarDimStorage<T>::operator()(Math1D::Vector<size_t>& pos) const {

  assert(pos.size() == dim_.size());

  uint data_pos = 0;

  for (size_t k=0; k < dim_.size(); k++) {

    assert(pos[k] < dim_[k]);

    if (k > 0)
      data_pos *= dim_[k];

    data_pos += pos[k];
  }

  return data_[data_pos];
}


template <typename T>
T& VarDimStorage<T>::operator()(Math1D::Vector<size_t>& pos) {

  assert(pos.size() == dim_.size());

  uint data_pos = 0;

  for (size_t k=0; k < dim_.size(); k++) {

    assert(pos[k] < dim_[k]);

    if (k > 0)
      data_pos *= dim_[k];

    data_pos += pos[k];
  }

  return data_[data_pos];
}



#endif
