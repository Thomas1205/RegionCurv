/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#include "makros.hh"

template<typename T>
class Storage3D {
public:

  Storage3D();

  //copy constructor
  Storage3D(const Storage3D& toCopy);

  Storage3D(size_t xDim, size_t yDim, size_t zDim);

  Storage3D(size_t xDim, size_t yDim, size_t zDim, T default_value);

  ~Storage3D();

  OPTINLINE T& operator()(size_t x, size_t y, size_t z) const;

  virtual const std::string& name() const;

  inline size_t size() const;

  inline size_t xDim() const;

  inline size_t yDim() const;

  inline size_t zDim() const;

  inline T* direct_access();

  inline const T* direct_access() const;

  inline T& direct_access(size_t i);

  inline T direct_access(size_t i) const;

  void operator=(const Storage3D<T>& toCopy);

  //existing positions are copied, new ones are uninitialized
  void resize(size_t newxDim, size_t newyDim, size_t newzDim);

  //existing positions are copied, new ones are uninitialized
  void resize(size_t newxDim, size_t newyDim, size_t newzDim, T default_value);
  
  //all elements are uninitialized after this operation
  void resize_dirty(size_t newxDim, size_t newyDim, size_t newzDim);

protected:
  size_t xDim_;
  size_t yDim_;
  size_t zDim_;
  size_t size_;

  T* data_;
  static const std::string stor3D_name_;
};


/******************************************** implementation **************************************************/
template<typename T>
/*static*/ const std::string Storage3D<T>::stor3D_name_ = "unnamed Storage3D";

template<typename T>
Storage3D<T>::Storage3D() : xDim_(0), yDim_(0), zDim_(0), size_(0), data_(0) {}

template<typename T>
Storage3D<T>::Storage3D(const Storage3D<T>& toCopy) {

  xDim_ = toCopy.xDim();
  yDim_ = toCopy.yDim();
  zDim_ = toCopy.zDim();
  size_ = toCopy.size();

  data_ = new T[size_];

  for (size_t i=0; i < size_; i++)
    data_[i] = toCopy.direct_access(i);
  
  //this is faster for basic types but it fails for complex types where e.g. arrays have to be copied
  //memcpy(data_,toCopy.direct_access(),size_*sizeof(T));  
}

template<typename T>
Storage3D<T>::Storage3D(size_t xDim, size_t yDim, size_t zDim) : xDim_(xDim), yDim_(yDim), zDim_(zDim) {

  size_ = xDim_*yDim_*zDim_;
  data_ = new T[size_];
}

template<typename T>
Storage3D<T>::Storage3D(size_t xDim, size_t yDim, size_t zDim, T default_value) :
  xDim_(xDim), yDim_(yDim), zDim_(zDim) {

  size_ = xDim_*yDim_*zDim_;
  data_ = new T[size_];
  for (size_t i=0; i < size_; i++) {
    data_[i] = default_value;
  }
}


template<typename T>
Storage3D<T>::~Storage3D() {

  if (data_ != 0)
    delete[] data_;
}

template<typename T>
OPTINLINE T& Storage3D<T>::operator()(size_t x, size_t y, size_t z) const {
#ifdef SAFE_MODE
  if (x >= xDim_ || y >= yDim_ || z >= zDim_) {
    INTERNAL_ERROR << "     invalid access on element (" << x << "," << y << "," << z << ") of 3D-storage \"" 
		   << this->name() << "\":" << std::endl;
    std::cerr << "     dimensions " << xDim_ << "x" << yDim_ << "x" << zDim_ << " exceeded. Exiting..." << std::endl;
    exit(1);
  }
#endif
  return data_[(y*xDim_+x)*zDim_+z];
}

template<typename T>
/*virtual*/ const std::string& Storage3D<T>::name() const {
  return stor3D_name_;
}

template<typename T>
inline size_t Storage3D<T>::size() const {
  return size_;
}

template<typename T>
inline size_t Storage3D<T>::xDim() const {
  return xDim_;
}

template<typename T>
inline size_t Storage3D<T>::yDim() const {
  return yDim_;
}

template<typename T>
inline size_t Storage3D<T>::zDim() const {
  return zDim_;
}

template<typename T>
inline T* Storage3D<T>::direct_access() {
  return data_;
}

template<typename T>
inline const T* Storage3D<T>::direct_access() const {
  return data_;
}

template<typename T>
inline T& Storage3D<T>::direct_access(size_t i) {
  return data_[i];
}

template<typename T>
inline T Storage3D<T>::direct_access(size_t i) const {
  return data_[i];
}

template<typename T>
void Storage3D<T>::operator=(const Storage3D<T>& toCopy) {

  if (size_ != toCopy.size()) {
    if (data_ != 0) {
      delete[] data_;
    }

    size_ = toCopy.size();
    data_ = new T[size_];
  }

  xDim_ = toCopy.xDim();
  yDim_ = toCopy.yDim();
  zDim_ = toCopy.zDim();
  assert(size_ == xDim_*yDim_*zDim_);

  for (size_t i=0; i < size_; i++)
    data_[i] = toCopy.direct_access(i);
  
  //this is faster for basic types but it fails for complex types where e.g. arrays have to be copied
  //memcpy(data_,toCopy.direct_access(),size_*sizeof(T));
}


//existing positions are copied, new ones are uninitialized
template<typename T>
void Storage3D<T>::resize(size_t newxDim, size_t newyDim, size_t newzDim) {

  size_t new_size = newxDim*newyDim*newzDim;
  T* new_data = new T[new_size];

  if (data_ != 0) {
    
    //copy existing elements
    for (size_t x=0; x < std::min(xDim_,newxDim); x++) {
      for (size_t y=0; y < std::min(yDim_,newyDim); y++) {
	for (size_t z=0; z < std::min(zDim_,newzDim); z++) {
	  new_data[(y*newxDim+x)*newzDim+z] = data_[(y*xDim_+x)*zDim_+z];
	}
      }
    }

    delete[] data_;
  }
  data_ = new_data;
  size_ = new_size;
  xDim_ = newxDim;
  yDim_ = newyDim;
  zDim_ = newzDim;
}

//existing positions are copied, new ones are uninitialized
template<typename T>
void Storage3D<T>::resize(size_t newxDim, size_t newyDim, size_t newzDim, T default_value) {

  size_t new_size = newxDim*newyDim*newzDim;
  T* new_data = new T[new_size];
  for (size_t i=0; i < new_size; i++)
    new_data[i] = default_value;

  if (data_ != 0) {
    
    //copy existing elements
    for (size_t x=0; x < std::min(xDim_,newxDim); x++) {
      for (size_t y=0; y < std::min(yDim_,newyDim); y++) {
	for (size_t z=0; z < std::min(zDim_,newzDim); z++) {
	  new_data[(y*newxDim+x)*newzDim+z] = data_[(y*xDim_+x)*zDim_+z];
	}
      }
    }

    delete[] data_;
  }
  data_ = new_data;
  size_ = new_size;
  xDim_ = newxDim;
  yDim_ = newyDim;
  zDim_ = newzDim;
}

//all elements are uninitialized after this operation
template<typename T>
void Storage3D<T>::resize_dirty(size_t newxDim, size_t newyDim, size_t newzDim) {

  if (data_ != 0)
    delete[] data_;

  xDim_ = newxDim;
  yDim_ = newyDim;
  zDim_ = newzDim;
  size_ = xDim_*yDim_*zDim_;

  data_ = new T[size_];
}
