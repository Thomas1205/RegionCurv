/*** written by Thomas Schoenemann as a private person without employment, July 2009 ***/
/*** if you desire the checked version, make sure your compiler defines the option SAFE_MODE on the command line ***/

#ifndef STORAGE2D_HH
#define STORAGE2D_HH

#include "makros.hh"
#include <algorithm>

//two-dimensional container class for objects of any type T
//(i.e. neither mathematical nor streaming operations need to be defined on T)
template<typename T>
class Storage2D {
public:

  //default constructor
  Storage2D();

  Storage2D(size_t xDim, size_t yDim);

  Storage2D(size_t xDim, size_t yDim, T default_value);

  //copy constructor
  Storage2D(const Storage2D<T>& toCopy);

  ~Storage2D();

  virtual const std::string& name() const;

  //saves all existing entries, new positions contain undefined data
  void resize(size_t newxDim, size_t newyDim);

  //saves all existing entries, new positions are filled with <code> fill_value </code>
  void resize(size_t newxDim, size_t newyDim, T fill_value);

  //all elements are uninitialized after this operation
  void resize_dirty(size_t newxDim, size_t newyDim);

  void set_constant(T new_constant);

  //access on an element
  OPTINLINE T& operator()(size_t x, size_t y) const;

  void operator=(const Storage2D<T>& toCopy);

  inline T* direct_access();

  inline const T* direct_access() const;

  inline T& direct_access(size_t i);

  inline T direct_access(size_t i) const;

  inline T value(size_t i) const;

  inline size_t xDim() const;

  inline size_t yDim() const;
  
  inline size_t size() const;

protected:

  T* data_;
  size_t xDim_;
  size_t yDim_;
  size_t size_;
  static const std::string stor2D_name_;
};

template<typename T>
class NamedStorage2D : public Storage2D<T> {
public:

  NamedStorage2D();

  NamedStorage2D(std::string name);
  
  NamedStorage2D(size_t xDim, size_t yDim, std::string name);
  
  NamedStorage2D(size_t xDim, size_t yDim, T default_value, std::string name);

  virtual const std::string& name() const;

  inline void operator=(const Storage2D<T>& toCopy);

  //NOTE: the name is NOT copied
  inline void operator=(const NamedStorage2D<T>& toCopy);

protected:
  std::string name_;
};


/**************************** implementation **************************************/

template<typename T>
/*static*/ const std::string Storage2D<T>::stor2D_name_ = "unnamed 2Dstorage";

//constructors
template<typename T>
Storage2D<T>::Storage2D() : data_(0), xDim_(0), yDim_(0), size_(0) {}

template<typename T>
Storage2D<T>::Storage2D(size_t xDim, size_t yDim) : xDim_(xDim), yDim_(yDim) {

  size_ = xDim_*yDim_;
  data_ = new T[size_];
}

template<typename T>
Storage2D<T>::Storage2D(size_t xDim, size_t yDim, T default_value) : xDim_(xDim), yDim_(yDim) {

  size_ = xDim_*yDim_;
  data_ = new T[size_];
  for (size_t i=0; i < size_; i++)
    data_[i] = default_value;
}

//copy constructor
template<typename T>
Storage2D<T>::Storage2D(const Storage2D<T>& toCopy) {

  xDim_ = toCopy.xDim();
  yDim_ = toCopy.yDim();
  size_ = toCopy.size();

  assert(size_ == xDim_*yDim_);

  if (size_ == 0)
    data_ = 0;
  else {
    data_ = new T[size_];

    for (size_t i = 0; i < size_; i++)
      data_[i] = toCopy.value(i);

    //this is faster for basic types but it fails for complex types where e.g. arrays have to be copied    
    //memcpy(data_,toCopy.direct_access(),size_*sizeof(T));
  }
}

//destructor
template <typename T>
Storage2D<T>::~Storage2D() {
  if (data_ != 0)
    delete[] data_;
}

template <typename T>
void Storage2D<T>::set_constant(T new_constant) {

  for (size_t i=0; i < size_; i++)
    data_[i] = new_constant;
}

template<typename T>
const std::string& Storage2D<T>::name() const {
  return Storage2D<T>::stor2D_name_;
}

template<typename T>
inline T* Storage2D<T>::direct_access() { return data_; }

template<typename T>
inline const T* Storage2D<T>::direct_access() const { return data_; }

template<typename T>
inline T& Storage2D<T>::direct_access(size_t i) { return data_[i]; }

template<typename T>
inline T Storage2D<T>::direct_access(size_t i) const {
  return data_[i];
}

template<typename T>
inline T Storage2D<T>::value(size_t i) const { return data_[i]; }

template<typename T>
inline size_t Storage2D<T>::xDim() const { return xDim_; }

template<typename T>
inline size_t Storage2D<T>::yDim() const { return yDim_; }

template<typename T>
inline size_t Storage2D<T>::size() const { return size_; }

template <typename T>
OPTINLINE T& Storage2D<T>::operator()(size_t x, size_t y) const {
#ifdef SAFE_MODE
  if (x >= xDim_ || y >= yDim_) {
    INTERNAL_ERROR << "    access on element(" << x << "," << y 
		   << ") exceeds storage dimensions of (" << xDim_ << "," << yDim_ << ")" << std::endl;
    std::cerr << "   in 2Dstorage \"" << this->name() << "\". exiting." << std::endl;  
    exit(1);
  }
#endif
  return data_[y*xDim_+x];
}

template <typename T>
void Storage2D<T>::operator=(const Storage2D<T>& toCopy) {

  if (size_ != toCopy.size()) {
    if (data_ != 0)
      delete[] data_;

    size_ = toCopy.size();
    data_ = new T[size_];
  }

  xDim_ = toCopy.xDim();
  yDim_ = toCopy.yDim();
  assert(size_ == xDim_*yDim_);

  const uint size = size_;

  for (size_t i = 0; i < size; i++)
    data_[i] = toCopy.value(i);
  
  //this is faster for basic types but it fails for complex types where e.g. arrays have to be copied
  //memcpy(data_,toCopy.direct_access(),size_*sizeof(T));
}

template<>
void Storage2D<int>::operator=(const Storage2D<int>& toCopy);

template<>
void Storage2D<uint>::operator=(const Storage2D<uint>& toCopy);

template<>
void Storage2D<float>::operator=(const Storage2D<float>& toCopy);

template<>
void Storage2D<double>::operator=(const Storage2D<double>& toCopy);

template<>
void Storage2D<long double>::operator=(const Storage2D<long double>& toCopy);


template <typename T>
void Storage2D<T>::resize(size_t newxDim, size_t newyDim) {

  if (data_ == 0) {
    data_ = new T[newxDim*newyDim];
  }
  else {

    T* new_data = new T[newxDim*newyDim];

    /* copy data */
    for (size_t y=0; y < std::min(yDim_,newyDim); y++)
      for (size_t x=0; x < std::min(xDim_,newxDim); x++)
	new_data[y*newxDim+x] = data_[y*xDim_+x];

    delete[] data_;
    data_ = new_data;
  }
    
  xDim_ = newxDim;
  yDim_ = newyDim;
  size_ = xDim_*yDim_;
}

template <typename T>
void Storage2D<T>::resize(size_t newxDim, size_t newyDim, T fill_value) {

  if (data_ == 0) {
    data_ = new T[newxDim*newyDim];

//     for (size_t i=0; i < newxDim*newyDim; i++)
//       data_[i] = fill_value;

    std::fill_n(data_,newxDim*newyDim,fill_value);
  }
  else {

    T* new_data = new T[newxDim*newyDim];
    for (size_t i=0; i < newxDim*newyDim; i++)
      new_data[i] = fill_value;

    /* copy data */
    for (size_t y=0; y < std::min(yDim_,newyDim); y++)
      for (size_t x=0; x < std::min(xDim_,newxDim); x++)
	new_data[y*newxDim+x] = data_[y*xDim_+x];

    delete[] data_;
    data_ = new_data;
  }
    
  xDim_ = newxDim;
  yDim_ = newyDim;
  size_ = xDim_*yDim_;
}

template<typename T>
void Storage2D<T>::resize_dirty(size_t newxDim, size_t newyDim) {

  if (data_ != 0) {
    delete[] data_;
  }

  xDim_ = newxDim;
  yDim_ = newyDim;
  size_ = xDim_*yDim_;

  data_ = new T[size_];
}

/***** implementation of NamedStorage2D ********/

template<typename T>
NamedStorage2D<T>::NamedStorage2D() : Storage2D<T>(), name_("yyy") {}

template<typename T>
NamedStorage2D<T>::NamedStorage2D(std::string name) : Storage2D<T>(), name_(name) {}

template<typename T>
NamedStorage2D<T>::NamedStorage2D(size_t xDim, size_t yDim, std::string name) : Storage2D<T>(xDim,yDim), name_(name) {}

template<typename T>
NamedStorage2D<T>::NamedStorage2D(size_t xDim, size_t yDim, T default_value, std::string name) 
  : Storage2D<T>(xDim,yDim,default_value), name_(name) {}

template<typename T>
/*virtual*/ const std::string& NamedStorage2D<T>::name() const {
  return name_;
}

template<typename T>
inline void NamedStorage2D<T>::operator=(const Storage2D<T>& toCopy) {
  Storage2D<T>::operator=(toCopy);
}

//NOTE: the name is NOT copied
template<typename T>
inline void NamedStorage2D<T>::operator=(const NamedStorage2D<T>& toCopy) {
  Storage2D<T>::operator=(static_cast<Storage2D<T> >(toCopy));
}




#endif
