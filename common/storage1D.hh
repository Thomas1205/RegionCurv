/*-*-c++-*-*/ 
/*** first version written by Thomas Schoenemann as a private person without employment, September 2009 ***/
/*** much refined by Thomas Schoenemann  at Lund University, Sweden and the University of Pisa, Italy, 2010 - 2011 **/
/*** if you desire the checked version, make sure your compiler defines the option SAFE_MODE on the command line ***/

#ifndef STORAGE1D_HH
#define STORAGE1D_HH

#include "makros.hh"
#include <string>
#include <cstring>
#include <algorithm>

template<typename T>
class Storage1D {
public:

  Storage1D();
    
  Storage1D(size_t size);
  
  Storage1D(size_t size, T default_value);

  //copy constructor
  Storage1D(const Storage1D<T>& toCopy);
    
  ~Storage1D();
  
  virtual const std::string& name() const;
  
  OPTINLINE T& operator[](size_t i) const;

  void operator=(const Storage1D<T>& toCopy);

  //maintains the values of existing positions, new ones are undefined 
  void resize(size_t new_size);

  //maintains the values of exisitng positions, new ones are filled with <code> fill_value </code>
  void resize(size_t new_size, T fill_value);

  //all elements are undefined after this operation
  void resize_dirty(size_t new_size);
  
  inline size_t size() const;
  
  inline T* direct_access();

  inline const T* direct_access() const;

  inline T& direct_access(size_t i);
  
  inline T direct_access(size_t i) const;

  void set_constant(T constant);
  
protected:
  
  T* data_;
  size_t size_;
  static const std::string stor1D_name_;
};

template<typename T>
class NamedStorage1D : public Storage1D<T> {
public:

  NamedStorage1D();

  NamedStorage1D(std::string name);
  
  NamedStorage1D(size_t size, std::string name);
  
  NamedStorage1D(size_t size, T default_value, std::string name);

  virtual const std::string& name() const;

  inline void operator=(const Storage1D<T>& toCopy);

  //NOTE: the name is NOT copied
  inline void operator=(const NamedStorage1D<T>& toCopy);

protected:
  std::string name_;
};

template<typename T>
std::ostream& operator<<(std::ostream& s, const Storage1D<T>& v);

template<typename T>
bool operator==(const Storage1D<T>& v1, const Storage1D<T>& v2);

template<typename T>
bool operator!=(const Storage1D<T>& v1, const Storage1D<T>& v2);


/***********************/

//this class is meant to replace std::vector with its push_back() functionality.
// It has slightly less functionality, though. E.g. erase() is not available.
template<typename T>
class FlexibleStorage1D {
public:

  FlexibleStorage1D();

  FlexibleStorage1D(size_t reserved_size);

  //copy constructor
  FlexibleStorage1D(const FlexibleStorage1D<T>& toCopy);

  ~FlexibleStorage1D();

  virtual const std::string& name() const;
  
  OPTINLINE T& operator[](size_t i) const;

  void resize(size_t size, bool exact_fit = false);

  void set_constant(T val);

  void operator=(const FlexibleStorage1D<T>& toCopy);

  size_t append(T val);

  void append(Storage1D<T>& toAppend);

  void append(FlexibleStorage1D<T>& toAppend);

  size_t size() const;

  size_t reserved_size() const;

  T* direct_access();

  const T* direct_access() const;

protected:

  T* data_;
  size_t size_;
  size_t reserved_size_;
  static const std::string flex_stor1D_name_;
};

template<typename T>
std::ostream& operator<<(std::ostream& s, const FlexibleStorage1D<T>& v);

template<typename T>
bool operator==(const FlexibleStorage1D<T>& v1, const FlexibleStorage1D<T>& v2);


template<typename T>
class NamedFlexibleStorage1D : public FlexibleStorage1D<T> {
public:

  NamedFlexibleStorage1D();

  NamedFlexibleStorage1D(const std::string& name);

  NamedFlexibleStorage1D(size_t reserved_size, const std::string& name);

  //copy constructors
  NamedFlexibleStorage1D(const NamedFlexibleStorage1D<T>& toCopy);

  NamedFlexibleStorage1D(const FlexibleStorage1D<T>& toCopy);

  virtual const std::string& name() const;

  //operators
  void operator=(const NamedFlexibleStorage1D<T>& toCopy);

  void operator=(const FlexibleStorage1D<T>& toCopy);

protected:
  std::string name_;
};

/********************************************** implementation ************************************/

/******* implementation of Storage1D *********/

template<typename T>
/*static*/ const std::string Storage1D<T>::stor1D_name_ = "unnamed 1Dstorage";

template<typename T>
Storage1D<T>::Storage1D(): data_(0), size_(0) {}

template<typename T>
Storage1D<T>::Storage1D(size_t size): size_(size) {
  data_ = new T[size];    
}

template<typename T>
Storage1D<T>::Storage1D(size_t size, T default_value): size_(size) {
  data_ = new T[size_];    

//    for (size_t i=0; i < size_; i++)
//      data_[i] = default_value;

  //std::fill_n(data_,size,default_value); //experimental result: fill_n is usually faster
  std::fill(data_, data_+size, default_value); //fill and fill_n are of equal speed
}

//copy constructor
template<typename T>
Storage1D<T>::Storage1D(const Storage1D<T>& toCopy) {

  size_ = toCopy.size();
  data_ = new T[size_];
 
  const std::size_t size = size_;
   
  for (size_t i=0; i < size; i++) {
    data_[i] = toCopy.direct_access(i);
  }

  //this is faster for basic types but it fails for complex types where e.g. arrays have to be copied
  //memcpy(data_,toCopy.direct_access(),size_*sizeof(T));
}

template<>
Storage1D<int>::Storage1D(const Storage1D<int>& toCopy); 

template<>
Storage1D<uint>::Storage1D(const Storage1D<uint>& toCopy); 

template<>
Storage1D<float>::Storage1D(const Storage1D<float>& toCopy); 

template<>
Storage1D<double>::Storage1D(const Storage1D<double>& toCopy); 


template<typename T>
void Storage1D<T>::set_constant(T constant) {
  
  //     for (size_t i=0; i < size_; i++) 
  //       data_[i] = constant;
  std::fill_n(data_,size_,constant); //experimental result: fill_n is usually faster
}


template<typename T>
Storage1D<T>::~Storage1D() {
  if (data_ != 0)
    delete[] data_;     
}

template<typename T>
/*virtual*/ const std::string& Storage1D<T>::name() const {
  return Storage1D::stor1D_name_;
} 

template<typename T>
inline size_t Storage1D<T>::size() const {
  return size_;
}

template<typename T>
inline T* Storage1D<T>::direct_access() {
  return data_;
}

template<typename T>
inline const T* Storage1D<T>::direct_access() const {
  return data_;
}

template<typename T>
inline T& Storage1D<T>::direct_access(size_t i) {
  return data_[i];
}

template<typename T>
inline T Storage1D<T>::direct_access(size_t i) const {
  return data_[i];
}

template<typename T>
OPTINLINE T& Storage1D<T>::operator[](size_t i) const {
#ifdef SAFE_MODE
  if (i >= size_) {
    INTERNAL_ERROR << "    invalid access on element " << i 
		   << " for Storage1D " <<  "\"" << this->name() << "\" of type " << typeid(T).name()
		   << " with " << size_ << " elements. exiting." << std::endl;
    exit(1);  
  }
#endif
  return data_[i];
}

template<typename T>
void Storage1D<T>::operator=(const Storage1D<T>& toCopy) {

  if (size_ != toCopy.size()) {

    if (data_ != 0)
      delete[] data_;

    size_ = toCopy.size();
    data_ = new T[size_];
  }
    
  for (size_t i=0; i < size_; i++) {
     data_[i] = toCopy.direct_access(i);
  }

  //this is faster for basic types but it fails for complex types where e.g. arrays have to be copied
  //memcpy(data_,toCopy.direct_access(),size_*sizeof(T));
}

template<>
void Storage1D<uint>::operator=(const Storage1D<uint>& toCopy);

template<>
void Storage1D<int>::operator=(const Storage1D<int>& toCopy);

template<>
void Storage1D<float>::operator=(const Storage1D<float>& toCopy);

template<>
void Storage1D<double>::operator=(const Storage1D<double>& toCopy);

template<>
void Storage1D<long double>::operator=(const Storage1D<long double>& toCopy);


//maintains the values of existing positions, new ones are undefined 
template<typename T>
void Storage1D<T>::resize(size_t new_size) {

  if (data_ == 0) {
    data_ = new T[new_size];
  }
  else if (size_ != new_size) {
    T* new_data = new T[new_size];

    for (size_t i=0; i < std::min(size_,new_size); i++)
      new_data[i] = data_[i];
    
    //this is faster for basic types but it fails for complex types where e.g. arrays have to be copied
    //memcpy(new_data,data_,std::min(size_,new_size)*sizeof(T));
    
    delete[] data_;
    data_ = new_data;
  }

  size_ = new_size;
}



//maintains the values of exisitng positions, new ones are filled with <code> fill_value </code>
template<typename T>
void Storage1D<T>::resize(size_t new_size, T fill_value) {

  if (data_ == 0) {
    data_ = new T[new_size];
    for (size_t i=0; i < new_size; i++)
      data_[i] = fill_value;
  }
  else if (size_ != new_size) {
    T* new_data = new T[new_size];
    for (size_t i=size_; i < new_size; i++)
      new_data[i] = fill_value;

    for (size_t i=0; i < std::min(size_,new_size); i++)
      new_data[i] = data_[i];
    //memcpy(new_data,data_,std::min(size_,new_size)*sizeof(T));
    
    delete[] data_;
    data_ = new_data;
  }

  size_ = new_size;
}

//all elements are undefined after this operation
template<typename T>
void Storage1D<T>::resize_dirty(size_t new_size) {

  if (size_ != new_size) {
    if (data_ != 0)
      delete[] data_;

    data_ = new T[new_size];
  }
  size_ = new_size;
}

/******** implementation of NamedStorage1D ***************/ 

template<typename T>
NamedStorage1D<T>::NamedStorage1D() : Storage1D<T>(), name_("yyy") {}

template<typename T>
NamedStorage1D<T>::NamedStorage1D(std::string name) : Storage1D<T>(), name_(name) {}
  
template<typename T>
NamedStorage1D<T>::NamedStorage1D(size_t size, std::string name) : Storage1D<T>(size), name_(name) {}
  
template<typename T>
NamedStorage1D<T>::NamedStorage1D(size_t size, T default_value, std::string name) : 
  Storage1D<T>(size,default_value), name_(name) {}

template<typename T>
/*virtual*/ const std::string& NamedStorage1D<T>::name() const {
  return name_;
}

template<typename T>
inline void NamedStorage1D<T>::operator=(const Storage1D<T>& toCopy) {
  Storage1D<T>::operator=(toCopy);
}

//NOTE: the name is NOT copied
template<typename T>
inline void NamedStorage1D<T>::operator=(const NamedStorage1D<T>& toCopy) {
  Storage1D<T>::operator=(static_cast<const Storage1D<T>&>(toCopy));
}


template<typename T>
std::ostream& operator<<(std::ostream& s, const Storage1D<T>& v) {

  s << "[ ";
  for (int i=0; i < ((int) v.size()) - 1; i++)
    s << v[i] << ",";
  if (v.size() > 0)
    s << v[v.size()-1];
  s << " ]";
  
  return s;
}

template<typename T>
bool operator==(const Storage1D<T>& v1, const Storage1D<T>& v2) {

  if (v1.size() != v2.size())
    return false;
  
  for (size_t k=0; k < v1.size(); k++) {
    if (v1[k] != v2[k])
      return false;
  }
  return true;
}

template<typename T>
bool operator!=(const Storage1D<T>& v1, const Storage1D<T>& v2) {
  return !operator==(v1,v2);
}


/******* implementation of FlexibleStorage1D *********/

template<typename T>
/*static*/ const std::string FlexibleStorage1D<T>::flex_stor1D_name_ = "unnamed flexible 1Dstorage";

template<typename T>
FlexibleStorage1D<T>::FlexibleStorage1D() : size_(0) {
  reserved_size_ = 4;
  data_ = new T[reserved_size_];
}

template<typename T>
FlexibleStorage1D<T>::FlexibleStorage1D(size_t reserved_size)  : size_(0), reserved_size_(reserved_size) {
  data_ = new T[reserved_size_];
}

//copy constructor
template<typename T>
FlexibleStorage1D<T>::FlexibleStorage1D(const FlexibleStorage1D<T>& toCopy) {

  size_ = toCopy.size();
  reserved_size_ = toCopy.reserved_size();
  
  data_ = new T[reserved_size_];
  for (uint k=0; k < toCopy.size(); k++)
    data_[k] = toCopy[k];
}

template<>
FlexibleStorage1D<uint>::FlexibleStorage1D(const FlexibleStorage1D<uint>& toCopy);

template<typename T>
void FlexibleStorage1D<T>::operator=(const FlexibleStorage1D<T>& toCopy) {

  uint new_res = toCopy.reserved_size();
  if (new_res != reserved_size_) {
    reserved_size_ = new_res;

    if (data_ != 0)
      delete[] data_;
    data_ = new T[reserved_size_];
  }

  size_ = toCopy.size();

  for (uint k=0; k < size_; k++)
    data_[k] = toCopy[k];
}

template<>
void FlexibleStorage1D<uint>::operator=(const FlexibleStorage1D<uint>& toCopy);

template<typename T>
/*virtual*/ const std::string& FlexibleStorage1D<T>::name() const {
  return flex_stor1D_name_;
}

template<typename T>
FlexibleStorage1D<T>::~FlexibleStorage1D()  {
  delete [] data_;
}

template<typename T>
void FlexibleStorage1D<T>::set_constant(T val) {

  for (size_t k=0; k < size_; k++)
    data_[k] = val;
}

template<typename T>
size_t FlexibleStorage1D<T>::size() const {
  return size_;
}

template<typename T>
size_t FlexibleStorage1D<T>::reserved_size() const {
  return reserved_size_;
}

template<typename T>
size_t FlexibleStorage1D<T>::append(T val) {

  if (size_ == reserved_size_) {

    reserved_size_ = size_t(1.2 * reserved_size_) + 4;

    T* new_data = new T[reserved_size_];
    for (uint k=0; k < size_; k++)
      new_data[k] = data_[k];

    delete[] data_;
    data_ = new_data;
  }

  const uint k = size_;
  data_[k] = val;

  size_++;

  return k;
}

template<typename T>
void FlexibleStorage1D<T>::append(Storage1D<T>& toAppend) {

  if (reserved_size_ < size_ + toAppend.size()) {

    reserved_size_ = size_ + toAppend.size() + 2;

    T* new_data = new T[reserved_size_];
    for (uint k=0; k < size_; k++)
      new_data[k] = data_[k];

    delete[] data_;
    data_ = new_data;
  }

  for (uint k=0; k < toAppend.size(); k++) {
    data_[size_] = toAppend[k];
    size_++;
  }
}

template<typename T>
void FlexibleStorage1D<T>::append(FlexibleStorage1D<T>& toAppend) {

  if (reserved_size_ < size_ + toAppend.size()) {

    reserved_size_ = size_ + toAppend.size() + 2;

    T* new_data = new T[reserved_size_];
    for (uint k=0; k < size_; k++)
      new_data[k] = data_[k];

    delete[] data_;
    data_ = new_data;
  }

  for (uint k=0; k < toAppend.size(); k++) {
    data_[size_] = toAppend[k];
    size_++;
  }
}

template<typename T>
void FlexibleStorage1D<T>::resize(size_t size, bool exact_fit) {

  if (size > reserved_size_ || size < (reserved_size_ / 3) ) {
    
    reserved_size_ = size;
    T* new_data = new T[reserved_size_];
    for (uint k=0; k < std::min(size_,size); k++)
      new_data[k] = data_[k];
    
    delete[] data_;
    data_ = new_data;
  }

  if (size < size_)
    size_ = size;
  
  if (exact_fit && size_ != reserved_size_) {

    reserved_size_ = size_;
    T* new_data = new T[reserved_size_];
    for (uint k=0; k < size_; k++)
      new_data[k] = data_[k];
    
    delete[] data_;
    data_ = new_data;
  }
}

template<typename T>
OPTINLINE T& FlexibleStorage1D<T>::operator[](size_t i) const {

#ifdef SAFE_MODE
  if (i >= size_) {
    INTERNAL_ERROR << "    invalid access on element " << i 
		   << " for FlexibleStorage1D " <<  "\"" << this->name() << "\" of type " << typeid(T).name()
		   << " with " << size_ << " (valid) elements. exiting." << std::endl;
    exit(1);  
  }
#endif
  return data_[i];
}

template<typename T>
T* FlexibleStorage1D<T>::direct_access() {

  return data_;
}

template<typename T>
const T* FlexibleStorage1D<T>::direct_access() const {

  return data_;
}

template<typename T>
std::ostream& operator<<(std::ostream& s, const FlexibleStorage1D<T>& v) {

  s << "[ ";
  for (int i=0; i < ((int) v.size()) - 1; i++)
    s << v[i] << ",";
  if (v.size() > 0)
    s << v[v.size()-1];
  s << " ]";
  
  return s;
}

template<typename T>
bool operator==(const FlexibleStorage1D<T>& v1, const FlexibleStorage1D<T>& v2) {

  if (v1.size() != v2.size())
    return false;

  for (size_t k=0; k < v1.size(); k++) {
    if (v1[k] != v2[k])
      return false;
  }
  return true;
}


/***********************************/

template<typename T>
NamedFlexibleStorage1D<T>::NamedFlexibleStorage1D() : name_("unfs1d") {}

template<typename T>
NamedFlexibleStorage1D<T>::NamedFlexibleStorage1D(const std::string& name) : name_(name) {
}

template<typename T>
NamedFlexibleStorage1D<T>::NamedFlexibleStorage1D(size_t reserved_size, const std::string& name) :
  FlexibleStorage1D<T>(reserved_size), name_(name) {}

//Note: the name is NOT copied
template<typename T>
NamedFlexibleStorage1D<T>::NamedFlexibleStorage1D(const NamedFlexibleStorage1D<T>& toCopy) : 
  FlexibleStorage1D<T>(toCopy), name_("unfs1d") {
}

template<typename T>
NamedFlexibleStorage1D<T>::NamedFlexibleStorage1D(const FlexibleStorage1D<T>& toCopy) : 
  FlexibleStorage1D<T>(toCopy), name_("unfs1d") {
}

template<typename T>
/*virtual*/ const std::string& NamedFlexibleStorage1D<T>::name() const {
  return name_;
}

template<typename T>
void NamedFlexibleStorage1D<T>::operator=(const NamedFlexibleStorage1D<T>& toCopy) {
  FlexibleStorage1D<T>::operator=(toCopy);
}

template<typename T>
void NamedFlexibleStorage1D<T>::operator=(const FlexibleStorage1D<T>& toCopy) {
  FlexibleStorage1D<T>::operator=(toCopy);
}

#endif
