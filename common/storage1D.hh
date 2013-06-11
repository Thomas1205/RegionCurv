/*-*-c++-*-*/ 
/*** first version written by Thomas Schoenemann as a private person without employment, September 2009 ***/
/*** much refined by Thomas Schoenemann  at Lund University, Sweden, the University of Pisa, Italy, ***
 *** and the University of Düsseldorf, Germany 2010 - 2012 **/
/*** if you desire the checked version, make sure your compiler defines the option SAFE_MODE on the command line ***/


#ifndef STORAGE1D_HH
#define STORAGE1D_HH

#include "makros.hh"
#include <string>
#include <cstring>
#include <algorithm>

template<typename T, typename ST=size_t>
class Storage1D {
public:

  Storage1D();
    
  Storage1D(ST size);
  
  Storage1D(ST size, T default_value);

  //copy constructor
  Storage1D(const Storage1D<T,ST>& toCopy);
    
  ~Storage1D();
  
  virtual const std::string& name() const;
  
  OPTINLINE const T& operator[](ST i) const;

  OPTINLINE T& operator[](ST i);

  void operator=(const Storage1D<T,ST>& toCopy);

#ifdef SAFE_MODE
  //for some reason g++ allows to assign an object of type T, but this does NOT produce the effect one would expect
  // => define this operator in safe mode, only to check that such an assignment is not made
  void operator=(const T& invalid_object);
#endif

  //maintains the values of existing positions, new ones are undefined 
  void resize(ST new_size);

  //maintains the values of exisitng positions, new ones are filled with <code> fill_value </code>
  void resize(ST new_size, T fill_value);

  //all elements are undefined after this operation
  void resize_dirty(ST new_size);
  
  inline ST size() const;
  
  inline T* direct_access();

  inline const T* direct_access() const;

  inline T& direct_access(ST i);
  
  inline T direct_access(ST i) const;

  void set_constant(T constant);
  
protected:
  
  T* data_;
  ST size_;
  static const std::string stor1D_name_;
};

template<typename T, typename ST=size_t>
class NamedStorage1D : public Storage1D<T,ST> {
public:

  NamedStorage1D();

  NamedStorage1D(std::string name);
  
  NamedStorage1D(ST size, std::string name);
  
  NamedStorage1D(ST size, T default_value, std::string name);

  virtual const std::string& name() const;

  inline void operator=(const Storage1D<T,ST>& toCopy);

  //NOTE: the name is NOT copied
  inline void operator=(const NamedStorage1D<T,ST>& toCopy);

protected:
  std::string name_;
};

template<typename T, typename ST>
std::ostream& operator<<(std::ostream& s, const Storage1D<T,ST>& v);

template<typename T, typename ST>
bool operator==(const Storage1D<T,ST>& v1, const Storage1D<T,ST>& v2);

template<typename T, typename ST>
bool operator!=(const Storage1D<T,ST>& v1, const Storage1D<T,ST>& v2);

template<typename T,typename ST>
bool operator<(const Storage1D<T,ST>& v1, const Storage1D<T,ST>& v2); 

template<typename T,typename ST>
bool operator<=(const Storage1D<T,ST>& v1, const Storage1D<T,ST>& v2); 

template<typename T,typename ST>
bool operator>(const Storage1D<T,ST>& v1, const Storage1D<T,ST>& v2); 

template<typename T,typename ST>
bool operator>=(const Storage1D<T,ST>& v1, const Storage1D<T,ST>& v2); 

namespace Makros {


  template<typename T, typename ST>
  class Typename<Storage1D<T,ST> > {
  public:

    std::string name() const {

      return "Storage1D<" + Makros::Typename<T>() + "," + Makros::Typename<ST>() + "> ";
    }
  };

  template<typename T>
  class Typename<Storage1D<T> > {
  public:

    std::string name() const {

      return "Storage1D<" + Makros::Typename<T>() + "> ";
    }
  };

  template<typename T, typename ST>
  class Typename<NamedStorage1D<T,ST> > {
  public:

    std::string name() const {

      return "NamedStorage1D<" + Makros::Typename<T>() + "," + Makros::Typename<ST>() + "> ";
    }
  };

  template<typename T>
  class Typename<NamedStorage1D<T> > {
  public:

    std::string name() const {

      return "NamedStorage1D<" + Makros::Typename<T>() + "> ";
    }
  };
  
}

/***********************/

//this class is meant to replace std::vector with its push_back() functionality.
// It has slightly less functionality, though. E.g. erase() is not available.
template<typename T, typename ST=size_t>
class FlexibleStorage1D {
public:

  FlexibleStorage1D();

  FlexibleStorage1D(ST reserved_size);

  //copy constructor
  FlexibleStorage1D(const FlexibleStorage1D<T,ST>& toCopy);

  ~FlexibleStorage1D();

  virtual const std::string& name() const;
  
  OPTINLINE T& operator[](ST i) const;

  void resize(ST size, bool exact_fit = false);

  void set_constant(T val);

  void operator=(const FlexibleStorage1D<T,ST>& toCopy);

  ST append(T val);

  void append(Storage1D<T,ST>& toAppend);

  void append(FlexibleStorage1D<T,ST>& toAppend);

  ST size() const;

  ST reserved_size() const;

  T* direct_access();

  const T* direct_access() const;

protected:

  T* data_;
  ST size_;
  ST reserved_size_;
  static const std::string flex_stor1D_name_;
};

template<typename T, typename ST>
std::ostream& operator<<(std::ostream& s, const FlexibleStorage1D<T,ST>& v);

template<typename T, typename ST>
bool operator==(const FlexibleStorage1D<T,ST>& v1, const FlexibleStorage1D<T,ST>& v2);


template<typename T, typename ST=size_t>
class NamedFlexibleStorage1D : public FlexibleStorage1D<T,ST> {
public:

  NamedFlexibleStorage1D();

  NamedFlexibleStorage1D(const std::string& name);

  NamedFlexibleStorage1D(ST reserved_size, const std::string& name);

  //copy constructors
  NamedFlexibleStorage1D(const NamedFlexibleStorage1D<T,ST>& toCopy);

  NamedFlexibleStorage1D(const FlexibleStorage1D<T,ST>& toCopy);

  virtual const std::string& name() const;

  //operators
  void operator=(const NamedFlexibleStorage1D<T,ST>& toCopy);

  void operator=(const FlexibleStorage1D<T,ST>& toCopy);

protected:
  std::string name_;
};

/********************************************** implementation ************************************/

/******* implementation of Storage1D *********/

template<typename T,typename ST>
/*static*/ const std::string Storage1D<T,ST>::stor1D_name_ = "unnamed 1Dstorage";

template<typename T,typename ST>
Storage1D<T,ST>::Storage1D(): data_(0), size_(0) {}

template<typename T,typename ST>
Storage1D<T,ST>::Storage1D(ST size): size_(size) {
  data_ = new T[size];    
}

template<typename T,typename ST>
Storage1D<T,ST>::Storage1D(ST size, T default_value): size_(size) {
  data_ = new T[size_];    

  std::fill(data_, data_+size, default_value); //fill and fill_n are of equal speed
}

//copy constructor
template<typename T,typename ST>
Storage1D<T,ST>::Storage1D(const Storage1D<T,ST>& toCopy) {

  size_ = toCopy.size();
  data_ = new T[size_];
 
  const ST size = size_;
   
  for (ST i=0; i < size; i++) {
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

template<>
Storage1D<ushort>::Storage1D(const Storage1D<ushort>& toCopy); 


template<typename T,typename ST>
void Storage1D<T,ST>::set_constant(T constant) {
  
  std::fill_n(data_,size_,constant); //experimental result: fill_n is usually faster
}


template<typename T,typename ST>
Storage1D<T,ST>::~Storage1D() {
  if (data_ != 0)
    delete[] data_;     
}

template<typename T,typename ST>
/*virtual*/ const std::string& Storage1D<T,ST>::name() const {
  return Storage1D::stor1D_name_;
} 

template<typename T,typename ST>
inline ST Storage1D<T,ST>::size() const {
  return size_;
}

template<typename T,typename ST>
inline T* Storage1D<T,ST>::direct_access() {
  return data_;
}

template<typename T,typename ST>
inline const T* Storage1D<T,ST>::direct_access() const {
  return data_;
}

template<typename T,typename ST>
inline T& Storage1D<T,ST>::direct_access(ST i) {
  return data_[i];
}

template<typename T,typename ST>
inline T Storage1D<T,ST>::direct_access(ST i) const {
  return data_[i];
}

template<typename T,typename ST>
OPTINLINE const T& Storage1D<T,ST>::operator[](ST i) const {
#ifdef SAFE_MODE
  if (i >= size_) {

    INTERNAL_ERROR << "    invalid access on element " << i 
		   << " for Storage1D " <<  "\"" << this->name() << "\" of type " 
		   << Makros::Typename<T>()
		   << " with " << size_ << " elements. exiting." << std::endl;
    exit(1);  
  }
#endif
  return data_[i];
}


template<typename T,typename ST>
OPTINLINE T& Storage1D<T,ST>::operator[](ST i) {
#ifdef SAFE_MODE
  if (i >= size_) {

    INTERNAL_ERROR << "    invalid access on element " << i 
		   << " for Storage1D \"" << this->name() << "\" of type " 
		   << Makros::Typename<T>()
		   << " with " << size_ << " elements. exiting." << std::endl;
    exit(1);  
  }
#endif
  return data_[i];
}


template<typename T,typename ST>
void Storage1D<T,ST>::operator=(const Storage1D<T,ST>& toCopy) {

  if (size_ != toCopy.size()) {

    if (data_ != 0)
      delete[] data_;

    size_ = toCopy.size();
    data_ = new T[size_];
  }
    
  for (ST i=0; i < size_; i++) {
     data_[i] = toCopy.direct_access(i);
  }

  //this is faster for basic types but it fails for complex types where e.g. arrays have to be copied
  //memcpy(data_,toCopy.direct_access(),size_*sizeof(T));
}

#ifdef SAFE_MODE
    //for some reason g++ allows to assign an object of type T, but this does NOT produce the effect one would expect
    // => define this operator in safe mode, only to check that such an assignment is not made
template<typename T,typename ST>
void Storage1D<T,ST>::operator=(const T& invalid_object) {
  INTERNAL_ERROR << "assignment of an atomic entity to Storage1D \"" << this->name() << "\" of type " 
		   << Makros::Typename<T>()
		   << " with " << size_ << " elements. exiting." << std::endl;
}
#endif


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

template<>
void Storage1D<ushort>::operator=(const Storage1D<ushort>& toCopy);


//maintains the values of existing positions, new ones are undefined 
template<typename T,typename ST>
void Storage1D<T,ST>::resize(ST new_size) {

  if (data_ == 0) {
    data_ = new T[new_size];
  }
  else if (size_ != new_size) {
    T* new_data = new T[new_size];

    for (ST i=0; i < std::min(size_,new_size); i++)
      new_data[i] = data_[i];
    
    //this is faster for basic types but it fails for complex types where e.g. arrays have to be copied
    //memcpy(new_data,data_,std::min(size_,new_size)*sizeof(T));
    
    delete[] data_;
    data_ = new_data;
  }

  size_ = new_size;
}

template<>
void Storage1D<float>::resize(size_t new_size);

template<>
void Storage1D<double>::resize(size_t new_size);

template<>
void Storage1D<long double>::resize(size_t new_size);

template<>
void Storage1D<int>::resize(size_t new_size);

template<>
void Storage1D<uint>::resize(size_t new_size);


//maintains the values of existing positions, new ones are filled with <code> fill_value </code>
template<typename T,typename ST>
void Storage1D<T,ST>::resize(ST new_size, T fill_value) {

  if (data_ == 0) {
    data_ = new T[new_size];
    for (size_t i=0; i < new_size; i++)
      data_[i] = fill_value;
  }
  else if (size_ != new_size) {
    T* new_data = new T[new_size];

    // for (size_t i=size_; i < new_size; i++)
    //   new_data[i] = fill_value;
    if (new_size > size_)
      std::fill_n(new_data+size_,new_size-size_,fill_value);

    for (size_t i=0; i < std::min(size_,new_size); i++)
      new_data[i] = data_[i];

    
    delete[] data_;
    data_ = new_data;
  }

  size_ = new_size;
}

template<>
void Storage1D<float>::resize(size_t new_size, float fill_value);

template<>
void Storage1D<double>::resize(size_t new_size, double fill_value);

template<>
void Storage1D<long double>::resize(size_t new_size, long double fill_value);

template<>
void Storage1D<int>::resize(size_t new_size, int fill_value);

template<>
void Storage1D<uint>::resize(size_t new_size, uint fill_value);



//all elements are undefined after this operation
template<typename T,typename ST>
void Storage1D<T,ST>::resize_dirty(ST new_size) {

  if (size_ != new_size) {
    if (data_ != 0)
      delete[] data_;

    data_ = new T[new_size];
  }
  size_ = new_size;
}

/******** implementation of NamedStorage1D ***************/ 

template<typename T,typename ST>
NamedStorage1D<T,ST>::NamedStorage1D() : Storage1D<T,ST>(), name_("yyy") {}

template<typename T,typename ST>
NamedStorage1D<T,ST>::NamedStorage1D(std::string name) : Storage1D<T,ST>(), name_(name) {}
  
template<typename T,typename ST>
NamedStorage1D<T,ST>::NamedStorage1D(ST size, std::string name) : Storage1D<T,ST>(size), name_(name) {}
  
template<typename T,typename ST>
NamedStorage1D<T,ST>::NamedStorage1D(ST size, T default_value, std::string name) : 
  Storage1D<T,ST>(size,default_value), name_(name) {}

template<typename T,typename ST>
/*virtual*/ const std::string& NamedStorage1D<T,ST>::name() const {
  return name_;
}

template<typename T,typename ST>
inline void NamedStorage1D<T,ST>::operator=(const Storage1D<T,ST>& toCopy) {
  Storage1D<T,ST>::operator=(toCopy);
}

//NOTE: the name is NOT copied
template<typename T,typename ST>
inline void NamedStorage1D<T,ST>::operator=(const NamedStorage1D<T,ST>& toCopy) {
  Storage1D<T,ST>::operator=(static_cast<const Storage1D<T,ST>&>(toCopy));
}


template<typename T,typename ST>
std::ostream& operator<<(std::ostream& s, const Storage1D<T,ST>& v) {

  s << "[ ";
  for (int i=0; i < ((int) v.size()) - 1; i++)
    s << v[i] << ",";
  if (v.size() > 0)
    s << v[v.size()-1];
  s << " ]";
  
  return s;
}

template<typename T,typename ST>
bool operator==(const Storage1D<T,ST>& v1, const Storage1D<T,ST>& v2) {

  if (v1.size() != v2.size())
    return false;
  
  for (ST k=0; k < v1.size(); k++) {
    if (v1[k] != v2[k])
      return false;
  }
  return true;
}

template<typename T,typename ST>
bool operator!=(const Storage1D<T,ST>& v1, const Storage1D<T,ST>& v2) {
  return !operator==(v1,v2);
}

template<typename T,typename ST>
bool operator<(const Storage1D<T,ST>& v1, const Storage1D<T,ST>& v2) {

  for (ST k=0; k < std::min(v1.size(),v2.size()); k++) {
    if (v1[k] != v2[k])
      return (v1[k] < v2[k]);
  }
  
  return (v1.size() < v2.size());
}

template<typename T,typename ST>
bool operator<=(const Storage1D<T,ST>& v1, const Storage1D<T,ST>& v2) {

  for (ST k=0; k < std::min(v1.size(),v2.size()); k++) {
    if (v1[k] != v2[k])
      return (v1[k] < v2[k]);
  }
  
  return (v1.size() <= v2.size());
}

template<typename T,typename ST>
bool operator>(const Storage1D<T,ST>& v1, const Storage1D<T,ST>& v2) {

  return !operator<=(v1,v2);
}

template<typename T,typename ST>
bool operator>=(const Storage1D<T,ST>& v1, const Storage1D<T,ST>& v2) {

  return !operator<(v1,v2);
}

/******* implementation of FlexibleStorage1D *********/

template<typename T, typename ST>
/*static*/ const std::string FlexibleStorage1D<T,ST>::flex_stor1D_name_ = "unnamed flexible 1Dstorage";

template<typename T, typename ST>
FlexibleStorage1D<T,ST>::FlexibleStorage1D() : size_(0) {
  reserved_size_ = 4;
  data_ = new T[reserved_size_];
}

template<typename T, typename ST>
FlexibleStorage1D<T,ST>::FlexibleStorage1D(ST reserved_size)  : size_(0), reserved_size_(reserved_size) {
  data_ = new T[reserved_size_];
}

//copy constructor
template<typename T, typename ST>
FlexibleStorage1D<T,ST>::FlexibleStorage1D(const FlexibleStorage1D<T,ST>& toCopy) {

  size_ = toCopy.size();
  reserved_size_ = toCopy.reserved_size();
  
  data_ = new T[reserved_size_];
  for (uint k=0; k < toCopy.size(); k++)
    data_[k] = toCopy[k];
}

template<>
FlexibleStorage1D<uint>::FlexibleStorage1D(const FlexibleStorage1D<uint>& toCopy);

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::operator=(const FlexibleStorage1D<T,ST>& toCopy) {

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

template<typename T, typename ST>
/*virtual*/ const std::string& FlexibleStorage1D<T,ST>::name() const {
  return flex_stor1D_name_;
}

template<typename T, typename ST>
FlexibleStorage1D<T,ST>::~FlexibleStorage1D()  {
  delete [] data_;
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::set_constant(T val) {

  for (ST k=0; k < size_; k++)
    data_[k] = val;
}

template<typename T, typename ST>
ST FlexibleStorage1D<T,ST>::size() const {
  return size_;
}

template<typename T, typename ST>
ST FlexibleStorage1D<T,ST>::reserved_size() const {
  return reserved_size_;
}

template<typename T, typename ST>
ST FlexibleStorage1D<T,ST>::append(T val) {

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

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::append(Storage1D<T,ST>& toAppend) {

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

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::append(FlexibleStorage1D<T,ST>& toAppend) {

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

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::resize(ST size, bool exact_fit) {

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

template<typename T, typename ST>
OPTINLINE T& FlexibleStorage1D<T,ST>::operator[](ST i) const {

#ifdef SAFE_MODE
  if (i >= size_) {
    INTERNAL_ERROR << "    invalid access on element " << i 
		   << " for FlexibleStorage1D " <<  "\"" << this->name() << "\" of type " 
		   << Makros::Typename<T>()
		   << " with " << size_ << " (valid) elements. exiting." << std::endl;
    exit(1);  
  }
#endif
  return data_[i];
}

template<typename T, typename ST>
T* FlexibleStorage1D<T,ST>::direct_access() {

  return data_;
}

template<typename T, typename ST>
const T* FlexibleStorage1D<T,ST>::direct_access() const {

  return data_;
}

template<typename T, typename ST>
std::ostream& operator<<(std::ostream& s, const FlexibleStorage1D<T,ST>& v) {

  s << "[ ";
  for (int i=0; i < ((int) v.size()) - 1; i++)
    s << v[i] << ",";
  if (v.size() > 0)
    s << v[v.size()-1];
  s << " ]";
  
  return s;
}

template<typename T, typename ST>
bool operator==(const FlexibleStorage1D<T,ST>& v1, const FlexibleStorage1D<T,ST>& v2) {

  if (v1.size() != v2.size())
    return false;

  for (size_t k=0; k < v1.size(); k++) {
    if (v1[k] != v2[k])
      return false;
  }
  return true;
}


/***********************************/

template<typename T, typename ST>
NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D() : name_("unfs1d") {}

template<typename T, typename ST>
NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D(const std::string& name) : name_(name) {
}

template<typename T, typename ST>
NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D(ST reserved_size, const std::string& name) :
  FlexibleStorage1D<T,ST>(reserved_size), name_(name) {}

//Note: the name is NOT copied
template<typename T, typename ST>
NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D(const NamedFlexibleStorage1D<T,ST>& toCopy) : 
  FlexibleStorage1D<T,ST>(toCopy), name_("unfs1d") {
}

template<typename T, typename ST>
NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D(const FlexibleStorage1D<T,ST>& toCopy) : 
  FlexibleStorage1D<T,ST>(toCopy), name_("unfs1d") {
}

template<typename T, typename ST>
/*virtual*/ const std::string& NamedFlexibleStorage1D<T,ST>::name() const {
  return name_;
}

template<typename T, typename ST>
void NamedFlexibleStorage1D<T,ST>::operator=(const NamedFlexibleStorage1D<T,ST>& toCopy) {
  FlexibleStorage1D<T,ST>::operator=(toCopy);
}

template<typename T, typename ST>
void NamedFlexibleStorage1D<T,ST>::operator=(const FlexibleStorage1D<T,ST>& toCopy) {
  FlexibleStorage1D<T,ST>::operator=(toCopy);
}

#endif
