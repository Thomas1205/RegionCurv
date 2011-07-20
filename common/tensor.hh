/**** written by Thomas Schoenemann as a private person without employment, October 2009 ***/

#ifndef TENSOR_HH
#define TENSOR_HH

#include <cmath>
#include <fstream>
#include <cstring>
#include <algorithm>

#include "storage3D.hh"
#include "fileio.hh"

namespace Math3D {

  
  /*** Tensor class, i.e. mathematical operations need to be defined on T ****/
  template<typename T>
  class Tensor : public Storage3D<T> {
  public:
    
    Tensor();
    
    Tensor(size_t xDim, size_t yDim, size_t zDim);

    Tensor(size_t xDim, size_t yDim, size_t zDim, T default_value);
    
    ~Tensor();

    virtual const std::string& name() const;

    void operator+=(const Tensor<T>& toAdd);

    void operator-=(const Tensor<T>& toSub);
    
    void operator*=(const T scalar); 
    
    double norm() const;

    double sqr_norm() const;

    inline double norm(size_t x, size_t y) const;

    inline double sqr_norm(size_t x, size_t y) const;

    inline T sum(size_t x, size_t y) const;

    inline T min(size_t x, size_t y) const;
    
    T sum() const;

    T max() const;
    
    T min() const;

    T max_abs() const;
    
    T max(size_t z) const;

    T min(size_t z) const;

    T max_abs(size_t z) const;

    double max_vector_norm() const;

    void set_constant(const T constant);

    //returns if the operation was successful
    bool savePPM(std::string filename, size_t max_intensity, bool fit_to_range = true) const;

  protected:
    static const std::string tensor_name_; 
  };

  /**** named Tensor class ****/
  template <typename T>
  class NamedTensor : public Tensor<T> {
  public:

    NamedTensor();
    
    NamedTensor(std::string name);

    NamedTensor(size_t xDim, size_t yDim, size_t zDim, std::string name);

    NamedTensor(size_t xDim, size_t yDim, size_t zDim, T default_value, std::string name);
    
    ~NamedTensor();

    virtual const std::string& name() const;

    void set_name(std::string name);

    inline void operator=(const Tensor<T>& toCopy);
    
    //NOTE: the name is NOT copied
    inline void operator=(const NamedTensor<T>& toCopy);

  protected:
    std::string name_;
  };


  /**** stand-alone operators ****/
  template<typename T>
  Tensor<T> operator+(const Tensor<T>& v1, const Tensor<T>& v2);

  template<typename T>
  Tensor<T> operator-(const Tensor<T>& v1, const Tensor<T>& v2);

/********************************** implementation **********************************/

  /*** implementation of (unnamed) Tensor ***/
  template<typename T>
  /*static*/ const std::string Tensor<T>::tensor_name_ = "unnamed tensor";

  template<typename T>
  Tensor<T>::Tensor() : Storage3D<T>() {}
  
  template<typename T>
  Tensor<T>::Tensor(size_t xDim, size_t yDim, size_t zDim) : Storage3D<T>(xDim,yDim,zDim) {}
  
  template<typename T>
  Tensor<T>::Tensor(size_t xDim, size_t yDim, size_t zDim, T default_value) :
    Storage3D<T>(xDim,yDim,zDim,default_value) {}
  
  template<typename T>
  Tensor<T>::~Tensor() {}

  template<typename T>
  /*virtual*/ const std::string& Tensor<T>::name() const {
    return tensor_name_;
  }

  template<typename T>
  void Tensor<T>::operator+=(const Tensor<T>& toAdd) {
    if (toAdd.xDim() != Storage3D<T>::xDim_ 
	|| toAdd.yDim() != Storage3D<T>::yDim_ 
	|| toAdd.zDim() != Storage3D<T>::zDim_) {
      INTERNAL_ERROR << "    illegal addition of tensor \"" << toAdd.name() << "\" to tensor \"" 
		     << this->name() << "\":" << std::endl
		     << "    sizes " << toAdd.xDim() << "x" << toAdd.yDim() << "x" << toAdd.zDim()
		     << " and " << Storage3D<T>::xDim_ << "x" << Storage3D<T>::yDim_ << "x" 
		     << Storage3D<T>::zDim_ << " are incompatible. Exiting..."
		     << std::endl;
      exit(1);
    }

    const size_t size = Storage3D<T>::size_;

    size_t i;
    for (i=0; i < size; i++)
      Storage3D<T>::data_[i] += toAdd.direct_access(i);
  }
  
  template<typename T>
  void Tensor<T>::operator-=(const Tensor<T>& toSub) {
    if (toSub.xDim() != Storage3D<T>::xDim_ 
	|| toSub.yDim() != Storage3D<T>::yDim_ 
	|| toSub.zDim() != Storage3D<T>::zDim_) {
      INTERNAL_ERROR << "    illegal subtraction of tensor \"" << toSub.name() << "\" from tensor \"" 
		     << this->name() << "\":" << std::endl
		     << "    sizes " << toSub.xDim() << "x" << toSub.yDim() << "x" << toSub.zDim()
		     << " and " << Storage3D<T>::xDim_ << "x" << Storage3D<T>::yDim_ << "x" 
		     << Storage3D<T>::zDim_ << " are incompatible. Exiting..."
		     << std::endl;
      exit(1);
    }

    const size_t size = Storage3D<T>::size_;

    size_t i;
    for (i=0; i < size; i++)
      Storage3D<T>::data_[i] -= toSub.direct_access(i);
  }

  template<typename T>
  void Tensor<T>::operator*=(const T scalar) {
    
    const size_t size = Storage3D<T>::size_;

    size_t i;
    for (i=0; i < size; i++)
      Storage3D<T>::data_[i] *= scalar;
  }

  template<>
  void Tensor<float>::operator*=(const float scalar);

  template<>
  void Tensor<double>::operator*=(const double scalar);

  template<typename T>
  double Tensor<T>::norm() const {

    double result = 0.0;
    T temp;
    for (size_t i=0; i < Storage3D<T>::size_; i++) {
      temp = Storage3D<T>::data_[i];
      result += temp*temp;
    }
    return sqrt(result);
  }

  template<typename T>
  double Tensor<T>::sqr_norm() const {

    double result = 0.0;
    T temp;
    for (size_t i=0; i < Storage3D<T>::size_; i++) {
      temp = Storage3D<T>::data_[i];
      result += temp*temp;
    }
    return result;
  }

  template<typename T>
  inline T Tensor<T>::sum(size_t x, size_t y) const {

    T result = (T) 0;
    size_t offs = (y*Storage3D<T>::xDim_+x)*Storage3D<T>::zDim_;

    for (size_t z=0; z < Storage3D<T>::zDim_; z++) {
      result += Storage3D<T>::data_[offs+z];
    }
      
    return result;
  }

  template<typename T>
  inline T Tensor<T>::min(size_t x, size_t y) const {

    return *std::min_element(Storage3D<T>::data_+(y*Storage3D<T>::xDim_+x)*Storage3D<T>::zDim_,
			     Storage3D<T>::data_+(y*Storage3D<T>::xDim_+x+1)*Storage3D<T>::zDim_);
  }

  
  template<typename T>
  double Tensor<T>::norm(size_t x, size_t y) const {

    double result = 0.0;
    T temp;
    size_t offs = (y*Storage3D<T>::xDim_+x)*Storage3D<T>::zDim_;

    for (size_t z=0; z < Storage3D<T>::zDim_; z++) {
      temp = Storage3D<T>::data_[offs+z];
      result += temp*temp;
    }
      
    return sqrt(result);
  }

  template<typename T>
  double Tensor<T>::sqr_norm(size_t x, size_t y) const {

    double result = 0.0;
    T temp;
    size_t offs = (y*Storage3D<T>::xDim_+x)*Storage3D<T>::zDim_;

    for (size_t z=0; z < Storage3D<T>::zDim_; z++) {
      temp = Storage3D<T>::data_[offs+z];
      result += temp*temp;
    }
      
    return result;
  }

  template<typename T>
  T Tensor<T>::sum() const {

    T result = 0.0;
    for (size_t i=0; i < Storage3D<T>::size(); i++) {
      result += Storage3D<T>::data_[i];
    }

    return result;
  }

  template<typename T>
  T Tensor<T>::max() const {
    
//     T max_el = std::numeric_limits<T>::min();
//     for (size_t i=0; i < Storage3D<T>::size_; i++)
//       max_el = std::max(Storage3D<T>::data_[i],max_el);
//     return max_el;
    
    return *std::max_element(Storage3D<T>::data_,Storage3D<T>::data_+Storage3D<T>::size_);
  }

  template<>
  float Tensor<float>::max() const;
    
  template<typename T>
  T Tensor<T>::min() const {
    
//     T min_el = std::numeric_limits<T>::max();
//     for (size_t i=0; i < Storage3D<T>::size_; i++)
//       min_el = std::min(Storage3D<T>::data_[i],min_el);
//     return min_el;
    
    return *std::min_element(Storage3D<T>::data_,Storage3D<T>::data_+Storage3D<T>::size_);
  }
  
  template<>
  float Tensor<float>::min() const;

  template<typename T>
  T Tensor<T>::max_abs() const {

    T max_el = std::numeric_limits<T>::min();
    for (size_t i=0; i < Storage3D<T>::size_; i++) {
      T cur = Storage3D<T>::data_[i];
      if (cur < 0)
	cur *= -1;
      max_el = std::max(cur,max_el);
    }
   
    return max_el;
  }

  template<typename T>
  T Tensor<T>::max(size_t z) const {

    T max_el = std::numeric_limits<T>::min();
    for (size_t i=z; i < Storage3D<T>::size_; i+=z)
      max_el = std::max(Storage3D<T>::data_[i],max_el);
      
    return max_el;
  }
    
  template<typename T>
  T Tensor<T>::min(size_t z) const {

    T min_el = std::numeric_limits<T>::max();
    for (size_t i=z; i < Storage3D<T>::size_; i+=z)
      min_el = std::min(Storage3D<T>::data_[i],min_el);

    return min_el;
  }

  template<typename T>
  T Tensor<T>::max_abs(size_t z) const {

    T max_el = std::numeric_limits<T>::min();
    for (size_t i=z; i < Storage3D<T>::size_; i+=z) {
      T cur = Storage3D<T>::data_[i];
      if (cur < 0)
	cur *= (T) -1;
      max_el = std::max(cur,max_el);
    }
   
    return max_el;
  }

  template<typename T>
  inline double Tensor<T>::max_vector_norm() const {

    double max_norm = 0.0;

    for (uint y=0; y < Storage3D<T>::yDim_; y++) {
      for (uint x=0; x < Storage3D<T>::xDim_; x++) {
	double cur_norm = 0.0;
	size_t base = (y*Storage3D<T>::xDim_+x)*Storage3D<T>::zDim_;
	for (uint z=0; z < Storage3D<T>::zDim_; z++) {
	  T cur_datum = Storage3D<T>::data_[base+z];
	  cur_norm += cur_datum*cur_datum;
	}
	cur_norm = sqrt(cur_norm);
	max_norm = std::max(max_norm,cur_norm);
      }
    }
    return max_norm;
  }

  template<typename T>
  void Tensor<T>::set_constant(const T constant) {

    const size_t size = Storage3D<T>::size_;

    size_t i;
    for (i=0; i < size; i++)
      Storage3D<T>::data_[i] = constant;
  }

  template<typename T>
  bool Tensor<T>::savePPM(std::string filename, size_t max_intensity, bool fit_to_range) const {
    
    size_t zDim = this->zDim();
    
    if (zDim != 1 && zDim != 3) {
      std::cerr << "WARNING: cannot save a tensor with " << zDim << " channels as either pgm or ppm. Operation aborted."
		<< std::endl;
      return false;
    }

    std::ofstream of(filename.c_str());

    if (!of.is_open()) {
      IO_ERROR << " while saving PGM: could not write file \"" << filename 
	       << "\". Please check if the path is correct." << std::endl;
      return false;
    }

    if (zDim == 1)
      of << "P5\n"; 
    else
      of << "P6\n";
    of << Storage3D<T>::xDim_ << " " << Storage3D<T>::yDim_ << "\n" << max_intensity; 

    //Reopen in binary mode to avoid silent conversion from '\n' to "\r\n" under Windows
    of.close();
    of.open(filename.c_str(), std::ios::binary | std::ios::app);
    of << '\n';

    for (size_t i=0; i < Storage3D<T>::size_; i++) {

      if (max_intensity < 256) {
	T cur_datum = Storage3D<T>::data_[i];
	if (fit_to_range) {
	  cur_datum = std::max(cur_datum,(T) 0);
	  cur_datum = std::min(cur_datum,(T) max_intensity);
	}
	uchar c = uchar(cur_datum);
	of << c;
      }
      else {
	TODO("handle sizes > 255 when saving PPMs (or PGMs)");
      }
    }

    return true;
  }
  
    
  /*** implementation of NamedTensor ***/

  template<typename T>
  NamedTensor<T>::NamedTensor() : Tensor<T>(), name_("yyy") {}
  
  template<typename T>
  NamedTensor<T>::NamedTensor(std::string name) : Tensor<T>(), name_(name) {}
  
  template<typename T>
  NamedTensor<T>::NamedTensor(size_t xDim, size_t yDim, size_t zDim, std::string name) :
    Tensor<T>(xDim,yDim,zDim), name_(name) {}
  
  template<typename T>
  NamedTensor<T>::NamedTensor(size_t xDim, size_t yDim, size_t zDim, T default_value, std::string name) :
    Tensor<T>(xDim,yDim,zDim,default_value), name_(name) {}
  
  template<typename T>
  NamedTensor<T>::~NamedTensor() {}
  
  template<typename T>
  /*virtual*/ const std::string& NamedTensor<T>::name() const {
    return name_;
  }

  template<typename T>
  void NamedTensor<T>::set_name(std::string name) {
    name_ = name;
  }
  
  template<typename T>
  inline void NamedTensor<T>::operator=(const Tensor<T>& toCopy) {
    Tensor<T>::operator=(toCopy);
  }
  
  //NOTE: the name is NOT copied
  template<typename T>
  inline void NamedTensor<T>::operator=(const NamedTensor<T>& toCopy) {
    Tensor<T>::operator=(toCopy);
  }

  /*** implementation of stand-alone operators and routines ****/
  template<typename T>
  Tensor<T> operator+(const Tensor<T>& v1, const Tensor<T>& v2) {

    if (v1.xDim() != v2.xDim() || v1.yDim() != v2.yDim() || v1.zDim() != v2.zDim() ) {

      INTERNAL_ERROR << "    cannot add vectors \"" << v1.name() << "\" and \"" << v2.name() 
		     << "\":" << std::endl
		     << "    sizes " <<  v1.xDim() << "x" << v1.yDim() << "x" << v1.zDim() << " and "  
		     << v2.xDim() << "x" << v2.yDim() << "x" << v2.zDim()  << " are incompatible. Exiting..." << std::endl;
      exit(1);
    }

    Tensor<T> result(v1.xDim(), v1.yDim(), v1.zDim());
    const size_t size = v1.size();

    size_t i;
    for (i=0; i < size; i++)
      result.direct_access(i) = v1.direct_access(i)+v2.direct_access(i);

    return result;
  }
  
  template<typename T>
  Tensor<T> operator-(const Tensor<T>& v1, const Tensor<T>& v2) {

    if (v1.xDim() != v2.xDim() || v1.yDim() != v2.yDim() || v1.zDim() != v2.zDim() ) {

      INTERNAL_ERROR << "    cannot subtract vector \"" << v2.name() << "\" from \"" << v1.name() 
		     << "\":" << std::endl
		     << "    sizes " <<  v2.xDim() << "x" << v2.yDim() << "x" << v2.zDim() << " and "  
		     << v1.xDim() << "x" << v1.yDim() << "x" << v1.zDim()  << " are incompatible. Exiting..." << std::endl;
      exit(1);
    }

    Tensor<T> result(v1.xDim(), v1.yDim(), v1.zDim());
    const size_t size = v1.size();
    
    size_t i;
    for (i=0; i < size; i++)
      result.direct_access(i) = v1.direct_access(i)-v2.direct_access(i);

    return result;
  }


} //end of namespace Math3D


#endif
