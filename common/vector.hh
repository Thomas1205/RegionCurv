/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef VECTOR_HH
#define VECTOR_HH

#include <cmath>
#include <limits>

#include "storage1D.hh"

namespace Math1D {

  /**********************************************/
  /*************** unnamed vector ***************/
  /**********************************************/
  template<typename T>
  class Vector : public ::Storage1D<T> {
  public:

    Vector();

    Vector(size_t size);
    
    Vector(size_t size, T default_value);

    Vector(const Vector<T>& toCopy);

    ~Vector();

    T sum() const;
        
    /*** maximal element ***/
    T max() const;
        
    /*** minimal element ***/
    T min() const;
        
    /*** maximal absolute element = l-infinity norm ***/
    T max_abs() const;
        
    /*** L2-norm of the vector ***/
    double norm() const;
 
    double sqr_norm() const;
        
    /*** L1-norm of the vector ***/
    double norm_l1() const;
        
    void operator+=(const Vector<T>& v);

    void operator-=(const Vector<T>& v);
        
    void operator*=(T constant);
        
    virtual const std::string& name() const;

    //void set_constant(T constant);
        
  protected:
    static const std::string vector_name_;
  };
    
  /**********************************************/
  /***************** named vector ****************/
  /**********************************************/
  template<typename T>
  class NamedVector : public Vector<T> {
  public:
    NamedVector();
        
    NamedVector(std::string name);
        
    NamedVector(size_t size, std::string name);
    
    NamedVector(size_t size, T default_value, std::string name);
    
    ~NamedVector();
        
    void set_name(std::string new_name);
    
    virtual const std::string& name() const;   

    inline void operator=(const Vector<T>& v);
    
    //NOTE: the name is NOT copied
    inline void operator=(const NamedVector<T>& v);


  protected:
    std::string name_;
  };
    
  /***********************************************/    
  /*************** operators *********************/
  /***********************************************/
  template<typename T>    
  Vector<T> operator+(const Vector<T>& v1, const Vector<T>& v2);

  template<typename T>    
  Vector<T> operator-(const Vector<T>& v1, const Vector<T>& v2);
    
  //scalar product of two vectors
  template<typename T>
  T operator%(const Vector<T>& v1, const Vector<T>& v2);

  template<typename T>
  Vector<T> cross(const Vector<T>& v1, const Vector<T>& v2); 

  //streaming
  template<typename T>
  std::ostream& operator<<(std::ostream& s, const Vector<T>& v);


  /******************************************** implementation *****************************************************/

  template<typename T>
  /*static*/ const std::string Vector<T>::vector_name_ = "unnamed vector";

  template<typename T>
  Vector<T>::Vector() : Storage1D<T>() {}

  template<typename T>
  Vector<T>::Vector(size_t size) : Storage1D<T>(size) {}

  template<typename T>
  Vector<T>::Vector(size_t size, T default_value) : Storage1D<T>(size,default_value) {}

  template<typename T>
  Vector<T>::Vector(const Vector<T>& toCopy) : Storage1D<T>(static_cast<const Storage1D<T>&>(toCopy)) {}

  template<typename T>
  Vector<T>::~Vector() {}

  template<typename T>
  T Vector<T>::sum() const {

    T result = 0.0;
    for (size_t i=0; i < Storage1D<T>::size(); i++) {
      result += Storage1D<T>::data_[i];
    }

    return result;
  }

  /*** maximal element ***/
  template<typename T>
  T Vector<T>::max() const {
    
//     T maxel = std::numeric_limits<T>::min();
//     for (size_t i=0; i < Storage1D<T>::size_; i++) {
//       if (Storage1D<T>::data_[i] > maxel)
// 	maxel = Storage1D<T>::data_[i];
//     }        
//    return maxel;
    return *std::max_element(Storage1D<T>::data_, Storage1D<T>::data_ + Storage1D<T>::size_);
  }

  template<>
  float Vector<float>::max() const;
    
  template<typename T>    
  T Vector<T>::min() const {
    
//     T minel = std::numeric_limits<T>::max();
//     for (size_t i=0; i < Storage1D<T>::size_; i++) {
//       if (Storage1D<T>::data_[i] < minel)
// 	minel = Storage1D<T>::data_[i];
//     }        
//     return minel;    

    return *std::min_element(Storage1D<T>::data_, Storage1D<T>::data_ + Storage1D<T>::size_);
  }

  template<>
  float Vector<float>::min() const;
        
  /*** maximal absolute element = l-infinity norm ***/
  template<typename T>    
  T Vector<T>::max_abs() const {
    
    T maxel = (T) 0;
    for (size_t i=0; i < Storage1D<T>::size_; i++) {
      T candidate = Storage1D<T>::data_[i];
      if (candidate < ((T) 0))
	candidate *= (T) -1;
      if (candidate > maxel)
	maxel = candidate;
    }
        
    return maxel;
  }
        
  /*** L2-norm of the vector ***/
  template<typename T>   
  double Vector<T>::norm() const {
    
    double result = 0.0;
    for (size_t i=0; i < Storage1D<T>::size_; i++) {
      const double cur = (double) Storage1D<T>::data_[i];
      result += cur*cur;
    }
        
    return sqrt(result);
  }
    
  template<typename T>   
  double Vector<T>::sqr_norm() const {
    
    double result = 0.0;
    for (size_t i=0; i < Storage1D<T>::size_; i++) {
      const double cur = (double) Storage1D<T>::data_[i];
      result += cur*cur;
    }
        
    return result;
  }
        
  /*** L1-norm of the vector ***/
  template<typename T>   
  double Vector<T>::norm_l1() const {
    
    double result = 0.0;
    for (size_t i=0; i < Storage1D<T>::size_; i++) {
      result += std::abs(Storage1D<T>::data_[i]);
    }
    
    return result;    
  }

  //template specialization (note that uints are never negative!)
  template<>
  double Vector<uint>::norm_l1() const;

  //template specialization (note that uchars are never negative!)
  template<>
  double Vector<uchar>::norm_l1() const;

  //template specialization (note that ushorts are never negative!)
  template<>
  double Vector<ushort>::norm_l1() const;

        
  template<typename T>       
  void Vector<T>::operator+=(const Vector<T>& v) {
    
    const size_t size = Storage1D<T>::size_;

    if (v.size() != size) {
      INTERNAL_ERROR << "   cannot add vector \"" << v.name() << "\" to vector \"" 
		     << this->name() << "\":" << std::endl
		     << "   sizes " << v.size() << " and " << this->size() << " mismatch. Exiting..."
		     << std::endl;   
      exit(0);
    }
    
    size_t i;
    for (i=0; i < size; i++)
      Storage1D<T>::data_[i] += v.direct_access(i);
  }

  template<typename T>
  void Vector<T>::operator-=(const Vector<T>& v) {
    
    const size_t size = Storage1D<T>::size_;

    if (v.size() != size) {
      INTERNAL_ERROR << "   cannot subtract vector \"" << v.name() << "\" from vector \"" 
		     << this->name() << "\":" << std::endl
		     << "   sizes " << v.size() << " and " << this->size() << " mismatch. Exiting..."
		     << std::endl;   
      exit(0);
    }
    
    size_t i;
    for (i=0; i < size; i++)
      Storage1D<T>::data_[i] -= v.direct_access(i);
  }
        
  template<typename T>       
  void Vector<T>::operator*=(T constant) {
    
    const size_t size = Storage1D<T>::size_;
    
    size_t i;
    for (i=0; i < size; i++)
      Storage1D<T>::data_[i] *= constant;
  }

  template<>
  void Vector<float>::operator*=(const float scalar);

  template<>
  void Vector<double>::operator*=(const double scalar);
 
  template<typename T>
  /*virtual*/ const std::string& Vector<T>::name() const {
    return Vector<T>::vector_name_;
  }   

//   template<typename T>
//   void Vector<T>::set_constant(T constant) {

//     std::fill_n(Storage1D<T>::data_,Storage1D<T>::size_,constant); //experimental result: fill_n is usually faster
//   }

  /************** implementation of NamedVector **********/

  template<typename T>
  NamedVector<T>::NamedVector() : Vector<T>(), name_("yyy") {}
       
  template<typename T>
  NamedVector<T>::NamedVector(std::string name) : Vector<T>(), name_(name) {}
        
  template<typename T>
  NamedVector<T>::NamedVector(size_t size, std::string name) : Vector<T>(size), name_(name) {}
    
  template<typename T>
  NamedVector<T>::NamedVector(size_t size, T default_value, std::string name) :
    Vector<T>(size,default_value), name_(name) {}
    
  template<typename T>
  NamedVector<T>::~NamedVector() {}
        
  template<typename T>
  void NamedVector<T>::set_name(std::string new_name) {
    name_ = new_name;
  }

  template<typename T>
  /*virtual*/ const std::string& NamedVector<T>::name() const {
    return name_;
  }   

  template<typename T>
  inline void NamedVector<T>::operator=(const Vector<T>& v) {    
    Storage1D<T>::operator=(v);
  }

  //NOTE: the name is NOT copied
  template<typename T>
  inline void NamedVector<T>::operator=(const NamedVector<T>& v) {    
    Storage1D<T>::operator=(v);
  }

  /************** implementation of stand-alone routines **********************/

  template<typename T>    
  Vector<T> operator+(const Vector<T>& v1, const Vector<T>& v2) {

    const size_t size = v1.size();

    if (size != v2.size()) {
      INTERNAL_ERROR << "     cannot add vectors \"" << v1.name() << "\" and \""
		     << v2.name() << "\":" << std::endl
		     << "    sizes " << v1.size() << " and " << v2.size() << " mismatch. Exiting..."
		     << std::endl;
      exit(1);
    }
        
    Vector<T> result(size);
    size_t i;
    for (i = 0; i < size; i++)
      result.direct_access(i) = v1.direct_access(i) + v2.direct_access(i);

    return result;
  }

  template<typename T>    
  Vector<T> operator-(const Vector<T>& v1, const Vector<T>& v2) {

    const size_t size = v1.size();

    if (size != v2.size()) {
      INTERNAL_ERROR << "     cannot subtract vector \"" << v2.name() << "\" from \""
		     << v1.name() << "\":" << std::endl
		     << "    sizes " << v2.size() << " and " << v1.size() << " mismatch. Exiting..."
		     << std::endl;
      exit(1);
    }
        
    Vector<T> result(size);
    size_t i;
    for (i = 0; i < size; i++)
      result.direct_access(i) = v1.direct_access(i) - v2.direct_access(i);

    return result;
  }

  template<typename T>
  T operator%(const Vector<T>& v1, const Vector<T>& v2) {
    
    const size_t size = v1.size();

    if (size != v2.size()) {
      INTERNAL_ERROR << "     cannot compute scalar product of vectors \"" 
		     << v1.name() << "\" and \"" << v2.name() << "\":" << std::endl
		     << "      sizes " << v1.size() << " and " << v2.size() << " mismatch. exiting." 
		     << std::endl;
      exit(1); 
    }
        
    T result = (T) 0;
    size_t i;
    for (i=0; i < size; i++)
      result += v1.direct_access(i)*v2.direct_access(i);
            
    return result;
  }


  template<typename T>
  std::ostream& operator<<(std::ostream& s, const Vector<T>& v) {

    s << "[ ";
    for (int i=0; i < ((int) v.size()) - 1; i++)
      s << v[i] << ",";
    if (v.size() > 0)
      s << v[v.size()-1];
    s << " ]";

    return s;
  }


  template<typename T>
  Vector<T> cross(const Vector<T>& v1, const Vector<T>& v2) {
    
    if (v1.size() != 3 || v2.size() != 3) {
      INTERNAL_ERROR << "      the cross product is only defined for vectors of size 3." << std::endl
		     << "                  here applied for vectors of size " << v1.size() << " and "
		     << v2.size() << ". exiting." << std::endl;
      exit(1);    
    } 
    
    Vector<T> result(3);
    result[0] = v1[1]*v2[2] - v1[2]*v2[1];
    result[1] = v1[2]*v2[0] - v1[0]*v2[2];
    result[2] = v1[0]*v2[1] - v1[1]*v2[0];
        
    return result;
  } 


}//end of namespace Math1D


#endif
