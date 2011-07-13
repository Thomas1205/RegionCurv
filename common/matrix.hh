/**** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef MATRIX_HH
#define MATRIX_HH

#include <cmath>
#include <fstream>

#include "storage2D.hh"
#include "vector.hh"

namespace Math2D {

  /**************** unnamed Matrix ***********/

  //matrix class, i.e. mathematical and streaming operations need to be defined on T
  template<typename T>
  class Matrix : public ::Storage2D<T> {
  public:

    /*---- constructors -----*/
    Matrix();

    Matrix(size_t xDim, size_t yDim);

    Matrix(size_t xDim, size_t yDim, T default_value);

    /*---- destructor ----*/
    ~Matrix();

    virtual const std::string& name() const;

    void set_constant(T constant);

    T sum() const;
    
    /*** maximal element ***/
    T max() const;
        
    /*** minimal element ***/
    T min() const;
        
    /*** maximal absolute element = l-infinity norm ***/
    T max_abs() const;
        
    /*** L2-norm of the matrix ***/
    double norm() const;
 
    /*** squared L2-norm ***/
    double sqr_norm() const;
        
    /*** L1-norm of the matrix ***/
    double norm_l1() const;
    
    //---- mathematical operators ----

    //addition of another matrix of equal dimensions
    void operator+=(const Matrix<T>& toAdd);

    void operator-=(const Matrix<T>& toSub);

    //multiplication with a scalar
    void operator*=(const T scalar);

    //returns if the operation was successful
    bool savePGM(const std::string& filename, size_t max_intensity, bool fit_to_range = true) const;

  protected:
    static const std::string matrix_name_;
  };

  /******************** Named Matrix ************************/

  template <typename T>
  class NamedMatrix : public Matrix<T> {
  public:
    NamedMatrix();

    NamedMatrix(std::string name);

    NamedMatrix(size_t xDim, size_t yDim, std::string name);

    NamedMatrix(size_t xDim, size_t yDim, T default_value, std::string name);

    ~NamedMatrix();

    inline void operator=(const Matrix<T>& toCopy);

    //NOTE: does NOT copy the name
    inline void operator=(const NamedMatrix<T>& toCopy);

    virtual const std::string& name() const;

    void set_name(std::string new_name);

  protected:
    std::string name_;
  };

  /***************** stand-alone operators and routines ********************/

  template<typename T>
  Matrix<T> operator+(const Matrix<T>& m1, const Matrix<T>& m2);

  template<typename T>
  Matrix<T> operator*(const Matrix<T>& m1, const Matrix<T>& m2); 

  //streaming
  template <typename T>
  std::ostream& operator<<(std::ostream& s, const Matrix<T>& m);   

  template<typename T>
  Matrix<T> transpose(const Matrix<T>& m);
  
  template<typename T>
  Math1D::Vector<T> operator*(const Matrix<T>& m, const Math1D::Vector<T>& v);


  /******************************** implementation ********************************/

  /****** implementation of (unnamed) Matrix ******/

  template<typename T>
  /*static*/ const std::string Matrix<T>::matrix_name_ = "unnamed matrix";

  template<typename T>
  Matrix<T>::Matrix() : Storage2D<T>() {}

  template<typename T>
  Matrix<T>::Matrix(size_t xDim, size_t yDim) : Storage2D<T>(xDim, yDim)  {}

  template<typename T>
  Matrix<T>::Matrix(size_t xDim, size_t yDim, T default_value) : Storage2D<T>(xDim, yDim) {
    for (size_t i=0; i < Storage2D<T>::size_; i++)
      Storage2D<T>::data_[i] = default_value;
  }
  
  template<typename T>
  Matrix<T>::~Matrix() {}

  template<typename T>
  /*virtual*/ const std::string& Matrix<T>::name() const {
    return Matrix<T>::matrix_name_;
  }

  template<typename T>
  void Matrix<T>::set_constant(T constant) {
    for (size_t i=0; i < Storage2D<T>::size_; i++)
      Storage2D<T>::data_[i] = constant;
  }

  template<typename T>
  T Matrix<T>::sum() const {

    T result = (T) 0;
    for (size_t i=0; i < Storage2D<T>::size_; i++)
      result += Storage2D<T>::data_[i];

    return result;
  }

  /*** maximal element ***/
  template<typename T>
  T Matrix<T>::max() const {
    
//     T maxel = std::numeric_limits<T>::min();
//     for (size_t i=0; i < Storage2D<T>::size_; i++) {
//       if (Storage2D<T>::data_[i] > maxel)
// 	maxel = Storage2D<T>::data_[i];
//     }
//     return maxel;

    return *std::max_element(Storage2D<T>::data_,Storage2D<T>::data_+Storage2D<T>::size_);
  }

  template<>
  float Matrix<float>::max() const;
    
  /*** minimal element ***/
  template<typename T>    
  T Matrix<T>::min() const {
    
//     T minel = std::numeric_limits<T>::max();
//     for (size_t i=0; i < Storage2D<T>::size_; i++) {
//       if (Storage2D<T>::data_[i] < minel)
// 	minel = Storage2D<T>::data_[i];
//     }
//     return minel;    

    return *std::min_element(Storage2D<T>::data_,Storage2D<T>::data_+Storage2D<T>::size_);
  }

  template<>
  float Matrix<float>::min() const;
        
  /*** maximal absolute element = l-infinity norm ***/
  template<typename T>    
  T Matrix<T>::max_abs() const {
    
    T maxel = (T) 0;
    for (size_t i=0; i < Storage2D<T>::size_; i++) {
      T candidate = Storage2D<T>::data_[i];
      if (candidate < ((T) 0))
	candidate *= (T) -1;
      if (candidate > maxel)
	maxel = candidate;
    }
        
    return maxel;
  }
        
  /*** L2-norm of the matrix ***/
  template<typename T>   
  double Matrix<T>::norm() const {
    
    double result = 0.0;
    for (size_t i=0; i < Storage2D<T>::size_; i++) {
      double cur = (double) Storage2D<T>::data_[i];
      result += cur*cur;
    }
        
    return sqrt(result);
  }
    
  template<typename T>   
  double Matrix<T>::sqr_norm() const {
    
    double result = 0.0;
    for (size_t i=0; i < Storage2D<T>::size_; i++) {
      double cur = (double) Storage2D<T>::data_[i];
      result += cur*cur;
    }
        
    return result;
  }
        
  /*** L1-norm of the matrix ***/
  template<typename T>   
  double Matrix<T>::norm_l1() const {
    
    double result = 0.0;
    for (size_t i=0; i < Storage2D<T>::size_; i++) {
      result += std::abs(Storage2D<T>::data_[i]);
    }
    
    return result;    
  }

  //template specialization (note that uints are never negative, so there is no need to call abs())
  template<>   
  double Matrix<uint>::norm_l1() const;

  //template specialization (note that ushorts are never negative, so there is no need to call abs())
  template<>   
  double Matrix<ushort>::norm_l1() const;

  //template specialization (note that uchars are never negative, so there is no need to call abs())
  template<>   
  double Matrix<uchar>::norm_l1() const;


  //addition of another matrix of equal dimensions
  template<typename T>
  void Matrix<T>::operator+=(const Matrix<T>& toAdd) {
    
    if (toAdd.xDim() != Storage2D<T>::xDim_ || toAdd.yDim() != Storage2D<T>::yDim_) {
      INTERNAL_ERROR << "    dimension mismatch in matrix addition(+=): (" 
		     << Storage2D<T>::xDim_ << "," << Storage2D<T>::yDim_ << ") vs. ("
		     << toAdd.xDim() << "," << toAdd.yDim() << ")." << std::endl;
      std::cerr << "     When adding matrix\"" << toAdd.name() << "\" to  matrix \""
		<< this->name() << "\". Exiting" << std::endl;
      exit(1);
    }
    else {
      assert(Storage2D<T>::size_ == Storage2D<T>::xDim_*Storage2D<T>::yDim_);
      for (size_t i=0; i < Storage2D<T>::size_; i++)
	Storage2D<T>::data_[i] += toAdd.direct_access(i);
    }
  }

  template<typename T>
  void Matrix<T>::operator-=(const Matrix<T>& toSub) {
    
    if (toSub.xDim() != Storage2D<T>::xDim_ || toSub.yDim() != Storage2D<T>::yDim_) {
      INTERNAL_ERROR << "    dimension mismatch in matrix subtraction(-=): (" 
		     << Storage2D<T>::xDim_ << "," << Storage2D<T>::yDim_ << ") vs. ("
		     << toSub.xDim() << "," << toSub.yDim() << ")." << std::endl;
      std::cerr << "     When subtracting matrix\"" << toSub.name() << "\" from  matrix \""
		<< this->name() << "\". Exiting" << std::endl;
      exit(1);
    }
    else {
      assert(Storage2D<T>::size_ == Storage2D<T>::xDim_*Storage2D<T>::yDim_);
      for (size_t i=0; i < Storage2D<T>::size_; i++)
	Storage2D<T>::data_[i] -= toSub.direct_access(i);
    }
  }

  //@returns if the operation was successful
  template<typename T>
  bool Matrix<T>::savePGM(const std::string& filename, size_t max_intensity, bool fit_to_range) const {

    std::ofstream of(filename.c_str());

    if (!of.is_open()) {
      IO_ERROR << " while saving PGM: could not write file \"" << filename 
	       << "\". Please check if the path is correct." << std::endl;
      return false;
    }

    of << "P5\n" << Storage2D<T>::xDim_ << " " << Storage2D<T>::yDim_ << "\n" << max_intensity << "\n";

    //Reopen in binary mode to avoid silent conversion from '\n' to "\r\n" under Windows
    of.close();
    of.open(filename.c_str(), std::ios::binary | std::ios::app);

    for (size_t i=0; i < Storage2D<T>::size_; i++) {
      
      if (max_intensity < 256) {
	T cur_datum = Storage2D<T>::data_[i];
	if (fit_to_range) {
	  cur_datum = std::max(cur_datum,(T) 0);
	  cur_datum = std::min(cur_datum,(T) max_intensity);
	}
	uchar c = cur_datum;
	of << c;
      }
      else {
	TODO("handle max_intensity > 255 when saving PGMs");
      }
    }

    return true;
  }


  /***************** implementation of Named Matrix ***********************/

  template<typename T>
  NamedMatrix<T>::NamedMatrix() : Matrix<T>(), name_("zzz") {}
  
  template<typename T>
  NamedMatrix<T>::NamedMatrix(std::string name) : Matrix<T>(), name_(name) {}
  
  template<typename T>
  NamedMatrix<T>::NamedMatrix(size_t xDim, size_t yDim, std::string name) : 
    Matrix<T>(xDim, yDim), name_(name) {}
  
  template<typename T>
  NamedMatrix<T>::NamedMatrix(size_t xDim, size_t yDim, T default_value, std::string name) :
    Matrix<T>(xDim,yDim,default_value), name_(name) {}
  
  template<typename T>
  NamedMatrix<T>::~NamedMatrix() {}

  template<typename T>
  inline void NamedMatrix<T>::operator=(const Matrix<T>& toCopy) {
    Matrix<T>::operator=(toCopy);
  }

  //NOTE: does NOT copy the name
  template<typename T>
  inline void NamedMatrix<T>::operator=(const NamedMatrix<T>& toCopy) {
    Matrix<T>::operator=(toCopy);
  }
  
  template<typename T>
  /*virtual*/ const std::string& NamedMatrix<T>::name() const {
    return name_;
  }

  template<typename T>
  void NamedMatrix<T>::set_name(std::string new_name) {
    name_ = new_name;
  }

  /***************** implementation of stand-alone operators **************/

  //multiplication with a scalar
  template<typename T>
  void Matrix<T>::operator*=(const T scalar) {

    assert(Storage2D<T>::size_ == Storage2D<T>::xDim_*Storage2D<T>::yDim_);
    size_t i;
    for (i=0; i < Storage2D<T>::size_; i++)
      Storage2D<T>::data_[i] *= scalar;    
  }

  template<>
  void Matrix<float>::operator*=(const float scalar);

  template<>
  void Matrix<double>::operator*=(const double scalar);

  //implementation of stand-alone operators
  template<typename T>
  Matrix<T> operator+(const Matrix<T>& m1, const Matrix<T>& m2) {

    if (m1.xDim != m2.xDim || m1.yDim != m2.yDim) {

      INTERNAL_ERROR << "     dimension mismatch in matrix addition(+): (" 
		     << m1.xDim() << "," << m1.yDim() << ") vs. ("
		     << m2.xDim() << "," << m2.yDim() << ")." << std::endl;
      std::cerr << "     When adding matrices \"" << m1.name() << "\" and\"" 
		<< m2.name() << "\". Exiting..." << std::endl;

      exit(1);
    }

    Matrix<T> result(m1.xDim(),m1.yDim());
    size_t i;
    for (i=0; i < m1.size(); i++)
      result.direct_access(i) = m1.value(i) + m2.value(i);

    return result;
  }

  template<typename T>
  Matrix<T> operator*(const Matrix<T>& m1, const Matrix<T>& m2) {

    if (m1.xDim() != m2.yDim()) {
      INTERNAL_ERROR << "     dimension mismatch in matrix multiplication(*): ("
		     << m1.xDim() << "," << m1.yDim() << ") vs. ("
		     << m2.xDim() << "," << m2.yDim() << ")." << std::endl;
      std::cerr << "     When multiplying matrices \"" << m1.name() << "\" and \"" 
		<< m2.name() << "\". Exiting..." << std::endl;
      exit(1);
    }
    
    const size_t xDim = m2.xDim();
    const size_t yDim = m1.yDim();
    const size_t zDim = m1.xDim();

    Matrix<T> result(xDim,yDim);
    size_t y,x,z;
    T sum;

    for (y=0; y < yDim; y++) {
      for (x=0; x < xDim; x++) {
	
	sum = (T) 0;
	for (z=0; z < zDim; z++) {
	  //sum += m1(z,y) * m2(x,z);
	  sum += m1.direct_access(y*zDim+z) * m2.direct_access(z*xDim+x);
	}

	//result(x,y) = sum;
	result.direct_access(y*xDim+x) = sum;
      }
    }

    return result;
  }

  // streaming
  template <typename T>
  std::ostream& operator<<(std::ostream& s, const Matrix<T>& m) {

    const size_t xDim = m.xDim();
    const size_t yDim = m.yDim();

    if (xDim == 0 || yDim == 0)
      s << "()";
    else if (yDim == 1) {
      s << "(";
      for (size_t x=0; x < xDim; x++) 
	s << " " << m.direct_access(x); 
      s << " )";
    }
    else {
      s << "( " << std::endl;
      for (size_t y=0; y < yDim; y++) {
	for (size_t x=0; x < xDim; x++) {
	  s << " " << m.direct_access(y*xDim+x);
	}
	s << std::endl;
      }
      s << ")";
    }
      
    return s;
  }


  template<typename T>
  Matrix<T> transpose(const Matrix<T>& m) {
    
    const size_t xDim = m.xDim();
    const size_t yDim = m.yDim();

    Matrix<T> result(yDim,xDim);
    size_t x,y;
    for (x=0; x < xDim; x++) {
      for (y=0; y < yDim; y++) {
	result.direct_access(x*yDim+y) = m.direct_access(y*xDim+x);
	//result(y,x) = m(x,y);
      }
    }

    return result;
  }


  template<typename T>
  Math1D::Vector<T> operator*(const Matrix<T>& m, const Math1D::Vector<T>& v) {

    const size_t xDim = m.xDim();
    const size_t yDim = m.yDim();

    if (xDim != v.size()) {
      INTERNAL_ERROR << "     cannot multiply matrix \"" << m.name() << "\" with vector \""
		     << v.name() << "\":" << std::endl
		     << "     dimensions " << xDim << "x" << yDim << " and " << v.size() 
		     << " mismatch. Exiting..." << std::endl;
      exit(1);
    }

    Math1D::Vector<T> result(yDim);
    size_t y,x;
    T sum;
    for (y=0; y < yDim; y++) {
      sum = (T) 0;
      for (x=0; x < xDim; x++)
	sum += m(x,y) * v[x];
      result[y] = sum;
    }

    return result;
  }

} //end of namespace Math2D

#endif
