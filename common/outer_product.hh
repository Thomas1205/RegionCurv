/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef OUTER_PRODUCT_HH
#define OUTER_PRODUCT_HH

#include "matrix.hh"
#include "vector.hh"

template<typename T>
Math2D::Matrix<T> outer_product(const Math1D::Vector<T>& v1, const Math1D::Vector<T>& v2);

template<typename T>
void add_outer_product(const Math1D::Vector<T>& v1, const Math1D::Vector<T>& v2,
		       Math2D::Matrix<T>& m);

template<typename T>
void set_outer_product(const Math1D::Vector<T>& v1, const Math1D::Vector<T>& v2,
		       Math2D::Matrix<T>& m);


/******************** implementation **********************/

template<typename T>
Math2D::Matrix<T> outer_product(const Math1D::Vector<T>& v1, const Math1D::Vector<T>& v2) {

  Math2D::Matrix<T> result(v2.size(),v1.size());
  
  for (uint y=0; y < v1.size(); y++)
    for (uint x=0; x < v2.size(); x++) 
      result(x,y) = v1[y]*v2[x];
  
  return result;
}

template<typename T>
void add_outer_product(const Math1D::Vector<T>& v1, const Math1D::Vector<T>& v2,
		       Math2D::Matrix<T>& m) {

  if (m.xDim() != v2.size() || m.yDim() != v1.size()) {

    INTERNAL_ERROR << " cannot add outer product of vectors \"" << v1.name() << "\" and \""
		   << v2.name() << "\" to matrix \"" << m.name() << "\": "
		   << std::endl;
    std::cerr << "    " << v2.size() << "x" << v1.size() << " does not equal " 
	      << m.xDim() << "x" << m.yDim() << ". Exiting..." << std::endl;
    exit(1);
  }

  for (uint y=0; y < v1.size(); y++)
    for (uint x=0; x < v2.size(); x++) 
      m(x,y) += v1[y]*v2[x];
}


template<typename T>
void set_outer_product(const Math1D::Vector<T>& v1, const Math1D::Vector<T>& v2,
		       Math2D::Matrix<T>& m) {

  if (m.xDim() != v2.size() || m.yDim() != v1.size()) {

    INTERNAL_ERROR << " cannot add outer product of vectors \"" << v1.name() << "\" and \""
		   << v2.name() << "\" to matrix \"" << m.name() << "\": "
		   << std::endl;
    std::cerr << "    " << v2.size() << "x" << v1.size() << " does not equal " 
	      << m.xDim() << "x" << m.yDim() << ". Exiting..." << std::endl;
    exit(1);
  }

  for (uint y=0; y < v1.size(); y++)
    for (uint x=0; x < v2.size(); x++) 
      m(x,y) = v1[y]*v2[x];
}


#endif
