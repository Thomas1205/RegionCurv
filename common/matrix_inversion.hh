/**** written by Thomas Schoenemann as a private person without employment, September 2009 ***/
/****  these routines are not members of the class Matrix<T> as inverses do not exist for e.g. integral types  ****/

#ifndef MATRIX_INVERSION_HH
#define MATRIX_INVERSION_HH

#include "matrix.hh"

class NotInvertibleException {};

template<typename T>
void invert_matrix(const Math2D::Matrix<T>& toInvert, Math2D::Matrix<T>& result) throw (NotInvertibleException);

//the matrix toInvert is modified in case the dimension of the matrix is greater than 2.
//toInvert then contains the identity matrix (up to roundoff errors) after the call.  
template<typename T>
void invert_and_destroy_matrix(Math2D::Matrix<T>& toInvert, Math2D::Matrix<T>& result) 
  throw (NotInvertibleException);


/********* implementation **********/
/** these routines are mostly intended for all those who need a quick solution, i.e. are too lazy to link to 
** one of the standard packages. Most likely there are faster and more robust solutions out there. **/

template<typename T>
void invert_matrix(const Math2D::Matrix<T>& toInvert, Math2D::Matrix<T>& result) throw (NotInvertibleException) {

  size_t dim = toInvert.xDim();
  if (toInvert.yDim() != dim) {

    INTERNAL_ERROR << " cannot invert non-square matrix \"" << toInvert.name() << "\" of size " 
		   << dim << "x" << toInvert.yDim() << ". Exiting..." << std::endl;  
  }

  if (result.xDim() != dim || result.yDim() != dim)
    result.resize_dirty(dim,dim);

  if (dim == 1) {
    result.direct_access(0) = 1.0 / toInvert.direct_access(0);
  }
  else if (dim == 2) {
    double det = toInvert.direct_access(0)*toInvert.direct_access(3) 
      - toInvert.direct_access(1)*toInvert.direct_access(2);
    
    if (det == 0) {
      throw NotInvertibleException();
    }

    result.direct_access(0) = toInvert.direct_access(3);
    result.direct_access(1) = -toInvert.direct_access(1);
    result.direct_access(2) = -toInvert.direct_access(2);
    result.direct_access(3) = toInvert.direct_access(0);

    result *= (1.0 / det);
  }
  else {

    //    std::cerr << "----------------------------------------" << std::endl;

    Math2D::Matrix<T> copy = toInvert;
    invert_and_destroy_matrix(copy,result);
  }

}


template<typename T>
void invert_and_destroy_matrix(Math2D::Matrix<T>& toInvert, Math2D::Matrix<T>& result) 
  throw (NotInvertibleException) {

  size_t dim = toInvert.xDim();
  if (toInvert.yDim() != dim) {

    INTERNAL_ERROR << " cannot invert non-square matrix \"" << toInvert.name() << "\" of size " 
		   << dim << "x" << toInvert.yDim() << ". Exiting..." << std::endl;  
  }


  if (result.xDim() != dim || result.yDim() != dim)
    result.resize_dirty(dim,dim);

  if (dim == 1) {
    result.direct_access(0) = 1.0 / toInvert.direct_access(0);
  }
  else if (dim == 2) {
    double det = toInvert.direct_access(0)*toInvert.direct_access(3) 
      - toInvert.direct_access(1)*toInvert.direct_access(2);
    
    if (det == 0) {
      throw NotInvertibleException();
    }

    result.direct_access(0) = toInvert.direct_access(3);
    result.direct_access(1) = -toInvert.direct_access(1);
    result.direct_access(2) = -toInvert.direct_access(2);
    result.direct_access(3) = toInvert.direct_access(0);

    result *= (1.0 / det);
  }
  else {

    //    std::cerr << "----------------------------------------" << std::endl;
    
    for (uint i=0; i < dim*dim; i++)
      result.direct_access(i) = 0.0;
    for (uint i=0; i < dim*dim; i+= dim+1)
      result.direct_access(i) = 1.0;

    for (uint y=0; y < dim; y++) {

//       std::cerr << "--y: " << y << std::endl;
//       std::cerr << "toInvert: " << toInvert << std::endl;
//       std::cerr << "result: " << result << std::endl;
      
      /** find next pivot row **/
      T max_val = toInvert(y,y);
      uint arg_max = y;

      for (uint yy=y+1; yy < dim; yy++) {
	T cur = toInvert(y,yy);
	if (max_val < cur) {
	  max_val = cur;
	  arg_max = yy;
	} 
      }

      //std::cerr << "max: " << max_val << ", arg_max: " << arg_max << std::endl;

      /** swap rows if necessary **/
      if (arg_max != y) {
	
// 	for (uint x=y; x <  dim; x++) {
// 	  swap(toInvert(x,y),toInvert(x,arg_max));
// 	}
	std::swap_ranges(toInvert.direct_access()+y*dim+y, toInvert.direct_access()+(y+1)*dim,
			 toInvert.direct_access()+arg_max*dim+y);
	
// 	for (uint x=0; x < dim; x++) {
// 	  swap(result(x,y),result(x,arg_max));
// 	}
	std::swap_ranges(result.direct_access()+y*dim,result.direct_access()+(y+1)*dim,
			 result.direct_access()+arg_max*dim);
	
      }
      
//       std::cerr << "toInvert after swap: " << toInvert << std::endl;
//       std::cerr << "result after swap: " << result << std::endl;

      if (std::abs(max_val) <= 1e-10)
	throw NotInvertibleException();

      /** scale row y **/
      if (max_val != 1.0) {
	T inv_max = 1.0 / max_val;
	for (uint x=y; x < dim; x++) {
	  toInvert(x,y) *= inv_max;
	}
	for (uint x=0; x < dim; x++)
	  result(x,y) *= inv_max;
      }

//       std::cerr << "toInvert after scale: " << toInvert << std::endl;
//       std::cerr << "result after scale: " << result << std::endl;

      /** perform gaussian elimination **/
      for (uint yy=0; yy < dim; yy++) {
	if (yy != y) {
	  T scale = -toInvert(y,yy);
	  toInvert(y,yy) = 0.0;
	  for (uint x=y+1; x < dim; x++)
	    toInvert(x,yy) += scale * toInvert(x,y);
	  
	  for (uint x=0; x < dim; x++)
	    result(x,yy) += scale * result(x,y);
	}
      }
    }

    //std::cerr << "----------------------------------------" << std::endl;
  }

}




#endif
