/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef MAKROS_HH
#define MAKROS_HH

#include <cassert>
#include <iostream>
#include <limits> //includes numeric_limits
#include <sstream>
#include <string>
#include <iomanip>
#include <cstdlib> //includes the exit-function

//#ifdef WIN32
namespace{
bool isnan(double x) {
  return (x != x);
}
}
#define M_PI 3.1415926535897931
//#endif


/******************** Data Macros *****************************/
typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned char uchar;

#define MIN_DOUBLE -1.0*std::numeric_limits<double>::max()
#define MAX_DOUBLE std::numeric_limits<double>::max()
#define HIGH_DOUBLE (0.1*MAX_DOUBLE)
#define EPS_DOUBLE std::numeric_limits<double>::epsilon()
#define MIN_FLOAT  -1.0f*std::numeric_limits<float>::max()
#define MAX_FLOAT  std::numeric_limits<float>::max()
#define HIGH_FLOAT (0.1f*MAX_FLOAT)
#define EPS_FLOAT  std::numeric_limits<float>::epsilon()
#define MAX_UINT std::numeric_limits<uint>::max()
#define MAX_USHORT std::numeric_limits<ushort>::max()

enum NormType {L1,L2};

/**** helpful routines ****/

template<typename T>
std::string toString(T obj, uint width=1) {
  
  std::ostringstream s;
  
  s << std::setw(width) << std::setfill('0') << obj;
  return s.str();
}

/***********************/
template <typename T>
T convert(const std::string s) {
  
  std::istringstream is(s);
  T result;
  
  is >> result;
  if (is.bad() || is.fail()) {
    std::cerr << "ERROR: conversion of \"" << s << "\" failed. exiting." << std::endl; 
    exit(1);
  }
  if (!is.eof()) {
  
    //check if the string contains additional characters that are not whitespace
    char c;
    while (is >> c) {
      if (c != ' ' && c != '\n' && c != 13 && c != 10) {
        std::cerr << "WARNING AFTER CONVERSION: string contains additional characters" << std::endl;
        break;
      }
    }
  }
  
  return result;
}

template<>
uint convert<uint>(const std::string s); 


/********************* Code Macros ****************************/
#define TODO(s) { std::cerr << "TODO ERROR[" << __FILE__ << ":" << __LINE__ << "]: feature \"" << (s) << "\" is currently not implemented. exiting..." << std::endl; exit(1); } 
#define EXIT(s) { std::cerr << s << std::endl; exit(1); }
#define MAKENAME(s) std::string(#s) + std::string("[") + std::string(__FILE__) + std::string(":") + toString(__LINE__) + std::string("]")

#ifdef SAFE_MODE
#define OPTINLINE
#else
#define OPTINLINE inline
#endif

#define INTERNAL_ERROR std::cerr << "INTERNAL ERROR[" << __FILE__ << ":" << __LINE__ << "]:" << std::endl
#define USER_ERROR std::cerr << "ERROR: "
#define IO_ERROR std::cerr << "I/O ERROR[" << __FILE__ << ":" << __LINE__ << "]:" << std::endl

template<typename T>
inline T sign(T arg) {

  if (arg < ((T) 0.0) )
    return ((T) -1.0);
  else if (arg == ((T) 0.0))
    return ((T) 0.0);
  else
    return ((T) 1.0);
}

template<typename T>
inline T robust_sign(T arg, T tolerance) {

  if (arg < ((T) -tolerance) )
    return ((T) -1.0);
  else if (arg > ((T) tolerance))
    return ((T) 1.0);
  else
    return ((T) 0.0);
}

//NOTE: prefetch will only work on an x86/x86_64 architecture with SSE1 (Pentium 3 or higher)
// if you are compiling on a different architecture simply remove the asm-statements
template<typename T>
inline void prefetcht0(const T* ptr) {
  asm __volatile__ ("prefetcht0 %[ptr]" : : [ptr] "m" (ptr[0]));
}

template<typename T>
inline void prefetcht1(const T* ptr) {
  asm ("prefetcht1 %[ptr]" : : [ptr] "m" (ptr[0]));
}

template<typename T>
inline void prefetcht2(const T* ptr) {
  asm ("prefetcht2 %[ptr]" : : [ptr] "m" (ptr[0]));
}

namespace Makros {

  inline float max(float* data, size_t nData) {
    float max_val=MIN_FLOAT;
    float cur_datum;
    size_t i;

#if !defined(USE_SSE) || USE_SSE < 2 

    for (i=0; i < nData; i++) {
      cur_datum = data[i];
      if (cur_datum > max_val)
	max_val = cur_datum;
    }
#else
    //movups is part of SSE2

    float tmp[4] = {MIN_FLOAT,MIN_FLOAT,MIN_FLOAT,MIN_FLOAT};
    float* fptr;

    asm __volatile__ ("movups %[tmp], %%xmm6" : : [tmp] "m" (tmp[0]) : "xmm6");
    for (i=0; (i+4) <= nData; i += 4) {
      fptr = data+i;
      asm __volatile__ ("movups %[fptr], %%xmm7" : : [fptr] "m" (fptr[0]) : "xmm7");
      asm __volatile__ ("maxps %%xmm7, %%xmm6" : : : "xmm6");
    }
    asm __volatile__ ("movups %%xmm6, %[tmp]" : [tmp] "=m" (tmp[0]) : : );
    for (i=0; i < 4; i++)
      max_val = std::max(max_val,tmp[i]);

    for (i= nData - (nData % 4); i < nData; i++) {
      cur_datum = data[i];
      if (cur_datum > max_val)
	max_val = cur_datum;
    }
#endif    
    
    return max_val;
  }

  inline float min(float* data, size_t nData) {
    float min_val=MAX_FLOAT;
    float cur_datum;
    size_t i;

#if !defined(USE_SSE) || USE_SSE < 2 
    for (i=0; i < nData; i++) {
      cur_datum = data[i];
      if (cur_datum < min_val)
	min_val = cur_datum;
    }
#else
    //movups is part of SSE2

    float tmp[4] = {MAX_FLOAT,MAX_FLOAT,MAX_FLOAT,MAX_FLOAT};
    float* fptr;

    asm __volatile__ ("movups %[tmp], %%xmm6" : : [tmp] "m" (tmp[0]) : "xmm6");
    for (i=0; (i+4) <= nData; i += 4) {
      fptr = data+i;
      asm __volatile__ ("movups %[fptr], %%xmm7" : : [fptr] "m" (fptr[0]) : "xmm7");
      asm __volatile__ ("minps %%xmm7, %%xmm6" : : );
    }
    asm __volatile__ ("movups %%xmm6, %[tmp]" : [tmp] "=m" (tmp[0]) :  : "xmm6");
    for (i=0; i < 4; i++)
      min_val = std::min(min_val,tmp[i]);

    for (i= nData - (nData % 4); i < nData; i++) {
      cur_datum = data[i];
      if (cur_datum < min_val)
	min_val = cur_datum;
    }
#endif    
    
    return min_val;
  }

  
  inline void find_max_and_argmax(float* data, size_t nData, float& max_val, size_t& arg_max) {

    max_val = MIN_FLOAT;
    arg_max = MAX_UINT; 
    size_t i;
    float cur_val;

#if !defined(USE_SSE) || USE_SSE < 4
    for (i=0; i < nData; i++) {
      cur_val = data[i];
      
      if (cur_val > max_val) {
	max_val = cur_val;
	arg_max = i;
      }
    }
#else
    //blendvps is part of SSE4

    assert(nData <= 17179869183);

    float tmp[4] = {MIN_FLOAT,MIN_FLOAT,MIN_FLOAT,MIN_FLOAT};
    float* fptr;
    
    wchar_t itemp[4] = {1,1,1,1}; 

    asm __volatile__ ("movups %[tmp], %%xmm6" : : [tmp] "m" (tmp[0]) : "xmm6" );
    asm __volatile__ ("xorps %%xmm5, %%xmm5" : : : "xmm5"); //sets xmm5 (= argmax) to zero
    asm __volatile__ ("movups %[itemp], %%xmm4" : : [itemp] "m" (itemp[0]) : "xmm4"); //vector of 1s
    asm __volatile__ ("xorps %%xmm3, %%xmm3" : : : "xmm3"); //contains candidate argmax
    for (i=0; (i+4) <= nData; i += 4) {
      fptr = data+i;
      asm __volatile__ ("movups %[fptr], %%xmm7" : : [fptr] "m" (fptr[0]) : "xmm7");
      asm __volatile__ ("movaps %%xmm7, %%xmm0" : : : "xmm0");
      asm __volatile__ ("cmpnleps %%xmm6, %%xmm0" : : : "xmm0"); //mask is stored in xmm0
      asm __volatile__ ("blendvps %%xmm7, %%xmm6" : : : "xmm6"); //xmm0 is implicit argument
      asm __volatile__ ("blendvps %%xmm3, %%xmm5" : : : "xmm5");
      asm __volatile__ ("paddd %%xmm4, %%xmm3" : : : "xmm3");
    }
    asm __volatile__ ("movups %%xmm6, %[tmp]" : [tmp] "=m" (tmp[0]) : : );
    asm __volatile__ ("movups %%xmm5, %[itemp]" : [itemp] "=m" (itemp[0]) : );
    
    for (i=0; i < 4; i++) {
      cur_val = tmp[i];
      if (cur_val > max_val) {
	max_val = cur_val;
	arg_max = 4*itemp[i] + i;
      }
    }

    for (i= nData - (nData % 4); i < nData; i++) {
      cur_val = data[i];
      if (cur_val > max_val) {
	max_val = cur_val;
	arg_max = i;
      }
    }    
#endif  
  }

  inline void find_min_and_argmin(float* data, size_t nData, float& min_val, size_t& arg_min) {

    min_val = MAX_FLOAT;
    arg_min = MAX_UINT; 
    size_t i;
    float cur_val;

#if !defined(USE_SSE) || USE_SSE < 4
    //#if 1
    for (i=0; i < nData; i++) {
      cur_val = data[i];
      
      if (cur_val < min_val) {
	min_val = cur_val;
	arg_min = i;
      }
    }
#else
    //blendvps is part of SSE4
    
    assert(nData <= 17179869183);

    volatile float tmp[4] = {MAX_FLOAT,MAX_FLOAT,MAX_FLOAT,MAX_FLOAT};
    float* fptr;
    
    volatile wchar_t itemp[4] = {1,1,1,1}; 

    asm __volatile__ ("movups %[tmp], %%xmm6" : : [tmp] "m" (tmp[0]) : "xmm6" );
    asm __volatile__ ("xorps %%xmm5, %%xmm5" : : : "xmm5" ); //sets xmm5 (= argmax) to zero
    asm __volatile__ ("movups %[itemp], %%xmm4" : : [itemp] "m" (itemp[0]) : "xmm4"); //vector of 1s
    asm __volatile__ ("xorps %%xmm3, %%xmm3" : : : "xmm3"); //contains candidate argmax
    for (i=0; (i+4) <= nData; i += 4) {
      fptr = data+i;
      asm __volatile__ ("movups %[fptr], %%xmm7" : : [fptr] "m" (fptr[0]) : "xmm7");
      asm __volatile__ ("movaps %%xmm7, %%xmm0" : : : "xmm7", "xmm0");
      asm __volatile__ ("cmpltps %%xmm6, %%xmm0" : : : "xmm0"); //mask is stored in xmm0
      asm __volatile__ ("blendvps %%xmm7, %%xmm6" : : : "xmm6"); //xmm0 is implicit argument
      asm __volatile__ ("blendvps %%xmm3, %%xmm5" : : : "xmm5");
      asm __volatile__ ("paddd %xmm4, %xmm3");
    }
    asm __volatile__ ("movups %%xmm6, %[tmp]" : [tmp] "=m" (tmp[0]) : : "xmm6" );
    asm __volatile__ ("movups %%xmm5, %[itemp]" : [itemp] "=m" (itemp[0]) : : "xmm5");

    //std::cerr << "intermediate minval: " << min_val << std::endl;

    for (i=0; i < 4; i++) {
      cur_val = tmp[i];
      //std::cerr << "cur val: " << cur_val << std::endl;
      if (cur_val < min_val) {
	min_val = cur_val;
	arg_min = 4*itemp[i] + i;
      }
    }

    //std::cerr << "minval: " << min_val << std::endl;

    for (i= nData - (nData % 4); i < nData; i++) {
      cur_val = data[i];
      if (cur_val < min_val) {
	min_val = cur_val;
	arg_min = i;
      }
    }
#endif  
  }

  inline void mul_array(float* data, size_t nData, const float constant) {

    size_t i;
#if !defined(USE_SSE) || USE_SSE < 2
    for (i=0; i < nData; i++) {
      data[i] *= constant;
    }
#else
    float temp[4];
    float* fptr;
    for (i=0; i < 4; i++)
      temp[i] = constant;
    asm volatile ("movups %[temp], %%xmm7" : : [temp] "m" (temp[0]) : "xmm7" );

    for (i=0; i+4 <= nData; i+=4) {
      fptr = data + i;
      asm volatile ("movups %[fptr], %%xmm6" : : [fptr] "m" (fptr[0]) : "xmm6" );
      asm volatile ("mulps %%xmm7, %%xmm6" : : : "xmm6");
      asm volatile ("movups %%xmm6, %[fptr]" : [fptr] "=m" (fptr[0]) : : );
    }
    
    for (i= nData - (nData % 4); i < nData; i++) {
      data[i] *= constant;
    }
#endif
  }
  
  inline void mul_array(double* data, size_t nData, const double constant) {

    size_t i;
#if !defined(USE_SSE) || USE_SSE < 2
    for (i=0; i < nData; i++) {
      data[i] *= constant;
    }
#else
    double temp[2];
    double* dptr;
    for (i=0; i < 2; i++)
      temp[i] = constant;
    asm volatile ("movupd %[temp], %%xmm7" : : [temp] "m" (temp[0]) : "xmm7" );

    for (i=0; i+2 <= nData; i+=2) {
      dptr = data + i;
      asm volatile ("movupd %[dptr], %%xmm6" : : [dptr] "m" (dptr[0]) : "xmm6" );
      asm volatile ("mulpd %%xmm7, %%xmm6" : : : "xmm6");
      asm volatile ("movupd %%xmm6, %[dptr]" : [dptr] "=m" (dptr[0]) : : );
    }
    
    for (i= nData - (nData % 2); i < nData; i++) {
      data[i] *= constant;
    }
#endif
  }

  //performs data[i] -= factor*data2[i] for each i
  //this is a frequent operation in the conjugate gradient algorithm
  inline void array_subtract_multiple(double* data, size_t nData, double factor, 
				      const double* data2) {

    size_t i;
#if !defined(USE_SSE) || USE_SSE < 2
    for (i=0; i < nData; i++)
      data[i] -= factor*data2[i];
#else
    double temp[2];
    double* dptr;
    const double* cdptr;
    for (i=0; i < 2; i++)
      temp[i] = factor;
    asm volatile ("movupd %[temp], %%xmm7" : : [temp] "m" (temp[0]) : "xmm7" );
    for (i=0; i+2 <= nData; i+=2) {
      cdptr = data2+i;
      asm volatile ("movupd %[cdptr], %%xmm6" : : [cdptr] "m" (cdptr[0]) : "xmm6" );
      asm volatile ("mulpd %%xmm7, %%xmm6" : : : "xmm6");
      dptr = data+i;
      asm volatile ("movupd %[dptr], %%xmm5" : : [dptr] "m" (dptr[0]) : "xmm5" );
      asm volatile ("subpd %%xmm6, %%xmm5" : : : "xmm5");
      asm volatile ("movupd %%xmm5, %[dptr]" : [dptr] "=m" (dptr[0]) : : );
    }
    
    for (i= nData - (nData % 2); i < nData; i++) 
      data[i] -= factor*data2[i];
#endif
  }
  


} //end of namespace Makros


/*** use std::swap instead!! ***/
// template<typename T>
// inline void swap(T& x, T& y) {
//   T z = x;
//   x = y;
//   y = z;
// }

#endif
