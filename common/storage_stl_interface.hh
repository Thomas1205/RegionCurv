/**** written by Thomas Schoenemann as a private person without employment, April 2013 ****/

#ifndef STORAGE_STL_INTERFACE
#define STORAGE_STL_INTERFACE

#include <vector>
#include "storage1D.hh"
#include "storage2D.hh"
#include "storage3D.hh"

template<typename T, typename ST>
void assign(Storage1D<T,ST>& target, const std::vector<T>& source); 

template<typename T, typename ST>
void assign(std::vector<T>& target, const Storage1D<T,ST>& source); 

template<typename T1, typename T2, typename ST>
void assign(Storage1D<T1,ST>& target, const Storage1D<T2,ST>& source);

template<typename T1, typename T2, typename ST>
void assign(Storage2D<T1,ST>& target, const Storage2D<T2,ST>& source);

template<typename T1, typename T2, typename ST>
void assign(Storage3D<T1,ST>& target, const Storage3D<T2,ST>& source);

/************** implementation *************/


template<typename T, typename ST>
void assign(Storage1D<T,ST>& target, const std::vector<T>& source) {
  
  //TODO: think about std::copy()
  
  target.resize_dirty(source.size());
  for (uint k=0; k < source.size(); k++)
    target[k] = source[k];
}

template<typename T, typename ST>
void assign(std::vector<T>& target, const Storage1D<T,ST>& source) {

  //TODO: think about std::copy()

  target.clear();
  target.reserve(source.size());

  for (uint k=0; k < source.size(); k++)
    target.push_back(source[k]);
}

template<typename T1, typename T2, typename ST>
void assign(Storage1D<T1,ST>& target, const Storage1D<T2,ST>& source) {

  target.resize_dirty(source.size());
  for (uint k=0; k < source.size(); k++)
    target[k] = (T1) source[k];  
}

template<typename T1, typename T2, typename ST>
void assign(Storage2D<T1,ST>& target, const Storage2D<T2,ST>& source) {

  target.resize_dirty(source.xDim(),source.yDim());

  for (uint k=0; k < source.size(); k++)
    target.direct_access(k) = (T1) source.direct_access(k);
}

template<typename T1, typename T2, typename ST>
void assign(Storage3D<T1,ST>& target, const Storage3D<T2,ST>& source) {

  target.resize_dirty(source.xDim(),source.yDim(),source.zDim());

  for (uint k=0; k < source.size(); k++)
    target.direct_access(k) = (T1) source.direct_access(k);
}

#endif
