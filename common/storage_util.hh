/**** written by Thomas Schoenemann as a private person without employment, March 2013 ****/

#ifndef STORAGE_UTIL_HH
#define STORAGE_UTIL_HH


#include "storage1D.hh"
#include <algorithm>

template<typename T>
inline void sort_storage1D(Storage1D<T>& stor) {

  std::sort(stor.direct_access(),stor.direct_access()+stor.size());
}

template<typename T, typename ST>
inline ST find_in_storage1D(const Storage1D<T,ST>& stor, T element) {

  return std::find(stor.direct_access(),stor.direct_access()+stor.size(),element) - stor.direct_access();
}

template<typename T, typename ST>
inline bool contains(const Storage1D<T,ST>& stor, T element) {

  T* end = stor.direct_access()+stor.size();

  return (std::find(stor.direct_access(),end,element) != end);
}
#endif
