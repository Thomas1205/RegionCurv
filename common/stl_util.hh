/*** written by Thomas Schoenemann as a private person without employment, March 2013 ****/

#ifndef STL_UTIL_HH
#define STL_UTIL_HH

#include "makros.hh"
#include <vector>
#include <set>
#include <map>
#include <algorithm>


template <typename T>
T vec_sum(const std::vector<T>& vec);

template <typename T>
T set_sum(const std::set<T>& s);

template <typename T>
inline typename std::vector<T>::iterator vec_find(std::vector<T>& vec, T element);

template <typename T>
inline typename std::vector<T>::const_iterator vec_find(const std::vector<T>& vec, T element);

template <typename T>
inline bool contains(const std::set<T> s, T element);


template <typename T>
inline void vec_sort(std::vector<T>& vec);


namespace Makros {

  template<typename T>
  class Typename<std::vector<T> > {
  public:

    std::string name() const {

      return "std::vector<" + Makros::Typename<T>() + "> ";
    }
  };

  template<typename T>
  class Typename<std::set<T> > {
  public:

    std::string name() const {

      return "std::set<" + Makros::Typename<T>() + "> ";
    }
  };

  template<typename T1, typename T2>
  class Typename<std::map<T1,T2> > {
  public:

    std::string name() const {

      return "std::map<" + Makros::Typename<T1>() + "," + Makros::Typename<T1>() + "> ";
    }
  };

  template<typename T1, typename T2>
  class Typename<std::pair<T1,T2> > {
  public:

    std::string name() const {

      return "std::pair<" + Makros::Typename<T1>() + "," + Makros::Typename<T1>() + "> ";
    }
  };


}


/*********** implementation *********/

template <typename T>
T vec_sum(const std::vector<T>& vec) {

  T sum = T();

  for (typename std::vector<T>::const_iterator it = vec.begin(); it != vec.end(); it++)
    sum += *it;

  return sum;
}

template <typename T>
T set_sum(const std::set<T>& s) {

  T sum = T();

  for (typename std::set<T>::const_iterator it = s.begin(); it != s.end(); it++)
    sum += *it;

  return sum;
}

template <typename T>
inline typename std::vector<T>::const_iterator vec_find(const std::vector<T>& vec, T element) {
  
  return std::find(vec.begin(),vec.end(),element);
}

template <typename T>
inline typename std::vector<T>::iterator vec_find(std::vector<T>& vec, T element) {
  
  return std::find(vec.begin(),vec.end(),element);
}

template <typename T>
inline bool contains(const std::set<T> s, T element) {

  return s.find(element) != s.end();
}

template <typename T>
inline void vec_sort(std::vector<T>& vec) {

  std::sort(vec.begin(),vec.end());
}

#endif
