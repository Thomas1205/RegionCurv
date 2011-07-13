/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef STL_OUT_HH
#define STL_OUT_HH

#include <map>
#include <vector>
#include <set>
#include <iostream>

/*************** declarations ******************/

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec);

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& s);

template<typename TK, typename TE>
std::ostream& operator<<(std::ostream& os, const std::map<TK,TE>& m);

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1,T2>& p);


/*************** implementation ****************/

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {

    os << "[ ";
    for (size_t i = 0; i < vec.size(); i++) {
        os << vec[i];
        if (i+1 != vec.size())
            os << ", ";
    }   
    os << " ]";

    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& s) {

    os << "{ ";
    for (typename std::set<T>::const_iterator it=s.begin(); it != s.end(); ) {
        os << (*it);
        it++;
        if (it != s.end())
            os << ", ";
    }   
    os << " }";

    return os;
}

template<typename TK, typename TE>
std::ostream& operator<<(std::ostream& os, const std::map<TK,TE>& m) {

    os << "[ ";
    for (typename std::map<TK,TE>::const_iterator it=m.begin(); it != m.end(); ) {
        os << it->first << "->" << it->second;
        it++;
        if (it != m.end())
            os << ", ";
    }       
    os << " ]";

    return os;
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1,T2>& p) {

  os << "(" << p.first << "," << p.second << ")";

  return os;
}





#endif
