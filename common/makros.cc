/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#include "makros.hh"
#include <map>

template<>
uint convert<uint>(const std::string s) {

  uint result = 0;
  char c;
  uint i=0;
  for (; i < s.size(); i++) {
    c = s[i];

    if (c < '0' || c > '9') {
      std::cerr << "ERROR: conversion of \"" << s << "\" to uint failed. Exiting." << std::endl; 
      exit(1);
    }
    result = 10*result + (c - '0');
  }

  return result;
}

namespace Makros {

  std::map<std::string,std::string> typename_map;

  void register_typename(const std::string& id, const std::string& fullname) {
    typename_map[id] = fullname;
  }
  
  std::string get_typename(const std::string& id) {
    
    std::map<std::string,std::string>::iterator it = typename_map.find(id);
    if (it == typename_map.end())
      return id;
    else
      return it->second;
  }
}
