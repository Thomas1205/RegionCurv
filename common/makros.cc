/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#include "makros.hh"

template<>
uint convert<uint>(const std::string s) {

  uint result = 0;
  char c;
  uint i=0;
  for (; i < s.size(); i++) {
    c = s[i];

    if (c < '0' || c > '9') {
      std::cerr << "ERROR: conversion of \"" << s << "\" failed. exiting." << std::endl; 
      exit(1);
    }
    result = 10*result + (c - '0');
  }

  return result;
}
