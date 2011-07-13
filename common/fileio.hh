/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef FILEIO_HH
#define FILEIO_HH

#include <cstdlib>
#include <cstdio>
#include "makros.hh"
#include "stringprocessing.hh"

class FileTruncatedException {};

class InvalidCharacterException {
public:

  InvalidCharacterException(char c);
  char c_;
};


void readCharacter(FILE* fptr, char& c) throw (FileTruncatedException);


//ignores leading whitespace
//returns last read character (which was not part of the number)
template<typename T>
char read_natural_number(FILE* fptr, T& number) 
  throw (FileTruncatedException,InvalidCharacterException) {
  
  number = 0;
  char c;
  //1. ignore whitespace until a number is found
  while (true) {
    size_t nRead = fread(&c,1,1,fptr);

    if (nRead != 1)
      throw FileTruncatedException();
    if (!is_whitespace(c)) {
      if (c >= '0' && c <= '9') {
	number = c - '0';
	break;
      }
      else
	throw InvalidCharacterException(c);
    }
  }
  while (true) {
    size_t nRead = fread(&c,1,1,fptr);

    if (nRead != 1)
      throw FileTruncatedException();
    if (c >= '0' && c <= '9')
      number = 10*number + ( c - '0');
    else
      break;
  }  
  return c;
}

char read_ws_until(FILE* fptr, char* allowed_chars, size_t nCharsListed) 
  throw (FileTruncatedException,InvalidCharacterException);

char read_until(FILE* fptr, char* allowed_chars, size_t nCharsListed) 
  throw (FileTruncatedException);






#endif
