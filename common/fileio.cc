/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#include "fileio.hh"

InvalidCharacterException::InvalidCharacterException(char c) : c_(c) {}


void readCharacter(FILE* fptr, char& c) throw (FileTruncatedException) { 

  size_t nRead = fread(&c,1,1,fptr);

  if (nRead != 1)
    throw FileTruncatedException();
}


char read_ws_until(FILE* fptr, char* allowed_chars, size_t nCharsListed) 
  throw (FileTruncatedException,InvalidCharacterException) {

  char c;
  
  while (true) {
    size_t nRead = fread(&c,1,1,fptr);
    if (nRead != 1)
      throw FileTruncatedException();

    bool found = false;

    for (size_t i=0; i < nCharsListed; i++) {
      if (c == allowed_chars[i]) {
	found = true;
	break;
      }
    }

    if (found)
      break;
    else if (!is_whitespace(c))
      throw InvalidCharacterException(c);
  }
  
  return c;
}

char read_until(FILE* fptr, char* allowed_chars, size_t nCharsListed) 
  throw (FileTruncatedException) {

  char c;
  
  while (true) {
    size_t nRead = fread(&c,1,1,fptr);
    if (nRead != 1)
      throw FileTruncatedException();

    bool found = false;

    for (size_t i=0; i < nCharsListed; i++) {
      if (c == allowed_chars[i]) {
	found = true;
	break;
      }
    }

    if (found)
      break;
  }
  
  return c;
}



