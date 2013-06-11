/*** first version written by Thomas Schoenemann as a private person without employment, November 2009 ***/
/*** refined at the University of Düsseldorf, Germany, 2012 ***/


#ifndef STRINGPROCESSING_HH
#define STRINGPROCESSING_HH

#include <string>
#include <vector>

bool is_whitespace(char c);

bool is_uppercase(char c);

char downcase(char c);

std::string downcase(std::string s);

bool is_natural_number(const std::string s);

//to avoid the awkward strcmp routine
bool strings_equal(std::string s1, std::string s2);

bool string_ends_with(std::string s, std::string suffix);

bool string_starts_with(std::string s, std::string prefix);

void tokenize(const std::string& s, std::vector<std::string>& tokens, char separator, bool empty_tokens=false);

void tokenize_with_stringsep(const std::string& s, std::vector<std::string>& tokens, std::string sep_string, bool empty_tokens=false);


#endif
