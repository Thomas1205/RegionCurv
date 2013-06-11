/*** first version written by Thomas Schoenemann as a private person without employment, November 2009 ***/
/*** refined at the University of Düsseldorf, Germany, 2012 ***/

#include "stringprocessing.hh"
#include "makros.hh"

bool is_whitespace(char c) {

  if (c == ' ' || c == '\n' || c == '\t' || c == 10 || c == '\r')
    return true;
  return false;
}

bool is_uppercase(char c) {

  if (c >= 'A' && c <= 'Z')
    return true;
//   if (c == 'Ä' || c == 'Ö' || c == 'Ü')
//      return true;
  
  return false;
}

char downcase(char c) {

  if (c >= 'A' && c <= 'Z')
    return c + ('a' - 'A');

  //outcommented since gcc claims these are several characters
//   if (c == 'Ä')
//     return 'ä';
//   if (c == 'Ö')
//     return 'ö';
//   if (c == "Ü")
//     return "ü";
//   if (c == 'É')
//     return 'é'

  return c;
}

static std::vector<std::pair<std::string,std::string> > replacements;

std::string downcase(std::string s) {

  if (replacements.size() == 0) {
    replacements.push_back(std::make_pair("Ä","ä"));
    replacements.push_back(std::make_pair("Å","å"));
    replacements.push_back(std::make_pair("À","à"));
    replacements.push_back(std::make_pair("Á","á"));
    replacements.push_back(std::make_pair("Ã","ã"));
    replacements.push_back(std::make_pair("Æ","æ"));
    replacements.push_back(std::make_pair("Ā","ā"));
    replacements.push_back(std::make_pair("Ă","ă"));
    replacements.push_back(std::make_pair("Ą","ą"));
    replacements.push_back(std::make_pair("É","é"));
    replacements.push_back(std::make_pair("È","è"));
    replacements.push_back(std::make_pair("Ê","ê"));
    replacements.push_back(std::make_pair("Ë","ë"));
    replacements.push_back(std::make_pair("Ē","ē"));
    replacements.push_back(std::make_pair("Ĕ","ĕ"));
    replacements.push_back(std::make_pair("Ė","ė"));
    replacements.push_back(std::make_pair("Ę","ę"));
    replacements.push_back(std::make_pair("Ě","ě"));
    replacements.push_back(std::make_pair("Î","î"));
    replacements.push_back(std::make_pair("Í","í"));
    replacements.push_back(std::make_pair("Ì","ì"));
    replacements.push_back(std::make_pair("Ï","ï"));
    replacements.push_back(std::make_pair("Ĩ","ĩ"));
    replacements.push_back(std::make_pair("Ĭ","ĭ"));
    replacements.push_back(std::make_pair("Į","į"));
    replacements.push_back(std::make_pair("Ö","ö"));
    replacements.push_back(std::make_pair("Ø","ø"));
    replacements.push_back(std::make_pair("Œ","œ"));
    replacements.push_back(std::make_pair("Ó","ó"));
    replacements.push_back(std::make_pair("Ô","ô"));
    replacements.push_back(std::make_pair("Ò","ò"));
    replacements.push_back(std::make_pair("Õ","õ"));
    replacements.push_back(std::make_pair("Ō","ō"));
    replacements.push_back(std::make_pair("Ŏ","ŏ"));
    replacements.push_back(std::make_pair("Ő","ő"));
    replacements.push_back(std::make_pair("Ü","ü"));
    replacements.push_back(std::make_pair("ú","Ú"));
    replacements.push_back(std::make_pair("Ù","ù"));
    replacements.push_back(std::make_pair("Û","û"));
    replacements.push_back(std::make_pair("Ũ","ũ"));
    replacements.push_back(std::make_pair("Ū","ū"));
    replacements.push_back(std::make_pair("Ŭ","ŭ"));
    replacements.push_back(std::make_pair("Ů","ů"));
    replacements.push_back(std::make_pair("Ű","ű"));
    replacements.push_back(std::make_pair("Ų","ų"));
    replacements.push_back(std::make_pair("Ç","ç"));
    replacements.push_back(std::make_pair("Đ","đ"));
    replacements.push_back(std::make_pair("Č","č"));
    replacements.push_back(std::make_pair("Ċ","ċ"));
    replacements.push_back(std::make_pair("Ĉ","ĉ"));
    replacements.push_back(std::make_pair("Ğ","ğ"));
    replacements.push_back(std::make_pair("Ġ","ġ"));
    replacements.push_back(std::make_pair("Ģ","ģ"));
    replacements.push_back(std::make_pair("Ĥ","ĥ"));
    replacements.push_back(std::make_pair("Ħ","ħ"));
    replacements.push_back(std::make_pair("Ĵ","ĵ"));
    replacements.push_back(std::make_pair("Ķ","ķ"));
    replacements.push_back(std::make_pair("Ł","ł"));
    replacements.push_back(std::make_pair("Ĺ","ĺ"));
    replacements.push_back(std::make_pair("Ļ","ļ"));
    replacements.push_back(std::make_pair("Ľ","ľ"));
    replacements.push_back(std::make_pair("Ŀ","ŀ"));
    replacements.push_back(std::make_pair("Ñ","ñ"));
    replacements.push_back(std::make_pair("Ń","ń"));
    replacements.push_back(std::make_pair("Ņ","ņ"));
    replacements.push_back(std::make_pair("Ň","ň"));
    replacements.push_back(std::make_pair("Ŋ","ŋ"));
    replacements.push_back(std::make_pair("Ŕ","ŕ"));
    replacements.push_back(std::make_pair("Ŗ","ŗ"));
    replacements.push_back(std::make_pair("Ř","ř"));
    replacements.push_back(std::make_pair("Š","š"));
    replacements.push_back(std::make_pair("Ś","ś"));
    replacements.push_back(std::make_pair("Ŝ","ŝ"));
    replacements.push_back(std::make_pair("Ş","ş"));
    replacements.push_back(std::make_pair("Š","š"));
    replacements.push_back(std::make_pair("Ţ","ţ"));
    replacements.push_back(std::make_pair("Ť","ť"));
    replacements.push_back(std::make_pair("Ŧ","ŧ"));
    replacements.push_back(std::make_pair("Ŵ","ŵ"));
    replacements.push_back(std::make_pair("Ÿ","ÿ"));
    replacements.push_back(std::make_pair("Ý","ý"));
    replacements.push_back(std::make_pair("Ź","ź"));
    replacements.push_back(std::make_pair("Ž","ž"));
    replacements.push_back(std::make_pair("Ż","ż"));
    replacements.push_back(std::make_pair("Þ","þ"));
  }
  std::string ls=s;
  for (uint k=0; k < ls.size(); k++) {
    ls[k] = downcase(ls[k]);

    for (uint j=0; j < replacements.size(); j++) {
      if (ls.substr(k,replacements[j].first.size()) == replacements[j].first)
        ls.replace(k,replacements[j].first.size(),replacements[j].second);
    }
  }

  return ls;
}

bool is_natural_number(const std::string s) {

  for (uint i=0; i < s.size(); i++) {
    if (s[i] < '0' || s[i] > '9')
      return false;
  }
  return true;
}


//to avoid the awkward strcmp routine
bool strings_equal(std::string s1, std::string s2) {
    return (s1 == s2);
}

void tokenize(const std::string& s, std::vector<std::string>& tokens, char separator, bool empty_tokens) {

  bool last_sep = true;
  tokens.clear();

  for (uint i=0; i < s.size(); i++) {

    if (s[i] == separator) {
      
      if (empty_tokens)
	tokens.push_back(std::string());

      last_sep = true;
    }
    else {
      if (last_sep)
	tokens.push_back(std::string());
      tokens.back() += s[i];
      last_sep = false;
    }
  }
}

void tokenize_with_stringsep(const std::string& org_s, std::vector<std::string>& tokens, std::string sep_string, bool empty_tokens) {

  tokens.clear();

  std::string s = org_s;

  while (true) {

    size_t pos = s.find(sep_string);

    if (pos >= s.size()) {
      if (s != "")
	tokens.push_back(s);
      break;
    }
    else {
      if (empty_tokens || pos != 0)
	tokens.push_back(s.substr(0,pos));

      s = s.substr(pos+sep_string.size());
    }
  }
}

bool string_ends_with(std::string s, std::string suffix) {

  if (s.size() < suffix.size())
    return false;
  else
    return s.substr(s.size() - suffix.size(), suffix.size()) == suffix;
}

// kaeshammer
bool string_starts_with(std::string s, std::string prefix) {
  return s.compare(0, prefix.length(), prefix) == 0;
}


