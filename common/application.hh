/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef APPLICATION_HH
#define APPLICATION_HH

#include <map>
#include <vector>
#include <set>
#include <fstream>
#include "makros.hh"
#include "stl_util.hh"

enum ParamType {
  flag,
  mandWithValue,
  optWithValue, 
  mandInFilename,
  optInFilename,
  mandOutFilename,
  optOutFilename
};


struct ParamDescr {
  std::string name_;  
  ParamType type_;
  bool default_given_;
  std::string  default_value_;
};

class Application {
public:

  Application(uint argc, char** argv, ParamDescr* param_list, uint nParams);

  ~Application();

  bool is_set(std::string name);
    
  std::string getParam(std::string name);

protected:
  std::map<std::string,std::string> param_value_;
  std::set<std::string> flags_;
    
  std::set<std::string> known_flags_;
  std::set<std::string> known_params_;    

  uint nParams_;
  ParamDescr* param_list_;
  std::string app_name_;
  std::vector<std::string> argv_;
};

#endif
