/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#include "application.hh"

bool Application::is_set(std::string name) {

  if (known_flags_.find(name) != known_flags_.end()) {
  
    return (flags_.find(name) != flags_.end());
  }
  else if (known_params_.find(name) != known_params_.end()) {

    return (param_value_.find(name) != param_value_.end());
  }
  else {
    INTERNAL_ERROR << "   \"" << name << "\" is not a known flag. Exiting..." << std::endl;
    exit(1); 
  }

//   if (known_flags_.find(name) == known_flags_.end()) {
//     INTERNAL_ERROR << "   \"" << name << "\" is not a known flag. Exiting..." << std::endl;
//     exit(1); 
//   } 
    
//   return (flags_.find(name) != flags_.end());
}

std::string Application::getParam(std::string name) {

  if (known_params_.find(name) == known_params_.end()) {
    USER_ERROR << "    \"" << name << "\" is not a known parameter. Exiting..." << std::endl;
    exit(1); 
  }

  if (param_value_.find(name) == param_value_.end()) {
    INTERNAL_ERROR << "    no value for the parameter\"" << name 
		   << "\" specified on command line. Exiting..." << std::endl;
    exit(1); 
  }
    
  return param_value_[name];
}

Application::Application(uint argc, char** argv, ParamDescr* param_list, uint nParams) {

  app_name_ = argv[0];
  nParams_ = nParams;
  param_list_ = new ParamDescr[nParams];
  for (uint p=0; p  < nParams; p++)
    param_list_[p] = param_list[p];

  for (uint i=0; i < argc; i++)
    argv_.push_back(argv[i]);

  /*** phase 0: build database of known options ***/
  for (uint p=0; p  < nParams; p++) {
        
    std::string cur_name = param_list[p].name_;
        
    if ((known_flags_.find(cur_name) != known_flags_.end()) ||
	(known_params_.find(cur_name) != known_params_.end())) {
      INTERNAL_ERROR << "    parameter \"" << cur_name << "\" is listed twice. Exiting."
		     << std::endl;
      exit(1);     
    }       
    
    if (param_list[p].type_ == flag) {
      known_flags_.insert(cur_name);    
    }
    else {
      known_params_.insert(cur_name);
    }
  }

  /*** phase 1: scan command line ***/
  for (uint i=1; i < argc; i++) {
    
    std::string cur_arg = argv[i];
        
    bool found = false;
    for (uint p=0; p  < nParams && !found; p++) {
        
      if (cur_arg == param_list[p].name_) {
            
	found = true;

	if (param_list[p].type_ != flag && i == argc-1) {
	  USER_ERROR << "    no value specified for parameter \"" << cur_arg << "\". Exiting..."
		     << std::endl;
	  exit(1);
	}

	switch (param_list[p].type_) {
	case flag : {
	  flags_.insert(cur_arg);
	  break;
	};
	case mandWithValue: {
	  param_value_[cur_arg] = argv[i+1];
	  break;
	}
	case optWithValue: {
	  param_value_[cur_arg] = argv[i+1];
	  break;
	}
	case mandInFilename: {
	  param_value_[cur_arg] = argv[i+1];
	  break;
	}
	case optInFilename: {
	  param_value_[cur_arg] = argv[i+1];
	  break;
	}
	case mandOutFilename: {
	  param_value_[cur_arg] = argv[i+1];
	  break;
	}
	case optOutFilename: {
	  param_value_[cur_arg] = argv[i+1];
	  break;
	}
	default: {
	  INTERNAL_ERROR << "    invalid type specified for option \"" << cur_arg 
		    << "\". Exiting..." << std::endl;
	  exit(1);
	}
	}
            
	if (param_list[p].type_ != flag)
	  i++;
      }
        
    } 
        
    if (!found) {
      USER_ERROR << " \"" << cur_arg << "\" is not a valid option. Exiting... " << std::endl;
      exit(1);
    }

  }   //end for(i...

  /*** phase 2: add default values and check if mandatory values are given ****/
  for (uint p=0; p  < nParams; p++) {

    ParamType type = param_list[p].type_;
    std::string name = param_list[p].name_;
    bool is_optional = (type == optWithValue) || (type == optInFilename) || (type == optOutFilename);

    if (is_optional && param_list[p].default_given_ && param_value_.find(name) == param_value_.end())
      param_value_[name] = param_list[p].default_value_;

    if (!is_optional && type != flag && param_value_.find(name) == param_value_.end()) {

      USER_ERROR << "  no value specified for mandatory parameter \"" << name << "\". Exiting..."
		 << std::endl;
      exit(1);
    }
  }
    
  /*** phase 3: check filenames for infiles and outfiles (outfiles might e.g. non-existing directories in the path). ***/
  for (uint p=0; p  < nParams; p++) {

    ParamType type = param_list[p].type_;
    std::string name = param_list[p].name_;

    if (type == mandInFilename || (type == optInFilename && param_value_.find(name) != param_value_.end())) {
      std::ifstream in(param_value_[name].c_str());

      if (!in.is_open()) {
	IO_ERROR << "   file \"" << param_value_[name] << "\" does not exist. Exiting..." 
		 << std::endl;
	exit(1);
      }
      in.close();
    }

    if (type == mandOutFilename || (type == optOutFilename && param_value_.find(name) != param_value_.end())) {

      /** a) check if file exists **/
      std::ifstream instream(param_value_[name].c_str());
      if (!instream.is_open()) {

	/** b) check if infofile can be written **/
	std::string infoname = param_value_[name] + ".info";
	std::ofstream of(infoname.c_str());
	if (!of.is_open()) {

	  IO_ERROR << "   outfile \"" << param_value_[name] << "\" cannot be created." 
		   << std::endl;
	  std::cerr << "   Please check if all directories in the specified path exist. Exiting...." 
		    << std::endl;
	  exit(1);
	}
	else
	  of.close();
      }
    }
  }
}


Application::~Application() {

  /*** NOTE: this code is only executed if the application terminates normally. NOT if the exit() function is used ****/

  /*** 1. write .info-files ****/
  for (uint p=0; p  < nParams_; p++) {

    ParamType type = param_list_[p].type_;
    std::string name = param_list_[p].name_;
    
    if (type == mandOutFilename || (type == optOutFilename && param_value_.find(name) != param_value_.end())) {

      std::string infoname = param_value_[name] + ".info";
      std::ofstream of(infoname.c_str());

      of << "**** created by " << app_name_ << std::endl;
      of << "flags: " << std::endl;
      for (std::set<std::string>::iterator it = flags_.begin(); it != flags_.end(); it++)
	of << "  " << (*it) << std::endl;
      
      of << "arguments: " << std::endl;
      for (std::map<std::string,std::string>::iterator it = param_value_.begin(); it != param_value_.end(); it++)
	of << "  " << it->first << " " << it->second << std::endl;
    }
  }  

  /*** 2. "last_call" and <applicationname>.last_call ***/

  for (uint i=0; i < 2; i++) {

    std::string filename;

    if (i == 0)
      filename = "last_call";
    else {
      filename = argv_[0];
      filename += ".last_call";
    }

    std::ofstream of(filename.c_str());

    for (size_t j=0; j < argv_.size(); j++)
      of << argv_[j] << " ";
    of << std::endl;
    of.close();
  }
  

  delete[] param_list_;
}
