/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef GRAYIMAGE_HH
#define GRAYIMAGE_HH

#include <cstdio>
#include <fstream>

#include "matrix.hh"
#include "fileio.hh"

namespace Math2D {

  template <typename T>
  class GrayImage : public Math2D::Matrix<T> {
  public:
    
    /**** constructors *****/
    
    GrayImage();
    
    //copy constructor
    GrayImage(const GrayImage<T>& toCopy);
    
    GrayImage(size_t xDim, size_t yDim);
        
    GrayImage(size_t xDim, size_t yDim, T default_value);

    GrayImage(size_t xDim, size_t yDim, T default_value, T max_intensity);
        
    GrayImage(const std::string& filename);
        
    void invert_image();

    virtual const std::string& name() const;
        
    //returns if the operation was successful
    bool loadPGM(const std::string& filename);

    //returns if the operation was successful
    bool savePGM(const std::string& filename, bool fit_to_range = true) const;
    
    T max_intensity() const; 

    inline void operator=(const GrayImage<T>& toCopy);
    
  protected:
    T max_intensity_;
    static const std::string grayim_name_;
  };


  template<typename T>
  class NamedGrayImage : public GrayImage<T> {
  public:

    NamedGrayImage();

    NamedGrayImage(const std::string& name);
    
    NamedGrayImage(size_t xDim, size_t yDim, const std::string& name);
        
    NamedGrayImage(size_t xDim, size_t yDim, T default_value, const std::string& name);

    NamedGrayImage(size_t xDim, size_t yDim, T default_value, T max_intensity, const std::string& name);
        
    NamedGrayImage(const std::string& filename, const std::string& name);

    virtual const std::string& name() const;

    //NOTE: the name is NOT copied
    inline void operator=(const NamedGrayImage<T>& toCopy);

    inline void operator=(const GrayImage<T>& toCopy);
    
  protected:
    std::string name_;
  };


/****************************** implementation **********************************/

  /*********** implementation of GrayImage **************/
  template<typename T>
  /*static*/ const std::string GrayImage<T>::grayim_name_ = "unnamed gray image";

  template<typename T>
  GrayImage<T>::GrayImage() : Matrix<T>(), max_intensity_((T) 255) {}


  //copy constructor
  template<typename T>
  GrayImage<T>::GrayImage(const GrayImage<T>& toCopy) : 
    Matrix<T>(static_cast<const Matrix<T>&>(toCopy)), max_intensity_(toCopy.max_intensity()) {
  }

  template<typename T>
  GrayImage<T>::GrayImage(size_t xDim, size_t yDim) : Matrix<T>(xDim,yDim), max_intensity_((T) 255) {}
    
  template<typename T>
  GrayImage<T>::GrayImage(size_t xDim, size_t yDim, T default_value) : 
    Matrix<T>(xDim,yDim,default_value), max_intensity_((T) 255) {}
    
  template<typename T>
  GrayImage<T>::GrayImage(size_t xDim, size_t yDim, T default_value, T max_intensity) : 
    Matrix<T>(xDim,yDim,default_value), max_intensity_(max_intensity) {}

  template<typename T>
  GrayImage<T>::GrayImage(const std::string& filename) : Matrix<T>(), max_intensity_((T) 255)  {
    loadPGM(filename);
  }
    
  template<typename T>
  const std::string& GrayImage<T>::name() const {
    return grayim_name_;
  }

  template<typename T>
  void GrayImage<T>::invert_image() {
    for (size_t i=0; i < Matrix<T>::size_; i++) {
      Matrix<T>::data_[i] = max_intensity - Matrix<T>::data_[i];
    }
  }
    
  template<typename T>
  bool GrayImage<T>::loadPGM(const std::string& filename) {
    
    FILE* fptr = fopen(filename.c_str(),"rb");

    if (fptr == 0) {
      IO_ERROR << "   unable to open file \"" << filename << "\"."
	       << std::endl;
      return false;
    }
    
    char firststop[3] = {'P','p','#'};
    char newline[1] = {'\n'};

    uchar c;
    short wc;
    size_t number;
    size_t nRead;

    try {

      //1. read preamble
      while (true) {
	c = read_ws_until(fptr,firststop,3);
      
	if (c == '#') {
	  //comment line, read until line end
	  read_until(fptr,newline,1);
	}
	else
	  break;
      }

      c = read_natural_number(fptr,number);
      if (number==6) {
	TODO("handle ppms in loadPGM. Please look at colorimage.hh");
      }
      else if (number!=5) {
	IO_ERROR << "   incompatible file header." << std::endl;
	return false;
      }

      if (!is_whitespace(c)) 
	throw InvalidCharacterException(c);
      
      //2. read x-dimension
      bool read_xdim = true;
      do {

	read_xdim = true;
	
	try {
	  c = read_natural_number(fptr,Matrix<T>::xDim_);
	  if (!is_whitespace(c)) 
	    throw InvalidCharacterException(c);
	}
	catch (InvalidCharacterException ice) {
	  read_xdim = false;
	  if (ice.c_ == '#') {
	    read_until(fptr,newline,1);
	  }
	  else
	    throw;
	}
      } while (!read_xdim);

      //std::cerr << "xDim = " << Matrix<T>::xDim_ << std::endl;

      //3. read y-dimension
      c = read_natural_number(fptr,Matrix<T>::yDim_);
      if (!is_whitespace(c)) 
	throw InvalidCharacterException(c);

      //std::cerr << "yDim = " << Matrix<T>::yDim_ << std::endl;

      //4. read maximal intensity
      c = read_natural_number(fptr,max_intensity_);
      if (!is_whitespace(c)) 
	throw InvalidCharacterException(c);

      //std::cerr << "max-intensity = " << max_intensity_ << std::endl;

      //5. read until line end, then read image data
      if (c != '\n')
	c = read_ws_until(fptr,newline,1);
      if (c != '\n')
	throw InvalidCharacterException(c);
      
      if (Matrix<T>::data_ != 0)
	delete[] Matrix<T>::data_;
      
      Matrix<T>::size_ = Matrix<T>::xDim_*Matrix<T>::yDim_;
      Matrix<T>::data_ = new T[Matrix<T>::size_];

      //NOTE: currently we assume ther is no line-feed character after the newline

      for (size_t i=0; i < Matrix<T>::size_; i++) {
	
	if (max_intensity_ < 256) {
	  nRead = fread(&c,1,1,fptr);
	  Matrix<T>::data_[i] = c;
	}
	else {
	  nRead = fread(&wc,2,1,fptr);
	  Matrix<T>::data_[i] = wc;
	}
	if (nRead != 1)
	  throw FileTruncatedException();
      }

      //test if there are additional bytes in the file
      nRead = fread(&c,1,1,fptr);
      if (nRead != 0)
	std::cerr << "WARNING: additional bytes in image file. Perhaps a line-feed was handled incorrectly?"
		  << std::endl;
    }
    catch (FileTruncatedException) {
      IO_ERROR << "  file appears to be truncated." << std::endl;
      Matrix<T>::xDim_ = 0;
      Matrix<T>::yDim_ = 0;
      Matrix<T>::size_ = Matrix<T>::xDim_*Matrix<T>::yDim_;
      return false;
    }
    catch (InvalidCharacterException ice) {
      IO_ERROR << "    read invalid character \"" << ice.c_ 
		<< "\"." << std::endl;
      Matrix<T>::xDim_ = 0;
      Matrix<T>::yDim_ = 0;
      Matrix<T>::size_ = Matrix<T>::xDim_*Matrix<T>::yDim_;
      return false;
    }

    return true;

    fclose(fptr);
  }

  template<typename T>
  bool GrayImage<T>::savePGM(const std::string& filename, bool fit_to_range) const {

    return Matrix<T>::savePGM(filename.c_str(), (size_t) max_intensity_, fit_to_range);

//     std::ofstream of(filename.c_str());

//     if (!of.is_open()) {
//       IO_ERROR << " while saving PGM: could not write file \"" << filename 
// 	       << "\". Please check if the path is correct." << std::endl;
//       return false;
//     }

//     of << "P5\n" << Matrix<T>::xDim_ << " " << Matrix<T>::yDim_ << "\n" << ((size_t) max_intensity_) << "\n";
//     for (size_t i=0; i < Matrix<T>::size_; i++) {

//       if (max_intensity_ < 256) {
// 	uchar c = Matrix<T>::data_[i];
// 	of << c;
//       }
//       else {
// 	TODO("handle sizes > 255 when saving PGMs");
//       }
//     }

//     return true;
  }

    
  template<typename T>
  T GrayImage<T>::max_intensity() const {
    return max_intensity_;
  }

  template<typename T>
  inline void GrayImage<T>::operator=(const GrayImage<T>& toCopy) {
    Matrix<T>::operator=(static_cast<const Matrix<T>&>(toCopy));
    max_intensity_ = toCopy.max_intensity();
  }

  
  /***** implementation of NamedGrayImage ******/

  template<typename T>
  NamedGrayImage<T>::NamedGrayImage() : GrayImage<T>(), name_("yyy") {}

  template<typename T>
  NamedGrayImage<T>::NamedGrayImage(const std::string& name) : GrayImage<T>(), name_(name) {}
  
  template<typename T>
  NamedGrayImage<T>::NamedGrayImage(size_t xDim, size_t yDim, const std::string& name) :
    GrayImage<T>(xDim,yDim), name_(name) {
  }
  
  template<typename T>
  NamedGrayImage<T>::NamedGrayImage(size_t xDim, size_t yDim, T default_value, const std::string& name) :
    GrayImage<T>(xDim,yDim,default_value), name_(name) {}
  
  template<typename T>
  NamedGrayImage<T>::NamedGrayImage(size_t xDim, size_t yDim, T default_value, T max_intensity, 
				    const std::string& name) :
    GrayImage<T>(xDim,yDim,default_value,max_intensity), name_(name) {}
  
  template<typename T>
  NamedGrayImage<T>::NamedGrayImage(const std::string& filename, const std::string& name) :
    GrayImage<T>(filename), name_(name) {}
  
  template<typename T>
  /*virtual*/ const std::string& NamedGrayImage<T>::name() const {
    return name_;
  }
  
  //NOTE: the name is NOT copied
  template<typename T>
  inline void NamedGrayImage<T>::operator=(const NamedGrayImage<T>& toCopy) {
    GrayImage<T>::operator=(toCopy);
  }
  
  template<typename T>
  inline void NamedGrayImage<T>::operator=(const GrayImage<T>& toCopy) {
    GrayImage<T>::operator=(toCopy);
  }

}



#endif



