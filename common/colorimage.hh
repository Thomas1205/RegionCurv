#ifndef COLORIMAGE_HH
#define COLORIMAGE_HH

#include "tensor.hh"

namespace Math3D {

  /*** ColorImage ***/
  template<typename T>
  class ColorImage : public Tensor<T> {
  public:

    ColorImage();
    
    ColorImage(std::string filename);

    ColorImage(size_t xDim, size_t yDim, size_t zDim);

    ColorImage(size_t xDim, size_t yDim, size_t zDim, T default_value);

    ColorImage(size_t xDim, size_t yDim, size_t zDim, T default_value, T max_intensity);

    //@returns if the operation was successful
    bool loadPPM(std::string filename);

    //@returns if the operation was successful   
    bool savePPM(std::string filename, bool fit_to_range = true) const;

    virtual const std::string& name() const;

    T max_intensity() const;

    void set_max_intensity(T max_intensity);

    inline void operator=(const Tensor<T>& toCopy);
 
  protected:
    T max_intensity_;
    static const std::string colorim_name_;
  };


  /*** NamedColorImage ***/
  template<typename T>
  class NamedColorImage : public ColorImage<T> {
  public:
    
    NamedColorImage();

    NamedColorImage(std::string name);

    NamedColorImage(std::string filename, std::string name);

    NamedColorImage(size_t xDim, size_t yDim, size_t zDim, std::string name);

    NamedColorImage(size_t xDim, size_t yDim, size_t zDim, T default_value, std::string name);

    NamedColorImage(size_t xDim, size_t yDim, size_t zDim, T default_value, T max_intensity, std::string name);

    ~NamedColorImage();

    inline void operator=(const Tensor<T>& toCopy);

    inline void operator=(const ColorImage<T>& toCopy);

    inline void operator=(const NamedColorImage<T>& toCopy);

    virtual const std::string& name() const;
    
  protected:
    std::string name_;
  };

/********************************* implementation ***************************************/

  template<typename T>
  /*static*/ const std::string ColorImage<T>::colorim_name_;

  template<typename T>
  ColorImage<T>::ColorImage() : Tensor<T>(), max_intensity_(255) {}
  
  template<typename T>
  ColorImage<T>::ColorImage(std::string filename) : Tensor<T>(), max_intensity_(255) {
    loadPPM(filename);
  }

  template<typename T>
  ColorImage<T>::ColorImage(size_t xDim, size_t yDim, size_t zDim) : 
    Tensor<T>(xDim,yDim,zDim), max_intensity_(255) {
  }

  template<typename T>
  ColorImage<T>::ColorImage(size_t xDim, size_t yDim, size_t zDim, T default_value) : 
    Tensor<T>(xDim,yDim,zDim,default_value), max_intensity_(255) {
  } 

  template<typename T>
  ColorImage<T>::ColorImage(size_t xDim, size_t yDim, size_t zDim, T default_value, T max_intensity) :
    Tensor<T>(xDim,yDim,zDim,default_value), max_intensity_(max_intensity) {
  }

  //@returns if the operation was successful
  template<typename T>
  bool ColorImage<T>::loadPPM(std::string filename) {
    
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
	Tensor<T>::zDim_ = 3;
      }
      else if (number == 5) {
	Tensor<T>::zDim_ = 1;
      }
      else {
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
	  c = read_natural_number(fptr,Tensor<T>::xDim_);
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

      //std::cerr << "xDim = " << Tensor<T>::xDim_ << std::endl;

      //3. read y-dimension
      c = read_natural_number(fptr,Tensor<T>::yDim_);
      if (!is_whitespace(c)) 
	throw InvalidCharacterException(c);

      //std::cerr << "yDim = " << Tensor<T>::yDim_ << std::endl;

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
      
      if (Tensor<T>::data_ != 0)
	delete[] Tensor<T>::data_;
      
      Tensor<T>::size_ = Tensor<T>::xDim_*Tensor<T>::yDim_*Tensor<T>::zDim_;
      Tensor<T>::data_ = new T[Tensor<T>::size_];

      //NOTE: currently we assume ther is no line-feed character after the newline

      for (size_t i=0; i < Tensor<T>::size_; i++) {
	
	if (max_intensity_ < 256) {
	  nRead = fread(&c,1,1,fptr);
	  Tensor<T>::data_[i] = c;
	}
	else {
	  nRead = fread(&wc,2,1,fptr);
	  Tensor<T>::data_[i] = wc;
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
      Tensor<T>::xDim_ = 0;
      Tensor<T>::yDim_ = 0;
      Tensor<T>::size_ = Tensor<T>::xDim_*Tensor<T>::yDim_;
      return false;
    }
    catch (InvalidCharacterException ice) {
      IO_ERROR << "    read invalid character \"" << ice.c_ 
		<< "\"." << std::endl;
      Tensor<T>::xDim_ = 0;
      Tensor<T>::yDim_ = 0;
      Tensor<T>::size_ = Tensor<T>::xDim_*Tensor<T>::yDim_;
      return false;
    }

    return true;

    fclose(fptr);
  }

  //@returns if the operation was successful
  template<typename T>
  bool ColorImage<T>::savePPM(std::string filename, bool fit_to_range) const {
    return Tensor<T>::savePPM(filename,max_intensity_,fit_to_range);
  }

  template<typename T>
  /*virtual*/ const std::string& ColorImage<T>::name() const {
    return colorim_name_;
  }

  template<typename T>
  T ColorImage<T>::max_intensity() const {
    return max_intensity_;
  }
  
  template<typename T>
  void ColorImage<T>::set_max_intensity(T max_intensity) {
    max_intensity_ = max_intensity;
  }

  template<typename T>
  inline void ColorImage<T>::operator=(const Tensor<T>& toCopy) {
    Tensor<T>::operator=(toCopy);
  }


  /****** implementation of NamedColorImage **************/
  template<typename T>
  NamedColorImage<T>::NamedColorImage() : ColorImage<T>(), name_("yyy") {}

  template<typename T>
  NamedColorImage<T>::NamedColorImage(std::string name) : ColorImage<T>(), name_(name) {}

  template<typename T>
  NamedColorImage<T>::NamedColorImage(std::string filename, std::string name) : ColorImage<T>(filename), name_(name) {}
  
  template<typename T>
  NamedColorImage<T>::NamedColorImage(size_t xDim, size_t yDim, size_t zDim, std::string name) :
    ColorImage<T>(xDim,yDim,zDim), name_(name) {}
  
  template<typename T>
  NamedColorImage<T>::NamedColorImage(size_t xDim, size_t yDim, size_t zDim, T default_value, std::string name) :
    ColorImage<T>(xDim,yDim,zDim,default_value), name_(name) {}
  
  template<typename T>
  NamedColorImage<T>::NamedColorImage(size_t xDim, size_t yDim, size_t zDim, T default_value, 
				      T max_intensity, std::string name) :
    ColorImage<T>(xDim,yDim,zDim,default_value,max_intensity), name_(name) {}
  
  template<typename T>
  NamedColorImage<T>::~NamedColorImage() {}

  template<typename T>
  inline void NamedColorImage<T>::operator=(const Tensor<T>& toCopy) {
    ColorImage<T>::operator=(toCopy);
  }
  
  template<typename T>
  inline void NamedColorImage<T>::operator=(const ColorImage<T>& toCopy) {
    ColorImage<T>::operator=(toCopy);
  }
  
  template<typename T>
  inline void NamedColorImage<T>::operator=(const NamedColorImage<T>& toCopy) {
    ColorImage<T>::operator=(static_cast<const ColorImage<T>&>(toCopy));
  }

  template<typename T>
  /*virtual*/ const std::string& NamedColorImage<T>::name() const {
    return name_;
  }

  
} //end of namespace Math3D


#endif
