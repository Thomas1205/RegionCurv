/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#include "conversion.hh"

void make_color_image(const Math2D::Matrix<float>& source, Math3D::ColorImage<float>& target) {

  uint xDim = uint( source.xDim() );
  uint yDim = uint( source.yDim() );

  target.resize_dirty(xDim,yDim,3);

  for (uint z=0; z < 3; z++)
    for (uint y=0; y < yDim; y++)
      for (uint x=0; x < xDim; x++)
	target(x,y,z) = source(x,y);
  
  if (source.max() <= 255)
    target.set_max_intensity(255);
  else 
    target.set_max_intensity(65535);
}

void make_color_image(const Math3D::Tensor<float>& source, Math3D::ColorImage<float>& target) {

  if (source.zDim() == 0) {

    INTERNAL_ERROR << " color image with 0 channels encountered. Exiting..." << std::endl;
    exit(0);
  }
  else if (source.zDim() == 3)
    target = source;
  else {

    if (source.zDim() != 1)
      std::cerr << " WARNING: conversion to color image from tensor with " << source.zDim() << " channels" << std::endl;
    uint xDim = uint( source.xDim() );
    uint yDim = uint( source.yDim() );
    target.resize_dirty(xDim,yDim,3);
    
    for (uint z=0; z < 3; z++)
      for (uint y=0; y < yDim; y++)
	for (uint x=0; x < xDim; x++)
	  target(x,y,z) = source(x,y,0);
  }
  
  
  if (source.max() <= 255)
    target.set_max_intensity(255);
  else 
    target.set_max_intensity(65535);
}
