/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef CONVERSION_HH
#define CONVERSION_HH

#include "grayimage.hh"
#include "colorimage.hh"

void make_color_image(const Math2D::Matrix<float>& source, Math3D::ColorImage<float>& target);

void make_color_image(const Math3D::Tensor<float>& source, Math3D::ColorImage<float>& target);



#endif
