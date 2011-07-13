/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef DRAW_SEGMENTATION_HH
#define DRAW_SEGMENTATION_HH

#include "matrix.hh"
#include "colorimage.hh"

void draw_segmentation(const Math2D::Matrix<uint>& segmentation, Math3D::Tensor<float>& image,
		       float r=255, float g=255, float b=255);


#endif
