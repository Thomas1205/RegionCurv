/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#include "draw_segmentation.hh"

void draw_segmentation(const Math2D::Matrix<uint>& segmentation, Math3D::Tensor<float>& image,
		       float r, float g, float b) {

  assert(segmentation.xDim() == image.xDim());
  assert(segmentation.yDim() == image.yDim());
  assert(image.zDim() == 3);
  
  uint xDim = uint( image.xDim() );
  uint yDim = uint( image.yDim() );

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {
      
      uint cur_seg = segmentation(x,y);

      bool edge = false;
      if (x > 0 && segmentation(x-1,y) != cur_seg)
	edge = true;
      if (x+1 < xDim  && segmentation(x+1,y) != cur_seg)
	edge = true;

      if (y > 0 && segmentation(x,y-1) != cur_seg)
	edge = true;
      if (y+1 < yDim && segmentation(x,y+1) != cur_seg)
	edge = true;

      if (edge) {
	image(x,y,0) = r;
	image(x,y,1) = g;
	image(x,y,2) = b;
      }
    }
  }

}
