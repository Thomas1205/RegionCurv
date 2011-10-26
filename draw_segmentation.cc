/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#include "draw_segmentation.hh"

int intabs(int i)
{
  return i>0?i:-i;
}

bool isedge(uint a, uint b)
{
  return intabs(int(a)-int(b)) > 100;
}

void draw_segmentation(const Math2D::Matrix<uint>& segmentation, Math3D::Tensor<float>& image,
		       float r, float g, float b) {

  assert(segmentation.xDim() == image.xDim());
  assert(segmentation.yDim() == image.yDim());
  assert(image.zDim() == 3);
  
  uint xDim = uint( image.xDim() );
  uint yDim = uint( image.yDim() );

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {
      
      int cur_seg = segmentation(x,y);

      bool edge = false;
      if (x > 0 && isedge(segmentation(x-1,y),cur_seg) )
	edge = true;
      if (x+1 < xDim  && isedge(segmentation(x+1,y),cur_seg) )
	edge = true;

      if (y > 0 && isedge(segmentation(x,y-1),cur_seg) )
	edge = true;
      if (y+1 < yDim && isedge(segmentation(x,y+1),cur_seg) )
	edge = true;

      if (edge) {
	image(x,y,0) = r;
	image(x,y,1) = g;
	image(x,y,2) = b;
      }
    }
  }

}
