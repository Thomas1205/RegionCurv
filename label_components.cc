/********** written by Thomas Schoenemann as an employee of Lund University, Sweden, 2010 ***************/

#include "label_components.hh"

void label(const Math2D::Matrix<float>& mask, Math2D::Matrix<uint>& components, uint x, uint y, uint cur_label) {

  components(x,y) = cur_label;

  if (x > 0 && mask(x-1,y) < 128.0 && components(x-1,y) == 0)
    label(mask,components,x-1,y,cur_label);
  
  if (x+1 < components.xDim() && mask(x+1,y) < 128.0 && components(x+1,y) == 0)
    label(mask,components,x+1,y,cur_label);
  
  if (y > 0 && mask(x,y-1) < 128.0 && components(x,y-1) == 0)
    label(mask,components,x,y-1,cur_label);
  
  if (y+1 < components.yDim() && mask(x,y+1) < 128.0 && components(x,y+1) == 0)
    label(mask,components,x,y+1,cur_label);  
}

uint label_components(const Math2D::Matrix<float>& mask, Math2D::Matrix<uint>& components) {

  uint xDim = mask.xDim();
  uint yDim = mask.yDim();

  components.resize(xDim,yDim);
  components.set_constant(0);
  
  uint next_component = 1;

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {

      if (mask(x,y) == 0 && components(x,y) == 0) {

	label(mask,components,x,y,next_component);

	next_component++;
      }
    }
  }

  return next_component-1;
}
