/***** written by Thomas Schoenemann as an employee of Lund University, Sweden, August 2010 ******/

#ifndef LINE_DRAWING_HH
#define LINE_DRAWING_HH

#include "vector.hh"
#include "matrix.hh"
#include "tensor.hh"

template<typename T>
void draw_line(Math3D::Tensor<T>& tensor, uint x1, uint y1, uint x2, uint y2,
	       const Math1D::Vector<T>& marker);

template<typename T>
void draw_line(Math2D::Matrix<T>& matrix, uint x1, uint y1, uint x2, uint y2,
	       T marker);


/******************* implementation ********************/

template<typename T>
void draw_line(Math3D::Tensor<T>& tensor, uint x1, uint y1, uint x2, uint y2,
	       const Math1D::Vector<T>& marker) {

  if (x1 >= tensor.xDim())
    x1 = tensor.xDim()-1;
  if (x2 >= tensor.xDim())
    x2 = tensor.xDim()-1;
  if (y1 >= tensor.yDim())
    y1 = tensor.yDim()-1;
  if (y2 >= tensor.yDim())
    y2 = tensor.yDim()-1;

  uint zDim = tensor.zDim();

  assert(zDim == marker.size());

  float dx = ((int) x2) - ((int) x1);
  float dy = ((int) y2) - ((int) y1);
  
  if (fabs(dx) >= fabs(dy)) {

    if (dx < 0) {
      dx *= -1.0; //note: this cancels out again
      dy *= -1.0;
      std::swap(x1,x2);
      std::swap(y1,y2);
    }

    float slope = dy/dx;
    float fy = y1;

    for (uint x=x1; x <= x2; x++) {

      if (x != x1)
	fy += slope;
      uint y = round(fy);

      for (uint z=0; z < zDim; z++)
	tensor(x,y,z) = marker[z];
    }

  }
  else {

    if (dy < 0) {
      dy *= -1.0; //note: this cancels out again
      dx *= -1.0;
      std::swap(x1,x2);
      std::swap(y1,y2);
    }

    float slope = dx/dy;
    float fx = x1;

    for (uint y=y1; y <= y2; y++) {

      if (y != y1)
	fx += slope;

      uint x = round(fx);
      
      for (uint z=0; z < zDim; z++)
	tensor(x,y,z) = marker[z];      
    }
  }
}


template<typename T>
void draw_line(Math2D::Matrix<T>& matrix, uint x1, uint y1, uint x2, uint y2,
	       T marker) {

  if (x1 >= matrix.xDim())
    x1 = uint( matrix.xDim()-1 );
  if (x2 >= matrix.xDim())
    x2 = uint( matrix.xDim()-1 );
  if (y1 >= matrix.yDim())
    y1 = uint( matrix.yDim()-1 );
  if (y2 >= matrix.yDim())
    y2 = uint( matrix.yDim()-1 );

  float dx = float(((int) x2) - ((int) x1));
  float dy = float(((int) y2) - ((int) y1));
  
  if (fabs(dx) >= fabs(dy)) {

    if (dx < 0) {
      dx *= -1.0; //note: this cancels out again
      dy *= -1.0;
      std::swap(x1,x2);
      std::swap(y1,y2);
    }

    float slope = dy/dx;
    float fy = float(y1);

    for (uint x=x1; x <= x2; x++) {

      if (x != x1)
	fy += slope;
      uint y = uint(fy+0.5);
      
      matrix(x,y) = marker;
    }

  }
  else {

    if (dy < 0) {
      dy *= -1.0; //note: this cancels out again
      dx *= -1.0;
      std::swap(x1,x2);
      std::swap(y1,y2);
    }

    float slope = dx/dy;
    float fx = float(x1);

    for (uint y=y1; y <= y2; y++) {

      if (y != y1)
	fx += slope;

      uint x = uint(fx+0.5);
      
      matrix(x,y) = marker;      
    }
  }

}



#endif
