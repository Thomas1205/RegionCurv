/******** written by Thomas Schoenemann as an employee of Lund University, Sweden, June 2010 ****/
/******** hybrid version writtem by Yubin Kuang as an employee of Lund University, Sweden, September 2010 ****/

#include"lp_inpainting.hh"

#include "mesh2D.hh"
#include "sparse_matrix_description.hh"
#include <coin/ClpSimplex.hpp>
#include <coin/ClpFactorization.hpp>

#include "tensor.hh"
#include "outer_product.hh"
#include "timing.hh"
#include "smoothing.hh"

template<typename T>
int round(T t) 
{
	return int(t + 0.5);
}

enum InpaintStatus {ToInpaint, Border, Irrelevant};
enum EdgeStatus {InteriorEdge, BorderEdge};

#ifdef HAS_GUROBI
#include "gurobi_c++.h"
#endif

#ifdef HAS_CPLEX
#include <ilcplex/cplex.h>
#endif

#ifdef HAS_XPRESS
#include "xprs.h" 
#endif

//@returns the number of faces per pixel
uint generate_mesh(uint neighborhood, const Math2D::Matrix<InpaintStatus>& status,  
		   Math2D::Matrix<uint>& face_offset, Math2D::Matrix<uint>& nFaces,
		   Mesh2D& mesh, bool legacy = false) {

  //DEBUG
  //legacy = true;
  //END_DEBUG

  const uint xDim = status.xDim();
  const uint yDim = status.yDim();

  uint nAreasPerPixel = 1;
  if (neighborhood == 4) {
    nAreasPerPixel = 1;
  }
  else if (neighborhood == 8) {
    nAreasPerPixel = 4;
  }
  else if (neighborhood == 16) {
    nAreasPerPixel = 32;
  }
  else {
    INTERNAL_ERROR << "invalid neighborhood \"" << neighborhood << "\". Exiting." << std::endl;
    exit(1);
  }

  /** generate mesh points **/
  for (uint y = 0; y <= yDim; y++) {
    for (uint x = 0; x <= xDim; x++) {
      mesh.add_point(Mesh2DPoint(x,y));
    }
  }

  const uint diagPointOffset = mesh.nPoints();

  Math2D::Matrix<uint> point_offset(xDim,yDim,MAX_UINT/2);

  if (neighborhood == 8) {
    for (uint y=0; y < yDim; y++)
      for (uint x=0; x < xDim; x++)
	mesh.add_point(Mesh2DPoint(x+0.5,y+0.5));
  }
  else if (neighborhood == 16) {

    for (uint y=0; y <= yDim; y++) {
      for (uint x=0; x < xDim; x++)
	mesh.add_point( Mesh2DPoint(x+0.5, y) );
    }

    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x <= xDim; x++)
	mesh.add_point( Mesh2DPoint(x, y+0.5) );
    }

    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x < xDim; x++) {

	point_offset(x,y) = mesh.nPoints();

	mesh.add_point(Mesh2DPoint(x+0.4,y+0.2));
	mesh.add_point(Mesh2DPoint(x+0.6,y+0.2));
	mesh.add_point(Mesh2DPoint(x+0.5,y+0.25)); 
	mesh.add_point(Mesh2DPoint(x+0.33333333,y+0.33333333));

	mesh.add_point(Mesh2DPoint(x+0.66666667,y+0.33333333));
	mesh.add_point(Mesh2DPoint(x+0.2,y+0.4));
	mesh.add_point(Mesh2DPoint(x+0.8,y+0.4)); 
	mesh.add_point(Mesh2DPoint(x+0.25,y+0.5));

	mesh.add_point(Mesh2DPoint(x+0.5,y+0.5));
	mesh.add_point(Mesh2DPoint(x+0.75,y+0.5)); 
	mesh.add_point(Mesh2DPoint(x+0.2,y+0.6));
	mesh.add_point(Mesh2DPoint(x+0.33333333,y+0.666666667));

	mesh.add_point(Mesh2DPoint(x+0.66666667,y+0.666666667));
	mesh.add_point(Mesh2DPoint(x+0.8,y+0.6)); 
	mesh.add_point(Mesh2DPoint(x+0.5,y+0.75));

	mesh.add_point(Mesh2DPoint(x+0.4,y+0.8));
	mesh.add_point(Mesh2DPoint(x+0.6,y+0.8)); 
      }
    }

  }

  std::cerr << "generating edges" << std::endl;

  /** generate mesh edges **/
  if (neighborhood <= 8) {
    //vertical edges
    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x <= xDim; x++) {

	if ((x < xDim && status(x,y)) < Irrelevant || (x > 0 && status(x-1,y) < Irrelevant))
	  mesh.add_edge(y*(xDim+1)+x, (y+1)*(xDim+1)+x);
      }
    }

    //horizontal edges
    for (uint y=0; y <= yDim; y++) {
      for (uint x=0; x < xDim; x++) {
	
	if ((y < yDim && status(x,y) < Irrelevant) || (y > 0 && status(x,y-1) < Irrelevant))
	    mesh.add_edge(y*(xDim+1)+x,y*(xDim+1)+x+1);
      }
    }

    //diagonal edges
    if (neighborhood == 8) {

      for (uint y=0; y < yDim; y++) {
	for (uint x=0; x < xDim; x++) {

	  if (status(x,y) < Irrelevant && (legacy || status(x,y) != Border)) {
	  
	    const uint diagPoint = diagPointOffset+y*xDim+x;
	    
	    mesh.add_edge(y*(xDim+1)+x,diagPoint);
	    mesh.add_edge(y*(xDim+1)+x+1,diagPoint);
	    mesh.add_edge((y+1)*(xDim+1)+x,diagPoint);
	    mesh.add_edge((y+1)*(xDim+1)+x+1,diagPoint);
	  }
	}
      }
    }
  }

  if (neighborhood == 16) {
    
    uint hmidpoint_offset = (yDim+1)*(xDim+1);

    for (uint y=0; y <= yDim; y++) {
      for (uint x=0; x < xDim; x++) {
      
	if ((y < yDim && status(x,y) < Irrelevant) || (y > 0 && status(x,y-1) < Irrelevant)) { 

	  mesh.add_edge( y*(xDim+1)+x, hmidpoint_offset + y*(xDim)+x );
	  mesh.add_edge( y*(xDim+1)+x+1, hmidpoint_offset + y*(xDim)+x );
	}
      }
    }

    uint vmidpoint_offset = hmidpoint_offset + (yDim+1)*xDim;

    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x <= xDim; x++) {

	if ((x < xDim && status(x,y)) < Irrelevant || (x > 0 && status(x-1,y) < Irrelevant)) {

	  mesh.add_edge( y*(xDim+1)+x, vmidpoint_offset + y*(xDim+1)+x );
	  mesh.add_edge( (y+1)*(xDim+1)+x, vmidpoint_offset + y*(xDim+1)+x );
	}
      }
    }

    //add interior edges

    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x < xDim; x++) {
	  
	if (status(x,y) < Irrelevant && (legacy || status(x,y) != Border)) {
	    
	  const uint po = point_offset(x,y);
	  
	  mesh.add_edge( y*(xDim+1)+x , po + 0);
	  mesh.add_edge(  hmidpoint_offset + y*xDim+x , po + 0);
	  mesh.add_edge(  hmidpoint_offset + y*xDim+x , po + 1);
	  mesh.add_edge( y*(xDim+1)+x+1 , po + 1); //edge 4
	  
	  mesh.add_edge( y*(xDim+1)+x , po + 5);
	  mesh.add_edge( y*(xDim+1)+x , po + 3);
	  mesh.add_edge( po + 0 , po + 3);
	  mesh.add_edge( po + 0 , po + 2); //edge 8
	  
	  mesh.add_edge( po + 1 , po + 2);
	  mesh.add_edge( po + 1 , po + 4);
	  mesh.add_edge( y*(xDim+1)+x+1 , po + 4 );
	  mesh.add_edge( y*(xDim+1)+x+1 , po + 6 ); //edge 12
	  
	  mesh.add_edge( vmidpoint_offset + y*(xDim+1)+x , po + 5 );
	  mesh.add_edge( po + 3 , po + 5);
	  mesh.add_edge( po + 2 , po + 3);
	  mesh.add_edge( po + 2 , po + 4);  //edge 16
	    
	  mesh.add_edge( po + 4 , po + 6);
	  mesh.add_edge( vmidpoint_offset + y*(xDim+1)+x+1 , po + 6);
	  mesh.add_edge( po + 5 , po + 7);
	  mesh.add_edge( po + 3 , po + 7); //edge 20
	  
	  mesh.add_edge( po + 3, po + 8);
	  mesh.add_edge( po + 4, po + 8);
	  mesh.add_edge( po + 4, po + 9);
	  mesh.add_edge( po + 6, po + 9); //edge 24
	  
	  mesh.add_edge( vmidpoint_offset + y*(xDim+1)+x , po + 10 );
	  mesh.add_edge( po + 7, po + 10 );
	  mesh.add_edge( po + 7, po + 11 );
	  mesh.add_edge( po + 8, po + 11 ); //edge 28
	  
	  mesh.add_edge( po + 8, po + 12 );
	  mesh.add_edge( po + 9, po + 12 );
	  mesh.add_edge( po + 9, po + 13 );
	  mesh.add_edge( vmidpoint_offset + y*(xDim+1)+x+1, po + 13 ); //edge 32
	    
	  mesh.add_edge( po + 10, po + 11 ); 
	  mesh.add_edge( po + 12, po + 13 ); 
	  mesh.add_edge( (y+1)*(xDim+1)+x, po + 10  );
	  mesh.add_edge( (y+1)*(xDim+1)+x, po + 11  );  //edge 36
	    
	  mesh.add_edge( po + 11, po + 15 );
	  mesh.add_edge( po + 11, po + 14 );
	  mesh.add_edge( po + 12, po + 14 );
	  mesh.add_edge( po + 12, po + 16 );  //edge 40
	  
	  mesh.add_edge( (y+1)*(xDim+1)+x+1, po + 12  );
	  mesh.add_edge( (y+1)*(xDim+1)+x+1, po + 13  );
	  mesh.add_edge( (y+1)*(xDim+1)+x, po + 15  );
	  mesh.add_edge( po + 14, po + 15 ); //edge 44
	  
	  mesh.add_edge( po + 14, po + 16 );
	  mesh.add_edge( (y+1)*(xDim+1)+x+1, po + 16 );
	  mesh.add_edge(  hmidpoint_offset + (y+1)*xDim+x , po + 15);
	  mesh.add_edge(  hmidpoint_offset + (y+1)*xDim+x , po + 16);
	}
      }    
    }
  }

  std::cerr << "generating faces" << std::endl;

  std::vector<uint> face;

  /** generate mesh faces **/
  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {

      if (status(x,y) < Irrelevant) {

	face_offset(x,y) = mesh.nFaces();

	if (neighborhood==4) {

	  nFaces(x,y) = 1;

	  face.clear();
	  face.push_back(mesh.find_edge(y*(xDim+1)+x,(y+1)*(xDim+1)+x).first);
	  face.push_back(mesh.find_edge((y+1)*(xDim+1)+x,(y+1)*(xDim+1)+x+1).first);
	  face.push_back(mesh.find_edge(y*(xDim+1)+x+1,(y+1)*(xDim+1)+x+1).first);
	  face.push_back(mesh.find_edge(y*(xDim+1)+x,y*(xDim+1)+x+1).first);
	
	  mesh.add_face(face);
	}
	else if (neighborhood==8) {

	  if (!legacy && status(x,y) == Border) {

	    nFaces(x,y) = 1;

	    face.clear();
	    face.push_back(mesh.find_edge(y*(xDim+1)+x,(y+1)*(xDim+1)+x).first);
	    face.push_back(mesh.find_edge((y+1)*(xDim+1)+x,(y+1)*(xDim+1)+x+1).first);
	    face.push_back(mesh.find_edge(y*(xDim+1)+x+1,(y+1)*(xDim+1)+x+1).first);
	    face.push_back(mesh.find_edge(y*(xDim+1)+x,y*(xDim+1)+x+1).first);
	    
	    mesh.add_face(face);
	  }
	  else {

	    nFaces(x,y) = 4;
	    
	    const uint diagPoint = diagPointOffset+y*xDim+x;
	    
	    face.clear();
	    face.push_back(mesh.find_edge(y*(xDim+1)+x,(y+1)*(xDim+1)+x).first);
	    face.push_back(mesh.find_edge((y+1)*(xDim+1)+x,diagPoint).first);
	    face.push_back(mesh.find_edge(y*(xDim+1)+x,diagPoint).first);
	    mesh.add_face(face);
	    
	    face.clear();
	    face.push_back(mesh.find_edge((y+1)*(xDim+1)+x,(y+1)*(xDim+1)+x+1).first);
	    face.push_back(mesh.find_edge((y+1)*(xDim+1)+x+1,diagPoint).first);
	    face.push_back(mesh.find_edge((y+1)*(xDim+1)+x,diagPoint).first);
	    mesh.add_face(face);
	    
	    face.clear();
	    face.push_back(mesh.find_edge((y+1)*(xDim+1)+x+1,y*(xDim+1)+x+1).first);
	    face.push_back(mesh.find_edge(y*(xDim+1)+x+1,diagPoint).first);
	    face.push_back(mesh.find_edge((y+1)*(xDim+1)+x+1,diagPoint).first);
	    mesh.add_face(face);
	    
	    face.clear();
	    face.push_back(mesh.find_edge(y*(xDim+1)+x,y*(xDim+1)+x+1).first);
	    face.push_back(mesh.find_edge(y*(xDim+1)+x,diagPoint).first);
	    face.push_back(mesh.find_edge(y*(xDim+1)+x+1,diagPoint).first);
	    mesh.add_face(face);
	  }
	}
	else if (neighborhood==16) {

	  if (!legacy && status(x,y) == Border) {

	    nFaces(x,y) = 1;
	    
	    uint hmidpoint_offset = (yDim+1)*(xDim+1);
	    uint vmidpoint_offset = hmidpoint_offset + (yDim+1)*xDim;

	    face.clear();
	    face.push_back(mesh.find_edge(y*(xDim+1)+x,vmidpoint_offset + y*(xDim+1)+x).first);
	    face.push_back(mesh.find_edge(vmidpoint_offset + y*(xDim+1)+x,(y+1)*(xDim+1)+x).first);
	    
	    face.push_back(mesh.find_edge((y+1)*(xDim+1)+x,hmidpoint_offset + (y+1)*(xDim)+x).first);
	    face.push_back(mesh.find_edge(hmidpoint_offset + (y+1)*(xDim)+x,(y+1)*(xDim+1)+x+1).first);

	    face.push_back(mesh.find_edge(vmidpoint_offset + y*(xDim+1)+x+1,(y+1)*(xDim+1)+x+1).first);	    
	    face.push_back(mesh.find_edge(y*(xDim+1)+x+1,vmidpoint_offset + y*(xDim+1)+x+1).first);
	    
	    face.push_back(mesh.find_edge(hmidpoint_offset + y*(xDim)+x,y*(xDim+1)+x+1).first);
	    face.push_back(mesh.find_edge(y*(xDim+1)+x,hmidpoint_offset + y*(xDim)+x).first);

	    mesh.add_face(face);
	  }
	  else {

	    nFaces(x,y) = 32;
	    
	    const uint po = point_offset(x,y);
	    const uint hmidpoint_offset = (yDim+1)*(xDim+1);
	    const uint vmidpoint_offset = hmidpoint_offset + (yDim+1)*xDim;
	  
	    face.clear(); //face 0
	    face.push_back( mesh.find_edge( y*(xDim+1)+x, po + 0 ).first );
	    face.push_back( mesh.find_edge( po + 0, hmidpoint_offset + y*xDim+x).first );
	    face.push_back( mesh.find_edge( y*(xDim+1)+x, hmidpoint_offset + y*xDim+x  ).first );
	    mesh.add_face(face);
	    
	    
	    face.clear(); //face 1
	    face.push_back( mesh.find_edge( hmidpoint_offset + y*xDim+x, po + 0).first );
	    face.push_back( mesh.find_edge( po + 0, po + 2).first );
	    face.push_back( mesh.find_edge( po + 2, po + 1).first );
	    face.push_back( mesh.find_edge( po + 1, hmidpoint_offset + y*xDim+x).first );
	    mesh.add_face(face);
	    
	    face.clear(); //face 2
	    face.push_back( mesh.find_edge( hmidpoint_offset + y*xDim+x, po + 1).first );
	    face.push_back( mesh.find_edge( po + 1, y*(xDim+1) + x+1).first );
	    face.push_back( mesh.find_edge( y*(xDim+1) + x+1, hmidpoint_offset + y*xDim+x).first );
	    mesh.add_face(face);
	    
	    face.clear(); //face 3
	    face.push_back( mesh.find_edge( y*(xDim+1)+x, vmidpoint_offset + y*(xDim+1)+x ).first );
	    face.push_back( mesh.find_edge( vmidpoint_offset + y*(xDim+1)+x, po + 5 ).first );
	    face.push_back( mesh.find_edge( po + 5, y*(xDim+1)+x ).first );
	    mesh.add_face(face);
	  
	    face.clear(); //face 4
	    face.push_back( mesh.find_edge( y*(xDim+1)+x, po + 5).first );
	    face.push_back( mesh.find_edge( po + 5, po + 3).first );
	    face.push_back( mesh.find_edge( po + 3, y*(xDim+1)+x).first );
	    mesh.add_face(face);
	    
	    face.clear(); //face 5
	    face.push_back( mesh.find_edge( y*(xDim+1)+x, po + 3).first );
	    face.push_back( mesh.find_edge( po + 3, po + 0).first );
	    face.push_back( mesh.find_edge( po + 0, y*(xDim+1)+x).first );
	    mesh.add_face(face);
	    
	    face.clear(); //face 6
	    face.push_back( mesh.find_edge( po + 0, po + 3).first );
	    face.push_back( mesh.find_edge( po + 3, po + 2).first );
	    face.push_back( mesh.find_edge( po + 2, po + 0).first );
	    mesh.add_face(face);
	  
	    face.clear(); //face 7
	    face.push_back( mesh.find_edge( po + 1, po + 2).first );
	    face.push_back( mesh.find_edge( po + 2, po + 4).first );
	    face.push_back( mesh.find_edge( po + 4, po + 1).first );
	    mesh.add_face(face);
	    
	    face.clear(); //face 8
	    face.push_back( mesh.find_edge( y*(xDim+1)+x+1, po + 1).first );
	    face.push_back( mesh.find_edge( po + 1, po + 4).first );
	    face.push_back( mesh.find_edge( po + 4, y*(xDim+1)+x+1).first );
	    mesh.add_face(face);
	    
	    face.clear(); //face 9
	    face.push_back( mesh.find_edge( y*(xDim+1)+x+1, po + 4).first );
	    face.push_back( mesh.find_edge( po + 4, po + 6).first );
	    face.push_back( mesh.find_edge( po + 6, y*(xDim+1)+x+1).first );
	    mesh.add_face(face);
	    
	    face.clear(); //face 10
	    face.push_back( mesh.find_edge( y*(xDim+1)+x+1, po + 6).first );
	    face.push_back( mesh.find_edge( po + 6, vmidpoint_offset + y*(xDim+1)+x+1).first );
	    face.push_back( mesh.find_edge( vmidpoint_offset + y*(xDim+1)+x+1, y*(xDim+1)+x+1 ).first );
	    mesh.add_face(face);
	    
	    face.clear(); //face 11
	    face.push_back( mesh.find_edge( vmidpoint_offset + y*(xDim+1)+x, po + 10 ).first );
	    face.push_back( mesh.find_edge( po + 10, po + 7).first );
	    face.push_back( mesh.find_edge( po + 7, po + 5).first );
	    face.push_back( mesh.find_edge( po + 5, vmidpoint_offset + y*(xDim+1)+x).first );
	    mesh.add_face(face);

	    face.clear(); //face 12
	    face.push_back( mesh.find_edge( po + 3, po + 5).first );
	    face.push_back( mesh.find_edge( po + 5, po + 7).first );
	    face.push_back( mesh.find_edge( po + 7, po + 3).first );
	    mesh.add_face(face);
	    
	    face.clear(); //face 13
	    face.push_back( mesh.find_edge( po + 2, po + 3).first );
	    face.push_back( mesh.find_edge( po + 3, po + 8).first );
	    face.push_back( mesh.find_edge( po + 8, po + 4).first );
	    face.push_back( mesh.find_edge( po + 4, po + 2).first );
	    mesh.add_face(face);
	    
	    face.clear(); //face 14
	    face.push_back( mesh.find_edge( po + 4, po + 9).first );
	    face.push_back( mesh.find_edge( po + 9, po + 6).first );
	    face.push_back( mesh.find_edge( po + 6, po + 4).first );
	    mesh.add_face(face);
	    
	    face.clear(); //face 15
	    face.push_back( mesh.find_edge( po + 6, po + 9).first );
	    face.push_back( mesh.find_edge( po + 9, po + 13).first );
	    face.push_back( mesh.find_edge( po + 13, vmidpoint_offset + y*(xDim+1)+x+1).first );
	    face.push_back( mesh.find_edge( vmidpoint_offset + y*(xDim+1)+x+1, po + 6).first );
	    mesh.add_face(face);

	    face.clear(); //face 16
	    face.push_back( mesh.find_edge( po + 3, po + 7).first );
	    face.push_back( mesh.find_edge( po + 7, po + 11).first );
	    face.push_back( mesh.find_edge( po + 11, po + 8).first );
	    face.push_back( mesh.find_edge( po + 8, po + 3).first );
	    mesh.add_face(face);

	    face.clear(); //face 17
	    face.push_back( mesh.find_edge( po + 4, po + 8).first );
	    face.push_back( mesh.find_edge( po + 8, po + 12).first );
	    face.push_back( mesh.find_edge( po + 12, po + 9).first );
	    face.push_back( mesh.find_edge( po + 9, po + 4).first );
	    mesh.add_face(face);
	  
	    face.clear(); //face 18
	    face.push_back( mesh.find_edge( po + 7, po + 10).first );
	    face.push_back( mesh.find_edge( po + 10, po + 11).first );
	    face.push_back( mesh.find_edge( po + 11, po + 7).first );
	    mesh.add_face(face);
	  
	    face.clear(); //face 19
	    face.push_back( mesh.find_edge( po + 8, po + 11).first );
	    face.push_back( mesh.find_edge( po + 11, po + 14).first );
	    face.push_back( mesh.find_edge( po + 14, po + 12).first );
	    face.push_back( mesh.find_edge( po + 12, po + 8).first );
	    mesh.add_face(face);

	    face.clear(); //face 20
	    face.push_back( mesh.find_edge( po + 9, po + 12).first );
	    face.push_back( mesh.find_edge( po + 12, po + 13).first );
	    face.push_back( mesh.find_edge( po + 13, po + 9).first );
	    mesh.add_face(face);

	    face.clear(); //face 21
	    face.push_back( mesh.find_edge( vmidpoint_offset + y*(xDim+1)+x, (y+1)*(xDim+1)+x).first );
	    face.push_back( mesh.find_edge( (y+1)*(xDim+1)+x, po + 10).first );
	    face.push_back( mesh.find_edge( po + 10, vmidpoint_offset + y*(xDim+1)+x).first );
	    mesh.add_face(face);

	    face.clear(); //face 22
	    face.push_back( mesh.find_edge( (y+1)*(xDim+1)+x, po + 11).first );
	    face.push_back( mesh.find_edge( po + 11, po + 10).first );
	    face.push_back( mesh.find_edge( po + 10, (y+1)*(xDim+1)+x).first );
	    mesh.add_face(face);

	    face.clear(); //face 23
	    face.push_back( mesh.find_edge( (y+1)*(xDim+1)+x, po + 15).first );
	    face.push_back( mesh.find_edge( po + 15, po + 11).first );
	    face.push_back( mesh.find_edge( po + 11, (y+1)*(xDim+1)+x).first );
	    mesh.add_face(face);

	    face.clear(); //face 24
	    face.push_back( mesh.find_edge( po + 11, po + 15).first );
	    face.push_back( mesh.find_edge( po + 15, po + 14).first );
	    face.push_back( mesh.find_edge( po + 14, po + 11).first );
	    mesh.add_face(face);

	    face.clear(); //face 25
	    face.push_back( mesh.find_edge( po + 12, po + 14).first );
	    face.push_back( mesh.find_edge( po + 14, po + 16).first );
	    face.push_back( mesh.find_edge( po + 16, po + 12).first );
	    mesh.add_face(face);

	    face.clear(); //face 26
	    face.push_back( mesh.find_edge( (y+1)*(xDim+1)+x+1, po + 12).first );
	    face.push_back( mesh.find_edge( po + 12, po + 16).first );
	    face.push_back( mesh.find_edge( po + 16, (y+1)*(xDim+1)+x+1).first );
	    mesh.add_face(face);

	    face.clear(); //face 27
	    face.push_back( mesh.find_edge( (y+1)*(xDim+1)+x+1, po + 13).first );
	    face.push_back( mesh.find_edge( po + 13, po + 12).first );
	    face.push_back( mesh.find_edge( po + 12, (y+1)*(xDim+1)+x+1).first );
	    mesh.add_face(face);
	
	    face.clear(); //face 28
	    face.push_back( mesh.find_edge( (y+1)*(xDim+1)+x+1, vmidpoint_offset + y*(xDim+1)+x+1 ).first );
	    face.push_back( mesh.find_edge( vmidpoint_offset + y*(xDim+1)+x+1, po + 13 ).first );
	    face.push_back( mesh.find_edge( po + 13, (y+1)*(xDim+1)+x+1 ).first );
	    mesh.add_face(face);

	    face.clear(); //face 29
	    face.push_back( mesh.find_edge( (y+1)*(xDim+1)+x, hmidpoint_offset + (y+1)*xDim+x ).first );
	    face.push_back( mesh.find_edge( hmidpoint_offset + (y+1)*xDim+x, po + 15 ).first );
	    face.push_back( mesh.find_edge( po + 15, (y+1)*(xDim+1)+x ).first );
	    mesh.add_face(face);

	    face.clear(); //face 30
	    face.push_back( mesh.find_edge( hmidpoint_offset + (y+1)*xDim+x, po + 16 ).first );
	    face.push_back( mesh.find_edge( po + 16, po + 14).first );
	    face.push_back( mesh.find_edge( po + 14, po + 15).first );
	    face.push_back( mesh.find_edge( po + 15, hmidpoint_offset + (y+1)*xDim+x ).first );
	    mesh.add_face(face);

	    face.clear(); //face 31
	    face.push_back( mesh.find_edge( (y+1)*(xDim+1)+x+1, po + 16 ).first );
	    face.push_back( mesh.find_edge( po + 16, hmidpoint_offset + (y+1)*xDim+x).first );
	    face.push_back( mesh.find_edge( hmidpoint_offset + (y+1)*xDim+x, (y+1)*(xDim+1)+x+1 ).first );
	    mesh.add_face(face);
	  }
	}
	else {
	  INTERNAL_ERROR << "connectivity \"" << neighborhood << "\" is not supported" << std::endl;
	}
      }
    }
  }

  std::cerr << "added " << mesh.nFaces() << " faces " << std::endl;

  if (xDim*yDim < 2000)
    mesh.draw("mesh.svg");
    
  return nAreasPerPixel;
}


//source : http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
void eig2x2(const Math2D::Matrix<double> M, Math1D::Vector<double>& v1, Math1D::Vector<double>& v2) {

  assert(M.xDim() == 2);
  assert(M.yDim() == 2);

  double offs = 0.5 * (M(0,0) + M(1,1));
  double det = M(0,0)*M(1,1) - M(0,1) * M(1,0);

  double lambda1 = offs + sqrt(offs*offs - det);
  double lambda2 = offs - sqrt(offs*offs - det);

  v1.resize(2);
  v2.resize(2);

  if (M(0,1) != 0.0) {
    v1[0] = lambda1 - M(1,1);
    v1[1] = M(0,1);
    v2[0] = lambda2 - M(1,1);
    v2[1] = M(0,1);    
  }
  else if (M(1,0) != 0.0) {
    v1[0] = M(1,0);
    v1[1] = lambda1 - M(0,0);
    v2[0] = M(1,0);
    v2[1] = lambda2 - M(0,0);
  }
  else if (M(0,0) >= M(1,1)) {
    v1[0] = 1.0;
    v1[1] = 0.0;
    v2[0] = 0.0;
    v2[1] = 1.0;
  }
  else {
    v1[0] = 0.0;
    v1[1] = 1.0;
    v2[0] = 1.0;
    v2[1] = 0.0;
  }

  v1 *= 1.0 / v1.norm();
  v2 *= 1.0 / v2.norm();
}


void compute_orientations(const Math2D::Matrix<float>& image, const Math2D::Matrix<InpaintStatus>& status,  
			  Math2D::Matrix<float>& horientation, Math2D::Matrix<float>& vorientation) {

  double rho = 4.0;
  double sigma = 1.5;

  uint xDim = image.xDim();
  uint yDim = image.yDim();

  int bound = 4;

  Math1D::Vector<double> ev1(2);
  Math1D::Vector<double> ev2(2);

  Math2D::Matrix<float> smooth_image = image;
  for (int x=0; x < (int) xDim; x++) {

    for (int y=0; y < (int) yDim; y++) {

      double sum = 0.0;
      double norm = 0.0;

      for (int yy = std::max(0,y-bound); yy <= std::min(y-bound,(int) yDim-1); yy++) {	
	for (int xx = std::max(0,x-bound); xx < std::min((int)xDim,x+bound+1); xx++) {

	  if (status(xx,yy) != ToInpaint) {

	    double sqr_dist_x = (yy-y) * (yy-y);
	    double sqr_dist_y = (xx-x) * (xx-x);
	    
	    double weight = exp(-2.0*sqr_dist_x / sigma) * exp(-2.0*sqr_dist_y / sigma);
	    sum += weight * image(xx,yy);
	    norm += weight;
	  }
	}
      }

      if (norm > 0.0)
	sum /= norm;

      smooth_image(x,y) = sum;
    }
  }  


  //NOTE: undefined regions are assigned a gradient of 0
  Math3D::Tensor<double> imgrad(xDim,yDim,2,0.0);
  for (uint x=0; x < xDim; x++) {

    for (uint y=0; y < yDim; y++) {

      if (status(x,y) != ToInpaint) {

	if (x > 0 && status(x-1,y) != ToInpaint) {
	  imgrad(x,y,0) = smooth_image(x,y) - smooth_image(x-1,y);
	}
	else if (x+1 < xDim && status(x+1,y) != ToInpaint) {
	  imgrad(x,y,0) = smooth_image(x+1,y) - smooth_image(x,y);	  
	}
	
	if (y > 0 && status(x,y-1) != ToInpaint) {
	  imgrad(x,y,1) =  smooth_image(x,y) - smooth_image(x,y-1);
	}
	else if (y+1 < yDim && status(x,y+1) != ToInpaint) {
	  imgrad(x,y,1) =  smooth_image(x,y+1) - smooth_image(x,y);	
	}
      }
    }
  }

  horientation.resize(xDim,yDim+1);
  horientation.set_constant(-100);

  for (int x=0; x < (int) xDim; x++) {

    //std::cerr << "x: " << x << std::endl;

    for (int y=0; y <= (int) yDim; y++) {

      Math2D::Matrix<double> structure_tensor(2,2,0.0);
      Math1D::Vector<double> grad(2);

      for (int xx = std::max(0,x-bound); xx < std::min((int)xDim,x+bound+1); xx++) {

	for (int yy = y-1; yy >= y-bound && yy >= 0; yy--) {
	  if (yy >= 0) {
	    grad[0] = imgrad(xx,yy,0);
	    grad[1] = imgrad(xx,yy,1);
	    Math2D::Matrix<double> temp = outer_product(grad,grad);

	    if (status(xx,yy) != ToInpaint) {
	     
	      double sqr_dist_x = (yy-y) * (yy-y);
	      double sqr_dist_y = (xx-x) * (xx-x);
	      
	      double weight = exp(-2.0*sqr_dist_x / rho) * exp(-2.0*sqr_dist_y / rho);
	      
	      temp *= weight;
	    }
	    else
	      temp *= 0.0;
	    structure_tensor += temp;
	  }
	}

	for (int yy = y; yy < y + bound && yy < (int) yDim; yy++) {
	  if (yy < (int) yDim) {
	    grad[0] = imgrad(xx,yy,0);
	    grad[1] = imgrad(xx,yy,1);
	    Math2D::Matrix<double> temp = outer_product(grad,grad);

	    if (status(xx,yy) != ToInpaint) {

	      double sqr_dist_x = (yy-y) * (yy-y);
	      double sqr_dist_y = (xx-x) * (xx-x);
	      
	      double weight = exp(-2.0*sqr_dist_x / rho) * exp(-2.0*sqr_dist_y / rho);
	      temp *= weight;
	    }
	    else
	      temp *= 0.0;
	    structure_tensor += temp;
	  }
	}
      }
	
      eig2x2(structure_tensor, ev1, ev2);

      horientation(x,y) = atan2(ev2[1],ev2[0]);

      //std::cerr << "orientation " << horientation(x,y) << std::endl;
    }
  }

  vorientation.resize(xDim+1,yDim);
  vorientation.set_constant(-100);
  
  for (int x=0; x < (int) xDim; x++) {
    for (int y=0; y < (int) yDim; y++) {

      Math2D::Matrix<double> structure_tensor(2,2,0.0);
      Math1D::Vector<double> grad(2);

      for (int yy = std::max(0,y-1); yy < std::min((int) yDim, y+bound+1); yy++) {

	for (int xx=x-1; xx >= x-bound; xx--) {
	  if (xx >= 0) {
	    grad[0] = imgrad(xx,yy,0);
	    grad[1] = imgrad(xx,yy,1);
	    structure_tensor += outer_product(grad,grad);
	  }
	}
	for(int xx=x; xx < x + bound; xx++) {
	  if (xx < (int) xDim) {
	    grad[0] = imgrad(xx,yy,0);
	    grad[1] = imgrad(xx,yy,1);
	    structure_tensor += outer_product(grad,grad);
	  }
	}
      }

      eig2x2(structure_tensor, ev1, ev2);

      vorientation(x,y) = atan2(ev2[1],ev2[0]);
    }
  }
}


double curv_weight(const Mesh2D& mesh, const Mesh2DEdgePair& pair, double curv_power=2.0,
		   double override_first = -100.0, double override_second = -100.0) {

  Mesh2DEdge e1 = mesh.edge(pair.first_edge_idx_);
  Mesh2DEdge e2 = mesh.edge(pair.second_edge_idx_);

  uint p2_idx = pair.common_point_idx_;
  uint p1_idx = (e1.from_idx_ == p2_idx) ? e1.to_idx_ : e1.from_idx_;
  uint p3_idx = (e2.from_idx_ == p2_idx) ? e2.to_idx_ : e2.from_idx_;

  double d1x = mesh.point(p2_idx).x_ - mesh.point(p1_idx).x_;
  double d1y = mesh.point(p2_idx).y_ - mesh.point(p1_idx).y_;

  double d2x = mesh.point(p3_idx).x_ - mesh.point(p2_idx).x_;
  double d2y = mesh.point(p3_idx).y_ - mesh.point(p2_idx).y_;

  double angle1 = atan2(d1y,d1x);
  double angle2 = atan2(d2y,d2x);

  //std::cerr << "org_angle1: " << angle1 << ", org_angle2: " << angle2 << std::endl;

  if (override_first != -100.0) {

    double new_angle = override_first;
    double dist = 1e300;

    for (int i=-1; i <= 1; i++) {

      double cand_angle = override_first + i*M_PI;
      double cand_dist = fabs(angle1 - cand_angle);

      //if (cand_dist < dist) {
      if (cand_angle >= -M_PI && cand_angle <= M_PI && cand_dist < dist) {
	dist = cand_dist;
	new_angle = cand_angle;
      }
    }
    angle1 = new_angle;
  }
  if (override_second != -100.0) {

    double new_angle = override_second;
    double dist = 1e300;

    for (int i=-1; i <= 1; i++) {

      double cand_angle = override_second + i*M_PI;
      double cand_dist = fabs(angle2 - cand_angle);

      if (cand_angle >= -M_PI && cand_angle <= M_PI && cand_dist < dist) {
	dist = cand_dist;
	new_angle = cand_angle;
      }
    }
    angle2 = new_angle;
  }

  double diff_angle = fabs(angle2-angle1);
  diff_angle = std::min(diff_angle, 2*M_PI - diff_angle);

  //std::cerr << "angle1: " << angle1 << ", angle2: " << angle2 << ", diff_angle: " << diff_angle << std::endl;

  return pow(diff_angle,curv_power);
}

double lp_inpaint(const Math2D::Matrix<float>& image, const Math2D::Matrix<float>& mask,
		  double lambda, double gamma, double curv_power, uint neighborhood, double energy_offset, std::string solver,
		  Math2D::Matrix<float>& inpainted_image, bool enforce_consistent_boundaries, 
		  bool enforce_regionedge, bool light_constraints, bool legacy) {

  //DEBUG
  //legacy = true;
  //END_DEBUG
  if (legacy)
    std::cerr << "WARNING: legacy mode!!!!!!!!!!!!!!!" << std::endl;

  uint light_factor = (light_constraints) ? 1 : 2;

  inpainted_image = image;

  uint xDim = image.xDim();
  uint yDim = image.yDim();

  Math2D::NamedMatrix<InpaintStatus> status(xDim,yDim,Irrelevant,MAKENAME(status));
  Math2D::NamedMatrix<uint> face_offset(xDim,yDim,MAX_UINT,MAKENAME(face_offset));
  Math2D::NamedMatrix<uint> nFaces(xDim,yDim,MAX_UINT,MAKENAME(nFaces));

  uint nToInpaint = 0;

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {
      if (mask(x,y) < 128.0) {
	status(x,y) = ToInpaint;
	nToInpaint++;
      }
    }
  }

  std::cerr << "to inpaint: " << nToInpaint << " pixels." << std::endl;

  float min_border = 1e36;
  float max_border = -1e36;

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {
      if (status(x,y) == Irrelevant) {
	
	bool border = false;

	if (x > 0 && status(x-1,y) == ToInpaint)
	  border = true;
	if (x+1 < xDim && status(x+1,y) == ToInpaint)
	  border = true;
	if (y > 0 && status(x,y-1) == ToInpaint)
	  border = true;
	if (y+1 < yDim && status(x,y+1) == ToInpaint)
	  border = true;
	if (x > 0 && y > 0 && status(x-1,y-1) == ToInpaint) 
	  border = true;
	if (x > 0 && y+1 < yDim && status(x-1,y+1) == ToInpaint) 
	  border = true;
       	if (x+1 < xDim && y > 0 && status(x+1,y-1) == ToInpaint) 
	  border = true;
	if (x+1 < xDim && y+1 < yDim && status(x+1,y+1) == ToInpaint) 
	  border = true;

	if (border) {
	  status(x,y) = Border;
	  min_border = std::min(min_border, image(x,y));
	  max_border = std::max(max_border, image(x,y));
	}
      }
    }
  }

  std::cerr << "range: [" << min_border << "," << max_border << "]" << std::endl;

  if (max_border == min_border) {

    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x < xDim; x++) {
	if (status(x,y) == ToInpaint) {
	  inpainted_image(x,y) = min_border;
	}
      }
    }
    
    return 0.0;
  }

  Math2D::Matrix<float> horientation; 
  Math2D::Matrix<float> vorientation;
  compute_orientations(image, status, horientation, vorientation);

  Mesh2D mesh;  
  generate_mesh(neighborhood,status,face_offset,nFaces,mesh);

  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  std::cerr << edge_pairs.size() << " edge pairs." << std::endl;

  uint nVars = mesh.nFaces() + 2*edge_pairs.size();
  uint nConstraints = 3*mesh.nEdges();

  const uint consistency_con_offs = nConstraints;
  if (enforce_consistent_boundaries) {  
    nConstraints += light_factor*mesh.nEdges();
  }
  const uint regionedge_constraints_offs = nConstraints;
  if (enforce_regionedge) {
    nConstraints += 2*light_factor*mesh.nEdges();
  }

  Math1D::NamedVector<double> cost(nVars,0.0,MAKENAME(cost));

  uint edge_pair_offset = mesh.nFaces();

  for (uint j=0; j < edge_pairs.size(); j++) {
    
    uint first = edge_pairs[j].first_edge_idx_;
    uint second = edge_pairs[j].second_edge_idx_;
    double weight = 0.5*lambda*(mesh.edge_length(first) + mesh.edge_length(second))
      + gamma * curv_weight(mesh,edge_pairs[j],curv_power);

    cost[edge_pair_offset+2*j] = weight;
    cost[edge_pair_offset+2*j+1] = weight;
  }

  Math1D::NamedVector<double> rhs_lower(nConstraints,0.0,MAKENAME(rhs_lower));
  Math1D::NamedVector<double> rhs_upper(nConstraints,0.0,MAKENAME(rhs_upper));

  Math1D::NamedVector<double> var_lb(nVars,0.0,MAKENAME(var_lb));
  Math1D::NamedVector<double> var_ub(nVars,max_border-min_border,MAKENAME(var_ub));

  for (uint v=edge_pair_offset; v < nVars; v++)
    var_lb[v] = 0.0;

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {

      if (status(x,y) == Border) {

	assert(face_offset(x,y) < nVars);

	for (uint f=face_offset(x,y); f < face_offset(x,y) + nFaces(x,y); f++) {
	 
	  var_lb[f] = image(x,y) - min_border;
	  var_ub[f] = image(x,y) - min_border;

	  assert(var_lb[f] >= 0.0);
	}
      }
    }
  }


  Math1D::NamedVector<EdgeStatus> edge_status(mesh.nEdges(),InteriorEdge,MAKENAME(edge_status));
  for (uint e=0; e < mesh.nEdges(); e++) {
      
    const std::vector<uint>& adjacent_regions = mesh.adjacent_faces(e); 
    
    bool all_border = true;
    
    for (size_t k=0; k < adjacent_regions.size(); k++) {
      
      uint region = adjacent_regions[k];
      
      if (var_lb[region] != var_ub[region])
	all_border = false;
    }

    if (all_border)
      edge_status[e] = BorderEdge;
  }

  if (!legacy) {
    
    for (uint p=0; p < edge_pairs.size(); p++) {
      
      assert(var_lb[edge_pair_offset+2*p] == 0.0);
      assert(var_lb[edge_pair_offset+2*p+1] == 0.0);

      Mesh2DEdgePair cur_pair = edge_pairs[p];
      uint first = cur_pair.first_edge_idx_;
      uint second = cur_pair.second_edge_idx_;
      
      if (edge_status[first] == InteriorEdge && edge_status[second] == InteriorEdge) {
	
	var_ub[edge_pair_offset+2*p]   = max_border - min_border;
	var_ub[edge_pair_offset+2*p+1] = max_border - min_border;
      }
      else {

	uint middle_point = edge_pairs[p].common_point_idx_;

	bool first_forward  = (mesh.edge(first).to_idx_ == middle_point);
	bool second_forward = (mesh.edge(second).to_idx_ == middle_point);
	
	if (edge_status[first] == BorderEdge) {
	
	  double drop = 0.0;
	  
	  for (uint k=0; k < mesh.adjacent_faces(first).size(); k++) {
	    
	    uint idx = mesh.adjacent_faces(first)[k];
	    assert(var_ub[idx] == var_lb[idx]);
	    
	    drop += var_ub[idx] * mesh.match(idx,first); 
	  }
	  
	  if (first_forward)
	    drop *= -1.0;
	  
	  if (drop >= 0.0) {
	    
	    var_ub[edge_pair_offset+2*p]   = drop;
	    var_ub[edge_pair_offset+2*p+1] = 0.0;
	  }
	  else {
	    
	    var_ub[edge_pair_offset+2*p]   = 0.0;
	    var_ub[edge_pair_offset+2*p+1] = - drop;
	  }
	  
	  //std::cerr << "drop: " << drop << std::endl;
	  //std::cerr << "lb: " << var_lb[edge_pair_offset+2*p] << ", ub: " << var_ub[edge_pair_offset+2*p] << std::endl;

	  assert(var_lb[edge_pair_offset+2*p]  <= var_ub[edge_pair_offset+2*p]);
	  assert(var_lb[edge_pair_offset+2*p+1]  <= var_ub[edge_pair_offset+2*p+1]);
	}

	if (edge_status[second] == BorderEdge) {

	  double drop = 0.0;
	  
	  for (uint k=0; k < mesh.adjacent_faces(second).size(); k++) {
	    
	    uint idx = mesh.adjacent_faces(second)[k];
	    assert(var_ub[idx] == var_lb[idx]);
	    
	    drop += var_ub[idx] * mesh.match(idx,second); 
	  }
	  
	  if (second_forward)
	    drop *= -1.0;
	  
	  if (drop >= 0.0) {
	    
	    var_ub[edge_pair_offset+2*p]   = 0.0;
	    var_ub[edge_pair_offset+2*p+1] = std::min(var_ub[edge_pair_offset+2*p+1] , drop);
	  }
	  else {
	    
	    var_ub[edge_pair_offset+2*p]   = std::min(var_ub[edge_pair_offset+2*p], -drop);
	    var_ub[edge_pair_offset+2*p+1] = 0.0;
	  }
	  
	}
      }
    }

    for (uint j=0; j < edge_pairs.size(); j++) {
      
      if (var_ub[edge_pair_offset+2*j] > 0.0 || var_ub[edge_pair_offset+2*j+1] > 0.0) {
      
	uint first  = edge_pairs[j].first_edge_idx_;
	uint second = edge_pairs[j].second_edge_idx_;
	
	Mesh2DEdge e1 = mesh.edge(first);
	Mesh2DEdge e2 = mesh.edge(second);
	
	double override_first  = -100.0;
	double override_second = -100.0;
	
	if (edge_status[first] == BorderEdge) {
	  int x1 = mesh.point(e1.from_idx_).x_;
	  int y1 = mesh.point(e1.from_idx_).y_;
	  int x2 = mesh.point(e1.to_idx_).x_;
	  int y2 = mesh.point(e1.to_idx_).y_;
	  
	  if (x1 == x2)
	    override_first = vorientation(x1,std::min(y1,y2));
	  else {
	    assert(y1 == y2);
	    override_first = horientation(std::min(x1,x2),y1);
	  }
	}
	if (edge_status[second] == BorderEdge) {
	  int x1 = mesh.point(e2.from_idx_).x_;
	  int y1 = mesh.point(e2.from_idx_).y_;
	  int x2 = mesh.point(e2.to_idx_).x_;
	  int y2 = mesh.point(e2.to_idx_).y_;
	  
	  if (x1 == x2)
	    override_second = vorientation(x1,std::min(y1,y2));
	  else {
	    assert(y1 == y2);
	    override_second = horientation(std::min(x1,x2),y1);
	  }
	}
	
	Mesh2DEdgePair rev_pair = edge_pairs[j];
	std::swap(rev_pair.first_edge_idx_,rev_pair.second_edge_idx_);
	
	double cw = std::min(curv_weight(mesh,edge_pairs[j],curv_power,override_first,override_second),
			     curv_weight(mesh,rev_pair,curv_power,override_second,override_first));
	
	double weight = 0.5*lambda*(mesh.edge_length(first) + mesh.edge_length(second))
	  + gamma * cw;
	
	cost[edge_pair_offset+2*j] = weight;
	cost[edge_pair_offset+2*j+1] = weight;
      }
    }
  }  
    
  
  std::cerr << "coding matrix" << std::endl;

  uint nEntries = 2*mesh.nEdges() + 8*edge_pairs.size(); 

  if (enforce_consistent_boundaries) {
    nEntries += light_factor*edge_pairs.size()
      + light_factor*mesh.nEdges(); //we also allocate space for the slack variables used with the convex solver
  }
  if (enforce_regionedge) {

    nEntries += 8*light_factor*edge_pairs.size() 
      + 2*light_factor*mesh.nEdges(); //Note: not exact number
  }
  

  SparseMatrixDescription<double> lp_descr(nEntries, nConstraints, nVars);

  /**** a) code surface continuation constraints *****/

  for (uint j=0; j < mesh.nEdges(); j++) {
    for (std::vector<uint>::const_iterator it = mesh.adjacent_faces(j).begin();
	 it != mesh.adjacent_faces(j).end(); it++) {

      lp_descr.add_entry(j,*it,mesh.match(*it,j));
    }
  }

  for (uint j=0; j < edge_pairs.size(); j++) {

    uint first_edge = edge_pairs[j].first_edge_idx_;
    uint second_edge = edge_pairs[j].second_edge_idx_;

    uint middle_point = edge_pairs[j].common_point_idx_;

    if (mesh.edge(first_edge).to_idx_ == middle_point) {
      lp_descr.add_entry(first_edge,edge_pair_offset+2*j,1);
    }
    else {
      lp_descr.add_entry(first_edge,edge_pair_offset+2*j,-1);
    }

    if (mesh.edge(second_edge).to_idx_ == middle_point) {
      lp_descr.add_entry(second_edge,edge_pair_offset+2*j+1,1);
    }
    else {
      lp_descr.add_entry(second_edge,edge_pair_offset+2*j+1,-1);
    }
  }

  /**** b) code boundary continuation constraints *****/
  uint boundary_con_offset = mesh.nEdges();

  for (uint j=0; j < edge_pairs.size(); j++) {

    uint first_edge = edge_pairs[j].first_edge_idx_;
    uint second_edge = edge_pairs[j].second_edge_idx_;

    uint middle_point = edge_pairs[j].common_point_idx_;

    if (mesh.edge(first_edge).to_idx_ == middle_point) {      
      lp_descr.add_entry(boundary_con_offset + 2*first_edge, edge_pair_offset+2*j, 1);
      lp_descr.add_entry(boundary_con_offset + 2*first_edge+1, edge_pair_offset+2*j+1, -1);
    }
    else {
      lp_descr.add_entry(boundary_con_offset + 2*first_edge+1, edge_pair_offset+2*j, 1);
      lp_descr.add_entry(boundary_con_offset + 2*first_edge, edge_pair_offset+2*j+1, -1);
    }

    if (mesh.edge(second_edge).from_idx_ == middle_point) {
      lp_descr.add_entry(boundary_con_offset + 2*second_edge, edge_pair_offset+2*j, -1);
      lp_descr.add_entry(boundary_con_offset + 2*second_edge+1, edge_pair_offset+2*j+1, 1);
    }
    else {
      lp_descr.add_entry(boundary_con_offset + 2*second_edge+1, edge_pair_offset+2*j, -1);
      lp_descr.add_entry(boundary_con_offset + 2*second_edge, edge_pair_offset+2*j+1, 1);
    }
  }

  if (enforce_consistent_boundaries) {

    //constraints in words: for each oriented edge, the pairs that start with this oriented edge 
    // and the pairs that end in the oppositely oriented edge may not sum to more than 1.0
    // (i.e. they are mutually exclusive)

    for (uint e=0; e < mesh.nEdges(); e++) {

      double bound = max_border - min_border;

      rhs_upper[consistency_con_offs+light_factor*e] =  bound;
      if (!light_constraints)
	rhs_upper[consistency_con_offs+light_factor*e+1] =  bound;
    }

    for (uint j=0; j < edge_pairs.size(); j++) {

      uint middle_point = edge_pairs[j].common_point_idx_;
      
      uint first_edge = edge_pairs[j].first_edge_idx_;
      uint second_edge = edge_pairs[j].second_edge_idx_;


      if (mesh.edge(first_edge).to_idx_ == middle_point) {      
	if (var_ub[edge_pair_offset+2*j] != 0.0)
	  lp_descr.add_entry(consistency_con_offs + light_factor*first_edge, edge_pair_offset+2*j, 1);
	if (var_ub[edge_pair_offset+2*j+1] != 0.0)
	  lp_descr.add_entry(consistency_con_offs + light_factor*first_edge, edge_pair_offset+2*j+1, 1);
      }
      else {
	if (!light_constraints) {
	  lp_descr.add_entry(consistency_con_offs + light_factor*first_edge+1, edge_pair_offset+2*j, 1);
	  lp_descr.add_entry(consistency_con_offs + light_factor*first_edge+1, edge_pair_offset+2*j+1, 1);
	}
      }
      
      if (mesh.edge(second_edge).to_idx_ == middle_point) {
	if (var_ub[edge_pair_offset+2*j] != 0.0)
	  lp_descr.add_entry(consistency_con_offs + light_factor*second_edge, edge_pair_offset+2*j, 1);
	if (var_ub[edge_pair_offset+2*j+1] != 0.0)
	  lp_descr.add_entry(consistency_con_offs + light_factor*second_edge, edge_pair_offset+2*j+1, 1);
      }
      else {
	if (!light_constraints) {
	  lp_descr.add_entry(consistency_con_offs + light_factor*second_edge+1, edge_pair_offset+2*j, 1);
	  lp_descr.add_entry(consistency_con_offs + light_factor*second_edge+1, edge_pair_offset+2*j+1, 1);
	}
      }
    }
  }

  if (enforce_regionedge) {

    uint rowoff = regionedge_constraints_offs;
    
    for (uint edge=0; edge < mesh.nEdges(); edge++) {
      
      //Get the two adjacent faces
      if (mesh.adjacent_faces(edge).size() != 2)
	{
	  //One of the edges is at the border of the image
	  
	  continue;
	}
      
      uint x1 = mesh.adjacent_faces(edge)[0];
      uint x2 = mesh.adjacent_faces(edge)[1];

      double bound = max_border - min_border;
      
      lp_descr.add_entry(rowoff+2*light_factor*edge  , x1, 1);
      lp_descr.add_entry(rowoff+2*light_factor*edge  , x2, 1);
      lp_descr.add_entry(rowoff+2*light_factor*edge+1, x1, -1);
      lp_descr.add_entry(rowoff+2*light_factor*edge+1, x2, -1);
      if (!light_constraints) {
	lp_descr.add_entry(rowoff+2*light_factor*edge+2, x1, 1);
	lp_descr.add_entry(rowoff+2*light_factor*edge+2, x2, 1);
	lp_descr.add_entry(rowoff+2*light_factor*edge+3, x1, -1);
	lp_descr.add_entry(rowoff+2*light_factor*edge+3, x2, -1);
      }

      rhs_lower[rowoff+2*light_factor*edge]   = 0;
      rhs_upper[rowoff+2*light_factor*edge]   = 2*bound;
      rhs_lower[rowoff+2*light_factor*edge+1] = -2*bound;
      rhs_upper[rowoff+2*light_factor*edge+1] = 0;
      if (!light_constraints) {
	rhs_lower[rowoff+2*light_factor*edge+2] = 0;
	rhs_upper[rowoff+2*light_factor*edge+2] = 2*bound;
	rhs_lower[rowoff+2*light_factor*edge+3] = -2*bound;
	rhs_upper[rowoff+2*light_factor*edge+3] = 0;
      }
    }
    
    for (uint j=0; j < edge_pairs.size(); j++) {
      uint first = edge_pairs[j].first_edge_idx_;
      uint second = edge_pairs[j].second_edge_idx_;	
      
      uint y = edge_pair_offset + 2*j;
      
      uint edge = first;
      if (mesh.adjacent_faces(edge).size() == 2) {
	lp_descr.add_entry(rowoff+2*light_factor*edge   ,y  , 1);
	lp_descr.add_entry(rowoff+2*light_factor*edge+1 ,y  , 1);
	if (!light_constraints) {
	  lp_descr.add_entry(rowoff+2*light_factor*edge+2 ,y+1, 1);
	  lp_descr.add_entry(rowoff+2*light_factor*edge+3 ,y+1, 1);
	}
      }
      
      edge = second;
      if (mesh.adjacent_faces(edge).size() == 2) {
	lp_descr.add_entry(rowoff+2*light_factor*edge   ,y+1, 1);
	lp_descr.add_entry(rowoff+2*light_factor*edge+1 ,y+1, 1);
	if (!light_constraints) {
	  lp_descr.add_entry(rowoff+2*light_factor*edge+2 ,y  , 1);
	  lp_descr.add_entry(rowoff+2*light_factor*edge+3 ,y  , 1);
	}
      }
    }
  }

  std::cerr << "matrix coded" << std::endl;

  Math1D::NamedVector<uint> row_start(nConstraints+1,MAKENAME(row_start));
  lp_descr.sort_by_row(row_start);

  bool solver_known = false;

  const double* lp_solution = 0;

  int error = 0;

#ifdef HAS_GUROBI

  Math1D::Vector<double> gurobi_solution;

  GRBenv   *grb_env   = 0;
  GRBmodel *grb_model = 0;

  if (solver == "gurobi") {

    solver_known = true;

    /* Create environment */

    error = GRBloadenv(&grb_env,NULL);
    GRBsetintparam(grb_env, GRB_INT_PAR_METHOD, GRB_METHOD_BARRIER);
    GRBsetdblparam(grb_env, "BarConvTol", 1e-10);
    GRBsetintparam(grb_env, "Crossover", 0);
    GRBsetintparam(grb_env, "CrossoverBasis", 1);
    GRBsetintparam(grb_env, "Presolve", 1);
    GRBsetintparam(grb_env, "PrePasses", 2);

    assert (!error && grb_env != NULL);

    /* Create an empty model */

    error = GRBnewmodel(grb_env, &grb_model, "curv-inpaint-lp", 0, NULL, NULL, NULL, NULL, NULL);
    assert(!error);
    
    Storage1D<char> vtype(nVars,GRB_CONTINUOUS);
    
    error = GRBaddvars(grb_model,nVars,0,NULL,NULL,NULL,cost.direct_access(),var_lb.direct_access(),
		       var_ub.direct_access(),vtype.direct_access(),NULL);
    assert(!error);
    
    error = GRBupdatemodel(grb_model);
    assert(!error);
    
    for (uint c=0; c < row_start.size()-1; c++) {
      
      if (rhs_lower[c] == rhs_upper[c]) {
	error = GRBaddconstr(grb_model, row_start[c+1]-row_start[c], ((int*) lp_descr.col_indices()) + row_start[c], 
			     lp_descr.value() + row_start[c], GRB_EQUAL, rhs_lower[c], NULL);
      }
      else {
	error = GRBaddrangeconstr(grb_model, row_start[c+1]-row_start[c], ((int*) lp_descr.col_indices()) + row_start[c], 
				  lp_descr.value() + row_start[c], rhs_lower[c], rhs_upper[c], NULL);
      }
      
      assert(!error);
    }

    /* Optimize model */
    error = GRBoptimize(grb_model);
    assert(!error);

    gurobi_solution.resize(nVars);
    
    for (uint v=0; v < nVars; v++)
      GRBgetdblattrelement(grb_model,"X",v, gurobi_solution.direct_access()+v);

    lp_solution = gurobi_solution.direct_access();
  }
#endif

#ifdef HAS_CPLEX

  Math1D::Vector<double> cplex_solution;
  CPXENVptr     cp_env = NULL;
  CPXLPptr      cp_lp = NULL;

  if (solver == "cplex") {

    solver_known = true;
  
    int cpx_status = 0;
    
    /* Initialize the CPLEX environment */
    
    cp_env = CPXopenCPLEX (&cpx_status);
    CPXsetintparam(cp_env, CPX_PARAM_BARCROSSALG, -1);
    //CPXsetintparam(cp_env, CPX_PARAM_STARTALG, CPX_ALG_BARRIER);
    //CPXsetintparam(cp_env, CPX_PARAM_MIPDISPLAY, 4);
    //CPXsetintparam(cp_env, CPX_PARAM_PREIND, CPX_OFF);
    //CPXsetintparam(cp_env, CPX_PARAM_PREPASS, 0);
    
    //CPXsetintparam(cp_env, CPX_PARAM_LPMETHOD, CPX_ALG_DUAL);
    
    /* If an error occurs, the status value indicates the reason for
       failure.  A call to CPXgeterrorstring will produce the text of
       the error message.  Note that CPXopenCPLEX produces no output,
       so the only way to see the cause of the error is to use
       CPXgeterrorstring.  For other CPLEX routines, the errors will
       be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON.  */
    
    if ( cp_env == NULL ) {
      char  errmsg[1024];
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (cp_env, cpx_status, errmsg);
      fprintf (stderr, "%s", errmsg);
      exit(1);
    }
    
    /* Turn on output to the screen */
    
    cpx_status = CPXsetintparam (cp_env, CPX_PARAM_SCRIND, CPX_ON);
    if ( cpx_status ) {
      fprintf (stderr,
	       "Failure to turn on screen indicator, error %d.\n", cpx_status);
      exit(1);
    }
    
    //necessary when using own cut generator (or heuristic??) with CPLEX
    //cpx_status = CPXsetintparam (env, CPX_PARAM_PREIND, CPX_OFF);
    
    //set problem data
    
    cp_lp = CPXcreateprob (cp_env, &cpx_status, "curv-inpaint-ilp");
    
    /* A returned pointer of NULL may mean that not enough memory
       was available or there was some other problem.  In the case of
       failure, an error message will have been written to the error
       channel from inside CPLEX.  In this example, the setting of
       the parameter CPX_PARAM_SCRIND causes the error message to
       appear on stdout.  */
    
    if ( cp_lp == NULL ) {
      fprintf (stderr, "Failed to create LP.\n");
      exit(1);
    }
    
    /* Now copy the problem data into the lp */
    
    char* row_sense = new char[nConstraints];
    for (uint c=0; c < nConstraints; c++) {
      
      if (rhs_lower[c] == rhs_upper[c]) {
	row_sense[c] = 'E';
      }
      else {
	assert(rhs_lower[c] == -MAX_DOUBLE || rhs_lower[c] == 0.0);
	row_sense[c] = 'L';
      }
    }
    
    int* row_count = new int[nConstraints];
    for (uint c=0; c < nConstraints; c++)
      row_count[c] = row_start[c+1] - row_start[c];
    
    cpx_status = CPXnewcols (cp_env, cp_lp, nVars, cost.direct_access(), var_lb.direct_access(), 
			     var_ub.direct_access(), NULL, NULL);
    if ( cpx_status )  
      exit(1);
    
    CPXaddrows(cp_env, cp_lp, 0, nConstraints, lp_descr.nEntries(), rhs_upper.direct_access(), row_sense, 
	       (int*) row_start.direct_access(), (int*) lp_descr.col_indices(), lp_descr.value(),
	       NULL, NULL);

    delete[] row_sense;
    delete[] row_count;

    std::cerr << "calling optimize" << std::endl;

    cpx_status = CPXbaropt(cp_env,cp_lp);
    //cpx_status = CPXlpopt(cp_env,cp_lp);

    if ( cpx_status ) {
      fprintf (stderr, "Failed to optimize MIP.\n");
      exit(1);
    }
  
    cplex_solution.resize_dirty(nVars);

    CPXsolution (cp_env, cp_lp, NULL, NULL, cplex_solution.direct_access(), NULL, NULL, NULL);
    
    lp_solution = cplex_solution.direct_access();
  }

#endif

#ifdef HAS_XPRESS

  Math1D::Vector<double> xpress_solution;

  XPRSprob xp_prob;

  if (solver == "xpress") {

    solver_known = true;

    int nReturn;

    nReturn=XPRSinit("/opt/xpressmp/");
    
    std::cerr << "nReturn: " << nReturn << std::endl;
    if (nReturn != 0) {
      
      char msg[512];
      XPRSgetlicerrmsg(msg,512);
      
      std::cerr << "error message: " << msg << std::endl;
    }
    
    assert(nReturn == 0);

    char banner[256];
    
    XPRSgetbanner(banner); printf("banner: %s \n",banner);
    
    nReturn=XPRScreateprob(&xp_prob);

    //XPRSsetintcontrol(xp_prob,XPRS_PRESOLVE, 0);
    XPRSsetintcontrol(xp_prob,XPRS_CROSSOVER,0);
    
    std::cerr << "nReturn: " << nReturn << std::endl;
  
    SparseMatrixDescription<double> lp_copy(lp_descr);
    
    Math1D::Vector<uint> col_start(nVars+1);
    lp_copy.sort_by_column(col_start);
    
    char* row_sense = new char[nConstraints];
    double* row_range = new double[nConstraints];
    
    for (uint c=0; c < nConstraints; c++) {
      
      if (rhs_lower[c] == rhs_upper[c]) {
	row_sense[c] = 'E';
	row_range[c] = 0.0;
      }
      else {
	row_sense[c] = 'R';
	row_range[c] = rhs_upper[c] - rhs_lower[c];
      }
    }
    
    nReturn = XPRSloadlp(xp_prob, "curv-inpaint-ilp", col_start.size()-1, row_start.size()-1, row_sense,
			 rhs_upper.direct_access(), row_range, cost.direct_access(), 
			 (int*) col_start.direct_access(), NULL, (int*) lp_copy.row_indices(), lp_copy.value(),
			 var_lb.direct_access(), var_ub.direct_access());
    
    delete[] row_sense;
    delete[] row_range;

    lp_copy.reset(0);
    
    // Math1D::Vector<char> col_type(nVars,'C');
    // Math1D::Vector<int> idx(nVars);
    // for (uint v=0; v < nVars; v++)
    //   idx[v] = v;
    
    // nReturn = XPRSchgcoltype(xp_prob, nVars, idx.direct_access(), col_type.direct_access());

    XPRSlpoptimize(xp_prob,"b");
    
    xpress_solution.resize_dirty(nVars);

    XPRSgetlpsol(xp_prob, xpress_solution.direct_access(), 0, 0, 0); 

    lp_solution = xpress_solution.direct_access();
  }
#endif

  if (!solver_known && solver != "clp") {

    std::cerr << "WARNING: solver\"" << solver << "\" unknown. Using Clp instead" << std::endl;
    solver = "clp";
  }

  ClpSimplex lpSolver;

  if (solver == "clp") {

    std::clock_t tStartCLP,tEndCLP;

    CoinPackedMatrix coinMatrix(false,(int*) lp_descr.row_indices(),(int*) lp_descr.col_indices(),
				lp_descr.value(),lp_descr.nEntries());
    
    lpSolver.loadProblem (coinMatrix, var_lb.direct_access(), var_ub.direct_access(),   
			  cost.direct_access(), rhs_lower.direct_access(), rhs_upper.direct_access());

    coinMatrix.cleanMatrix();
    
    tStartCLP = std::clock();
    
    //lpSolver.dual();

    ClpSolve solve_options;
    solve_options.setSolveType(ClpSolve::useDual);
    //solve_options.setSolveType(ClpSolve::useBarrier);
    solve_options.setPresolveType(ClpSolve::presolveNumber,5);
    lpSolver.initialSolve(solve_options);
  
    error = 1 - lpSolver.isProvenOptimal();

    tEndCLP = std::clock();

    if (error != 0)
      std::cerr << "!!!!!!!!!!!!!!LP-solver failed!!!!!!!!!!!!!!!!!!!" << std::endl;

    std::cerr << "CLP-time: " << diff_seconds(tEndCLP,tStartCLP) << " seconds. " << std::endl;

    lp_solution = lpSolver.primalColumnSolution();
  }



  //check for possible crossing points
  Storage1D<std::vector<uint> > point_pairs(mesh.nPoints());
  
  for (uint j=0; j < edge_pairs.size(); j++) {
    
    uint common_point = edge_pairs[j].common_point_idx_;
    point_pairs[common_point].push_back(j);
  }
  
  while (!legacy) {

    uint nAddedConstraints = 0;

    for (uint p=0; p < mesh.nPoints(); p++) {
      
      double sum = 0.0;
      for (std::vector<uint>::iterator it = point_pairs[p].begin(); it != point_pairs[p].end(); it++)
	sum += lp_solution[edge_pair_offset +  2*(*it)] + lp_solution[edge_pair_offset +  2*(*it) + 1];
      
      if (sum > max_border-min_border) {
	std::cerr << "checking if  point #" << p << " is a crossing point" << std::endl;
	  
	bool message_output = false;

	for (uint k1=0; k1 < point_pairs[p].size()-1; k1++) {
	  
	  uint pair1 = point_pairs[p][k1];
	  
	  for (uint k2=k1+1; k2 < point_pairs[p].size(); k2++) {
	    
	    uint pair2 = point_pairs[p][k2];
	    
	    double sum = lp_solution[edge_pair_offset +  2*pair1] + lp_solution[edge_pair_offset +  2*pair1 + 1]
	      + lp_solution[edge_pair_offset +  2*pair2] + lp_solution[edge_pair_offset +  2*pair2 + 1];
	    
	    if (sum >= 1.01*(max_border-min_border)) {

	      if (!message_output)
		std::cerr << "---> yes." << std::endl;
	      message_output = true;

	      int cols[4];
	      double coeffs[4] = {1.0,1.0,1.0,1.0};
	      cols[0] = edge_pair_offset +  2*pair1;
	      cols[1] = edge_pair_offset +  2*pair1 + 1;
	      cols[2] = edge_pair_offset +  2*pair2;
	      cols[3] = edge_pair_offset +  2*pair2 + 1;
	      
	      //note: adding constraints separately is VERY inefficient
	      if (solver == "clp")
		lpSolver.addRow(4, cols, coeffs, 0.0, max_border-min_border);
#ifdef HAS_GUROBI
	      if (solver == "gurobi")
		GRBaddconstr(grb_model,4,cols,coeffs,'L',max_border-min_border,NULL);
#endif
#ifdef HAS_XPRESS
	      if (solver == "xpress") {
		
		double new_rhs[1] = {max_border-min_border};
		double new_range[1] = {0.0};
		int new_start[2] = {0,4};
		XPRSaddrows(xp_prob, 1, 4, "L", new_rhs,new_range,new_start,cols,coeffs);
	      }
#endif
#ifdef HAS_CPLEX
	      if (solver == "cplex") {
		
		double new_rhs[1] = {max_border-min_border};
		int new_start[2] = {0,4};

		CPXaddrows(cp_env, cp_lp, 0, 1, 4, new_rhs, "L", new_start, cols, coeffs,  NULL, NULL);		
	      }
#endif
	      nAddedConstraints++;
	    }
	  }
	}
      }
    }

    if (nAddedConstraints == 0)
      break;

    if (solver == "clp") {
      lpSolver.dual();
      lp_solution = lpSolver.primalColumnSolution();
    }
#ifdef HAS_GUROBI
    if (solver == "gurobi") {
      
      int error = GRBupdatemodel(grb_model);
      error = GRBoptimize(grb_model);
      assert(!error);

      for (uint v=0; v < nVars; v++)
	GRBgetdblattrelement(grb_model,"X",v, gurobi_solution.direct_access()+v);
    }
#endif
#ifdef HAS_XPRESS
    if (solver == "xpress") {

      XPRSlpoptimize(xp_prob,"b");      
      XPRSgetlpsol(xp_prob, xpress_solution.direct_access(), 0, 0, 0); 

      lp_solution = xpress_solution.direct_access();
    }
#endif
#ifdef HAS_CPLEX
    if (solver == "cplex") {

      int cpx_status = CPXbaropt(cp_env,cp_lp);
      //int cpx_status = CPXlpopt(cp_env,cp_lp);

      if ( cpx_status ) {
	fprintf (stderr, "Failed to optimize MIP.\n");
	exit(1);
      }
  
      CPXsolution (cp_env, cp_lp, NULL, NULL, cplex_solution.direct_access(), NULL, NULL, NULL);
    
      lp_solution = cplex_solution.direct_access();
    }
#endif
  }

  uint nFrac=0;
    
  for (uint i=0; i < nVars; i++) {
      
    double val = lp_solution[i];
    
    if (fabs(round(val) - val) > 0.02)
      nFrac++;
  }
  std::cerr << nFrac << " fractional variables." << std::endl;
  
  double energy = energy_offset;
  for (uint i=0; i < nVars; i++)
    energy += cost[i] * lp_solution[i];
  
  std::cerr << "energy: " << energy << std::endl;
  
  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {

      //if (face_offset(x,y) != MAX_UINT) {
      if (status(x,y) < Irrelevant) {

	double sum = min_border;

	for (uint f=face_offset(x,y); f < face_offset(x,y) + nFaces(x,y); f++)
	  sum += lp_solution[f] * mesh.convex_area(f);

	assert (status(x,y) != Border || fabs(sum - image(x,y) < 0.01));
	  

	inpainted_image(x,y) = sum;
      }
    }
  }

  /**** clean-up ****/
#ifdef HAS_GUROBI

  if (grb_model != 0)
    GRBfreemodel(grb_model);
  if (grb_env != 0)
    GRBfreeenv(grb_env);
#endif
#ifdef HAS_XPRESS
  if (solver == "xpress")
    XPRSdestroyprob(xp_prob);
#endif
#ifdef HAS_CPLEX
  if (cp_lp != NULL)
    CPXfreeprob (cp_env, &cp_lp);
  if (cp_env != NULL)
    CPXcloseCPLEX (&cp_env);
#endif

  return energy;

}


double lp_inpaint_hybrid(const Math2D::Matrix<float>& image, const Math2D::Matrix<float>& mask,
			 double lambda, double gamma, double curv_power, uint neighborhood, uint nBin, double energy_offset, std::string solver,
			 Math2D::Matrix<float>& inpainted_image, bool enforce_consistent_boundaries, 
			 bool enforce_regionedge, bool enforce_level_consistency,  bool light_constraints, bool legacy) {

  uint light_factor = (light_constraints) ? 1 : 2;

  inpainted_image = image;
  Math2D::Matrix<float> level_image = image;

  uint xDim = image.xDim();
  uint yDim = image.yDim();

  Math2D::NamedMatrix<InpaintStatus> status(xDim,yDim,Irrelevant,MAKENAME(status));
  Math2D::NamedMatrix<uint> face_offset(xDim,yDim,MAX_UINT,MAKENAME(face_offset));
  Math2D::NamedMatrix<uint> nFaces(xDim,yDim,MAX_UINT,MAKENAME(nFaces));

  uint nToInpaint = 0;

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {
      if (mask(x,y) < 128.0) {
	status(x,y) = ToInpaint;
	nToInpaint++;
      }
    }
  }

  std::cerr << "to inpaint: " << nToInpaint << " pixels." << std::endl;

  float min_border = 1e36;
  float max_border = -1e36;
  float im_min_border = 1e36;
  float im_max_border = -1e36;

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {
      if (status(x,y) == Irrelevant) {

	bool border = false;

	if (x > 0 && status(x-1,y) == ToInpaint)
	  border = true;
	if (x+1 < xDim && status(x+1,y) == ToInpaint)
	  border = true;
	if (y > 0 && status(x,y-1) == ToInpaint)
	  border = true;
	if (y+1 < yDim && status(x,y+1) == ToInpaint)
	  border = true;

	if (x > 0 && y > 0 && status(x-1,y-1) == ToInpaint) 
	  border = true;
	if (x > 0 && y+1 < yDim && status(x-1,y+1) == ToInpaint) 
	  border = true;
	if (x+1 < xDim && y > 0 && status(x+1,y-1) == ToInpaint) 
	  border = true;
	if (x+1 < xDim && y+1 < yDim && status(x+1,y+1) == ToInpaint) 
	  border = true;

	if (border) {
	  status(x,y) = Border;
	  im_min_border = std::min(im_min_border, image(x,y));
	  im_max_border = std::max(im_max_border, image(x,y));
	}
      }
    }
  }

  int range = ceil(im_max_border - im_min_border);
  while ((range % nBin) != 0)
    range++;

  float bin_size = ((float) range)/nBin;


  Math2D::Matrix<float> horientation; 
  Math2D::Matrix<float> vorientation;
  compute_orientations(image, status, horientation, vorientation);

  Mesh2D mesh;  
  generate_mesh(neighborhood,status,face_offset,nFaces,mesh);

  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  std::cerr << edge_pairs.size() << " edge pairs." << std::endl;

  uint nVarsPerLayer = mesh.nFaces() + 2*edge_pairs.size();
  uint nConstraintsPerLayer = 3*mesh.nEdges();

  const uint consistency_con_offs = nConstraintsPerLayer;
  if (enforce_consistent_boundaries) {  
    nConstraintsPerLayer += light_factor*mesh.nEdges();
  }
  const uint regionedge_constraints_offs = nConstraintsPerLayer;
  if (enforce_regionedge) {
    nConstraintsPerLayer += 2*light_factor*mesh.nEdges();
  }
		
  uint nVars = nBin*nVarsPerLayer;
  uint nConstraints = nBin*nConstraintsPerLayer;
  if (enforce_level_consistency)
    nConstraints += (nBin-1)*nVarsPerLayer;

  Math1D::NamedVector<double> cost(nVars,0.0,MAKENAME(cost));

  uint edge_pair_offset = mesh.nFaces();

  for (uint j=0; j < edge_pairs.size(); j++) {

    uint first = edge_pairs[j].first_edge_idx_;
    uint second = edge_pairs[j].second_edge_idx_;
    double weight = 0.5*lambda*(mesh.edge_length(first) + mesh.edge_length(second))
      + gamma * curv_weight(mesh,edge_pairs[j],curv_power);
    
    for (uint kk=0; kk < nBin; kk++){
      cost[kk*nVarsPerLayer+edge_pair_offset+2*j] = weight;
      cost[kk*nVarsPerLayer+edge_pair_offset+2*j+1] = weight;
    }
  }
  Math1D::NamedVector<double> rhs_lower(nConstraints,0.0,MAKENAME(rhs_lower));
  Math1D::NamedVector<double> rhs_upper(nConstraints,0.0,MAKENAME(rhs_upper));

  //the initialization should be more specific
  Math1D::NamedVector<double> var_lb(nVars,0.0,MAKENAME(var_lb));
  Math1D::NamedVector<double> var_ub(nVars,bin_size,MAKENAME(var_ub));

  uint nEntriesPerLayer = 2*mesh.nEdges() + 6*edge_pairs.size(); 

  if (enforce_consistent_boundaries) {
    nEntriesPerLayer += light_factor*(2*edge_pairs.size()
			      + mesh.nEdges()); //we also allocate space for the slack variables used with the convex solver
  }
  if (enforce_regionedge) {
    nEntriesPerLayer += 8*light_factor*edge_pairs.size()
      + light_factor*2*mesh.nEdges();//we also allocate space for the slack variables used with the convex solver
  }

  uint nEntries = nBin*nEntriesPerLayer;
  if (enforce_level_consistency) {
    nEntries += (nBin-1)*mesh.nFaces();
  }

  std::cerr << "nEntries: " << nEntries << std::endl;
  SparseMatrixDescription<double> lp_descr(nEntries, nConstraints, nVars);

  uint bin_var_offset = 0;
  uint bin_cons_offset = 0;

  std::cerr << "bin_size: " << bin_size << std::endl;

  for (uint kk=0; kk < nBin; kk++){

    std::cerr << "loop K : " << kk << std::endl;

    bin_var_offset = kk*nVarsPerLayer;
    bin_cons_offset = kk*nConstraintsPerLayer;

    double bin_ub = im_min_border + (kk+1)*bin_size;
    double bin_lb = im_min_border + (kk)*bin_size;

    min_border = 0.0;
    max_border = -1e100;

    //setting up quantized images for level line recovery 
    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x < xDim; x++) {
	if (image(x,y) <= bin_lb)
	  level_image(x,y) = 0.0;
	else if (image(x,y) <= bin_ub && image(x,y) > bin_lb)
	  level_image(x,y) = image(x,y) - bin_lb;
	else { 
	  
 	  if (kk < nBin - 1)
 	    level_image(x,y) = bin_size;
 	  else
 	    level_image(x,y) = im_max_border - bin_lb; 
	}
	if (status(x,y) == Border) {
	  //min_border = std::min(min_border, level_image(x,y));
	  max_border = std::max(max_border, level_image(x,y));
	}	
      }
    }
    assert(min_border >= 0);
    //std::cerr << "max_border: " << max_border << std::endl;

    for (uint v=0; v < nVarsPerLayer; v++){
      var_lb[bin_var_offset+v] = min_border;
      var_ub[bin_var_offset+v] = max_border;
    } 

    for (uint v=edge_pair_offset; v < nVarsPerLayer; v++)
      var_lb[bin_var_offset+v] = 0.0;

    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x < xDim; x++) {

	if (status(x,y) == Border) {

	  assert(face_offset(x,y) < nVarsPerLayer);

	  for (uint f=face_offset(x,y); f < face_offset(x,y) + nFaces(x,y); f++) {

	    var_lb[bin_var_offset+f] = level_image(x,y);
	    var_ub[bin_var_offset+f] = level_image(x,y);
	  }
	}
      }
    }


    Math1D::NamedVector<EdgeStatus> edge_status(mesh.nEdges(),InteriorEdge,MAKENAME(edge_status));
    for (uint e=0; e < mesh.nEdges(); e++) {

      const std::vector<uint>& adjacent_regions = mesh.adjacent_faces(e); 

      bool all_border = true;

      for (size_t k=0; k < adjacent_regions.size(); k++) {

	uint region = adjacent_regions[k];

	if (var_lb[bin_var_offset+region] != var_ub[bin_var_offset+region])
	  all_border = false;
      }

      if (all_border)
	edge_status[e] = BorderEdge;
    }

    if (!legacy) {

      for (uint p=0; p < edge_pairs.size(); p++) {

	assert(var_lb[bin_var_offset+edge_pair_offset+2*p] == 0.0);
	assert(var_lb[bin_var_offset+edge_pair_offset+2*p+1] == 0.0);

	Mesh2DEdgePair cur_pair = edge_pairs[p];
	uint first = cur_pair.first_edge_idx_;
	uint second = cur_pair.second_edge_idx_;

	if (edge_status[first] == InteriorEdge && edge_status[second] == InteriorEdge) {

	  var_ub[bin_var_offset+edge_pair_offset+2*p]   = max_border - min_border;
	  var_ub[bin_var_offset+edge_pair_offset+2*p+1] = max_border - min_border;
	}
	else {

	  uint middle_point = edge_pairs[p].common_point_idx_;

	  bool first_forward  = (mesh.edge(first).to_idx_ == middle_point);
	  bool second_forward = (mesh.edge(second).to_idx_ == middle_point);

	  if (edge_status[first] == BorderEdge) {

	    double drop = 0.0;

	    for (uint k=0; k < mesh.adjacent_faces(first).size(); k++) {

	      uint idx = mesh.adjacent_faces(first)[k];
	      assert(var_ub[bin_var_offset+idx] == var_lb[bin_var_offset+idx]);

	      drop += var_ub[bin_var_offset+idx] * mesh.match(idx,first); 
	    }

	    if (first_forward)
	      drop *= -1.0;

	    if (drop >= 0.0) {

	      var_ub[bin_var_offset+edge_pair_offset+2*p]   = drop;
	      var_ub[bin_var_offset+edge_pair_offset+2*p+1] = 0.0;
	    }
	    else {

	      var_ub[bin_var_offset+edge_pair_offset+2*p]   = 0.0;
	      var_ub[bin_var_offset+edge_pair_offset+2*p+1] = - drop;
	    }

	    //std::cerr << "drop: " << drop << std::endl;
	    //std::cerr << "lb: " << var_lb[edge_pair_offset+2*p] << ", ub: " << var_ub[edge_pair_offset+2*p] << std::endl;

	    assert(var_lb[bin_var_offset+edge_pair_offset+2*p]  <= var_ub[bin_var_offset+edge_pair_offset+2*p]);
	    assert(var_lb[bin_var_offset+edge_pair_offset+2*p+1]  <= var_ub[bin_var_offset+edge_pair_offset+2*p+1]);
	  }

	  if (edge_status[second] == BorderEdge) {

	    double drop = 0.0;

	    for (uint k=0; k < mesh.adjacent_faces(second).size(); k++) {

	      uint idx = mesh.adjacent_faces(second)[k];
	      assert(var_ub[bin_var_offset+idx] == var_lb[bin_var_offset+idx]);

	      drop += var_ub[bin_var_offset+idx] * mesh.match(idx,second); 
	    }

	    if (second_forward)
	      drop *= -1.0;

	    if (drop >= 0.0) {

	      var_ub[bin_var_offset+edge_pair_offset+2*p]   = 0.0;
	      var_ub[bin_var_offset+edge_pair_offset+2*p+1] = std::min(var_ub[bin_var_offset+edge_pair_offset+2*p+1] , drop);
	    }
	    else {

	      var_ub[bin_var_offset+edge_pair_offset+2*p]   = std::min(var_ub[bin_var_offset+edge_pair_offset+2*p], -drop);
	      var_ub[bin_var_offset+edge_pair_offset+2*p+1] = 0.0;
	    }

	  }
	}
      }

      for (uint j=0; j < edge_pairs.size(); j++) {

	if (var_ub[bin_var_offset+edge_pair_offset+2*j] > 0.0 || var_ub[bin_var_offset+edge_pair_offset+2*j+1] > 0.0) {

	  uint first  = edge_pairs[j].first_edge_idx_;
	  uint second = edge_pairs[j].second_edge_idx_;

	  Mesh2DEdge e1 = mesh.edge(first);
	  Mesh2DEdge e2 = mesh.edge(second);

	  double override_first  = -100.0;
	  double override_second = -100.0;

	  if (edge_status[first] == BorderEdge) {
	    int x1 = mesh.point(e1.from_idx_).x_;
	    int y1 = mesh.point(e1.from_idx_).y_;
	    int x2 = mesh.point(e1.to_idx_).x_;
	    int y2 = mesh.point(e1.to_idx_).y_;

	    if (x1 == x2)
	      override_first = vorientation(x1,std::min(y1,y2));
	    else {
	      assert(y1 == y2);
	      override_first = horientation(std::min(x1,x2),y1);
	    }
	  }
	  if (edge_status[second] == BorderEdge) {
	    int x1 = mesh.point(e2.from_idx_).x_;
	    int y1 = mesh.point(e2.from_idx_).y_;
	    int x2 = mesh.point(e2.to_idx_).x_;
	    int y2 = mesh.point(e2.to_idx_).y_;

	    if (x1 == x2)
	      override_second = vorientation(x1,std::min(y1,y2));
	    else {
	      assert(y1 == y2);
	      override_second = horientation(std::min(x1,x2),y1);
	    }
	  }

	  Mesh2DEdgePair rev_pair = edge_pairs[j];
	  std::swap(rev_pair.first_edge_idx_,rev_pair.second_edge_idx_);

	  double cw = std::min(curv_weight(mesh,edge_pairs[j],curv_power,override_first,override_second),
			       curv_weight(mesh,rev_pair,curv_power,override_second,override_first));

	  double weight = 0.5*lambda*(mesh.edge_length(first) + mesh.edge_length(second))
	    + gamma * cw;

	  cost[bin_var_offset+edge_pair_offset+2*j] = weight;
	  cost[bin_var_offset+edge_pair_offset+2*j+1] = weight;
	}
      }
    }  

    std::cerr << "coding matrix" << std::endl;

    /**** a) code surface continuation constraints *****/

    for (uint j=0; j < mesh.nEdges(); j++) {
      for (std::vector<uint>::const_iterator it = mesh.adjacent_faces(j).begin();
	   it != mesh.adjacent_faces(j).end(); it++) {

	lp_descr.add_entry(bin_cons_offset+j,bin_var_offset+*it,mesh.match(*it,j));
      }
    }

    for (uint j=0; j < edge_pairs.size(); j++) {

      uint first_edge = edge_pairs[j].first_edge_idx_;
      uint second_edge = edge_pairs[j].second_edge_idx_;

      uint middle_point = edge_pairs[j].common_point_idx_;

      if (mesh.edge(first_edge).to_idx_ == middle_point) {
	lp_descr.add_entry(bin_cons_offset+first_edge,bin_var_offset+edge_pair_offset+2*j,1);
      }
      else {
	lp_descr.add_entry(bin_cons_offset+first_edge,bin_var_offset+edge_pair_offset+2*j,-1);
      }

      if (mesh.edge(second_edge).to_idx_ == middle_point) {
	lp_descr.add_entry(bin_cons_offset+second_edge,bin_var_offset+edge_pair_offset+2*j+1,1);
      }
      else {
	lp_descr.add_entry(bin_cons_offset+second_edge,bin_var_offset+edge_pair_offset+2*j+1,-1);
      }
    }

    /**** b) code boundary continuation constraints *****/
    uint boundary_con_offset = mesh.nEdges();

    for (uint j=0; j < edge_pairs.size(); j++) {

      uint first_edge = edge_pairs[j].first_edge_idx_;
      uint second_edge = edge_pairs[j].second_edge_idx_;
      
      uint middle_point = edge_pairs[j].common_point_idx_;

      if (mesh.edge(first_edge).to_idx_ == middle_point) {      
	lp_descr.add_entry(bin_cons_offset+boundary_con_offset + 2*first_edge,   bin_var_offset+edge_pair_offset+2*j, 1);
	lp_descr.add_entry(bin_cons_offset+boundary_con_offset + 2*first_edge+1, bin_var_offset+edge_pair_offset+2*j+1, -1);
      }
      else {
	lp_descr.add_entry(bin_cons_offset+boundary_con_offset + 2*first_edge+1, bin_var_offset+edge_pair_offset+2*j, 1);
	lp_descr.add_entry(bin_cons_offset+boundary_con_offset + 2*first_edge,   bin_var_offset+edge_pair_offset+2*j+1, -1);
      }

      if (mesh.edge(second_edge).from_idx_ == middle_point) {
	lp_descr.add_entry(bin_cons_offset+boundary_con_offset + 2*second_edge,   bin_var_offset+edge_pair_offset+2*j, -1);
	lp_descr.add_entry(bin_cons_offset+boundary_con_offset + 2*second_edge+1, bin_var_offset+edge_pair_offset+2*j+1, 1);
      }
      else {
	lp_descr.add_entry(bin_cons_offset+boundary_con_offset + 2*second_edge+1, bin_var_offset+edge_pair_offset+2*j, -1);
	lp_descr.add_entry(bin_cons_offset+boundary_con_offset + 2*second_edge,   bin_var_offset+edge_pair_offset+2*j+1, 1);
      }
    }

    if (enforce_consistent_boundaries) {

      //constraints in words: for each oriented edge, the pairs that start with this oriented edge 
      // and the pairs that end in the oppositely oriented edge may not sum to more than 1.0
      // (i.e. they are mutually exclusive)

      for (uint e=0; e < mesh.nEdges(); e++) {

 	double bound = max_border;

// 	if (edge_status[e] != BorderEdge) {
// 	  bound -= min_border;
// 	}

	rhs_upper[bin_cons_offset+consistency_con_offs+light_factor*e] =  bound;
	if (!light_constraints)
	  rhs_upper[bin_cons_offset+consistency_con_offs+light_factor*e+1] =  bound;
      }

      for (uint j=0; j < edge_pairs.size(); j++) {

	uint middle_point = edge_pairs[j].common_point_idx_;

	uint first_edge = edge_pairs[j].first_edge_idx_;
	uint second_edge = edge_pairs[j].second_edge_idx_;

	if (mesh.edge(first_edge).to_idx_ == middle_point) {      
	  lp_descr.add_entry(bin_cons_offset+consistency_con_offs + light_factor*first_edge, bin_var_offset+edge_pair_offset+2*j, 1);
	  lp_descr.add_entry(bin_cons_offset+consistency_con_offs + light_factor*first_edge, bin_var_offset+edge_pair_offset+2*j+1, 1);
	}
	else {
	  if (!light_constraints) {
	    lp_descr.add_entry(bin_cons_offset+consistency_con_offs + light_factor*first_edge+1, bin_var_offset+edge_pair_offset+2*j, 1);
	    lp_descr.add_entry(bin_cons_offset+consistency_con_offs + light_factor*first_edge+1, bin_var_offset+edge_pair_offset+2*j+1, 1);
	  }
	}

	if (mesh.edge(second_edge).to_idx_ == middle_point) {
	  lp_descr.add_entry(bin_cons_offset+consistency_con_offs + light_factor*second_edge, bin_var_offset+edge_pair_offset+2*j, 1);
	  lp_descr.add_entry(bin_cons_offset+consistency_con_offs + light_factor*second_edge, bin_var_offset+edge_pair_offset+2*j+1, 1);
	}
	else {
	  if (!light_constraints) {
	    lp_descr.add_entry(bin_cons_offset+consistency_con_offs + light_factor*second_edge+1, bin_var_offset+edge_pair_offset+2*j, 1);
	    lp_descr.add_entry(bin_cons_offset+consistency_con_offs + light_factor*second_edge+1, bin_var_offset+edge_pair_offset+2*j+1, 1);
	  }
	}
      }
    }

    if (enforce_regionedge) {
      uint rowoff = regionedge_constraints_offs;

      for (uint edge=0; edge < mesh.nEdges(); edge++) {

	//Get the two adjacent faces
	if (mesh.adjacent_faces(edge).size() != 2)
	  {
	    //One of the edges is at the border of the image

	    continue;
	  }

	uint x1 = bin_var_offset+mesh.adjacent_faces(edge)[0];
	uint x2 = bin_var_offset+mesh.adjacent_faces(edge)[1];

	double bound = max_border;

	lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge  , x1, 1);
	lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge  , x2, 1);
	lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge+1, x1, -1);
	lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge+1, x2, -1);
	if (!light_constraints) {
	  lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge+2, x1, 1);
	  lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge+2, x2, 1);
	  lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge+3, x1, -1);
	  lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge+3, x2, -1);
	}

	rhs_lower[bin_cons_offset+rowoff+2*light_factor*edge]   = 0;
	rhs_upper[bin_cons_offset+rowoff+2*light_factor*edge]   = 2*bound;
	rhs_lower[bin_cons_offset+rowoff+2*light_factor*edge+1] = -2*bound;
	rhs_upper[bin_cons_offset+rowoff+2*light_factor*edge+1] = 0;
	if (!light_constraints) {
	  rhs_lower[bin_cons_offset+rowoff+2*light_factor*edge+2] = 0;
	  rhs_upper[bin_cons_offset+rowoff+2*light_factor*edge+2] = 2*bound;
	  rhs_lower[bin_cons_offset+rowoff+2*light_factor*edge+3] = -2*bound;
	  rhs_upper[bin_cons_offset+rowoff+2*light_factor*edge+3] = 0;
	}
      }

      for (uint j=0; j < edge_pairs.size(); j++) {
	uint first = edge_pairs[j].first_edge_idx_;
	uint second = edge_pairs[j].second_edge_idx_;	

	uint y = bin_var_offset+mesh.nFaces() + 2*j;

	uint edge = first;
	if (mesh.adjacent_faces(edge).size() == 2) {
	  lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge   ,y  , 1);
	  lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge+1 ,y  , 1);
	  if (!light_constraints) {
	    lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge+2 ,y+1, 1);
	    lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge+3 ,y+1, 1);
	  }
	}

	edge = second;
	if (mesh.adjacent_faces(edge).size() == 2) {
	  lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge   ,y+1, 1);
	  lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge+1 ,y+1, 1);
	  if (!light_constraints) {
	    lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge+2 ,y  , 1);
	    lp_descr.add_entry(bin_cons_offset+rowoff+2*light_factor*edge+3 ,y  , 1);
	  }
	}
      }
    }
  }

  // coding level constraints xi>xi+1 on regions 
  uint ll_cons_offset = nBin*nConstraintsPerLayer;
  if (enforce_level_consistency){
    for (uint kk=0 ;kk < nBin - 1;kk++){
      uint level_cons_offset = ll_cons_offset+kk*nVarsPerLayer;
      uint level_var_offset = kk*nVarsPerLayer;
      for (uint jj=0 ;jj < mesh.nFaces(); jj++){

	if (var_ub[level_var_offset+jj] !=  var_lb[level_var_offset+jj]) {
	  //var not fixed
	  lp_descr.add_entry(level_cons_offset+jj,level_var_offset+jj,1);
	  lp_descr.add_entry(level_cons_offset+jj,level_var_offset+nVarsPerLayer+jj,-1);
	  rhs_lower[level_cons_offset+jj] = 0;
	  rhs_upper[level_cons_offset+jj] = max_border;
	}
      }
    }
  }

  Math1D::NamedVector<uint> row_start(nConstraints+1,MAKENAME(row_start));
  lp_descr.sort_by_row(row_start);

  nConstraints = row_start.size()-1; //since level constraints are not included for fixed pixels this
                                           //value may have decreased

  //std::cerr << "row_start.size(): " << row_start.size() << std::endl;
  //std::cerr << "estimated #cons: " << nConstraints << std::endl;

  //assert(row_start.size() == nConstraints+1);

  bool solver_known = false;

  const double* lp_solution = 0;

  int error = 0;

#ifdef HAS_GUROBI

  Math1D::Vector<double> gurobi_solution;

  if (solver == "gurobi") {

    solver_known = true;

    std::cerr << "using gurobi" << std::endl;
    
    GRBenv   *grb_env   = NULL;
    GRBmodel *grb_model = NULL;
    
    /* Create environment */
    
    error = GRBloadenv(&grb_env,NULL);
    GRBsetintparam(grb_env, GRB_INT_PAR_METHOD, GRB_METHOD_BARRIER);
    GRBsetdblparam(grb_env, "BarConvTol", 1e-10);
    GRBsetintparam(grb_env, "Crossover", 0);
    GRBsetintparam(grb_env, "CrossoverBasis", 1);
    GRBsetintparam(grb_env, "Presolve", 1);
    GRBsetintparam(grb_env, "PrePasses", 2);

    assert (!error && grb_env != NULL);

    /* Create an empty model */
    
    error = GRBnewmodel(grb_env, &grb_model, "curv-inpaint-lp", 0, NULL, NULL, NULL, NULL, NULL);
    assert(!error);
    
    NamedStorage1D<char> vtype(nVars,GRB_CONTINUOUS,MAKENAME(vtype));
    
    error = GRBaddvars(grb_model,nVars,0,NULL,NULL,NULL,cost.direct_access(),var_lb.direct_access(),
		       var_ub.direct_access(),vtype.direct_access(),NULL);
    assert(!error);
    
    error = GRBupdatemodel(grb_model);
    assert(!error);
    
    for (uint c=0; c < nConstraints; c++) {
      
      if (rhs_lower[c] == rhs_upper[c]) {
	error = GRBaddconstr(grb_model, row_start[c+1]-row_start[c], ((int*) lp_descr.col_indices()) + row_start[c], 
			     lp_descr.value() + row_start[c], GRB_EQUAL, rhs_lower[c], NULL);
      }
      else {
	error = GRBaddrangeconstr(grb_model, row_start[c+1]-row_start[c], ((int*) lp_descr.col_indices()) + row_start[c], 
				  lp_descr.value() + row_start[c], rhs_lower[c], rhs_upper[c], NULL);
      }
      
      assert(!error);
    }
    
    /* Optimize model */
    error = GRBoptimize(grb_model);
    assert(!error);
    
    gurobi_solution.resize(nVars);
    
    for (uint v=0; v < nVars; v++) {
      GRBgetdblattrelement(grb_model,"X",v, gurobi_solution.direct_access()+v);
    }
    lp_solution = gurobi_solution.direct_access();
    
    GRBfreemodel(grb_model);
    GRBfreeenv(grb_env);
  }
#endif

#ifdef HAS_CPLEX

  Math1D::Vector<double> cplex_solution;

  if (solver == "cplex") {

    solver_known = true;

    CPXENVptr     env = NULL;
    CPXLPptr      lp = NULL;
    int cpx_status = 0;
    
    /* Initialize the CPLEX environment */
    
    env = CPXopenCPLEX (&cpx_status);
    //CPXsetintparam(env, CPX_PARAM_STARTALG, CPX_ALG_BARRIER);
    //CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);
    //CPXsetintparam(env, CPX_PARAM_PREIND, CPX_OFF);
    //CPXsetintparam(env, CPX_PARAM_PREPASS, 0);
    
    //CPXsetintparam(env, CPX_PARAM_LPMETHOD, CPX_ALG_DUAL);
    
    /* If an error occurs, the status value indicates the reason for
       failure.  A call to CPXgeterrorstring will produce the text of
       the error message.  Note that CPXopenCPLEX produces no output,
       so the only way to see the cause of the error is to use
       CPXgeterrorstring.  For other CPLEX routines, the errors will
       be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON.  */
    
    if ( env == NULL ) {
      char  errmsg[1024];
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, cpx_status, errmsg);
      fprintf (stderr, "%s", errmsg);
      exit(1);
    }
    
    /* Turn on output to the screen */
    
    cpx_status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
    if ( cpx_status ) {
      fprintf (stderr,
	       "Failure to turn on screen indicator, error %d.\n", cpx_status);
      exit(1);
    }
    
    //necessary when using own cut generator (or heuristic??) with CPLEX
    //cpx_status = CPXsetintparam (env, CPX_PARAM_PREIND, CPX_OFF);
    
    //set problem data
    
    lp = CPXcreateprob (env, &cpx_status, "curv-inpaint-ilp");

    /* A returned pointer of NULL may mean that not enough memory
       was available or there was some other problem.  In the case of
       failure, an error message will have been written to the error
       channel from inside CPLEX.  In this example, the setting of
       the parameter CPX_PARAM_SCRIND causes the error message to
       appear on stdout.  */
    
    if ( lp == NULL ) {
      fprintf (stderr, "Failed to create LP.\n");
      exit(1);
    }
    
    /* Now copy the problem data into the lp */
    

    char* row_sense = new char[nConstraints];
    for (uint c=0; c < nConstraints; c++) {
      
      if (rhs_lower[c] == rhs_upper[c]) {
	row_sense[c] = 'E';
      }
      else {
	assert(rhs_lower[c] == -MAX_DOUBLE || rhs_lower[c] == 0.0);
	row_sense[c] = 'L';
      }
    }

    int* row_count = new int[nConstraints];
    for (uint c=0; c < nConstraints; c++)
      row_count[c] = row_start[c+1] - row_start[c];
    
    cpx_status = CPXnewcols (env, lp, nVars, cost.direct_access(), var_lb.direct_access(), 
			   var_ub.direct_access(), NULL, NULL);
    if ( cpx_status )  
      exit(1);
    
    CPXaddrows(env, lp, 0, nConstraints, lp_descr.nEntries(), rhs_upper.direct_access(), row_sense, 
	       (int*) row_start.direct_access(), (int*) lp_descr.col_indices(), lp_descr.value(),
	       NULL, NULL);
    
    delete[] row_sense;
    delete[] row_count;
  
    std::cerr << "calling optimize" << std::endl;
    
    cpx_status = CPXbaropt(env,lp);
    //cpx_status = CPXlpopt(env,lp);

    if ( cpx_status ) {
      fprintf (stderr, "Failed to optimize MIP.\n");
      exit(1);
    }
    
    cplex_solution.resize(nVars);

    CPXsolution (env, lp, NULL, NULL, cplex_solution.direct_access(), NULL, NULL, NULL);

    lp_solution = cplex_solution.direct_access();

    CPXfreeprob (env, &lp);
    CPXcloseCPLEX (&env);
  }

#endif

#ifdef HAS_XPRESS
  Math1D::Vector<double> xpress_solution;
  
  if (solver == "xpress") {

    solver_known = true;
    int nReturn;
    
    nReturn=XPRSinit("/opt/xpressmp/");
    
    std::cerr << "nReturn: " << nReturn << std::endl;
    if (nReturn != 0) {
      
      char msg[512];
      XPRSgetlicerrmsg(msg,512);
      
      std::cerr << "error message: " << msg << std::endl;
    }
    
    assert(nReturn == 0);
    
    XPRSprob xp_prob;
    char banner[256];
    
    XPRSgetbanner(banner); printf("banner: %s \n",banner);
    
    nReturn=XPRScreateprob(&xp_prob);
    
    //XPRSsetintcontrol(xp_prob,XPRS_PRESOLVE, 0);
    XPRSsetintcontrol(xp_prob,XPRS_CROSSOVER,0); 

    std::cerr << "nReturn: " << nReturn << std::endl;
    
    SparseMatrixDescription<double> lp_copy(lp_descr);
    
    Math1D::Vector<uint> col_start(nVars+1);
    lp_copy.sort_by_column(col_start);
    
    char* row_sense = new char[nConstraints];
    double* row_range = new double[nConstraints];
    
    for (uint c=0; c < nConstraints; c++) {

      if (rhs_lower[c] == rhs_upper[c]) {
	row_sense[c] = 'E';
	row_range[c] = 0.0;
      }
      else {
	row_sense[c] = 'R';
	row_range[c] = rhs_upper[c] - rhs_lower[c];
      }
    }
    
    nReturn = XPRSloadlp(xp_prob, "curv-inpaint-lp", col_start.size()-1, row_start.size()-1, row_sense,
			 rhs_upper.direct_access(), row_range, cost.direct_access(), 
			 (int*) col_start.direct_access(), NULL, (int*) lp_copy.row_indices(), lp_copy.value(),
			 var_lb.direct_access(), var_ub.direct_access());
    
    delete[] row_sense;
    delete[] row_range;
    
    Math1D::Vector<char> col_type(nVars,'C');
    Math1D::Vector<int> idx(nVars);
    for (uint v=0; v < nVars; v++)
      idx[v] = v;
    
    nReturn = XPRSchgcoltype(xp_prob, col_start.size(), idx.direct_access(), col_type.direct_access());
    
    XPRSlpoptimize(xp_prob,"b");
   
    xpress_solution.resize(nVars);

    XPRSgetlpsol(xp_prob, xpress_solution.direct_access(), 0, 0, 0); 

    lp_solution = xpress_solution.direct_access();
 
    nReturn=XPRSdestroyprob(xp_prob);
  }
#endif

  if (!solver_known && solver != "clp") {

    std::cerr << "WARNING: solver\"" << solver << "\" unknown. Using Clp instead" << std::endl;
    solver = "clp";
  }

  ClpSimplex lpSolver;

  if (solver == "clp") {

    std::clock_t tStartCLP,tEndCLP;

    CoinPackedMatrix coinMatrix(false,(int*) lp_descr.row_indices(),(int*) lp_descr.col_indices(),
				lp_descr.value(),lp_descr.nEntries());
    
    lpSolver.loadProblem (coinMatrix, var_lb.direct_access(), var_ub.direct_access(),   
			  cost.direct_access(), rhs_lower.direct_access(), rhs_upper.direct_access());

    coinMatrix.cleanMatrix();
    
    tStartCLP = std::clock();
    
    //lpSolver.dual();
    
    ClpSolve solve_options;
    solve_options.setSolveType(ClpSolve::useDual);
    //solve_options.setSolveType(ClpSolve::useBarrier);
    solve_options.setPresolveType(ClpSolve::presolveNumber,5);
    lpSolver.initialSolve(solve_options);
    
    error = 1 - lpSolver.isProvenOptimal();
  
    tEndCLP = std::clock();

    if (error != 0)
      std::cerr << "!!!!!!!!!!!!!!LP-solver failed!!!!!!!!!!!!!!!!!!!" << std::endl;

    std::cerr << "CLP-time: " << diff_seconds(tEndCLP,tStartCLP) << " seconds. " << std::endl;

    lp_solution = lpSolver.primalColumnSolution();

    //check for possible crossing points
    Storage1D<std::vector<uint> > point_pairs(mesh.nPoints());
      
    for (uint j=0; j < edge_pairs.size(); j++) {
	
      uint common_point = edge_pairs[j].common_point_idx_;
      point_pairs[common_point].push_back(j);
    }

    for (uint p=0; p < mesh.nPoints(); p++) {

      for (uint kk=0; kk < nBin; kk++){

	uint cur_offs = kk*nVarsPerLayer + edge_pair_offset;
      
	double sum = 0.0;
	for (std::vector<uint>::iterator it = point_pairs[p].begin(); it != point_pairs[p].end(); it++)
	  sum += lp_solution[cur_offs +  2*(*it)] + lp_solution[cur_offs +  2*(*it) + 1];
	
	if (sum > max_border) {
	  std::cerr << "checking if point #" << p << " in level " << kk << " is a crossing point" << std::endl;
	  
	  for (uint k1=0; k1 < point_pairs[p].size()-1; k1++) {
	    
	    uint pair1 = point_pairs[p][k1];
	    
	    for (uint k2=k1+1; k2 < point_pairs[p].size(); k2++) {
	      
	      uint pair2 = point_pairs[p][k2];
	      
	      double sum = lp_solution[cur_offs +  2*pair1] + lp_solution[cur_offs +  2*pair1 + 1]
		+ lp_solution[cur_offs +  2*pair2] + lp_solution[cur_offs +  2*pair2 + 1];
	      
	      if (sum >= 1.01*max_border) {
		std::cerr << " ----> yes" << std::endl;
	      }
	    }
	  }
	}
      }
    }
  }

  uint nFrac=0;

  for (uint i=0; i < nVars; i++) {

    double val = lp_solution[i];

    if (fabs(round(val) - val) > 0.02)
      nFrac++;
  }
  std::cerr << nFrac << "/" << nVars << " fractional variables." << std::endl;

  double energy = energy_offset;
  for (uint i=0; i < nVars; i++)
    energy += cost[i] * lp_solution[i];

  std::cerr << "energy: " << energy << std::endl;

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {
		
      if (status(x,y) == ToInpaint) {
	inpainted_image(x,y) = im_min_border;
      }
    }
  }
		   
  for (uint y=0; y < yDim; y++)   {
    for (uint x=0; x < xDim; x++)   {
      if (status(x,y) == ToInpaint) {
	for (uint kk=0; kk < nBin; kk++){
	  double sum = 0.0;
		
	  for (uint f=face_offset(x,y); f < face_offset(x,y) + nFaces(x,y); f++) {
	    sum += lp_solution[kk*nVarsPerLayer + f] * mesh.convex_area(f);
	  }
	  inpainted_image(x,y) += sum;		    	  
	}
	assert (status(x,y) != Border || fabs(inpainted_image(x,y) - image(x,y) < 0.01));
      }
    }
  } 


  return energy;
}
