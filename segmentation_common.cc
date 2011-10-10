/*** First version written by Thomas Schoenemann as a private person without employment, October 2009 ***/
/*** continued by Thomas Schoenemann as an employee of Lund University, Sweden, January 2010 ***/
/*** further extensions by Petter Strandmark as an employee of Lund University, Sweden, summer 2010 ***/ 

#include <iostream>
#include <set>

#include "segmentation_common.hh"
#include "lp_segmentation.hh"
#include "mesh2D.hh"
#include "timing.hh"

#include "gpcpetter.hh"

#include "Petter-Color.hh"


LPSegOptions::LPSegOptions() {
  neighborhood_ = 16;
  enforce_consistent_boundaries_ = false;
  enforce_consistent_points_ = false;
  enforce_regionedge_ = false;
  prevent_crossings_ = false;
  light_constraints_ = false;
  griddim_xDim_ = 16;
  griddim_yDim_ = 16;
  adaptive_mesh_n_ = -1;
  solver_ = "clp";
  lambda_ = 0;
  gamma_ = 0;
  fix_regions_ = false;
  gridtype_ = Square;
  save_matrix_ = false;
  bruckstein_ = true;
  output_factor_ = 5;
  debug_svg_ = false;
  refine_ = false,
  convex_prior_ = false;
}

void add_grid_output(uint region_index, double label, const Mesh2D& mesh, const Math2D::Matrix<double>& output)
{
  //using namespace std;

  std::vector<Mesh2DPoint> points;
  mesh.get_polygon_points(region_index, points);

  //Calculate min/max of x and y to figure out
  //which pixels we have to intersect with
  double minx = 1e20, maxx = -1e20, miny = 1e20, maxy = -1e20;
  for (uint i=0; i<points.size(); ++i) {
    double x = points[i].x_;
    double y = points[i].y_;
    if (x < minx) minx = x;
    if (x > maxx) maxx = x;
    if (y < miny) miny = y;
    if (y > maxy) maxy = y;
  }
  //Round to integer
  minx = floor(minx);
  maxx = ceil(maxx);
  miny = floor(miny);
  maxy = ceil(maxy);

  //Polygon which represents a pixel
  std::vector<Mesh2DPoint> pixel(4);

  for (int y=int(miny); y < int(maxy); ++y) {
    for (int x=int(minx); x < int(maxx); ++x) {
      if (x<0 || x >= ((int) output.xDim()) || y<0 || y >= ((int) output.yDim())) continue;

      //Intersect the region with this pixel
      pixel[0].x_ = x;
      pixel[0].y_ = y;
      pixel[1].x_ = x+1;
      pixel[1].y_ = y;
      pixel[2].x_ = x+1;
      pixel[2].y_ = y+1;
      pixel[3].x_ = x;
      pixel[3].y_ = y+1;

      std::vector<Mesh2DPoint> intersection = poly_intersection(pixel, points);

      double area = 0;
      if (intersection.size() > 2) {
        for (uint i=1; i<intersection.size()-1; ++i) {
          area += triangle_area(intersection[0], intersection[i], intersection[i+1]);
        }
      }

      if (label < 0) {
        output(x,y) += area * 0.5;
      }
      else {
        output(x,y) += area*double(label);
      }

    }
  }
}

void add_grid_output_mreg(uint face_index, const double* label, uint nRegions, const Mesh2D& mesh, const Math3D::Tensor<double>& output) {

  std::vector<Mesh2DPoint> points;
  mesh.get_polygon_points(face_index, points);

  //Calculate min/max of x and y to figure out
  //which pixels we have to intersect with
  double minx = 1e20, maxx = -1e20, miny = 1e20, maxy = -1e20;
  for (uint i=0; i<points.size(); ++i) {
    double x = points[i].x_;
    double y = points[i].y_;
    if (x < minx) minx = x;
    if (x > maxx) maxx = x;
    if (y < miny) miny = y;
    if (y > maxy) maxy = y;
  }
  //Round to integer
  minx = floor(minx);
  maxx = ceil(maxx);
  miny = floor(miny);
  maxy = ceil(maxy);

  //Polygon which represents a pixel
  std::vector<Mesh2DPoint> pixel(4);

  for (int y=int(miny); y < int(maxy); ++y) {
    for (int x=int(minx); x < int(maxx); ++x) {
      if (x<0 || x >= ((int) output.xDim()) || y<0 || y >= ((int) output.yDim())) continue;

      //Intersect the region with this pixel
      pixel[0].x_ = x;
      pixel[0].y_ = y;
      pixel[1].x_ = x+1;
      pixel[1].y_ = y;
      pixel[2].x_ = x+1;
      pixel[2].y_ = y+1;
      pixel[3].x_ = x;
      pixel[3].y_ = y+1;

      std::vector<Mesh2DPoint> intersection = poly_intersection(pixel, points);

      double area = 0;
      if (intersection.size() > 2) {
        for (uint i=1; i<intersection.size()-1; ++i) {
          area += triangle_area(intersection[0], intersection[i], intersection[i+1]);
        }
      }

      for (uint r = 0; r < nRegions; r++) {
	
	if (label[r] < 0.0) {
	  output(x,y,r) += area * 0.5;
	}
	else {
	  output(x,y,r) += area*double(label[r]);
	}
      }
    }
  }

}


double calculate_data_term(uint region_index, Mesh2D& mesh, const Math2D::Matrix<float>& data_term)
{
  using namespace std;

  std::vector<Mesh2DPoint> points;
  mesh.get_polygon_points(region_index, points);

  //Calculate min/max of x and y to figure out
  //which pixels we have to intersect with
  double minx = 1e20, maxx = -1e20, miny = 1e20, maxy = -1e20;
  for (uint i=0; i<points.size(); ++i) {
    double x = points[i].x_;
    double y = points[i].y_;
    if (x < minx) minx = x;
    if (x > maxx) maxx = x;
    if (y < miny) miny = y;
    if (y > maxy) maxy = y;
  }
  //Round to integer
  minx = floor(minx);
  maxx = ceil(maxx);
  miny = floor(miny);
  maxy = ceil(maxy);

  //Polygon which represents a pixel
  std::vector<Mesh2DPoint> pixel(4);

  //Go through every pixel that possibly intersects
  //the region
  double total = 0;
  double total_area = 0;
  for (int y=int(miny); y < int(maxy); ++y) {
    for (int x=int(minx); x < int(maxx); ++x) {
      if (x<0 || x>= ((int) data_term.xDim()) || y<0 || y >= ((int) data_term.yDim())) continue;

      //cerr << "test " << x << "," << y << "  ";
      //Intersect the region with this pixel
      pixel[0].x_ = x;
      pixel[0].y_ = y;
      pixel[1].x_ = x+1;
      pixel[1].y_ = y;
      pixel[2].x_ = x+1;
      pixel[2].y_ = y+1;
      pixel[3].x_ = x;
      pixel[3].y_ = y+1;

      std::vector<Mesh2DPoint> intersection = poly_intersection(pixel, points);

      double area = 0;
      if (intersection.size() > 2) {
        for (uint i=1; i<intersection.size()-1; ++i) {
          area += triangle_area(intersection[0], intersection[i], intersection[i+1]);
        }
      }

      total += area*data_term(x,y);
      total_area += area;
    }
  }

  return total;
}

double calculate_data_term(uint region_index, uint r, Mesh2D& mesh, const Math3D::Tensor<float>& data_term) {

  std::vector<Mesh2DPoint> points;
  mesh.get_polygon_points(region_index, points);

  //Calculate min/max of x and y to figure out
  //which pixels we have to intersect with
  double minx = 1e20, maxx = -1e20, miny = 1e20, maxy = -1e20;
  for (uint i=0; i<points.size(); ++i) {
    double x = points[i].x_;
    double y = points[i].y_;
    if (x < minx) minx = x;
    if (x > maxx) maxx = x;
    if (y < miny) miny = y;
    if (y > maxy) maxy = y;
  }
  //Round to integer
  minx = floor(minx);
  maxx = ceil(maxx);
  miny = floor(miny);
  maxy = ceil(maxy);

  //Polygon which represents a pixel
  std::vector<Mesh2DPoint> pixel(4);

  //Go through every pixel that possibly intersects
  //the region
  double total = 0;
  double total_area = 0;
  for (int y=int(miny); y < int(maxy); ++y) {
    for (int x=int(minx); x < int(maxx); ++x) {
      if (x<0 || x>= ((int) data_term.xDim()) || y<0 || y >= ((int) data_term.yDim())) continue;

      //cerr << "test " << x << "," << y << "  ";
      //Intersect the region with this pixel
      pixel[0].x_ = x;
      pixel[0].y_ = y;
      pixel[1].x_ = x+1;
      pixel[1].y_ = y;
      pixel[2].x_ = x+1;
      pixel[2].y_ = y+1;
      pixel[3].x_ = x;
      pixel[3].y_ = y+1;

      std::vector<Mesh2DPoint> intersection = poly_intersection(pixel, points);

      double area = 0;
      if (intersection.size() > 2) {
        for (uint i=1; i<intersection.size()-1; ++i) {
          area += triangle_area(intersection[0], intersection[i], intersection[i+1]);
        }
      }

      total += area*data_term(x,y,r);
      total_area += area;
    }
  }

  return total;
}


void generate_mesh(uint xDim, uint yDim, uint neighborhood, Mesh2D& mesh, bool silent, const Math2D::Matrix<int>* seeds) 
{
  if (!silent) {
    Petter::statusTry("Generating mesh...");
  }

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
      for (uint x=0; x < xDim; x++) {
        mesh.add_point( Mesh2DPoint(x+0.5, y) );
      }
    }

    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x <= xDim; x++) {
        mesh.add_point( Mesh2DPoint(x, y+0.5) );
      }
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

  //std::cerr << "generating edges" << std::endl;

  /** generate mesh edges **/
  if (neighborhood <= 8) {
    //horizontal edges
    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x <= xDim; x++) {
        mesh.add_edge(y*(xDim+1)+x, (y+1)*(xDim+1)+x);
      }
    }

    //vertical edges
    for (uint y=0; y <= yDim; y++) {
      for (uint x=0; x < xDim; x++) {
        mesh.add_edge(y*(xDim+1)+x,y*(xDim+1)+x+1);
      }
    }

    //diagonal edges
    if (neighborhood == 8) {
      for (uint y=0; y < yDim; y++) {
        for (uint x=0; x < xDim; x++) {

          if (seeds == 0 || (*seeds)(x,y) < 0) {
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

        mesh.add_edge( y*(xDim+1)+x, hmidpoint_offset + y*(xDim)+x );
        mesh.add_edge( y*(xDim+1)+x+1, hmidpoint_offset + y*(xDim)+x );
      }
    }

    uint vmidpoint_offset = hmidpoint_offset + (yDim+1)*xDim;

    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x <= xDim; x++) {

        mesh.add_edge( y*(xDim+1)+x  , vmidpoint_offset + y*(xDim+1)+x );
        mesh.add_edge( (y+1)*(xDim+1)+x  , vmidpoint_offset + y*(xDim+1)+x );
      }
    }

    //add interior edges

    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x < xDim; x++) {

        if (seeds == 0 || (*seeds)(x,y) < 0) {

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

  //std::cerr << "generating faces" << std::endl;

  std::vector<uint> face;

  /** generate mesh faces **/
  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {

      if (neighborhood==4) {
        face.clear();
        face.push_back(mesh.find_edge(y*(xDim+1)+x,(y+1)*(xDim+1)+x).first);
        face.push_back(mesh.find_edge((y+1)*(xDim+1)+x,(y+1)*(xDim+1)+x+1).first);
        face.push_back(mesh.find_edge(y*(xDim+1)+x+1,(y+1)*(xDim+1)+x+1).first);
        face.push_back(mesh.find_edge(y*(xDim+1)+x,y*(xDim+1)+x+1).first);

        mesh.add_face(face);
      }
      else if (neighborhood==8) {

        if (seeds == 0 || (*seeds)(x,y) < 0) {

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
        else {

          face.clear();
          face.push_back(mesh.find_edge(y*(xDim+1)+x,(y+1)*(xDim+1)+x).first);
          face.push_back(mesh.find_edge((y+1)*(xDim+1)+x,(y+1)*(xDim+1)+x+1).first);
          face.push_back(mesh.find_edge(y*(xDim+1)+x+1,(y+1)*(xDim+1)+x+1).first);
          face.push_back(mesh.find_edge(y*(xDim+1)+x,y*(xDim+1)+x+1).first);

          mesh.add_face(face);
        }
      }
      else if (neighborhood==16) {

        if (seeds == 0 || (*seeds)(x,y) < 0) {

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
        else {

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
      }
      else {
        INTERNAL_ERROR << "connectivity \"" << neighborhood << "\" is not supported" << std::endl;
      }
    }
  }

  if (!silent) {
    Petter::statusOK();
    std::cerr << mesh.nPoints() << " points, " << mesh.nFaces() << " faces, " << mesh.nEdges() << " edges." << std::endl;
  }

  //return nAreasPerPixel;
}

void add_hex_to_mesh(double x, double y, double w, Mesh2D& mesh) 
{
  double dw = 1.154700538379252*w; //diagonal length
  dw /= 2; //radius
  w /= 2;  //radius
  uint pc  = mesh.find_or_add_point( Mesh2DPoint(x,y) );
  std::vector<uint> p(12);
  p[0]  = mesh.find_or_add_point( Mesh2DPoint(x,y+dw) );
  p[1]  = mesh.find_or_add_point( Mesh2DPoint(x+0.5*w,y+0.866025403784439*w) );
  p[2]  = mesh.find_or_add_point( Mesh2DPoint(x+0.866025403784439*dw,y+0.5*dw) );
  p[3]  = mesh.find_or_add_point( Mesh2DPoint(x+w,y) );
  p[4]  = mesh.find_or_add_point( Mesh2DPoint(x+0.866025403784439*dw,y-0.5*dw) );
  p[5]  = mesh.find_or_add_point( Mesh2DPoint(x+0.5*w,y-0.866025403784439*w) );
  p[6]  = mesh.find_or_add_point( Mesh2DPoint(x,y-dw) );
  p[7]  = mesh.find_or_add_point( Mesh2DPoint(x-0.5*w,y-0.866025403784439*w) );
  p[8]  = mesh.find_or_add_point( Mesh2DPoint(x-0.866025403784439*dw,y-0.5*dw) );
  p[9]  = mesh.find_or_add_point( Mesh2DPoint(x-w,y) );
  p[10] = mesh.find_or_add_point( Mesh2DPoint(x-0.866025403784439*dw,y+0.5*dw) );
  p[11] = mesh.find_or_add_point( Mesh2DPoint(x-0.5*w,y+0.866025403784439*w) );

  std::vector<uint> l(24);
  for (int li=0;li<12;li++) {
    l[li] = mesh.find_or_add_edge(p[li],p[(li+1)%12]);
  }
  for (int li=0;li<12;li++) {
    l[12+li] = mesh.find_or_add_edge(pc,p[li]);
  }

  std::vector<uint> face;

  for (int fi=0;fi<12;fi++) {
    face.clear();
    face.push_back(l[12+fi]);
    face.push_back(l[fi]);
    face.push_back(l[12+(fi+1)%12]);
    mesh.add_face(face);
  }
}


void generate_hexagonal_mesh(uint xDim, uint yDim, double w, uint neighborhood, Mesh2D& mesh) 
{	
  using namespace Petter;

  statusTry("Generating hexagonal mesh...");
  const double dw = 1.154700538379252*w; //diagonal length
  const double s = dw/2; //Side length

  for (double x=0; x<=xDim+0.001; x+=w) {
    for (double y=0; y<yDim+dw; y+=s+dw) {
      add_hex_to_mesh(x,y,w, mesh);
    }
    for (double y=(s+dw)/2; y<yDim; y+=s+dw) {
      add_hex_to_mesh(x+0.5*w,y, w, mesh);
    }
  }
  statusOK();

  statusTry("Merging parallel edges...");
  mesh.merge_parallel_edges();
  statusOK();

  statusTry("Drawing...");
  if (mesh.nFaces() > 20000) {
    statusFailed();
  }
  else {
    mesh.draw("hex_mesh.svg", false);
    statusOK();
  }

  std::cerr << mesh.nPoints() << " points, " << mesh.nFaces() << " faces, " << mesh.nEdges() << " edges." << std::endl;
}


PixelFaceRelation::PixelFaceRelation() : face_idx_(MAX_UINT), share_(-1.0) {}

PixelFaceRelation::PixelFaceRelation(uint face_idx, float share) : face_idx_(face_idx), share_(share) {}


void compute_pixel_shares(const Mesh2D& mesh, uint xDim, uint yDim, Storage1D<PixelFaceRelation>& shares,
			  Math1D::Vector<uint>& share_start) {
  
  Storage1D< std::vector<PixelFaceRelation> > share_map(xDim*yDim);

  uint nShares = 0;

  for (uint f=0; f < mesh.nFaces(); f++) {

    std::vector<Mesh2DPoint> points;
    mesh.get_polygon_points(f, points);

    //Calculate min/max of x and y to figure out
    //which pixels we have to intersect with
    double minx = 1e20, maxx = -1e20, miny = 1e20, maxy = -1e20;
    for (uint i=0; i<points.size(); ++i) {
      double x = points[i].x_;
      double y = points[i].y_;
      if (x < minx) minx = x;
      if (x > maxx) maxx = x;
      if (y < miny) miny = y;
      if (y > maxy) maxy = y;
    }
    //Round to integer
    minx = floor(minx);
    maxx = ceil(maxx);
    miny = floor(miny);
    maxy = ceil(maxy);

    //Polygon which represents a pixel
    std::vector<Mesh2DPoint> pixel(4);

    for (int y=int(miny); y < int(maxy); ++y) {
      for (int x=int(minx); x < int(maxx); ++x) {
        if (x<0 || x >= ((int) xDim) || y<0 || y >= ((int) yDim)) continue;

        //Intersect the region with this pixel
        pixel[0].x_ = x;
        pixel[0].y_ = y;
        pixel[1].x_ = x+1;
        pixel[1].y_ = y;
        pixel[2].x_ = x+1;
        pixel[2].y_ = y+1;
        pixel[3].x_ = x;
        pixel[3].y_ = y+1;

        std::vector<Mesh2DPoint> intersection = poly_intersection(pixel, points);

        double area = 0;
        if (intersection.size() > 2) {
          for (uint i=1; i<intersection.size()-1; ++i) {
            area += triangle_area(intersection[0], intersection[i], intersection[i+1]);
          }
        }

        PixelFaceRelation new_rel(f,float(area / mesh.convex_area(f)));
        share_map[y*xDim+x].push_back(new_rel);
        nShares++;
      }
    }
  }

  shares.resize_dirty(nShares);  
  share_start.resize_dirty(xDim*yDim+1);

  uint nCur = 0;
  for (uint i=0; i < xDim*yDim; i++) {

    share_start[i] = nCur;
    for (size_t j=0; j < share_map[i].size(); j++) {

      shares[nCur+j] = share_map[i][j];
    }
    nCur += uint( share_map[i].size() );
  }
  share_start[xDim*yDim] = nCur;
}

