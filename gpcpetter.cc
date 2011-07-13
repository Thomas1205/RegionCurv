/**** written by Petter Strandmark as an employee of Lund University, Sweden, summer 2010 ***/

#include "gpcpetter.hh"
#include "mesh2D.hh"

#include "gpc.h"

//
// Computes the intersection of two polygons. 
// The polygons are specified as vectors of points
//
std::vector<Mesh2DPoint> poly_intersection(const std::vector<Mesh2DPoint>& vec1, const std::vector<Mesh2DPoint>& vec2)
{
  gpc_polygon poly1, poly2;
  poly1.num_contours = 1;
  poly2.num_contours = 1;
  int hole = 0;
  poly1.hole = &hole;
  poly2.hole = &hole;
  gpc_vertex_list list1,list2;
  poly1.contour = &list1;
  poly2.contour = &list2;

  list1.num_vertices = (int)vec1.size();
  list2.num_vertices = (int)vec2.size();
  gpc_vertex* vertices1 = new gpc_vertex[vec1.size()];
  gpc_vertex* vertices2 = new gpc_vertex[vec2.size()];
  list1.vertex = vertices1;
  list2.vertex = vertices2;
  for (uint i=0; i<vec1.size(); ++i) {
    vertices1[i].x = vec1[i].x_;
    vertices1[i].y = vec1[i].y_;
  }
  for (uint i=0; i<vec2.size(); ++i) {
    vertices2[i].x = vec2[i].x_;
    vertices2[i].y = vec2[i].y_;
  }

  gpc_polygon result;
  gpc_polygon_clip(GPC_INT, &poly1, &poly2, &result);
	
  std::vector<Mesh2DPoint> vecr;
  if (result.num_contours > 0) {
    for (int i=0; i<result.contour->num_vertices; ++i) {
      Mesh2DPoint coord;
      coord.x_ = result.contour->vertex[i].x;
      coord.y_ = result.contour->vertex[i].y;
      vecr.push_back(coord);
    }
  }
  delete[] vertices1;
  delete[] vertices2;
  gpc_free_polygon(&result);
  return vecr;
}
