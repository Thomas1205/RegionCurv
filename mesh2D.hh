/*** first version written by Thomas Schoenemann as a private person without employment, September 2009 ***/
/*** extended by Thomas Schoenemann and Petter Strandmark as employees of Lund University, Sweden, September 2010 ***/

#ifndef MESH2D_HH
#define MESH2D_HH

#include <vector>

#include "makros.hh"


struct Mesh2DPoint {

  Mesh2DPoint(double x, double y);
  Mesh2DPoint();

  bool operator==(const Mesh2DPoint& cmp) const;

  double x_;
  double y_;

  std::pair<double,double> coord_pair() const;

  std::ostream& operator<<(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, const Mesh2DPoint& point);

bool lines_cross(Mesh2DPoint p1, Mesh2DPoint p2, Mesh2DPoint q1, Mesh2DPoint q2, 
		 std::pair<double,double>& crossing_point);

bool line_pairs_with_meeting_point_cross(const Mesh2DPoint& p1, const Mesh2DPoint& p2, 
					 const Mesh2DPoint& q1, const Mesh2DPoint& q2,
					 const Mesh2DPoint& meeting_point);

double triangle_area(const Mesh2DPoint& p1, const Mesh2DPoint& p2, const Mesh2DPoint& p3);

struct Mesh2DEdge {

  Mesh2DEdge(uint from_idx, uint to_idx); 

  uint from_idx_;
  uint to_idx_;

  std::ostream& operator<<(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, const Mesh2DEdge& edge);

struct Mesh2DEdgePair {

  Mesh2DEdgePair(uint first_edge, uint second_edge, uint common_point);

  uint first_edge_idx_;
  uint second_edge_idx_;
  uint common_point_idx_;

  std::ostream& operator<<(std::ostream& os) const;
};

std::ostream& operator<<(std::ostream& os, const Mesh2DEdgePair& pair);

struct Mesh2DFace {

  std::vector<uint> edge_idx_;

};


class Mesh2D {
public:

  Mesh2D();

  void add_point(const Mesh2DPoint& new_point);

  void add_edge(uint from_idx, uint to_idx);

  void add_face(const std::vector<uint>& edge_indices);

  double edge_length(uint edge_idx) const;

  double convex_area(uint face_idx) const;

  bool contains_point(uint face_idx, double x, double y) const;

  void get_polygon_points(uint face_idx, std::vector<uint>& point_indices) const;

  void get_polygon_points(uint face_idx, std::vector<Mesh2DPoint>& points) const;

  int match(uint face_idx, uint edge_idx) const;

  //draw mesh in svg format
  void draw(std::string filename, bool drawfaces=true) const;

  void draw_labels(std::string filename, const int* labels) const;

  void draw_labels_with_pairs(std::string filename, const double* labels, const std::vector<Mesh2DEdgePair>& edge_pairs, 
                              double xDim=-1, double yDim=-1) const;

  void generate_edge_pair_list(std::vector<Mesh2DEdgePair>& vec) const;

  const std::vector<uint>& adjacent_faces(uint edge_idx) const;

  std::vector<Mesh2DPoint>::const_iterator point_start() const;

  std::vector<Mesh2DPoint>::const_iterator point_end() const;

  std::vector<Mesh2DEdge>::const_iterator edge_start() const;

  std::vector<Mesh2DEdge>::const_iterator edge_end() const;

  std::vector<Mesh2DFace>::const_iterator face_start() const;

  std::vector<Mesh2DFace>::const_iterator face_end() const;

  Mesh2DPoint point(uint point_idx) const;

  Mesh2DEdge edge(uint edge_idx) const;

  uint nPoints() const;

  uint nEdges() const;

  uint nFaces() const;

  uint find_point(Mesh2DPoint point) const;

  uint find_or_add_point(Mesh2DPoint point);

  uint find_or_add_edge(uint to, uint from);

  //Post-processing of generated irregular meshes to remove 
  //parallel edges
  void replace_edge(uint e1,uint e2); //Removes e1 and replaces it with e2
  void move_edge(uint p, uint p1, uint p2, uint e1, uint e2); //Changes e1 from (p,p1) to (p1,p2)
  void merge_parallel_edges(); //Fixes all parallel edges

  //returns the index of the edge and a flag which is true if the edge was found and false
  //if the reverse edge was found
  std::pair<uint,bool> find_edge(Mesh2DPoint from, Mesh2DPoint to) const;

  //returns the index of the edge and a flag which is true if the edge was found and false
  //if the reverse edge was found
  std::pair<uint,bool> find_edge(uint from, uint to) const;

  void enlarge(double x_factor, double y_factor);

  Mesh2DPoint face_center(uint idx);

protected:
  std::vector<Mesh2DPoint> point_;
  std::vector<Mesh2DEdge> edge_;
  std::vector<Mesh2DFace> face_;

  std::vector<std::vector<uint> > point_adjacent_edges_;
  std::vector<std::vector<uint> > edge_adjacent_faces_;

  double min_x_;
  double min_y_;
  double max_x_;
  double max_y_;
};


#endif
