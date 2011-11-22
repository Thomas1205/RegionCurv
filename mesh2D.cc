/*** first version written by Thomas Schoenemann as a private person without employment, September 2009 ***/
/*** extended by Thomas Schoenemann and Petter Strandmark as employees of Lund University, Sweden, September 2010 ***/

#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <set>
#include <map>
#include "mesh2D.hh"
#include "vector.hh"
#include "svg.hh"

//DEBUG
//#include "stl_out.hh"
//END_DEBUG

bool operator==(std::pair<double, double> q, Mesh2DPoint p) {

  return (p.x_ == q.first && p.y_ == q.second);
}

bool lines_cross(Mesh2DPoint p1, Mesh2DPoint p2, Mesh2DPoint q1, Mesh2DPoint q2, 
		 std::pair<double,double>& crossing_point) {

    crossing_point.first =  -100000;
    crossing_point.second = -100000;

//     if ( ((p1.first-p2.first) * (q1.first-q2.first) + (p1.second-p2.second) * (q1.second - q2.second)) == 0.0)
// 	//the lines are collinnear
// 	return false;


    Math1D::Vector<double> pn(3);
    pn[0] = p1.y_ - p2.y_;
    pn[1] = p2.x_ - p1.x_;

    double pnorm = sqrt(pn[0]*pn[0] + pn[1]*pn[1]);
    if (pnorm == 0.0) //if any of the two lines reduces to a point, we define that the lines do not cross;
	return false;

    pn[0] /= pnorm;
    pn[1] /= pnorm;
    pn[2] = -p1.x_*pn[0] - p1.y_ * pn[1];

    Math1D::Vector<double> qn(3);
    qn[0] = q1.y_ - q2.y_;
    qn[1] = q2.x_ - q1.x_;

    double qnorm = sqrt(qn[0]*qn[0] + qn[1]*qn[1]);
    if (qnorm == 0)
	return false;

    qn[0] /= qnorm;
    qn[1] /= qnorm;
    qn[2] = -q1.x_*qn[0] - q1.y_*qn[1];

    Math1D::Vector<double> cross_vec = Math1D::cross(pn, qn); //pn / qn;
    
//     std::cerr << "pn: " << pn << std::endl;
//     std::cerr << "qn: " << qn << std::endl; 
//     std::cerr << "cross product: " << cross << std::endl;

    if (cross_vec[2] == 0.0) {
	//std::cerr << "parallel" << std::endl;
	return false;
    }
    
    cross_vec[0] /= cross_vec[2];
    cross_vec[1] /= cross_vec[2];
    cross_vec[2] = 1.0;

    //std::cerr << "crossing point: " << cross << std::endl;

    //if (cross[0] == p1.first || cross[0] == p2.first || cross[0] == q1.first || cross[0] == q2.first) {
    std::pair<double,double> cp = std::make_pair(cross_vec[0],cross_vec[1]);
    if (cp == p1 || cp == p2 || cp == q1 || cp ==q2) {
	//std::cerr << "touch" << std::endl;
	//the lines are allowed to touch (but not cross) at a single point
	return false;
    }

    double px = (p1.x_ - cross_vec[0]) * (cross_vec[0] - p2.x_);
    double py = (p1.y_ - cross_vec[1]) * (cross_vec[1] - p2.y_);
    
    //std::cerr << "px: " << px << std::endl;
    //std::cerr << "py: " << py << std::endl;

    if (px <= 0.001 && py <= 0.001)
	return false; //either the crossing point is outside the segment or coincides with one of the endpoints

    double qx = (q1.x_ - cross_vec[0]) * (cross_vec[0] - q2.x_);
    double qy = (q1.y_ - cross_vec[1]) * (cross_vec[1] - q2.y_);

    if (qx <= 0.001 && qy <= 0.001)
	return false; //either the crossing point is outside the segment or coincides with one of the endpoints
    
    crossing_point.first = cross_vec[0];
    crossing_point.second = cross_vec[1];
    return true;
}

bool line_pairs_with_meeting_point_cross(const Mesh2D& mesh, const Mesh2DEdgePair& pair1, const Mesh2DEdgePair& pair2) {

  if (pair1.common_point_idx_ != pair2.common_point_idx_) {
    INTERNAL_ERROR << " line pairs do not meet in a mesh point. Exiting..." << std::endl;
    exit(1);
  }
  
  uint point = pair1.common_point_idx_;

  uint pair1_edge1 = pair1.first_edge_idx_;
  uint pair1_edge2 = pair1.second_edge_idx_;
  
  uint p1_idx = (mesh.edge(pair1_edge1).from_idx_ != point) ? 
    mesh.edge(pair1_edge1).from_idx_ : mesh.edge(pair1_edge1).to_idx_;
  
  uint p2_idx = (mesh.edge(pair1_edge2).from_idx_ != point) ? 
    mesh.edge(pair1_edge2).from_idx_ : mesh.edge(pair1_edge2).to_idx_;
  
  uint pair2_edge1 = pair2.first_edge_idx_;
  uint pair2_edge2 = pair2.second_edge_idx_;
  
  uint p3_idx = (mesh.edge(pair2_edge1).from_idx_ != point) ? 
    mesh.edge(pair2_edge1).from_idx_ : mesh.edge(pair2_edge1).to_idx_;
  
  uint p4_idx = (mesh.edge(pair2_edge2).from_idx_ != point) ? 
    mesh.edge(pair2_edge2).from_idx_ : mesh.edge(pair2_edge2).to_idx_;
  
  
  return ( line_pairs_with_meeting_point_cross(mesh.point(p1_idx), mesh.point(p2_idx), mesh.point(p3_idx), 
					       mesh.point(p4_idx), mesh.point(point)) );
}


bool line_pairs_with_meeting_point_cross(const Mesh2DPoint& p1, const Mesh2DPoint& p2, 
					 const Mesh2DPoint& q1, const Mesh2DPoint& q2,
					 const Mesh2DPoint& meeting_point) {

  Math1D::Vector<double> diff_p1(2);
  diff_p1[0] = meeting_point.x_ - p1.x_;
  diff_p1[1] = meeting_point.y_ - p1.y_;

  Math1D::Vector<double> diff_p2(2);
  diff_p2[0] = p2.x_ - meeting_point.x_;
  diff_p2[1] = p2.y_ - meeting_point.y_;

  Math1D::Vector<double> diff_q1(2);
  diff_q1[0] = meeting_point.x_ - q1.x_;
  diff_q1[1] = meeting_point.y_ - q1.y_;

  Math1D::Vector<double> diff_q2(2);
  diff_q2[0] = q2.x_ - meeting_point.x_;
  diff_q2[1] = q2.y_ - meeting_point.y_;
  
  double cos_p = (diff_p1 % diff_p2) / (diff_p1.norm() * diff_p2.norm());
  double cos_q = (diff_q1 % diff_q2) / (diff_q1.norm() * diff_q2.norm());

  std::pair<double,double> crossing_point;

  if (cos_p > 0.98 && cos_q > 0.98) {

    //two non-identical straight lines with a meeting point will always cross
    return true;
  }

  if (cos_p < cos_q) {

    //the p-line has the higher curvature
    
    //make sure the q-lines are long enough to detect crossings
    Mesh2DPoint aux1 = q1;
    aux1.x_ -= 100.0 * diff_q1[0];
    aux1.y_ -= 100.0 * diff_q1[1];

    Mesh2DPoint aux2 = q2;
    aux2.x_ += 100.0 * diff_q2[0];
    aux2.y_ += 100.0 * diff_q2[1];

    return (lines_cross(aux1,meeting_point,p1,p2,crossing_point)
	    || lines_cross(aux2,meeting_point,p1,p2,crossing_point) );
  }
  else {

    //the q-line has the higher curvature

    //make sure the p-lines are long enough to detect crossings
    Mesh2DPoint aux1 = p1;
    aux1.x_ -= 100.0 * diff_p1[0];
    aux1.y_ -= 100.0 * diff_p1[1];

    Mesh2DPoint aux2 = p2;
    aux2.x_ += 100.0 * diff_p2[0];
    aux2.y_ += 100.0 * diff_p2[1];

    return (lines_cross(aux1,meeting_point,q1,q2,crossing_point)
	    || lines_cross(aux2,meeting_point,q1,q2,crossing_point) );
  }

}

double triangle_area(const Mesh2DPoint& p1, const Mesh2DPoint& p2, const Mesh2DPoint& p3) {

  Math1D::Vector<double> d1(3,0.0);
  Math1D::Vector<double> d2(3,0.0);

  d1[0] = p2.x_ - p1.x_;
  d1[1] = p2.y_ - p1.y_;

  d2[0] = p3.x_ - p1.x_;
  d2[1] = p3.y_ - p1.y_;

  return 0.5 * cross(d1,d2).norm();
}

Mesh2DPoint::Mesh2DPoint(double x, double y) : x_(x), y_(y) {}
Mesh2DPoint::Mesh2DPoint() : x_(0), y_(0) {}

bool Mesh2DPoint::operator==(const Mesh2DPoint& cmp) const {

  return ( (std::abs(x_-cmp.x_) + std::abs(y_-cmp.y_)) <= 1e-4);
}

std::pair<double,double> Mesh2DPoint::coord_pair() const {

  return std::make_pair(x_,y_);
}

std::ostream& Mesh2DPoint::operator<<(std::ostream& os) const {

  os << "(" << x_ << "," << y_ << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const Mesh2DPoint& point) {

  os << "(" << point.x_ << "," << point.y_ << ")";
  return os;
}


Mesh2DEdge::Mesh2DEdge(uint from_idx, uint to_idx) : from_idx_(from_idx), to_idx_(to_idx) {}

std::ostream& Mesh2DEdge::operator<<(std::ostream& os) const {

  os << from_idx_ << " -> " << to_idx_;
  return os;
}

std::ostream& operator<<(std::ostream& os, const Mesh2DEdge& edge) {

  os << edge.from_idx_ << " -> " << edge.to_idx_;
  return os;

}

Mesh2DEdgePair::Mesh2DEdgePair(uint first_edge, uint second_edge, uint common_point) 
: first_edge_idx_(first_edge), second_edge_idx_(second_edge), common_point_idx_(common_point) 
{}

std::ostream& Mesh2DEdgePair::operator<<(std::ostream& os) const {

  os << "[" << first_edge_idx_ << "," << second_edge_idx_ << ", shared point=" <<  common_point_idx_ << "]";
  return os;
}

std::ostream& operator<<(std::ostream& os, const Mesh2DEdgePair& pair) {

  os << "[" << pair.first_edge_idx_ << "," << pair.second_edge_idx_ << ", shared point=" << pair.common_point_idx_ << "]";
  return os;
}

Mesh2D::Mesh2D() : min_x_(0.0), min_y_(0.0), max_x_(0.0), max_y_(0.0) {}

void Mesh2D::add_point(const Mesh2DPoint& new_point) {

  point_.push_back(new_point);
  point_adjacent_edges_.push_back(std::vector<uint>());

  max_x_ = std::max(max_x_,new_point.x_);
  max_y_ = std::max(max_x_,new_point.y_); //Bug
  min_x_ = std::min(min_x_,new_point.x_);
  min_y_ = std::min(min_y_,new_point.y_);  
}

void Mesh2D::add_edge(uint from_idx, uint to_idx) {

  assert(from_idx < point_.size());
  assert(to_idx < point_.size());

  edge_.push_back(Mesh2DEdge(from_idx,to_idx));
  edge_adjacent_faces_.push_back(std::vector<uint>());

  uint idx = uint( edge_.size()-1 );
  point_adjacent_edges_[from_idx].push_back(idx);
  point_adjacent_edges_[to_idx].push_back(idx);
}

void Mesh2D::add_face(const std::vector<uint>& edge_indices) {

  face_.push_back(Mesh2DFace());
  face_.back().edge_idx_ = edge_indices;

  uint idx = uint( face_.size()-1 );
  std::set<uint> point_indices;
  for (uint i=0; i < edge_indices.size(); i++) {
    uint edge_idx = edge_indices[i];
    assert(edge_idx < edge_.size());
    edge_adjacent_faces_[edge_indices[i]].push_back(idx);    

    uint next_edge_idx = edge_indices[(i+1) % edge_indices.size()]; 

    point_indices.clear();
    point_indices.insert(edge_[edge_idx].to_idx_);
    point_indices.insert(edge_[edge_idx].from_idx_);
    point_indices.insert(edge_[next_edge_idx].to_idx_);
    point_indices.insert(edge_[next_edge_idx].from_idx_);

    if (point_indices.size() != 3) {
      INTERNAL_ERROR << " not a proper face. Exiting..." << std::endl;
      exit(1);
    }
  }
}

double Mesh2D::edge_length(uint edge_idx) const {

  uint from_idx = edge_[edge_idx].from_idx_;
  uint to_idx = edge_[edge_idx].to_idx_;

  double diff_x = point_[from_idx].x_ - point_[to_idx].x_;
  double diff_y = point_[from_idx].y_ - point_[to_idx].y_;

  return sqrt(diff_x*diff_x + diff_y*diff_y);
}

double Mesh2D::convex_area(uint face_idx) const {

  assert(face_[face_idx].edge_idx_.size() > 2);

  std::vector<uint> points;
  get_polygon_points(face_idx,points);

  double sum = 0.0;
  for (uint j=2; j < points.size(); j++) {
    sum += triangle_area(point_[points[0]],point_[points[j-1]],point_[points[j]]);
  }

  return sum;
}

bool Mesh2D::contains_point(uint face_idx, double x, double y) const {

  std::vector<uint> points;
  get_polygon_points(face_idx,points);

  double last_dist = 0.0;

  for (uint i=0; i < points.size(); i++) {

    double cx = point_[points[i]].x_;
    double cy = point_[points[i]].y_;

    uint ni = (i+1) % points.size();

    double nx = point_[points[ni]].x_;
    double ny = point_[points[ni]].y_;

    double dist = x*(ny-cy) - y*(nx-cx) - cx*ny + cy*nx;

    if (dist*last_dist < 0)
      return false;

    last_dist = dist;
  }

  return true;
}

void Mesh2D::get_polygon_points(uint face_idx, std::vector<uint>& point_indices) const {

  point_indices.clear();

  std::map<uint,uint> point_count;

  uint nEdges = uint( face_[face_idx].edge_idx_.size() );

  for (uint j=0; j < nEdges; j++) {
    uint edge_idx = face_[face_idx].edge_idx_[j];
    uint next_edge_idx = face_[face_idx].edge_idx_[(j+1) % nEdges];

    point_count.clear();
    point_count[edge_[edge_idx].to_idx_]++;
    point_count[edge_[edge_idx].from_idx_]++;
    point_count[edge_[next_edge_idx].to_idx_]++;
    point_count[edge_[next_edge_idx].from_idx_]++;

    assert(point_count.size() == 3);

    for (std::map<uint,uint>::iterator it = point_count.begin(); it != point_count.end(); it++) {
      if (it->second == 2) {
        point_indices.push_back(it->first);
        break;
      }
    }
  }
}

void Mesh2D::get_polygon_points(uint face_idx, std::vector<Mesh2DPoint>& points) const 
{
  std::vector<uint> points_idx;
  get_polygon_points(face_idx,points_idx);
  points.clear();
  points.reserve(points_idx.size());
  for (uint i=0;i<points_idx.size(); ++i) {
    points.push_back( point_[points_idx[i]] );
  }
}


const std::vector<uint>& Mesh2D::adjacent_faces(uint edge_idx) const {
  return edge_adjacent_faces_[edge_idx];
}

//draw mesh in svg format
void Mesh2D::draw(std::string filename, bool drawfaces) const {
  using namespace std;

  std::ofstream of(filename.c_str());

  init_svg_file(of,uint(max_x_),uint(max_y_));

  std::string face_style = "fill:#eeeeee;stroke-width:0.03;stroke:black";
  std::string line_style = "stroke-width:0.1;stroke:black";

  if (drawfaces) {
    for (uint i=0; i < face_.size(); i++) {

      uint nEdges = uint( face_[i].edge_idx_.size() );

      std::vector<std::pair<double, double> > points;
      for (uint j=0; j < nEdges; j++) {

        std::map<uint,uint> point_count;

        uint edge_idx = face_[i].edge_idx_[j];
        uint next_edge_idx = face_[i].edge_idx_[(j+1) % nEdges];

        point_count[edge_[edge_idx].to_idx_]++;
        point_count[edge_[edge_idx].from_idx_]++;
        point_count[edge_[next_edge_idx].to_idx_]++;
        point_count[edge_[next_edge_idx].from_idx_]++;

        //std::cerr << "point count: " << point_count << std::endl;

        assert(point_count.size() == 3);
        uint end_point = MAX_UINT;
        uint start_point = MAX_UINT;

        for (std::map<uint,uint>::iterator it = point_count.begin(); it != point_count.end(); it++) {
          if (it->second == 2) {
            end_point = it->first;
            break;
          }
        }

        start_point = (edge_[edge_idx].to_idx_ == end_point) ? 
          edge_[edge_idx].from_idx_ : edge_[edge_idx].to_idx_;

        points.push_back(std::make_pair(point_[end_point].x_,point_[end_point].y_));
      }

      svg_draw_polygon(of,face_style,points);
    }
  }

  for (uint i=0; i < edge_.size(); i++) {
    Mesh2DPoint p1 = point_[edge_[i].from_idx_];
    Mesh2DPoint p2 = point_[edge_[i].to_idx_];
    if (!(p1 == p2)) {
      //Edge still used
      of << "<line style=\"" << line_style << "\"  x1=\"" << p1.x_ << "\" y1=\"" << p1.y_ 
        << "\" x2=\"" << p2.x_ << "\" y2=\"" << p2.y_ << "\" />" << std::endl;
    }
  }
  finish_svg_file(of);
}


//draw mesh in svg format
void Mesh2D::draw_labels(std::string filename, const int* labels) const {

  std::ofstream of(filename.c_str());

  init_svg_file(of,uint(max_x_),uint(max_y_));

  std::string fg_style = "fill:#eeeeee;stroke-width:0.01;stroke:black";
  std::string bg_style = "fill:#333333;stroke-width:0.01;stroke:black";
  std::string uk_style = "fill:#ff0000;stroke-width:0.01;stroke:black";

  for (uint i=0; i < face_.size(); i++) {

    uint nEdges = uint( face_[i].edge_idx_.size() );

    std::vector<std::pair<double, double> > points;

    //Should we draw this face?
    bool should_draw = false;

    for (uint j=0; j < nEdges; j++) {

      std::map<uint,uint> point_count;

      uint edge_idx = face_[i].edge_idx_[j];
      uint next_edge_idx = face_[i].edge_idx_[(j+1) % nEdges];

      std::vector<uint> adj_faces = adjacent_faces(edge_idx);
      for (uint k=0; k<adj_faces.size(); ++k) {
        if (labels[adj_faces[k]] != labels[i]) {
          should_draw = true;
        }
      }

      point_count[edge_[edge_idx].to_idx_]++;
      point_count[edge_[edge_idx].from_idx_]++;
      point_count[edge_[next_edge_idx].to_idx_]++;
      point_count[edge_[next_edge_idx].from_idx_]++;

      //std::cerr << "point count: " << point_count << std::endl;

      assert(point_count.size() == 3);
      uint end_point = MAX_UINT;
      uint start_point = MAX_UINT;

      for (std::map<uint,uint>::iterator it = point_count.begin(); it != point_count.end(); it++) {
        if (it->second == 2) {
          end_point = it->first;
          break;
        }
      }

      start_point = (edge_[edge_idx].to_idx_ == end_point) ? 
        edge_[edge_idx].from_idx_ : edge_[edge_idx].to_idx_;

      points.push_back(std::make_pair(point_[end_point].x_,point_[end_point].y_));
    }

    should_draw = true;
    if (should_draw || labels[i] <= 0) {
      if (labels[i] == 0) {
        svg_draw_polygon(of,bg_style,points);
      }
      else if (labels[i] < 0) {
        svg_draw_polygon(of,uk_style,points);
      }
      else {
        svg_draw_polygon(of,fg_style,points);
      }
    }
  }

  finish_svg_file(of);
}


//draw mesh in svg format
void Mesh2D::draw_labels_with_pairs(std::string filename, const double* labels, const std::vector<Mesh2DEdgePair>& edge_pairs, 
                                    double xDim, double yDim) const
{
  using namespace std;

  vector<string> graylevels(256);
  for (int i=0;i<256;++i) {
    stringstream sout;
    sout << setw(2) << setfill('0') << hex << i;
    graylevels[i] = sout.str();
  }

  ofstream of(filename.c_str());
  if (xDim < 0 || yDim < 0) {
    init_svg_file(of,uint(max_x_),uint(max_y_));
  }
  else {
    init_svg_file(of,uint(xDim),uint(yDim));
  }

  std::string line_style = "stroke-width:0.05;stroke:#4f4f4f;";

  for (uint i=0; i < face_.size(); i++) {

    uint nEdges = uint( face_[i].edge_idx_.size() );

    std::vector<std::pair<double, double> > points;

    //Should we draw this face?
    bool should_draw = false;

    for (uint j=0; j < nEdges; j++) {

      std::map<uint,uint> point_count;

      uint edge_idx = face_[i].edge_idx_[j];
      uint next_edge_idx = face_[i].edge_idx_[(j+1) % nEdges];

      std::vector<uint> adj_faces = adjacent_faces(edge_idx);
      for (uint k=0; k<adj_faces.size(); ++k) {
        if (labels[adj_faces[k]] != labels[i]) {
          should_draw = true;
        }
      }


      point_count[edge_[edge_idx].to_idx_]++;
      point_count[edge_[edge_idx].from_idx_]++;
      point_count[edge_[next_edge_idx].to_idx_]++;
      point_count[edge_[next_edge_idx].from_idx_]++;

      //std::cerr << "point count: " << point_count << std::endl;

      assert(point_count.size() == 3);
      uint end_point = MAX_UINT;
      uint start_point = MAX_UINT;

      for (std::map<uint,uint>::iterator it = point_count.begin(); it != point_count.end(); it++) {
        if (it->second == 2) {
          end_point = it->first;
          break;
        }
      }

      start_point = (edge_[edge_idx].to_idx_ == end_point) ? 
        edge_[edge_idx].from_idx_ : edge_[edge_idx].to_idx_;

      points.push_back(std::make_pair(point_[end_point].x_,point_[end_point].y_));
    }

    int color = int(labels[i]*255);
    if (color<0) color=0;
    if (color>255) color=255;

    //if (color < 255) {
    stringstream style;
    style << "fill:#" << graylevels[color] << graylevels[color] << graylevels[color] <<
      ";stroke-width:0.0;stroke:#7f7f7f7";
    svg_draw_polygon(of,style.str(),points);
    //}
  }

  for (uint i=0; i < edge_.size(); i++) {
    Mesh2DPoint p1 = point_[edge_[i].from_idx_];
    Mesh2DPoint p2 = point_[edge_[i].to_idx_];
    if (!(p1 == p2)) {
      //Edge still used
      of << "<line style=\"" << line_style << "\"  x1=\"" << p1.x_ << "\" y1=\"" << p1.y_ 
        << "\" x2=\"" << p2.x_ << "\" y2=\"" << p2.y_ << "\" />" << std::endl;
    }
  }

  for (size_t ind=0; ind < edge_pairs.size(); ind++) {
    for (int i=0; i<=1; ++i) {
      double val =labels[face_.size() + 2*ind+i];

      if (val > 0.1) {
        int color = int((1-val)*255);
        if (color<0) color=0;
        if (color>255) color=255;

        uint l1 = edge_pairs[ind].first_edge_idx_;
        Mesh2DPoint p1 = point_[edge_[l1].from_idx_];
        Mesh2DPoint p2 = point_[edge_[l1].to_idx_];
        of << "<g><line style=\"stroke-width:0.1;stroke:#ff" << graylevels[color] << graylevels[color] <<
          "\" x1=\"" << p1.x_ << "\" y1=\"" << p1.y_ <<
          "\" x2=\"" << p2.x_ << "\" y2=\"" << p2.y_ << "\" />" << std::endl;

        uint l2 = edge_pairs[ind].second_edge_idx_;
        p1 = point_[edge_[l2].from_idx_];
        p2 = point_[edge_[l2].to_idx_];
        of << "<line style=\"stroke-width:0.1;stroke:#ff" << graylevels[color] << graylevels[color] <<
          "\" x1=\"" << p1.x_ << "\" y1=\"" << p1.y_ <<
          "\" x2=\"" << p2.x_ << "\" y2=\"" << p2.y_ << "\" /></g>" << std::endl;
      }
    }
  }

  finish_svg_file(of);
}



int Mesh2D::match(uint face_idx, uint edge_idx) const {

  uint nEdges = uint( face_[face_idx].edge_idx_.size() );
  for (uint i=0; i < nEdges; i++) {

    if (face_[face_idx].edge_idx_[i] == edge_idx) {

      uint next_edge_idx = face_[face_idx].edge_idx_[(i+1) % nEdges];
      std::map<uint,uint> point_count;

      point_count[edge_[edge_idx].to_idx_]++;
      point_count[edge_[edge_idx].from_idx_]++;
      point_count[edge_[next_edge_idx].to_idx_]++;
      point_count[edge_[next_edge_idx].from_idx_]++;

      if (point_count[edge_[edge_idx].to_idx_] == 2)
        return 1;
      else if (point_count[edge_[edge_idx].from_idx_] == 2)
        return -1;
      else 
        return 0;
    }
  }

  return 0;
}


void Mesh2D::generate_edge_pair_list(std::vector<Mesh2DEdgePair>& vec) const {

  for (uint i=0; i < point_.size(); i++) {

    uint nEdges = uint( point_adjacent_edges_[i].size() );

    for (uint j=0; j < nEdges; j++) {
      for (uint j2 = j+1; j2 < nEdges; j2++)
        vec.push_back(Mesh2DEdgePair(point_adjacent_edges_[i][j],point_adjacent_edges_[i][j2],i));
    }
  }
}


std::vector<Mesh2DPoint>::const_iterator Mesh2D::point_start() const {
  return point_.begin();
}

std::vector<Mesh2DPoint>::const_iterator Mesh2D::point_end() const {
  return point_.end();
}

std::vector<Mesh2DEdge>::const_iterator Mesh2D::edge_start() const {
  return edge_.begin();
}

std::vector<Mesh2DEdge>::const_iterator Mesh2D::edge_end() const {
  return edge_.end();
}

std::vector<Mesh2DFace>::const_iterator Mesh2D::face_start() const {
  return face_.begin();
}

std::vector<Mesh2DFace>::const_iterator Mesh2D::face_end() const {
  return face_.end();
}

uint Mesh2D::nPoints() const {
  return uint( point_.size() );
}

uint Mesh2D::nEdges() const {
  return uint( edge_.size() );
}

uint Mesh2D::nFaces() const {
  return uint( face_.size() );
}

uint Mesh2D::find_point(Mesh2DPoint point) const {
  for (uint i=0; i < point_.size(); i++) {
    if (point_[i] == point)
      return i;
  }

  return MAX_UINT;
}


uint Mesh2D::find_or_add_point(Mesh2DPoint point) 
{
  uint ind = find_point(point);
  if (ind < MAX_UINT) {
    return ind;
  }
  else {
    add_point(point);
    return uint( point_.size()-1 );
  }
}


uint Mesh2D::find_or_add_edge(uint to, uint from)
{
  std::vector<uint>::const_iterator it;

  //
  // Check if this edge already exists
  //
  for (it = point_adjacent_edges_[to].begin(); it != point_adjacent_edges_[to].end(); it++) {
    uint i = *it;
    const Mesh2DEdge& edge = edge_[i];

    if (edge.from_idx_ == from && edge.to_idx_ == to) 
      return i;
    else if (edge.from_idx_ == to && edge.to_idx_ == from)
      return i;
  }

  //If not, add it
  add_edge(to, from);
  return uint( edge_.size() - 1 );
}

void Mesh2D::replace_edge(uint e1,uint e2)
{
  //Replace e1 with e2
  using namespace std;

  //Faces are now built up with e2 instead of e1
  vector<uint>::const_iterator it;
  for (it = edge_adjacent_faces_[e1].begin(); it != edge_adjacent_faces_[e1].end(); ++it) {
    uint f = *it;
    for (uint i=0;i<face_[f].edge_idx_.size();++i) {
      if (face_[f].edge_idx_[i] == e1) {
        face_[f].edge_idx_[i] = e2;
      }
    }
    edge_adjacent_faces_[e2].push_back(f);
  }

  //Remove from point_adjacent_edges
  uint p;
  vector<uint>::iterator invalid;
  p = edge_[e1].from_idx_;
  invalid = remove(point_adjacent_edges_[p].begin(), point_adjacent_edges_[p].end(), e1);
  point_adjacent_edges_[p].erase(invalid,point_adjacent_edges_[p].end());
  p = edge_[e1].to_idx_;
  invalid = remove(point_adjacent_edges_[p].begin(), point_adjacent_edges_[p].end(), e1);
  point_adjacent_edges_[p].erase(invalid,point_adjacent_edges_[p].end());

  //Remove from edge vector
  //edge_.erase(edge_.begin() + e1);
  edge_[e1].from_idx_ = 0;
  edge_[e1].to_idx_   = 0;
  edge_adjacent_faces_[e1].clear();
}

void Mesh2D::move_edge(uint p, uint p1, uint p2, uint e1, uint e2)
{
  // e1 : p->p1
  // e2 : p->p2
  // e1 is changed from p->p1 to p2->p1
  using namespace std;

  //Faces are now built up with both e1 and e2 instead of e1
  vector<uint>::const_iterator it;
  for (it = edge_adjacent_faces_[e1].begin(); it != edge_adjacent_faces_[e1].end(); ++it) {
    uint f = *it;
    vector<uint> newface;
    uint nEdges = uint( face_[f].edge_idx_.size() );

    for (uint j=0; j < nEdges; j++) {
      uint e = face_[f].edge_idx_[j];
      uint nexte = face_[f].edge_idx_[(j+1) % nEdges];
      if (e == e1) {
        if (edge_[nexte].from_idx_ == p1 || edge_[nexte].to_idx_ == p1) {
          newface.push_back(e2);
          newface.push_back(e1);
        }
        else {
          newface.push_back(e1);
          newface.push_back(e2);
        }
      }
      else {
        newface.push_back(e);
      }
    }
    face_[f].edge_idx_ = newface;

    //Add to edge_adjacent_faces
    edge_adjacent_faces_[e2].push_back(f);
  }

  //Change e1
  if (edge_[e1].from_idx_ == p) {
    edge_[e1].from_idx_ = p2;
    edge_[e1].to_idx_   = p1;
  }
  else {
    edge_[e1].from_idx_ = p1;
    edge_[e1].to_idx_   = p2;
  }

  //Remove from point_adjacent_edges
  vector<uint>::iterator invalid;
  invalid = remove(point_adjacent_edges_[p].begin(), point_adjacent_edges_[p].end(), e1);
  point_adjacent_edges_[p].erase(invalid,point_adjacent_edges_[p].end());

  //Add to point_adjacent_edges
  point_adjacent_edges_[p2].push_back(e1);
}

void Mesh2D::merge_parallel_edges()
{
  using namespace std;

  bool done;
  do {
    done = true;

    for (uint pi=0; pi < point_.size(); pi++) {
      Mesh2DPoint p = point_[pi];

      std::vector<uint>::const_iterator it1,it2;
      for (it1 = point_adjacent_edges_[pi].begin(); it1 != point_adjacent_edges_[pi].end(); it1++) {
        uint pi1 = edge_[*it1].from_idx_;
        if (pi1 == pi) pi1 = edge_[*it1].to_idx_;
        Mesh2DPoint p1 = point_[pi1];
        double dx1 = p.x_ - p1.x_;
        double dy1 = p.y_ - p1.y_;
        for (it2 = it1+1; it2 != point_adjacent_edges_[pi].end(); it2++) {
          uint pi2 = edge_[*it2].from_idx_;
          if (pi2 == pi) pi2 = edge_[*it2].to_idx_;
          Mesh2DPoint p2 = point_[pi2];
          double dx2 = p.x_ - p2.x_;
          double dy2 = p.y_ - p2.y_;
          //Do these two edges leaving p have the same slope?
          if ( abs(dx1*dy2 - dy1*dx2) < 1e-4 && dx1*dx2 > -1e-4 && dy1*dy2 > -1e-4 ) {
            //Are they the same edge?
            if (p1 == p2) {

              //cout << "(" << p.x_ << "," << p.y_ << ")->(" << p1.x_ << "," << p1.y_ << ") and " 
              //	 << "(" << p.x_ << "," << p.y_ << ")->(" << p2.x_ << "," << p2.y_ << ") are duplicates" << endl;

              //Remove the second edge
              replace_edge(*it2,*it1);
              done = false;
              goto endpoint;
            }
            //cout << "(" << p.x_ << "," << p.y_ << ")->(" << p1.x_ << "," << p1.y_ << ") and " 
            //		 << "(" << p.x_ << "," << p.y_ << ")->(" << p2.x_ << "," << p2.y_ << ") are parallel" << endl;

            //The longest of these edges is to be moved
            if ( dx1*dx1 + dy1*dy1 > dx2*dx2 + dy2*dy2 ) {
              move_edge(pi,pi1,pi2,*it1,*it2);
              done = false;
              goto endpoint;
            }
            else {
              move_edge(pi,pi2,pi1,*it2,*it1);
              done = false;
              goto endpoint;
            }
          }
        }
      }

endpoint: 
      ;
    }

  } while (!done);
}


//returns the index of the edge and a flag which is true if the edge was found and false
//if the reverse edge was found
std::pair<uint,bool> Mesh2D::find_edge(Mesh2DPoint to, Mesh2DPoint from) const {

  for (uint i=0; i < edge_.size(); i++) {
    if (point_[edge_[i].from_idx_] == from && point_[edge_[i].to_idx_] == to)
      return std::make_pair(i,true);
    else if (point_[edge_[i].from_idx_] == to && point_[edge_[i].to_idx_] == from)
      return std::make_pair(i,false);
  }

  return std::make_pair(MAX_UINT,false);
}

//returns the index of the edge and a flag which is true if the edge was found and false
//if the reverse edge was found
std::pair<uint,bool> Mesh2D::find_edge(uint to, uint from) const {

  uint i;
  std::vector<uint>::const_iterator it;

  for (it = point_adjacent_edges_[to].begin(); it != point_adjacent_edges_[to].end(); it++) {

    i = *it;
    const Mesh2DEdge& edge = edge_[i];

    if (edge.from_idx_ == from && edge.to_idx_ == to)
      return std::make_pair(*it,true);
    else if (edge.from_idx_ == to && edge.to_idx_ == from)
      return std::make_pair(*it,false);
  }

  std::cerr << "edge not found" << std::endl;
  return std::make_pair(MAX_UINT,false);
}


Mesh2DPoint Mesh2D::point(uint point_idx) const {
  return point_[point_idx];
}

Mesh2DEdge Mesh2D::edge(uint edge_idx) const {
  return edge_[edge_idx];
}


void Mesh2D::enlarge(double x_factor, double y_factor)
{
  for (uint i=0; i < point_.size(); ++i) {
    point_[i].x_ *= x_factor;
    point_[i].y_ *= y_factor;
  }
  min_x_ *= x_factor;
  max_x_ *= x_factor;
  min_y_ *= y_factor;
  max_y_ *= y_factor;
}


Mesh2DPoint Mesh2D::face_center(uint idx)
{
  std::vector<uint>::const_iterator it;
  Mesh2DPoint p(0,0);
  int n = 0;
  for (it = face_[idx].edge_idx_.begin(); it != face_[idx].edge_idx_.end(); it++) {
    Mesh2DEdge e = edge(*it);
    Mesh2DPoint p1 = point(e.from_idx_);
    Mesh2DPoint p2 = point(e.to_idx_);
    p.x_ += p1.x_;
    p.y_ += p1.y_;
    p.x_ += p2.x_;
    p.y_ += p2.y_;
    n++;
  }
  p.x_ /= (2*n);
  p.y_ /= (2*n);

  return p;
}
