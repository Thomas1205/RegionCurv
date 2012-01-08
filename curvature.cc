/*** written by Thomas Schoenemann as an employee of Lund University, Sweden, August 2010 ***/

#include "curvature.hh"

#include <cmath>

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

double pair_diff_angle(const Mesh2D& mesh, const Mesh2DEdgePair& pair) {
  const Mesh2DEdge e1 = mesh.edge(pair.first_edge_idx_);
  const Mesh2DEdge e2 = mesh.edge(pair.second_edge_idx_);

  const uint p2_idx = pair.common_point_idx_;
  const uint p1_idx = (e1.from_idx_ == p2_idx) ? e1.to_idx_ : e1.from_idx_;
  const uint p3_idx = (e2.from_idx_ == p2_idx) ? e2.to_idx_ : e2.from_idx_;

  const double d1x = mesh.point(p2_idx).x_ - mesh.point(p1_idx).x_;
  const double d1y = mesh.point(p2_idx).y_ - mesh.point(p1_idx).y_;

  const double d2x = mesh.point(p3_idx).x_ - mesh.point(p2_idx).x_;
  const double d2y = mesh.point(p3_idx).y_ - mesh.point(p2_idx).y_;

  const double angle1 = atan2(d1y,d1x);
  const double angle2 = atan2(d2y,d2x);

  double diff_angle = fabs(angle2-angle1);
  diff_angle = std::min(diff_angle, 2*M_PI - diff_angle);

  return diff_angle;
}

double curv_weight(const Mesh2D& mesh, const Mesh2DEdgePair& pair, double curv_power, bool bruckstein) {

  const Mesh2DEdge e1 = mesh.edge(pair.first_edge_idx_);
  const Mesh2DEdge e2 = mesh.edge(pair.second_edge_idx_);

  const uint p2_idx = pair.common_point_idx_;
  const uint p1_idx = (e1.from_idx_ == p2_idx) ? e1.to_idx_ : e1.from_idx_;
  const uint p3_idx = (e2.from_idx_ == p2_idx) ? e2.to_idx_ : e2.from_idx_;

  const double d1x = mesh.point(p2_idx).x_ - mesh.point(p1_idx).x_;
  const double d1y = mesh.point(p2_idx).y_ - mesh.point(p1_idx).y_;

  const double d2x = mesh.point(p3_idx).x_ - mesh.point(p2_idx).x_;
  const double d2y = mesh.point(p3_idx).y_ - mesh.point(p2_idx).y_;

  const double angle1 = atan2(d1y,d1x);
  const double angle2 = atan2(d2y,d2x);

  double diff_angle = fabs(angle2-angle1);
  diff_angle = std::min(diff_angle, 2*M_PI - diff_angle);

  double result = 0.0;

  if (bruckstein) {    
    const double l1 = sqrt(d1x*d1x + d1y*d1y);
    const double l2 = sqrt(d2x*d2x + d2y*d2y);
    diff_angle /= std::min(l1,l2);

    result = std::min(l1,l2) * pow(diff_angle,curv_power);
  }
  else
    result = pow(diff_angle,curv_power);

  return result;
}


double curv_weight(double x1, double y1, double x2, double y2, double x3, double y3, double curv_power, bool bruckstein) {

  const double d1x = x2 - x1;
  const double d1y = y2 - y1;

  const double d2x = x3 - x2;
  const double d2y = y3 - y2;

  const double angle1 = atan2(d1y,d1x);
  const double angle2 = atan2(d2y,d2x);

  double diff_angle = fabs(angle2-angle1);
  diff_angle = std::min(diff_angle, 2*M_PI - diff_angle);

  double result = 0.0;

  if (bruckstein) {    
    const double l1 = sqrt(d1x*d1x + d1y*d1y);
    const double l2 = sqrt(d2x*d2x + d2y*d2y);
    diff_angle /= std::min(l1,l2);

    result = std::min(l1,l2) * pow(diff_angle,curv_power);
  }
  else
    result = pow(diff_angle,curv_power);

  return result;
}
