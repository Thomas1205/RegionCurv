/*** written by Thomas Schoenemann as an employee of Lund University, Sweden, August 2010 ***/

#ifndef CURVATURE_HH
#define CURVATURE_HH

#include "mesh2D.hh"

double curv_weight(const Mesh2D& mesh, const Mesh2DEdgePair& pair, double curv_power=2.0, bool bruckstein = false);

double curv_weight(double x1, double y1, double x2, double y2, double x3, double y3, double curv_power = 2.0, bool bruckstein = false);


#endif

