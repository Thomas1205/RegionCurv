/**** written by Petter Strandmark as an employee of Lund University, Sweden, 2010 ****/

#ifndef QPBO_SEGMENTATION_H
#define QPBO_SEGMENTATION_H

#include "matrix.hh"
#include "lp_segmentation.hh"

//solves a segmentation problem with length and curvature regularity via an LP
//@param lambda: the weight for the length term
//@param beta  : the weight for the curvature term  
void err_function(char * err);
double qpbo_segment_curvreg(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset, Math2D::Matrix<uint>& segmentation);
double qpbo_segment_curvreg_new(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset, Math2D::Matrix<uint>& segmentation);


#endif
