/*** first version written by Thomas Schoenemann as a private person without employment, September 2009 ***/
/*** extended by Thomas Schoenemann and Petter Strandmark as employees of Lund University, Sweden, September 2010 - 2011 ***/

#ifndef LP_SEGMENTATION_HH
#define LP_SEGMENTATION_HH

#include "matrix.hh"
#include "tensor.hh"
#include "segmentation_common.hh"

//solves a segmentation problem with length regularity via an LP
double lp_segment_lenreg(const Math2D::Matrix<float>& data_term, const LPSegOptions& options,
                         double energy_offset, Math2D::Matrix<uint>& segmentation, const Math2D::Matrix<int>* fixed_labels = 0);

//solves a segmentation problem with length and curvature regularity via an LP
double lp_segment_curvreg(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset, 
                          Math2D::Matrix<uint>& segmentation, const Math2D::Matrix<int>* fixed_labels = 0);


double lp_segment_curvreg_message_passing(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset, 
                                          Math2D::Matrix<uint>& segmentation, std::string method = "bp",
                                          const Math2D::Matrix<int>* fixed_labels = 0);


double curv_icm(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset,
                Math2D::Matrix<uint>& segmentation);

#endif
