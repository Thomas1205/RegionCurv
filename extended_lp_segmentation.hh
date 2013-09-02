/*** written by Thomas Schoenemann as an employee of Lund University, Sweden, 2010 - 2011 ***/

#ifndef EXTENDED_LP_SEGMENTATION_HH
#define EXTENDED_LP_SEGMENTATION_HH

#include "matrix.hh"
#include "tensor.hh"
#include "segmentation_common.hh"

double lp_segment_pottscurvreg(const Math3D::Tensor<float>& data_term, const LPSegOptions& options, Math2D::Matrix<uint>& segmentation);

double multi_curv_icm(const Math3D::Tensor<float>& data_term, const LPSegOptions& options,
                      Math2D::Matrix<uint>& segmentation);

double factor_lp_segment_curvreg(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset, 
                                 Math2D::Matrix<uint>& segmentation, const Math2D::Matrix<int>* fixed_labels = 0);

double factor_lp_segment_pottscurvreg(const Math3D::Tensor<float>& data_term, const LPSegOptions& options, Math2D::Matrix<uint>& segmentation);

double factor_lp_segment_pottscurvreg_layered(const Math3D::Tensor<float>& data_term, const LPSegOptions& options, 
                                              Math2D::Matrix<uint>& segmentation);

double factor_lp_segment_curvreg_minsum_diffusion(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset, 
                                                  Math2D::Matrix<uint>& segmentation, uint nIter = 12000,
                                                  const Math2D::Matrix<int>* fixed_labels = 0);

double factor_lp_segment_curvreg_minsum_diffusion_memsave(const Math3D::Tensor<float>& data_term, 
                                                          const LPSegOptions& options, double energy_offset, 
                                                          Math2D::Matrix<uint>& segmentation, uint nIter = 12000,
                                                          const Math2D::Matrix<int>* fixed_labels = 0);

double factor_lp_segment_curvreg_message_passing(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset, 
                                                 Math2D::Matrix<uint>& segmentation, std::string method = "bp", uint nIter = 500, 
						 bool quiet = false, const Math2D::Matrix<int>* fixed_labels = 0, bool rescale=false);


double factor_curv_icm(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset,
                       Math2D::Matrix<uint>& segmentation);

double factor_curv_icm(const Math3D::Tensor<float>& data_term, const LPSegOptions& options,
                       Math2D::Matrix<uint>& segmentation);

#endif
