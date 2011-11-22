/**** written by Thomas Schoenemann as an employee of Lund University, Sweden, 2010 *****/

#ifndef LP_SEGMENTER_HH
#define LP_SEGMENTER_HH

#include "segmentation_common.hh"
#include "matrix.hh"
#include "tensor.hh"
#include "sparse_matrix_description.hh"

#include <coin/ClpSimplex.hpp>
#include <coin/ClpPlusMinusOneMatrix.hpp>
#include <coin/CbcModel.hpp>
#include <coin/OsiClpSolverInterface.hpp>

#ifdef HAS_GUROBI
#include "gurobi_c++.h"
#endif
#ifdef HAS_XPRESS
#include "xprs.h" 
#endif


class LpSegmenter {
public:

  LpSegmenter(const Math2D::Matrix<float>& image, const  LPSegOptions& options, uint nRegions,
	      bool shared_boundary_vars = false);

  ~LpSegmenter();

  double segment(uint nIter=1);

  const Math2D::Matrix<uint>& segmentation();

  double curv_energy();

  double curv_icm();

protected:

  double point_energy(uint point, const std::vector<Mesh2DEdgePair>& edge_pairs, 
		      const NamedStorage1D<std::vector<uint> >& point_pair,
		      const NamedStorage1D<std::vector<uint> >& point_edge);

  const Math2D::Matrix<float>& image_;
  Math3D::Tensor<float> data_term_;

  Math1D::Vector<double> mean_;

  Mesh2D mesh_;  

  LPSegOptions options_;

  uint nRegions_;
  uint xDim_;
  uint yDim_;

  uint nVars_; 
  uint nSlacks_;
  uint nConstraints_;
  uint nEntries_;

  bool shared_boundary_vars_;

  Math1D::Vector<double> cost_;

  Math1D::Vector<double> var_ub_;
  Math1D::Vector<double> var_lb_;

  Math1D::Vector<double> rhs_upper_;
  Math1D::Vector<double> rhs_lower_;

  Storage1D<PixelFaceRelation> shares_;
  Math1D::Vector<uint> share_start_;

  Math1D::Vector<uint> simplex_start_;

  std::vector<Mesh2DEdgePair> edge_pairs_;

  CoinPackedMatrix coinMatrix_;

  OsiClpSolverInterface lpSolver_;
  SparseMatrixDescription<char> lp_descr_;

  Math1D::Vector<uint> integral_solution_;

  std::string solver_;

#ifdef HAS_GUROBI
  GRBenv*   grb_env_;
  GRBmodel* grb_model_;
#endif
#ifdef HAS_XPRESS
  XPRSprob xp_prob_;
#endif

  Math2D::Matrix<uint> segmentation_;
  Math3D::Tensor<double> output_;

};


#endif
