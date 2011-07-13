/******** written by Thomas Schoenemann as an employee of Lund University, Sweden, June 2010 ****/
/******** extended by Yubin Kuang as an employee of Lund University, Sweden, September 2010 ****/

#ifndef LP_INPAINTING_HH
#define LP_INPAINTING_HH

#include "matrix.hh"

double lp_inpaint(const Math2D::Matrix<float>& image, const Math2D::Matrix<float>& mask,
		  double lambda, double gamma, double curv_power, uint neighborhood, double energy_offset, std::string solver,
		  Math2D::Matrix<float>& inpainted_image, bool enforce_boundary_consistency,
		  bool enforce_region_edge_consistency, bool light_constraints, bool legacy = false);


double lp_inpaint_hybrid(const Math2D::Matrix<float>& image, const Math2D::Matrix<float>& mask,
			 double lambda, double gamma, double curv_power, uint neighborhood, uint nBin, double energy_offset, std::string solver,
			 Math2D::Matrix<float>& inpainted_image, bool enforce_boundary_consistency,
			 bool enforce_region_edge_consistency,bool enforce_level_consistency,  bool light_constraints, bool legacy = false);


#endif
