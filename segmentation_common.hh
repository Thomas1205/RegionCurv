/**** written by Petter Strandmark as an employee of Lund University, Sweden, summer 2010 ***/

#ifndef SEGMENTATION_COMMON_HH
#define SEGMENTATION_COMMON_HH

#include "curvature.hh"
#include "mesh2D.hh"
#include "timing.hh"
#include "Petter-Color.hh"
#include "tensor.hh"
#include "matrix.hh"


struct LPSegOptions
{
  LPSegOptions();

  bool enforce_consistent_boundaries_;
  bool enforce_consistent_points_;
  bool enforce_regionedge_;
  bool prevent_crossings_;
  bool light_constraints_;
  
  int neighborhood_;
  int griddim_xDim_;
  int griddim_yDim_;
  int adaptive_mesh_n_;

  std::string solver_;
  
  std::string base_filename_;

  double lambda_;
  double gamma_;

  uint output_factor_;

  bool fix_regions_;
  
  bool save_matrix_;
  
  enum {Hex,Square} gridtype_;
  
  bool bruckstein_;

  bool debug_svg_;

  bool refine_;
};

//Generates a mesh with a given neighborhood structure
void generate_mesh(uint xDim, uint yDim, uint neighborhood, Mesh2D& mesh, bool silent=false,
		   const Math2D::Matrix<int>* seeds = 0);
//
// Generates a hexagonal mesh. Neighborhood is currently ignored.
//
void generate_hexagonal_mesh(uint xDim, uint yDim, double w, uint neighborhood, Mesh2D& mesh);
void add_hex_to_mesh(double x, double y, double w, Mesh2D& mesh);

//
// Adds the contribution from a region to an output image
//
void add_grid_output(uint region_index, double label, const Mesh2D& mesh, const Math2D::Matrix<double>& output);

void add_grid_output_mreg(uint face_index, const double* label, uint nRegions, const Mesh2D& mesh, const Math3D::Tensor<double>& output);

//
// Calculates the data term for a region given the data term in pixels
//
double calculate_data_term(uint region_index, Mesh2D& mesh, const Math2D::Matrix<float>& data_term);
double calculate_data_term(uint region_index, uint r, Mesh2D& mesh, const Math3D::Tensor<float>& data_term);

//
// Generates an adaptive square grid
//  Typical values of the neighborhood is 4,8 or 16
//
void generate_adaptive_mesh(const Math2D::Matrix<float>& data_term, Mesh2D& mesh, uint neighborhood, uint limit);
//
// Generates an adaptive hexagonal grid. Neighborhood is currently ignored
//
void generate_adaptive_hexagonal_mesh(const Math2D::Matrix<float>& data_term, Mesh2D& mesh, uint neighborhood, uint limit);

struct PixelFaceRelation {

  PixelFaceRelation();

  PixelFaceRelation(uint face_idx, float share);

  //the id of the pixel is stored elsewhere (to keep memory usage low)
  uint face_idx_;
  float share_; //how much of the face intersects with the pixel (a number between 0.0 and 1.0)
};

void compute_pixel_shares(const Mesh2D& mesh, uint xDim, uint yDim, Storage1D<PixelFaceRelation>& shares,
			  Math1D::Vector<uint>& share_start);


#endif
