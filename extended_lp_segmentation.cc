/*** written by Thomas Schoenemann as an employee of Lund University, Sweden, 2010 - 2011 ***/

#include "extended_lp_segmentation.hh"
#include "mesh2D.hh"
#include "grid.hh"
#include "sparse_matrix_description.hh"
#include "line_drawing.hh"
#include "timing.hh"
#include "curvature.hh"


#include "segmentation_common.hh"
#include "stl_out.hh"

#include "ClpSimplex.hpp"
#include "ClpPlusMinusOneMatrix.hpp"
#include "ClpFactorization.hpp"
#include "OsiClpSolverInterface.hpp"


#include "CbcModel.hpp"
#include "CglGomory/CglGomory.hpp"
#include "CglProbing/CglProbing.hpp"
#include "CglRedSplit/CglRedSplit.hpp"
#include "CglTwomir/CglTwomir.hpp"
#include "CglMixedIntegerRounding/CglMixedIntegerRounding.hpp"
#include "CglMixedIntegerRounding2/CglMixedIntegerRounding2.hpp"
#include "CglOddHole/CglOddHole.hpp"
#include "CglLandP/CglLandP.hpp"
#include "CglClique/CglClique.hpp"
#include "CglStored.hpp"

#include "conv_lp_solving.hh"

#define USE_PM_ONE_MATRIX

#ifdef HAS_VNK_GRAPH
#include "graph.h"
#endif

#ifdef HAS_GUROBI
#include "gurobi_c++.h"
#endif

#ifdef HAS_CPLEX
#include <ilcplex/cplex.h>
#endif

#ifdef HAS_XPRESS
#include "xprs.h" 
#endif


double lp_segment_pottscurvreg(const Math3D::Tensor<float>& data_term, const LPSegOptions& options, Math2D::Matrix<uint>& segmentation) {

  std::cerr << "CURVATURE POTTS MODEL" << std::endl;

  double lambda = options.lambda_;
  double gamma = options.gamma_;
  int neighborhood = options.neighborhood_;
  bool enforce_consistent_boundaries = options.enforce_consistent_boundaries_;
  bool enforce_regionedge = options.enforce_regionedge_;
  bool bruckstein = options.bruckstein_;
  
  if (options.light_constraints_) {
    std::cerr << "WARNING: light constraints currently not implemented" << std::endl;
  }

  uint light_factor = 2;

  assert(neighborhood <= 16); 

  uint xDim = uint( data_term.xDim() );
  uint yDim = uint( data_term.yDim() );
  uint nRegions = uint( data_term.zDim() );

  Mesh2D mesh;  
  if (options.gridtype_ == options.Square) {
    if (options.adaptive_mesh_n_ < 0) {
      double xfac = double(xDim) / double(options.griddim_xDim_);
      double yfac = double(yDim) / double(options.griddim_yDim_);
      generate_mesh( options.griddim_xDim_, options.griddim_yDim_, neighborhood, mesh, false, 0);
      mesh.enlarge(xfac,yfac);
    }
    else {
      //Adaptive mesh
      TODO("combination of curvature potts model and adaptive meshes");
      //generate_adaptive_mesh(data_term, mesh, neighborhood, options.adaptive_mesh_n_);
    }
  }
  else {

    if (options.adaptive_mesh_n_ < 0) {
      double xfac = double(xDim) / double(options.griddim_xDim_);
      double yfac = double(yDim) / double(options.griddim_yDim_);
      
      generate_hexagonal_mesh( xDim, yDim, 0.5*(xfac+yfac), neighborhood,mesh); //TODO: proper handling of the factor
    }
    else {
      //Adaptive mesh
      TODO("combination of curvature potts model and adaptive hexagonal meshes");
      //generate_adaptive_hexagonal_mesh(data_term, mesh, neighborhood, options.adaptive_mesh_n_);
    }
  }

  Storage1D<PixelFaceRelation> shares;
  Math1D::Vector<uint> share_start;
  compute_pixel_shares(mesh, xDim, yDim, shares, share_start);

  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  std::cerr << edge_pairs.size() << " edge pairs." << std::endl;

  uint nVars = uint( nRegions*mesh.nFaces() + 2*nRegions*edge_pairs.size() );
  uint nConstraints = mesh.nFaces() //simplex constraints
    + nRegions*3*mesh.nEdges();

  uint surface_con_offs = mesh.nFaces();
  uint boundary_con_offset = surface_con_offs + nRegions*mesh.nEdges();

  const uint consistency_con_offs = nConstraints;
  if (enforce_consistent_boundaries) {
    nConstraints += nRegions*2*mesh.nEdges();
  }
  const uint regionedge_constraints_offs = nConstraints;
  if (enforce_regionedge) {
    nConstraints += nRegions*4*mesh.nEdges();
  }

  Math1D::NamedVector<double> var_lb(nVars,0.0,MAKENAME(var_lb));
  Math1D::NamedVector<double> var_ub(nVars,1.0,MAKENAME(var_ub));
  Math1D::NamedVector<double> cost(nVars,0.0,MAKENAME(cost));
  
  for (uint i=0; i<mesh.nFaces(); ++i) {
    for (uint r = 0; r < nRegions; r++) 
      cost[i*nRegions+r] = calculate_data_term(i, r, mesh, data_term);
  }

  uint edge_pair_var_offs = nRegions*mesh.nFaces();

  for (uint j=0; j < edge_pairs.size(); j++) {

    uint first = edge_pairs[j].first_edge_idx_;
    uint nFirstAdjacent = uint( mesh.adjacent_faces(first).size() );
    
    uint second = edge_pairs[j].second_edge_idx_;
    uint nSecondAdjacent = uint( mesh.adjacent_faces(second).size() );

    double weight = 0.0;
    if (nFirstAdjacent > 1)
      weight += 0.5*lambda*mesh.edge_length(first);
    if (nSecondAdjacent > 1)
      weight += 0.5*lambda*mesh.edge_length(second);

    //do not penalize the image corners for their curvature
    if (nFirstAdjacent > 1 || nSecondAdjacent > 1)
      weight += gamma * curv_weight(mesh,edge_pairs[j],2.0,bruckstein);

    for (uint r = 0; r < nRegions; r++) {

      cost[edge_pair_var_offs+2*r*edge_pairs.size() + 2*j] = 0.5*weight;
      cost[edge_pair_var_offs+2*r*edge_pairs.size() + 2*j+1] = 0.5*weight;

      /*** check if (at the image border) one of the edge pairs is impossible. if so, set its upper bound to 0 ***/
      uint edge = first;
      if (mesh.adjacent_faces(edge).size() == 1) {
        int match = mesh.match(mesh.adjacent_faces(edge)[0],edge);

        uint y = uint( edge_pair_var_offs+2*r*edge_pairs.size()+2*j );

        if (edge_pairs[j].common_point_idx_ == mesh.edge(edge).from_idx_)
          match *= -1;

        if (match == -1)
          var_ub[y+1] = 0.0;
        else if (match == 1)
          var_ub[y] = 0.0;
      }

      edge = second;
      if (mesh.adjacent_faces(edge).size() == 1) {
        int match = mesh.match(mesh.adjacent_faces(edge)[0],edge);

        uint y = uint( edge_pair_var_offs+2*r*edge_pairs.size()+2*j );

        if (edge_pairs[j].common_point_idx_ == mesh.edge(edge).from_idx_)
          match *= -1;

        if (match == -1)
          var_ub[y] = 0.0;
        else if (match == 1)
          var_ub[y+1] = 0.0;
      }
    }
  }
  
  Math1D::NamedVector<double> rhs_lower(nConstraints,0.0,MAKENAME(rhs_lower));
  Math1D::NamedVector<double> rhs_upper(nConstraints,0.0,MAKENAME(rhs_upper));

  for (uint c=0; c < mesh.nFaces(); c++) {
    rhs_lower[c] = 1.0;
    rhs_upper[c] = 1.0;
  }

  std::cerr << "coding matrix" << std::endl;

  uint nEntries = uint( nRegions*(mesh.nFaces() + 2*mesh.nEdges() + 6*edge_pairs.size()) ); 
  if (enforce_consistent_boundaries) {
    nEntries += uint( 5*nRegions*edge_pairs.size() ); //we also allocate space for the slack variables used with the convex solver
  }
  if (enforce_regionedge) {
    nEntries += uint( nRegions*(16*edge_pairs.size() 
			  + 4*mesh.nEdges()) ); //we also allocate space for the slack variables used with the convex solver
  }

  SparseMatrixDescription<double> lp_descr(nEntries, nConstraints, nVars);

  /*** a) region simplex constraints ****/
  for (uint i=0; i < mesh.nFaces(); i++) {

    for (uint r=0; r < nRegions; r++)
      lp_descr.add_entry(i, i*nRegions + r, 1.0);
  }

  /*** b) surface continuation constraints ****/
  for (uint r=0; r < nRegions; r++) {

    uint cur_con_offs = surface_con_offs+r*mesh.nEdges();

    for (uint j=0; j < mesh.nEdges(); j++) {

      const std::vector<uint>& adjacent_faces = mesh.adjacent_faces(j);
      for (std::vector<uint>::const_iterator it = adjacent_faces.begin();
        it != adjacent_faces.end(); it++) {

          lp_descr.add_entry(cur_con_offs + j, nRegions*(*it)+r, mesh.match(*it,j));
      }
    }

    uint cur_edge_offs = uint( edge_pair_var_offs + 2*r*edge_pairs.size() );

    for (uint j=0; j < edge_pairs.size(); j++) {

      uint first_edge = edge_pairs[j].first_edge_idx_;
      uint second_edge = edge_pairs[j].second_edge_idx_;

      uint middle_point = edge_pairs[j].common_point_idx_;

      if (mesh.edge(first_edge).to_idx_ == middle_point) {
        lp_descr.add_entry(cur_con_offs+first_edge,cur_edge_offs+2*j,1);
      }
      else {
        lp_descr.add_entry(cur_con_offs+first_edge,cur_edge_offs+2*j,-1);
      }

      if (mesh.edge(second_edge).to_idx_ == middle_point) {
        lp_descr.add_entry(cur_con_offs+second_edge,cur_edge_offs+2*j+1,1);
      }
      else {
        lp_descr.add_entry(cur_con_offs+second_edge,cur_edge_offs+2*j+1,-1);
      }
    }

  }

  /*** c) boundary continuation constraints ****/
  for (uint j=0; j < edge_pairs.size(); j++) {

    uint first_edge = edge_pairs[j].first_edge_idx_;
    uint second_edge = edge_pairs[j].second_edge_idx_;

    uint middle_point = edge_pairs[j].common_point_idx_;

    for (uint r=0; r < nRegions; r++) {

      uint cur_edge_offs = uint( edge_pair_var_offs + 2*r*edge_pairs.size() );

      if (mesh.edge(first_edge).to_idx_ == middle_point) {      
        lp_descr.add_entry(boundary_con_offset + 2*(nRegions*first_edge+r), cur_edge_offs+2*j, 1);
        lp_descr.add_entry(boundary_con_offset + 2*(nRegions*first_edge+r)+1, cur_edge_offs+2*j+1, -1);
      }
      else {
        lp_descr.add_entry(boundary_con_offset + 2*(nRegions*first_edge+r)+1, cur_edge_offs+2*j, 1);
        lp_descr.add_entry(boundary_con_offset + 2*(nRegions*first_edge+r), cur_edge_offs+2*j+1, -1);
      }

      if (mesh.edge(second_edge).from_idx_ == middle_point) {
        lp_descr.add_entry(boundary_con_offset + 2*(nRegions*second_edge+r), cur_edge_offs+2*j, -1);
        lp_descr.add_entry(boundary_con_offset + 2*(nRegions*second_edge+r)+1, cur_edge_offs+2*j+1, 1);
      }
      else {
        lp_descr.add_entry(boundary_con_offset + 2*(nRegions*second_edge+r)+1, cur_edge_offs+2*j, -1);
        lp_descr.add_entry(boundary_con_offset + 2*(nRegions*second_edge+r), cur_edge_offs+2*j+1, 1);
      }
    }
  }

  uint nStandardEntries = lp_descr.nEntries();

  if (enforce_consistent_boundaries) {

    //constraints in words: for each oriented edge, the pairs that start with this oriented edge 
    // and the pairs that end in the oppositely oriented edge may not sum to more than 1.0
    // (i.e. they are mutually exclusive)

    for (uint c=consistency_con_offs; c < consistency_con_offs + 2*nRegions*mesh.nEdges(); c++) 
      rhs_upper[c] = 1.0;

    for (uint j=0; j < edge_pairs.size(); j++) {

      uint middle_point = edge_pairs[j].common_point_idx_;
      
      uint first_edge = edge_pairs[j].first_edge_idx_;
      uint second_edge = edge_pairs[j].second_edge_idx_;

      for (uint r = 0; r < nRegions; r++) {

	uint cur_con_offs = consistency_con_offs + 2*r*mesh.nEdges();
	uint cur_edge_offs = uint( edge_pair_var_offs + 2*r*edge_pairs.size() );

	if (mesh.edge(first_edge).to_idx_ == middle_point) {      
	  lp_descr.add_entry(cur_con_offs + 2*first_edge, cur_edge_offs+2*j, 1);
	  lp_descr.add_entry(cur_con_offs + 2*first_edge, cur_edge_offs+2*j+1, 1);
	}
	else {
	  lp_descr.add_entry(cur_con_offs + 2*first_edge+1, cur_edge_offs+2*j, 1);
	  lp_descr.add_entry(cur_con_offs + 2*first_edge+1, cur_edge_offs+2*j+1, 1);
	}
	
	if (mesh.edge(second_edge).to_idx_ == middle_point) {
	  lp_descr.add_entry(cur_con_offs + 2*second_edge, cur_edge_offs+2*j, 1);
	  lp_descr.add_entry(cur_con_offs + 2*second_edge, cur_edge_offs+2*j+1, 1);
	}
	else {
	  lp_descr.add_entry(cur_con_offs + 2*second_edge+1, cur_edge_offs+2*j, 1);
	  lp_descr.add_entry(cur_con_offs + 2*second_edge+1, cur_edge_offs+2*j+1, 1);
	}
      }
    }
  }
  if (enforce_regionedge) {

    uint rowoff = regionedge_constraints_offs;

    for (uint edge=0; edge < mesh.nEdges(); edge++) {
      
      //Get the two adjacent faces
      if (mesh.adjacent_faces(edge).size() != 2) {
	  //One of the edges is at the border of the image
      
	continue;
      }
      
      uint x1 = mesh.adjacent_faces(edge)[0];
      uint x2 = mesh.adjacent_faces(edge)[1];

      for (uint r=0; r < nRegions; r++) {

	uint cur_row = rowoff + 2*light_factor*r*mesh.nEdges() + 2*light_factor*edge; 
	
	lp_descr.add_entry(cur_row , nRegions*x1+r, 1);
	lp_descr.add_entry(cur_row , nRegions*x2+r, 1);
	lp_descr.add_entry(cur_row+1 , nRegions*x1+r, -1);
	lp_descr.add_entry(cur_row+1, nRegions*x2+r, -1);
	lp_descr.add_entry(cur_row+2, nRegions*x1+r, 1);
	lp_descr.add_entry(cur_row+2, nRegions*x2+r, 1);
	lp_descr.add_entry(cur_row+3, nRegions*x1+r, -1);
	lp_descr.add_entry(cur_row+3, nRegions*x2+r, -1);

	//NOTE (important for the construction with slacks: the binding constraints are in both cases the upper bounds)
	rhs_lower[cur_row]   = 0;
	rhs_upper[cur_row]   = 2;
	rhs_lower[cur_row+1] = -2;
	rhs_upper[cur_row+1] = 0;
	
	rhs_lower[cur_row+2] = 0;
	rhs_upper[cur_row+2] = 2;
	rhs_lower[cur_row+3] = -2;
	rhs_upper[cur_row+3] = 0;
      }
    }

    for (uint r = 0; r < nRegions; r++) {

      uint cur_edge_offs = uint( edge_pair_var_offs + 2*r*edge_pairs.size() );

      for (uint j=0; j < edge_pairs.size(); j++) {
	uint first = edge_pairs[j].first_edge_idx_;
	uint second = edge_pairs[j].second_edge_idx_;	
	
	uint y = cur_edge_offs + 2*j;
	
	uint edge = first;
	if (mesh.adjacent_faces(edge).size() == 2) {
	  
	  uint cur_row = rowoff + 2*light_factor*r*mesh.nEdges() + 2*light_factor*edge; 

	  lp_descr.add_entry(cur_row   ,y  , 1);
	  lp_descr.add_entry(cur_row+1 ,y  , 1);
	  lp_descr.add_entry(cur_row+2 ,y+1, 1);
	  lp_descr.add_entry(cur_row+3 ,y+1, 1);
	}
      
	edge = second;
	if (mesh.adjacent_faces(edge).size() == 2) {
	
	  uint cur_row = rowoff + 2*light_factor*r*mesh.nEdges() + 2*light_factor*edge; 

	  lp_descr.add_entry(cur_row   ,y+1, 1);
	  lp_descr.add_entry(cur_row+1 ,y+1, 1);
	  lp_descr.add_entry(cur_row+2 ,y  , 1);
	  lp_descr.add_entry(cur_row+3 ,y  , 1);
	}
      }
    }
  }

  Math1D::Vector<uint> row_start(nConstraints+1);
  lp_descr.sort_by_row(row_start);

  CoinPackedMatrix coinMatrix(false,(int*) lp_descr.row_indices(),(int*) lp_descr.col_indices(),
			      lp_descr.value(),lp_descr.nEntries());

  int error = 0;
  
  std::clock_t tStartCLP,tEndCLP;

  ClpSimplex lpSolver;
  lpSolver.loadProblem (coinMatrix, var_lb.direct_access(), var_ub.direct_access(),   
			cost.direct_access(), rhs_lower.direct_access(), rhs_upper.direct_access());

  coinMatrix.cleanMatrix();

  //lpSolver.writeMps("curv.mps");

  tStartCLP = std::clock();

  lpSolver.dual();

  const double* lp_solution = lpSolver.primalColumnSolution();


  //ClpSolve solve_options;
  //solve_options.setSolveType(ClpSolve::useDual);
  //solve_options.setSolveType(ClpSolve::useBarrier);
  //solve_options.setPresolveType(ClpSolve::presolveNumber,5);
  //lpSolver.initialSolve(solve_options);

  error = 1 - lpSolver.isProvenOptimal();

  tEndCLP = std::clock();

  if (error != 0)
    std::cerr << "!!!!!!!!!!!!!!LP-solver failed!!!!!!!!!!!!!!!!!!!" << std::endl;

  std::cerr << "CLP-time: " << diff_seconds(tEndCLP,tStartCLP) << " seconds. " << std::endl;

  uint nFrac = 0;
  double energy = 0.0;
  for (uint i=0; i < nVars; i++) {
    double val = lp_solution[i];
    energy += cost[i] * val;

    if (val > 0.01 && val < 0.99)
      nFrac++;
  }
  
  std::cerr << nFrac << "/" << nVars << " variables are fractional" << std::endl;

  //TODO: threshold and re-run the solver.

  /**** generate output ****/

  Math1D::Vector<int> labels(mesh.nFaces(),-1);
  
  uint out_factor = options.output_factor_;

  segmentation.resize(xDim*out_factor,yDim*out_factor);

  Math3D::Tensor<double> output(xDim*out_factor,yDim*out_factor,nRegions,0.0);
  mesh.enlarge(out_factor,out_factor);

  //re-compute pixel shares for the now larger mesh
  compute_pixel_shares(mesh, out_factor*xDim, out_factor*yDim, shares, share_start);

  uint seg_factor = 255 / nRegions;

  for (uint y=0; y < yDim*out_factor; y++) {
    for (uint x=0; x < xDim*out_factor; x++) {

      Math1D::Vector<double> votes(nRegions);

      for (uint r=0; r < nRegions; r++) {
	double sum = 0.0;
	for (uint k= share_start[y*(xDim*out_factor)+x]; k < share_start[y*(xDim*out_factor)+x+1]; k++) {
	  uint face = shares[k].face_idx_;
	  double cur_share = shares[k].share_;
	  sum += mesh.convex_area(face) * cur_share * lp_solution[face*nRegions + r];
	}

	votes[r] = sum;
      }

      double max_val = 0.0;
      uint arg_max = MAX_UINT;

      for (uint r=0; r < nRegions; r++) {

	double val = votes[r];
	if (val > max_val) {
	  max_val = val;
	  arg_max = r;
	}
      }

      segmentation(x,y) = arg_max * seg_factor;
    }
  }

  if (mesh.nFaces() > 20000) {
    //statusFailed();
  }
  else {
    mesh.draw_labels_with_pairs("mesh_lp.svg",lp_solution,edge_pairs,xDim,yDim);
  }

  return energy; 
}


double edge_ilp_segment_curvreg(const Math2D::Matrix<float>& data_term, double lambda, double gamma, 
				uint neighborhood, double energy_offset, Math2D::Matrix<uint>& segmentation,
				bool bruckstein) {

  bool node_constraints = true;
  bool edge_constraints = true;

  const uint xDim = uint( data_term.xDim() );
  const uint yDim = uint( data_term.yDim() );

  Math2D::NamedMatrix<double> hsum(xDim,yDim,0.0,MAKENAME(hsum));

  for (uint y=0; y < yDim; y++) {

    double sum = 0.0;

    for (uint x=0; x < xDim; x++) {

      sum += data_term(x,y);
      hsum(x,y) = sum;
    }
  }

  assert(neighborhood <= 8); // Stockes' theorem is currently only correctly applied for an 8-connectivity or less

  /******* generate grid ********/

  Grid grid(xDim+1,yDim+1);
 
  for (uint y=0; y <= yDim; y++) {
    for (uint x=0; x <= xDim; x++) {

      if (x > 0)
	grid.add_line(x,y,x-1,y);
      if (x+1 <= xDim)
	grid.add_line(x,y,x+1,y);
      if (y > 0)
	grid.add_line(x,y,x,y-1);
      if (y+1 <= yDim)
	grid.add_line(x,y,x,y+1);

      if (neighborhood >= 8) {

	if (x > 0 && y > 0) 
	  grid.add_line(x,y,x-1,y-1);
	if (x > 0 && y+1 <= yDim)
	  grid.add_line(x,y,x-1,y+1);
	if (x+1 <= xDim && y > 0)
	  grid.add_line(x,y,x+1,y-1);
	if (x+1 <= xDim && y+1 <= yDim)
	  grid.add_line(x,y,x+1,y+1);
      }
      if (neighborhood >= 16) {
	
	if (x > 1 && y > 0)
	  grid.add_line(x,y,x-2,y-1);
	if (x > 1 && y+1 < yDim)
	  grid.add_line(x,y,x-2,y+1);

	if (x+2 < xDim && y > 0)
	  grid.add_line(x,y,x+2,y-1);
	if (x+2 < xDim && y+1 < yDim)
	  grid.add_line(x,y,x+2,y+1);

	if (x > 0 && y > 1)
	  grid.add_line(x,y,x-1,y-2);
	if (x > 0 && y+2 < yDim)
	  grid.add_line(x,y,x-1,y+2);

	if (x+1 < xDim && y > 1)
	  grid.add_line(x,y,x+1,y-2);
	if (x+1 < xDim && y+2 < yDim)
	  grid.add_line(x,y,x+1,y+2);
      }
    }
  }

  uint nLines = grid.nLines();

  grid.generate_grid_line_pairs();
  
  uint nLinePairs = grid.nLinePairs();

  Math1D::NamedVector<float> line_data_term(nLines,0.0,MAKENAME(line_data_term));

  for (uint l=0; l < nLines; l++) {
    
    const GridLine& cur_line = grid.get_line(l);

    float dx = float(((int) cur_line.x2_) - ((int) cur_line.x1_) );
    float dy = float(((int) cur_line.y2_) - ((int) cur_line.y1_) );

    //std::cerr << "dx: " << dx << std::endl;
    //std::cerr << "dy: " << dy << std::endl;

    //float edge_length = sqrt(dx*dx+dy*dy);
    //float nx = -dy;
    //float ny = dx;

    uint min_x = std::min(cur_line.x2_, cur_line.x1_);
    uint min_y = std::min(cur_line.y2_, cur_line.y1_);

    float cost = 0.0;
    
    if (fabs(dx) == 1.0 && fabs(dy) == 1.0) {
      //diagonal line -> take the gradient of the pixel center

      if (min_x > 0)
	cost += float(sign(dy) * (hsum(min_x-1,min_y) + 0.5 * data_term(min_x,min_y)));
    }
    else if (dx == 0.0) {

      if (min_x > 0)
	cost += float(sign(dy) * hsum(min_x-1,min_y));
    }
    else if (dy == 0.0) {
      //no costs according to the current application of Stokes' theorem
    }

    line_data_term[l] = cost;
  }
  
  uint nVars = nLinePairs;
  uint nConstraints = nLines; // flow conservation

  const uint node_con_offs = nConstraints;

  if (node_constraints)
    nConstraints += (xDim+1)*(yDim+1);

  const uint edge_intersec_con_offs = nConstraints;

  if (edge_constraints) {
    if (neighborhood >= 8)
      nConstraints += xDim*yDim;
    assert(neighborhood <= 8);
    //TODO: 16-connect.
  }


  Math1D::NamedVector<double> var_lb(nVars,0.0,MAKENAME(var_lb));
  Math1D::NamedVector<double> var_ub(nVars,1.0,MAKENAME(var_ub));
  
  Math1D::NamedVector<double> cost(nVars,0.0,MAKENAME(cost)); //to be recomputed in every iteration
  
  Math1D::NamedVector<double> rhs_lower(nConstraints,0.0,MAKENAME(rhs_lower));
  Math1D::NamedVector<double> rhs_upper(nConstraints,0.0,MAKENAME(rhs_upper));

  uint nEntries = 2*nLinePairs; // for the flow conservation

  if (node_constraints)
    nEntries += nLinePairs;

  if (edge_constraints) {
    if (neighborhood >= 8)
      nEntries += 2*nVars;
    //TODO: 16-connect
    assert(neighborhood <= 8);
  }
  
  SparseMatrixDescription<double> lp_descr(nEntries, nConstraints, nVars);

  for (uint k=0; k < nLinePairs; k++) {

    const GridLinePair& pair = grid.get_line_pair(k);

    uint l1 = pair.line1_;
    uint l2 = pair.line2_;

    double sum = 0.5 * (line_data_term[l1] + line_data_term[l2]);

    const GridLine& line1 = grid.get_line(l1);
    const GridLine& line2 = grid.get_line(l2);

    double x1 = line1.x1_;
    double y1 = line1.y1_;
    double x2 = line1.x2_;
    double y2 = line1.y2_;
    double x3 = line2.x2_;
    double y3 = line2.y2_;
    
    assert(x2 == line2.x1_);
    assert(y2 == line2.y1_);

    double cw = curv_weight(x1, y1, x2, y2, x3, y3, 2.0, bruckstein);

    sum += gamma*cw + 0.5 * lambda * (line1.length() + line2.length());
    
    cost[k] = sum;
  }

  for (uint l=0; l < nLines; l++) {

    const GridLine& cur_line = grid.get_line(l);

    const std::vector<uint>& inpairs = cur_line.ending_line_pairs_;

    for (uint k=0; k < inpairs.size(); k++)
      lp_descr.add_entry(l,inpairs[k], 1.0);

    const std::vector<uint>& outpairs = cur_line.starting_line_pairs_;

    for (uint k=0; k < outpairs.size(); k++)
      lp_descr.add_entry(l,outpairs[k], -1.0);
  }

  if (node_constraints) {

    std::cerr << "adding node constraints" << std::endl;

    for (uint y = 0; y <= yDim; y++) {
      for (uint x = 0; x <= xDim; x++) {

	const uint row = node_con_offs + y*(xDim+1)+x;

	rhs_lower[row] = 0.0;
	rhs_upper[row] = 1.0;

	const GridNode& cur_node = grid.get_node(x,y);
	
	for (uint l_in = 0; l_in < cur_node.incoming_lines_.size(); l_in++) {
	  
	  const GridLine& cur_line = grid.get_line(cur_node.incoming_lines_[l_in]);

	  for (uint p=0; p < cur_line.ending_line_pairs_.size(); p++) {

	    lp_descr.add_entry(row, cur_line.ending_line_pairs_[p], 1.0);
	  }
	}
      }
    }
  }

  if (edge_constraints) {

    if (neighborhood >= 8) {

      for (uint y=0; y < yDim; y++) {
	for (uint x=0; x < xDim; x++) {
	  
	  uint row = edge_intersec_con_offs + y*xDim+x;

	  rhs_lower[row] = 0.0;
	  rhs_upper[row] = 2.0;
	  
	  for (uint i=0; i < 4; i++) {

	    uint cur_line_num = 0;
	    
	    if (i==0) 
	      cur_line_num = grid.find_line(x,y,x+1,y+1);
	    else if (i == 1)
	      cur_line_num = grid.find_line(x+1,y+1,x,y);
	    else if (i == 2)
	    cur_line_num = grid.find_line(x+1,y,x,y+1);
	    else if (i == 3)
	      cur_line_num = grid.find_line(x,y+1,x+1,y);
	    
	    const GridLine& cur_line = grid.get_line(cur_line_num);
	    
	    for (uint p=0; p < cur_line.starting_line_pairs_.size(); p++) 
	      lp_descr.add_entry(row, cur_line.starting_line_pairs_[p], 1.0);
	    for (uint p=0; p < cur_line.ending_line_pairs_.size(); p++) 
	      lp_descr.add_entry(row, cur_line.ending_line_pairs_[p], 1.0);
	  }
	}
      }
    }

  }

  Math1D::Vector<uint> row_start(nConstraints+1);
  lp_descr.sort_by_row(row_start);

#ifdef USE_SCIP

  SCIP* scip_env;
  SCIP_CALL(SCIPcreate(&scip_env) );
  SCIP_CALL(SCIPincludeDefaultPlugins(scip_env));
  SCIP_CALL(SCIPcreateProb(scip_env,"JointSeg",0,0,0,0,0,0));

  SCIP_CALL( SCIPsetIntParam(scip_env,"separating/poolfreq",1)  );
  SCIP_CALL( SCIPsetIntParam(scip_env,"separating/cmir/freq",1)  );

  Storage1D<SCIP_VAR*> scip_var(nVars,0);
  for (uint v=0; v < nVars; v++) {
    SCIP_CALL(SCIPcreateVar(scip_env,scip_var.direct_access()+v, "", var_lb[v], var_ub[v], cost[v],
			    SCIP_VARTYPE_BINARY, true, false, 0,0,0,0));

    SCIP_CALL(SCIPaddVar(scip_env,scip_var[v]));
  }

  Storage1D<SCIP_CONS*> scip_con(nConstraints,0);
  for (uint c=0; c < nConstraints; c++) {

    uint row_size = row_start[c+1] - row_start[c];

    Storage1D<SCIP_VAR*> cur_vars(row_size);
    for (uint v=0; v < row_size; v++)
      cur_vars[v] = scip_var[lp_descr.col_indices()[row_start[c] + v] ];
    
    SCIP_CALL( SCIPcreateConsLinear(scip_env, scip_con.direct_access() + c, "", row_size, 
				    cur_vars.direct_access(), lp_descr.value() + row_start[c],
				    rhs_lower[c],rhs_upper[c], true, true, true, true, true, false, false, false, false, false)  );

    SCIP_CALL( SCIPaddCons(scip_env, scip_con[c]) );
  }


  SCIP_CALL_ABORT( SCIPtransformProb(scip_env) );

  SCIP_CALL_ABORT( SCIPsetIntParam(scip_env,"display/verblevel",5) );

  SCIP_CALL( SCIPsolve(scip_env) );


  double best_energy = -100.0; //TODO: get correct value

  for (uint v=0; v < nVars; v++)
    SCIP_CALL( SCIPreleaseVar(scip_env,scip_var.direct_access()+v) );

  for (uint c=0; c < nConstraints; c++) 
    SCIP_CALL( SCIPreleaseCons(scip_env,scip_con.direct_access()+c) );

  SCIPfree(&scip_env);

#else
#ifdef USE_GUROBI

  std::cerr << "trying the gurobi solver" << std::endl;

  Math1D::Vector<double> grb_solution(nVars,0.0);

  GRBenv   *grb_env   = NULL;
  GRBmodel *grb_model = NULL;

  /* Create environment */
  int error;

  error = GRBloadenv(&grb_env,NULL);
  assert (!error && grb_env != NULL);

  //GRBsetintparam(grb_env, GRB_INT_PAR_LPMETHOD, GRB_METHOD_BARRIER);
  //GRBsetdblparam(grb_env, "BarConvTol", 1e-10);
  //GRBsetintparam(grb_env, "Crossover", 2);
  //GRBsetintparam(grb_env, "CrossoverBasis", 1);
  //GRBsetintparam(grb_env, "Presolve", 1);
  GRBsetintparam(grb_env, "PrePasses", 1);
  
  error = GRBsetdblparam(grb_env,"NodeLimit",100.0);
  assert (!error);

  //error = GRBsetintparam(grb_env,"PreCrush",1);
  //assert(!error);

  /* Create an empty model */

  error = GRBnewmodel(grb_env, &grb_model, "edge-curv-ilp", 0, NULL, NULL, NULL, NULL, NULL);
  assert(!error);

  Storage1D<char> vtype(nVars,GRB_BINARY);
  
  error = GRBaddvars(grb_model,nVars,0,NULL,NULL,NULL,cost.direct_access(),var_lb.direct_access(),
		     var_ub.direct_access(),vtype.direct_access(),NULL);
  assert(!error);

  error = GRBupdatemodel(grb_model);
  assert(!error);

  for (uint c=0; c < nConstraints; c++) {

    if (rhs_lower[c] == rhs_upper[c]) {

      error = GRBaddconstr(grb_model, row_start[c+1]-row_start[c], ((int*) lp_descr.col_indices()) + row_start[c], 
			   lp_descr.value() + row_start[c], GRB_EQUAL, rhs_lower[c], NULL);
    }
    else {
      
      error = GRBaddrangeconstr(grb_model, row_start[c+1]-row_start[c], ((int*) lp_descr.col_indices()) + row_start[c], 
				lp_descr.value() + row_start[c], rhs_lower[c], rhs_upper[c], NULL);      
    }

    if (error) {

      std::cerr << "error for con " << c << ": " << error << " " << GRBgeterrormsg(grb_env) << std::endl;
      error = GRBupdatemodel(grb_model);
      assert(!error);

      for (uint k=row_start[c]; k < row_start[c+1]; k++) {

	std::cerr << lp_descr.col_indices()[k] << ", ";
      }
      std::cerr << std::endl;
      assert(!error);
    }
  }

  error = GRBupdatemodel(grb_model);
  assert(!error);

  /* Optimize model */
  //GRBsetcallbackfunc( grb_model, grb_jointseg_callback, &c_info);

  error = GRBoptimize(grb_model);
  assert(!error);

  for (uint v=0; v < nVars; v++)
    GRBgetdblattrelement(grb_model,"X",v, grb_solution.direct_access()+v);

  const double* best_solution = grb_solution.direct_access();

  double best_energy = 0.0; 
  for (uint v=0; v < nVars; v++)
    best_energy += grb_solution[v] * cost[v];

  GRBfreemodel(grb_model);
  GRBfreeenv(grb_env);

#else

  CoinPackedMatrix coinMatrix(false,(int*) lp_descr.row_indices(),(int*) lp_descr.col_indices(),
			      lp_descr.value(),lp_descr.nEntries());

  OsiClpSolverInterface lpSolver;
  lpSolver.loadProblem(coinMatrix, var_lb.direct_access(), var_ub.direct_access(),   
		       cost.direct_access(), rhs_lower.direct_access(), rhs_upper.direct_access());
  
  coinMatrix.cleanMatrix();

  for (uint v=0; v < nVars; v++) {
    lpSolver.setInteger(v);
  }
  
  lpSolver.resolve();

  const double* lp_solution = lpSolver.getColSolution();

  double lp_energy = lpSolver.getObjValue();

  /******* check for fractional solutions ****/
  uint nNonInt = 0;
  for (uint v=0; v < nVars; v++) {

    double val = lp_solution[v];

    if (val > 0.01 && val < 0.99)
      nNonInt++;
  }
  
  double best_energy = lp_energy;
  const double* best_solution = lp_solution;

  std::cerr << nNonInt << " of " << nVars << " variables are fractional" << std::endl;

  lpSolver.messageHandler()->setLogLevel(0);
  
  CbcModel cbc_model(lpSolver);

  cbc_model.solver()->setIntParam(OsiMaxNumIterationHotStart,100);

  if (nNonInt > 0) {
    //if (false) {

    lpSolver.findIntegersAndSOS(false);
    lpSolver.setupForRepeatedUse();
    
    //cbc_model.messageHandler()->setLogLevel(0);
    
    CglGomory gomory_cut;
    gomory_cut.setLimit(50);
    gomory_cut.setLimitAtRoot(50);    
    gomory_cut.setAway(0.01);
    gomory_cut.setAwayAtRoot(0.01);
    //cbc_model.addCutGenerator(&gomory_cut,0,"Gomory Cut");

    CglProbing probing_cut;
    probing_cut.setUsingObjective(true);

//     probing_cut.setMaxPass(10);
//     probing_cut.setMaxPassRoot(50);
//     probing_cut.setMaxProbe(100);
//     probing_cut.setMaxProbeRoot(500);
//     probing_cut.setMaxLook(150);
//     probing_cut.setMaxLookRoot(1500);

    probing_cut.setMaxPass(5);
    probing_cut.setMaxPassRoot(15);
    probing_cut.setMaxProbe(40);
    //probing_cut.setMaxProbe(10);
    probing_cut.setMaxProbeRoot(100);
    //probing_cut.setMaxProbeRoot(20);
    probing_cut.setMaxLook(50);
    probing_cut.setMaxLookRoot(500);
    //probing_cut.setMaxLookRoot(50);

    //cbc_model.addCutGenerator(&probing_cut,0,"Probing Cut");
    
    CglRedSplit redsplit_cut;
    redsplit_cut.setLimit(1500);
    cbc_model.addCutGenerator(&redsplit_cut,0,"RedSplit Cut");

    CglMixedIntegerRounding2 mi2_cut;
    //cbc_model.addCutGenerator(&mi2_cut,0,"Mixed Integer 2");

    CglTwomir twomir_cut;
    //cbc_model.addCutGenerator(&twomir_cut,0,"Twomir Cut");

    CglLandP landp_cut;
    //cbc_model.addCutGenerator(&landp_cut,0,"LandP Cut");

    CglOddHole oddhole_cut;
    //cbc_model.addCutGenerator(&oddhole_cut,0,"OddHole Cut");

    CglClique clique_cut;
    //cbc_model.addCutGenerator(&clique_cut,0,"Clique Cut");

    CglStored stored_cut;
    //cbc_model.addCutGenerator(&stored_cut,0,"Stored Cut");

    //DO NOT USE THIS! Apparently incorrect!!!
    //cbc_model.findCliques(false,2,200);

    
    cbc_model.setMaximumNodes(100);

    cbc_model.branchAndBound();

    
    const double* cbc_solution = cbc_model.bestSolution();

    double energy = 0.0;

    for (uint v=0; v < nVars; v++) {
      
      double var_val = cbc_solution[v];
      
      energy += cost[v] * var_val;
    }

    best_solution = cbc_solution;
    best_energy = energy;
  }


#endif
#endif

  segmentation.set_constant(0);

  uint nFrac = 0;
  uint nActive = 0;

  /****  draw the resulting lines ****/
  for (uint v=0; v < nVars; v++) {

    double val = best_solution[v];

    if (val > 0.01 && val < 0.99)
      nFrac++;

    if (val > 0.01) {

      nActive++;

      const GridLinePair& pair = grid.get_line_pair(v);

      uint l1 = pair.line1_;
      uint l2 = pair.line2_;

      const GridLine& line1 = grid.get_line(l1);
      const GridLine& line2 = grid.get_line(l2);
      
      double x1 = line1.x1_;
      double y1 = line1.y1_;
      double x2 = line1.x2_;
      double y2 = line1.y2_;
      double x3 = line2.x2_;
      double y3 = line2.y2_;

      draw_line(segmentation,uint(x1),uint(y1),uint(x2),uint(y2),(uint) 1);
      draw_line(segmentation,uint(x2),uint(y2),uint(x3),uint(y3),(uint) 1);
    }
  }

  std::cerr << nActive << " variables are active" << std::endl;

  if (nFrac > 0)
    std::cerr << "WARNING: integer solution is actually fractional" << std::endl;

  std::cerr << "total energy: " << (best_energy + energy_offset) << std::endl;

  return best_energy + energy_offset; 
}

double clique_curv_energy(Mesh2D& mesh, const std::vector<Mesh2DEdgePair>& edge_pairs,
			  const Math2D::Matrix<float>& data_term, const LPSegOptions& options,
			  const Math1D::Vector<uint>& face_label) {

  double energy = 0.0;

  for (uint k=0; k < face_label.size(); k++) {

    if (face_label[k] == 1)
      energy += calculate_data_term(k, mesh, data_term);
  }


  double lambda = options.lambda_;
  double gamma = options.gamma_;
  int neighborhood = options.neighborhood_;
  bool bruckstein = options.bruckstein_;

  // double bin_terms = 0.0;
  // double tern_terms = 0.0;
  // double tern_a_terms = 0.0;
  // double tern_b_terms = 0.0;
  // double tern_c_terms = 0.0;
  // double fo_terms = 0.0;

  uint nTernFactors = 0;

  // process factors
  for (uint p=0; p < edge_pairs.size(); p++) {
    
    const uint first_edge = edge_pairs[p].first_edge_idx_;
    const uint second_edge = edge_pairs[p].second_edge_idx_;    
    
    double mul1 = (mesh.edge(first_edge).to_idx_ == edge_pairs[p].common_point_idx_) ? 1.0 : -1.0;
    double mul2 = (mesh.edge(second_edge).from_idx_ == edge_pairs[p].common_point_idx_) ? 1.0 : -1.0;
    
    const std::vector<uint>& adj1 = mesh.adjacent_faces(first_edge);
    const std::vector<uint>& adj2 = mesh.adjacent_faces(second_edge);

    double curvature_weight = gamma * curv_weight(mesh,edge_pairs[p],2.0,bruckstein);

    if (adj1.size() == 1 && adj2.size() == 1)
      curvature_weight = 0.0; //don't penalize the image corners

    //TODO: treat length separately as two-cliques

    double boundary_cost = curvature_weight;

    // no length-cost for the image border
    if (adj1.size() > 1)
      boundary_cost += 0.5*lambda*mesh.edge_length(first_edge);
    if (adj2.size() > 1)
      boundary_cost += 0.5*lambda*mesh.edge_length(second_edge);
    
    double pair_cost = boundary_cost;

    if (pair_cost == 0)
      continue;
    
    std::set<uint> adjacent_regions_set;
    std::map<uint,double> match1;
    std::map<uint,double> match2;
    
    for (std::vector<uint>::const_iterator it = adj1.begin(); it != adj1.end(); it++) {
      adjacent_regions_set.insert(*it);
      match1[*it] = mesh.match(*it,first_edge);
    }
    
    for (std::vector<uint>::const_iterator it = adj2.begin(); it != adj2.end(); it++) {
      adjacent_regions_set.insert(*it);
      match2[*it] = mesh.match(*it,second_edge);
    }
    
    Math1D::Vector<uint> adjacent_regions(adjacent_regions_set.size());
    uint k=0;
    for (std::set<uint>::iterator it = adjacent_regions_set.begin(); it != adjacent_regions_set.end(); it++) {
      adjacent_regions[k] = *it;
      k++;
    }      
    
    //add cost for the factor
    if (adjacent_regions.size() == 4) {

      assert(adj1.size() == 2);
      assert(adj2.size() == 2);

      Storage1D< Math3D::Tensor<float> > cost(2);
      for (uint i=0; i < 2; i++)
	cost[i].resize(2,2,2,0.0);
      
      cost[0](1,0,1) = pair_cost;
      cost[1](0,1,0) = pair_cost;

      uint v1 = adj1[0];
      uint v2 = adj1[1];
      
      if (match1[v1]*mul1 != 1) {
	assert(match1[v1]*mul1 == -1);
	std::swap(v1,v2);
      }

      uint v3 = adj2[0];
      uint v4 = adj2[1];
      
      if (match2[v3]*mul2 != 1) {
	assert(match2[v3]*mul2 == -1);
	std::swap(v3,v4);
      }
      
      energy += cost[face_label[v1]](face_label[v2],face_label[v3],face_label[v4]);
      //fo_terms += cost[face_label[v1]](face_label[v2],face_label[v3],face_label[v4]);
    }
    else if (adjacent_regions.size() == 3) {

      nTernFactors++;

      Math3D::Tensor<float> cost(2,2,2,0.0);
      
      if (adj1.size() == 1) {
	//boundary edge pair
	
	uint v1 = adj1[0];
	double m1 = match1[v1]*mul1;

	uint v2 = adj2[0];
	uint v3 = adj2[1];

	if (match2[v2]*mul2 != m1)
	  std::swap(v2,v3);

	cost(1,1,0) = pair_cost;

	energy += cost(face_label[v1],face_label[v2],face_label[v3]);
	//tern_terms += cost(face_label[v1],face_label[v2],face_label[v3]);
	//tern_a_terms += cost(face_label[v1],face_label[v2],face_label[v3]);
      }
      else if (adj2.size() == 1) {

	//boundary edge pair
	
	uint v1 = adj2[0];
	double m2 = match2[v1]*mul2;
	
	uint v2 = adj1[0];
	uint v3 = adj1[1];

	if (match1[v2]*mul1 != m2)
	    std::swap(v2,v3);

	cost(1,1,0) = pair_cost;

	energy += cost(face_label[v1],face_label[v2],face_label[v3]);
	//tern_terms += cost(face_label[v1],face_label[v2],face_label[v3]);
	//tern_b_terms += cost(face_label[v1],face_label[v2],face_label[v3]);
      }
      else {
	//interior edge pair with a region that borders both edges

	//std::cerr << "C" << std::endl;

	//find shared regions
	uint shared = MAX_UINT;
	uint v1 = MAX_UINT;
	uint v2 = MAX_UINT;
	
	for (uint k=0; k < 3; k++) {
	  uint v = adjacent_regions[k];
	  if (std::find(adj1.begin(), adj1.end(), v) != adj1.end() 
	      && std::find(adj2.begin(), adj2.end(), v) != adj2.end())
	    shared = v;
	  else if (v1 == MAX_UINT)
	    v1 = v;
	  else
	    v2 = v;
	}
	
	cost(0,1,1) = pair_cost;
	cost(1,0,0) = pair_cost;
	
	assert(shared < MAX_UINT);
	assert(v1 < MAX_UINT);
	assert(v2 < MAX_UINT);	  
      	
	energy += cost(face_label[shared],face_label[v1],face_label[v2]);
	//tern_terms += cost(face_label[shared],face_label[v1],face_label[v2]);
	//tern_c_terms += cost(face_label[shared],face_label[v1],face_label[v2]);
      }
    }
    else if (adjacent_regions.size() == 2) {

      Math2D::Matrix<float> cost(2,2,0.0);

      if (adj1.size() == 1 && adj2.size() == 1) {

	double val1 = match1[adj1[0]] * mul1;
	double val2 = match2[adj2[0]] * mul2;
	
	if (val1 == val2) { //both regions are on the same side of the line pair
	  cost(1,1) = pair_cost;
	}
	else {
	  std::cerr << "should this happen???" << std::endl;
	  cost(0,1) = pair_cost;
	  cost(1,0) = pair_cost;
	}
      }
      else {

	//there should only be one shared region
	assert(adj1.size() == 1 || adj2.size() == 1);

	if (find(adj1.begin(),adj1.end(),adjacent_regions[0]) != adj1.end()
	    && find(adj2.begin(),adj2.end(),adjacent_regions[0]) != adj2.end()) {
	  cost(1,0) = pair_cost;
	}
	else
	  cost(0,1) = pair_cost;
      }

      energy += cost(face_label[adjacent_regions[0]],face_label[adjacent_regions[1]]);
      //bin_terms += cost(face_label[adjacent_regions[0]],face_label[adjacent_regions[1]]);
    }
  }

  // std::cerr << "binary terms: " << bin_terms << std::endl;
  // std::cerr << "ternary terms: " << tern_terms << "  (" << tern_a_terms 
  // 	    << "," << tern_b_terms << "," << tern_c_terms << ")"<< std::endl;
  // std::cerr << "4th order terms: " << fo_terms << std::endl;

  // std::cerr << nTernFactors << " factors" << std::endl;

  return energy;
}


double clique_lp_segment_curvreg(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset, 
				 Math2D::Matrix<uint>& segmentation, const Math2D::Matrix<int>* fixed_labels) {
  
  //Get the options
  double lambda = options.lambda_;
  double gamma = options.gamma_;
  int neighborhood = options.neighborhood_;
  bool bruckstein = options.bruckstein_;
  
  assert(neighborhood <= 16); 

  uint xDim = uint( data_term.xDim() );
  uint yDim = uint( data_term.yDim() );

  Mesh2D mesh;  
  if (options.gridtype_ == options.Square) {
    if (options.adaptive_mesh_n_ < 0) {
      double xfac = double(xDim) / double(options.griddim_xDim_); 
      double yfac = double(yDim) / double(options.griddim_yDim_); 
      generate_mesh( options.griddim_xDim_, options.griddim_yDim_, neighborhood, mesh);
      mesh.enlarge(xfac,yfac);
    }
    else {
      std::cerr << "adaptive square mesh" << std::endl;

      //Adaptive mesh
      generate_adaptive_mesh(data_term, mesh, neighborhood, options.adaptive_mesh_n_);
    }
  }
  else {
    if (options.adaptive_mesh_n_ < 0) {
      double xfac = double(xDim) / double(options.griddim_xDim_); 
      double yfac = double(yDim) / double(options.griddim_yDim_); 
      generate_hexagonal_mesh( xDim, yDim, 0.5*(xfac+yfac), neighborhood,mesh); //TODO: proper handling of non-square images
    }
    else {
      std::cerr << "adaptive hex. mesh" << std::endl;

      //Adaptive mesh
      generate_adaptive_hexagonal_mesh(data_term, mesh, neighborhood, options.adaptive_mesh_n_);
    }
  }
  
  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  std::cerr << edge_pairs.size() << " edge pairs." << std::endl;

  uint nVars = mesh.nFaces() + 2*mesh.nEdges();  //region vars and length cliques
  uint nConstraints = mesh.nEdges();
  uint nEntries = 4*mesh.nEdges();

  Math1D::Vector<uint> pair_var_base(edge_pairs.size());
  Math1D::Vector<uint> pair_con_base(edge_pairs.size());

  double curv_threshold = 0.001;

  uint nQuadCliques = 0;
  uint nTripleCliques = 0;
  uint nDoubleCliques = 0;
  for (uint p=0; p < edge_pairs.size(); p++) {

    pair_var_base[p] = nVars;
    pair_con_base[p] = nConstraints;

    uint first_edge = edge_pairs[p].first_edge_idx_;
    uint second_edge = edge_pairs[p].second_edge_idx_;    

    const std::vector<uint>& adj1 = mesh.adjacent_faces(first_edge);
    const std::vector<uint>& adj2 = mesh.adjacent_faces(second_edge);

    double curv_cost = curv_weight(mesh,edge_pairs[p],2.0,bruckstein);

    if (adj1.size() == 1 && adj2.size() == 1)
      curv_cost = 0.0; // don't penalize the image corners

    if (curv_cost < curv_threshold)
      continue;

    std::set<uint> adjacent_regions;

    for (std::vector<uint>::const_iterator it = adj1.begin(); it != adj1.end(); it++)
      adjacent_regions.insert(*it);

    for (std::vector<uint>::const_iterator it = adj2.begin(); it != adj2.end(); it++)
      adjacent_regions.insert(*it);

    if (adjacent_regions.size() == 4) {
      nVars += 16;
      nConstraints += 4*2;
      nEntries += 16*4 + 4*2;

      nQuadCliques++;
    }
    else if (adjacent_regions.size() == 3) {
      nVars += 8;
      nConstraints += 3*2;
      nEntries += 8*3 + 3*2;
      
      nTripleCliques++;
    }
    else if (adjacent_regions.size() == 2) {
      nVars += 2; //here we use the construction via absolutes
      nConstraints += 1; 
      nEntries += 4;
      nDoubleCliques++;
    }
    
    //NOTE: even a count of 1 is possible, at the image corners. 
    // The cost can then be directly added to the cost of the face, so there is no need for extra variables
  }

  Math1D::Vector<double> var_lb(nVars,0.0);
  Math1D::Vector<double> var_ub(nVars,1.0);
  
  Math1D::Vector<double> cost(nVars,0.0);
  
  Math1D::Vector<double> rhs(nConstraints,0.0);

  for (uint i=0; i<mesh.nFaces(); ++i) {
    cost[i] = calculate_data_term(i, mesh, data_term);
    
    if (options.fix_regions_) {
      //Fix region variables
      if (cost[i] < 0) {
	var_lb[i] = 1.0;
      }
      else if (cost[i] > 0) {
	var_ub[i] = 0.0;
      }
      cost[i] = 0.0;
    }
  }

  if (fixed_labels != 0) {
    Storage1D<PixelFaceRelation> shares;
    Math1D::Vector<uint> share_start;
    compute_pixel_shares(mesh, xDim, yDim, shares, share_start);

    bool has_warned = false;

    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x < xDim; x++) {

	int fixed = (*fixed_labels)(x,y);
	if (fixed == 0) {
	  for (uint v=share_start[y*xDim+x]; v < share_start[y*xDim+x+1]; v++) {
	    uint f = shares[v].face_idx_;
	    float area_share = std::min(1.0f, shares[v].share_);
	    if (area_share >= 0.95) {
	      var_ub[f] = 0.0;
	    }
	    else if (!has_warned) {
	      std::cerr << Petter::RED << "WARNING: ignored a partial seed region" << Petter::NORMAL << std::endl;
	      has_warned = true;
	    }
	  }
	}
	else if (fixed == 1) {
	  for (uint v=share_start[y*xDim+x]; v < share_start[y*xDim+x+1]; v++) {
	    uint f = shares[v].face_idx_;
	    float area_share = std::min(1.0f, shares[v].share_);
	    if (area_share >= 0.95)
	      var_lb[f] = 1.0;
	    else if (!has_warned) {
	      std::cerr << Petter::RED << "WARNING: ignored a partial seed region" << Petter::NORMAL << std::endl;
	      has_warned = true;
	    }
	  }
	}
      }
    }
  }


  for (uint v=mesh.nFaces(); v < mesh.nFaces()+2*mesh.nEdges(); v++) {
    if (mesh.adjacent_faces((v-mesh.nFaces())/2).size() == 2)
      cost[v] = lambda;
  }
  
  SparseMatrixDescription<double> lp_descr(nEntries, nConstraints, nVars);

  for (uint e = 0; e < mesh.nEdges(); e++) {

    if (mesh.adjacent_faces(e).size() == 2) { 
      for (uint f=0; f < mesh.adjacent_faces(e).size(); f++) {

	uint var = mesh.adjacent_faces(e)[f];
	lp_descr.add_entry(e, var, mesh.match(var,e));
      }
      
      lp_descr.add_entry(e,mesh.nFaces()+2*e,1.0);
      lp_descr.add_entry(e,mesh.nFaces()+2*e+1,-1.0);
    }
    else {
      //don't penalize the image border
    }
  }
  
  /*** joint calculation of cost and constraint system ***/
  for (uint p=0; p < edge_pairs.size(); p++) {

    const uint first_edge = edge_pairs[p].first_edge_idx_;
    const uint second_edge = edge_pairs[p].second_edge_idx_;    

    double mul1 = (mesh.edge(first_edge).to_idx_ == edge_pairs[p].common_point_idx_) ? 1.0 : -1.0;
    double mul2 = (mesh.edge(second_edge).from_idx_ == edge_pairs[p].common_point_idx_) ? 1.0 : -1.0;

    const std::vector<uint>& adj1 = mesh.adjacent_faces(first_edge);
    const std::vector<uint>& adj2 = mesh.adjacent_faces(second_edge);

    double curvature_weight = gamma * curv_weight(mesh,edge_pairs[p],2.0,bruckstein);

    if (adj1.size() == 1 && adj2.size() == 1)
      curvature_weight = 0.0;

    if ((curvature_weight / gamma) < curv_threshold)
      continue; //length-based cost are handled via two-cliques

    double boundary_cost = curvature_weight;

    std::set<uint> adjacent_regions_set;
    std::map<uint,double> match1;
    std::map<uint,double> match2;

    for (std::vector<uint>::const_iterator it = adj1.begin(); it != adj1.end(); it++) {
      adjacent_regions_set.insert(*it);
      match1[*it] = mesh.match(*it,first_edge);
    }

    for (std::vector<uint>::const_iterator it = adj2.begin(); it != adj2.end(); it++) {
      adjacent_regions_set.insert(*it);
      match2[*it] = mesh.match(*it,second_edge);
    }

    const uint cur_var_base = pair_var_base[p];
    const uint cur_con_base = pair_con_base[p];

    Math1D::Vector<uint> adjacent_regions(adjacent_regions_set.size());
    uint k=0;
    for (std::set<uint>::iterator it = adjacent_regions_set.begin(); it != adjacent_regions_set.end(); it++) {
      adjacent_regions[k] = *it;
      k++;
    }

    if (adjacent_regions.size() == 4) {

      //note: in the absence of a direction change it is probably advisable to drop this clique
      // into introduce a two-clique for every edge, reflecting its length cost

      for (uint f=0; f < 4; f++) {
	
	lp_descr.add_entry(cur_con_base + 2*f, adjacent_regions[f], -1.0);
	lp_descr.add_entry(cur_con_base + 2*f+1, adjacent_regions[f], 1.0);
	rhs[cur_con_base + 2*f+1] = 1.0;
      }
      
      for (uint v_addon = 0; v_addon < 16; v_addon++) {

	double sum_e1 = 0.0;
	double sum_e2 = 0.0;

	for (uint f=0; f < 4; f++) {

	  uint mask = 1 << f;
	  assert(mask != 0);

	  //uint offs = (v_addon & mask) >> f;
	  uint offs = 1 - ((v_addon & mask) >> f);
	  lp_descr.add_entry(cur_con_base + 2*f+offs, cur_var_base + v_addon, 1.0);

	  sum_e1 += (1-offs) * match1[adjacent_regions[f]];
	  sum_e2 += (1-offs) * match2[adjacent_regions[f]];	  
	}

	sum_e1 *= mul1;
	sum_e2 *= mul2;

	if (fabs(sum_e1) > 0.5 && fabs(sum_e2) > 0.5 && sum_e1*sum_e2 > 0.0) {
	  cost[cur_var_base + v_addon] = boundary_cost;
	}
      }
    }
    else if (adjacent_regions.size() == 3) {
      
      for (uint f=0; f < 3; f++) {
	
	lp_descr.add_entry(cur_con_base + 2*f, adjacent_regions[f], -1.0);
	lp_descr.add_entry(cur_con_base + 2*f+1, adjacent_regions[f], 1.0);
	rhs[cur_con_base + 2*f+1] = 1.0;
      }
      
      for (uint v_addon = 0; v_addon < 8; v_addon++) {

	double sum_e1 = 0.0;
	double sum_e2 = 0.0;

	for (uint f=0; f < 3; f++) {
	
	  uint mask = 1 << f;
	  assert(mask != 0);

	  //uint offs = (v_addon & mask) >> f;
	  uint offs = 1 - ((v_addon & mask) >> f);
	  lp_descr.add_entry(cur_con_base + 2*f+offs, cur_var_base + v_addon, 1.0);
	  
	  sum_e1 += (1-offs) * match1[adjacent_regions[f]];
	  sum_e2 += (1-offs) * match2[adjacent_regions[f]];
	}

	sum_e1 *= mul1;
	sum_e2 *= mul2;

	if (fabs(sum_e1) > 0.5 && fabs(sum_e2) > 0.5 && sum_e1*sum_e2 > 0.0) {

	  //at the image boundary there can be triple cliques that do not share a region
	  if (adj1.size() > 1 && adj2.size() > 1) {
	    if (! (sum_e1*sum_e2 > 0.0)) {

	      const Mesh2DEdgePair& cur_edge_pair = edge_pairs[p];

	      std::cerr << "edge pair: " << cur_edge_pair << std::endl;
	      std::cerr << "first edge: " << mesh.edge(cur_edge_pair.first_edge_idx_) << std::endl;
	      std::cerr << "second edge: " << mesh.edge(cur_edge_pair.second_edge_idx_) << std::endl;
	      
	      std::cerr << "match1: " << match1 << std::endl;
	      std::cerr << "match2: " << match2 << std::endl;

	      std::cerr << "sum_e1: " << sum_e1 << std::endl;
	      std::cerr << "sum_e2: " << sum_e2 << std::endl;
	      
	      std::cerr << "mul1: " << mul1 << std::endl;
	      std::cerr << "mul2: " << mul2 << std::endl;
	      
	      std::cerr << "v_addon: " << v_addon << std::endl;
	    }

	    assert(sum_e1*sum_e2 > 0.0);
	  }
	  
	  cost[cur_var_base + v_addon] = boundary_cost;
	}
      }
    }
    else if (adjacent_regions.size() == 2) {

      double c01 = 0.0;
      double c10 = 0.0;

      if (adj1.size() == 1 && adj2.size() == 1) {

	double val1 = match1[adj1[0]] * mul1;
	double val2 = match2[adj2[0]] * mul2;
	
	if (val1 == val2) { //both regions are on the same side of the line pair
	  assert(boundary_cost == 0.0);
	}
	else {
	  std::cerr << "should this happen???" << std::endl;
	  c01 = boundary_cost;
	  c10 = boundary_cost;
	}
      }
      else {

	//there should only be one shared region
	assert(adj1.size() == 1 || adj2.size() == 1);

	if (find(adj1.begin(),adj1.end(),adjacent_regions[0]) != adj1.end()
	    && find(adj2.begin(),adj2.end(),adjacent_regions[0]) != adj2.end()) {
	  c10 = boundary_cost;
	}
	else
	  c01 = boundary_cost;
      }

      if (c01 != 0.0 || c10 != 0.0) {
	//here we use the construction via the absolute differences
	lp_descr.add_entry(cur_con_base, adjacent_regions[0], 1.0);
	lp_descr.add_entry(cur_con_base, adjacent_regions[1], -1.0);
	
	lp_descr.add_entry(cur_con_base, cur_var_base, 1.0);
	lp_descr.add_entry(cur_con_base, cur_var_base+1, -1.0);

	cost[cur_var_base] = c01;
	cost[cur_var_base+1] = c10;
      }
    }
  }

  Math1D::Vector<uint> row_start(nConstraints+1);
  lp_descr.sort_by_row(row_start,true);

#ifdef USE_OWN_CONV

  Math1D::Vector<double> conv_solution = var_lb;

  for (uint f=0; f < mesh.nFaces(); f++) {

    if (var_lb[f] != var_ub[f]) {

      if (cost[f] < 0.0)
	conv_solution[f] = 1.0;
    }
  }

  for (uint p=0; p < edge_pairs.size(); p++) {

    const uint first_edge = edge_pairs[p].first_edge_idx_;
    const uint second_edge = edge_pairs[p].second_edge_idx_;    

    std::set<uint> adjacent_regions_set;

    const std::vector<uint>& adj1 = mesh.adjacent_faces(first_edge);
    for (std::vector<uint>::const_iterator it = adj1.begin(); it != adj1.end(); it++) {
      adjacent_regions_set.insert(*it);
    }

    const std::vector<uint>& adj2 = mesh.adjacent_faces(second_edge);
    for (std::vector<uint>::const_iterator it = adj2.begin(); it != adj2.end(); it++) {
      adjacent_regions_set.insert(*it);
    }

    const uint cur_var_base = pair_var_base[p];

    Math1D::Vector<uint> adjacent_regions(adjacent_regions_set.size());
    uint k=0;
    for (std::set<uint>::iterator it = adjacent_regions_set.begin(); it != adjacent_regions_set.end(); it++) {
      adjacent_regions[k] = *it;
      k++;
    }

    //set induced clique order variables
    if (adjacent_regions.size() >= 3) {
      TODO("set solution");
    }
    else if (adjacent_regions.size() == 2) {
      double diff = conv_solution[adjacent_regions[0]] - conv_solution[adjacent_regions[1]];
      if (diff < 0.0)
	conv_solution[cur_var_base] = -diff;
      else
	conv_solution[cur_var_base+1] = diff;
    }
  }
  
  //std::cerr << "WARNING: the inequality constraints are currently not properly handled by the convex solver" << std::endl;
  uint nSlacks = 0;

  double stepsize_coeff = 0.5;

  std::cerr << "calling main routine" << std::endl;

  std::clock_t tStartConv,tEndConv;
  tStartConv = std::clock();

  eq_constrained_lp_solving_auglagrange_nesterov(nVars, nConstraints, cost.direct_access(), var_lb.direct_access(), var_ub.direct_access(),
 						 lp_descr, rhs.direct_access(), conv_solution.direct_access(), 
						 gamma*100.0, stepsize_coeff, 1500, 15, 1.25);
  tEndConv = std::clock();

  std::cerr << "convex solver needed " << diff_seconds(tEndConv,tStartConv) << " seconds. " << std::endl;
#endif

  CoinPackedMatrix coinMatrix(false,(int*) lp_descr.row_indices(),(int*) lp_descr.col_indices(),
			      lp_descr.value(),lp_descr.nEntries());

  ClpSimplex lpSolver;
  lpSolver.loadProblem(coinMatrix, var_lb.direct_access(), var_ub.direct_access(),   
		       cost.direct_access(), rhs.direct_access(), rhs.direct_access());

  coinMatrix.cleanMatrix();

  for (uint v=0; v < nVars; v++) {
    lpSolver.setInteger(v);
  }
  
  lpSolver.dual();

  const double* lp_solution = lpSolver.getColSolution();

  double lp_energy = lpSolver.getObjValue();

  std::cerr.precision(10);
  std::cerr << "original relaxed energy: " << (lp_energy + energy_offset) << std::endl;

  Math1D::Vector<uint> face_uint_solution(mesh.nFaces(),0);
  
  uint nNonInt = 0;
  for (uint v=0; v < mesh.nFaces(); v++) {
    
    double val = lp_solution[v];

    if (val > 0.01 && val < 0.99)
      nNonInt++;
    
    if (val < 0.5) 
      var_ub[v] = 0.0;
    else {
      var_lb[v] = 1.0;
      face_uint_solution[v] = 1;
    }
  }
  
  double energy = lp_energy;

  double disc_energy = clique_curv_energy(mesh, edge_pairs, data_term, options, face_uint_solution)
    + energy_offset;

  std::cerr << "solution energy found by discrete routine: " << disc_energy << std::endl;
  
  if (nNonInt > 0) {
    
    std::cerr << nNonInt << " non-integral region variables. Computing thresholded solution" << std::endl;
    
    //fix the region variables (according to the above thresholding) and re-run the solver
    // to get the boundary variables
    
    for (uint i=0; i < mesh.nFaces(); i++)
      lpSolver.setColumnBounds(i,var_lb[i],var_ub[i]);
    
    lpSolver.dual();
    
    double thresh_energy = lpSolver.getObjValue();
    uint nNonIntThresh = 0;
    for (uint i=0; i < nVars; i++) {
      double val = lpSolver.primalColumnSolution()[i]; 
      if (val > 0.01 && val < 0.99)
	nNonIntThresh++;
    }
    
    std::cerr << "energy of thresholded solution: " << (thresh_energy + energy_offset)
	      << "  (" << nNonIntThresh << " non-integral variables)" << std::endl;

    double bin_terms = 0.0;
    double tern_terms = 0.0;
    double fo_terms = 0.0;

    uint nTernFactors = 0;

    for (uint k=0; k < 2*mesh.nEdges(); k++)
      bin_terms += cost[mesh.nFaces() + k] * lpSolver.primalColumnSolution()[mesh.nFaces() + k];

    for (uint p=0; p < edge_pairs.size(); p++) {

      uint var_offs = pair_var_base[p];
      uint next_var_offs = (p+1 < edge_pairs.size()) ? pair_var_base[p+1] : nVars;
      
      uint diff = next_var_offs - var_offs;
      
      double sum = 0.0;
      for (uint k=var_offs; k < next_var_offs; k++) {
	sum += cost[k] * lpSolver.primalColumnSolution()[k];
      }

      if (diff == 2 || diff == 4)
	bin_terms += sum;
      else if (diff == 8) {
	tern_terms += sum;
	nTernFactors++;
      }
      else
	fo_terms += sum;
    }

    std::cerr << "binary terms: " << bin_terms << std::endl;
    std::cerr << "ternary terms: " << tern_terms << std::endl;
    std::cerr << "4th order terms: " << fo_terms << std::endl;
    
    std::cerr << nTernFactors << " ternary factors" << std::endl;

    energy = thresh_energy;
  }
  else
    std::cerr << "global optimum found!!" << std::endl;

  uint out_factor = options.output_factor_;

  segmentation.resize(xDim*out_factor,yDim*out_factor);

  mesh.enlarge(out_factor,out_factor);

  Math1D::Vector<int> labels(mesh.nFaces());
  
  Math2D::Matrix<double> output(xDim*out_factor,yDim*out_factor,0);
  for (uint i=0;i<mesh.nFaces();++i) {
    add_grid_output(i,lpSolver.primalColumnSolution()[i],mesh,output);	
    if (lpSolver.primalColumnSolution()[i] > 0.99) {
      labels[i] = 1;
    }
    else if (lpSolver.primalColumnSolution()[i] < 0.01) {
      labels[i] = 0;
    }
    else {
      labels[i] = -1;
    }
  }
  for (uint y=0; y < yDim*out_factor; y++) {
    for (uint x=0; x < xDim*out_factor; x++) {
      segmentation(x,y) = uint(output(x,y)*255.0);
    }
  }

  if (mesh.nFaces() > 20000) {
    //statusFailed();
  }
  else {
    mesh.draw_labels_with_pairs("mesh_lp.svg",lpSolver.primalColumnSolution(),edge_pairs,xDim,yDim);
  }

  return energy + energy_offset; 
}


double clique_lp_segment_curvreg_minsum_diffusion(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset, 
						  Math2D::Matrix<uint>& segmentation, const Math2D::Matrix<int>* fixed_labels) {

  //Get the options
  double lambda = options.lambda_;
  double gamma = options.gamma_;
  int neighborhood = options.neighborhood_;
  bool bruckstein = options.bruckstein_;
  
  assert(neighborhood <= 16); 

  uint xDim = uint( data_term.xDim() );
  uint yDim = uint( data_term.yDim() );

  Mesh2D mesh;  
  if (options.gridtype_ == options.Square) {
    if (options.adaptive_mesh_n_ < 0) {
      double xfac = double(xDim) / double(options.griddim_xDim_); 
      double yfac = double(yDim) / double(options.griddim_yDim_); 
      generate_mesh( options.griddim_xDim_, options.griddim_yDim_, neighborhood, mesh);
      mesh.enlarge(xfac,yfac);
    }
    else {
      std::cerr << "adaptive square mesh" << std::endl;

      //Adaptive mesh
      generate_adaptive_mesh(data_term, mesh, neighborhood, options.adaptive_mesh_n_);
    }
  }
  else {
    if (options.adaptive_mesh_n_ < 0) {
      double xfac = double(xDim) / double(options.griddim_xDim_); 
      double yfac = double(yDim) / double(options.griddim_yDim_); 
      generate_hexagonal_mesh( xDim, yDim, 0.5*(xfac+yfac), neighborhood,mesh); //TODO: proper handling of non-square images
    }
    else {
      std::cerr << "adaptive hex. mesh" << std::endl;

      //Adaptive mesh
      generate_adaptive_hexagonal_mesh(data_term, mesh, neighborhood, options.adaptive_mesh_n_);
    }
  }
  
  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  std::cerr << edge_pairs.size() << " edge pairs." << std::endl;

  uint nVars = 2*mesh.nFaces(); // for the diffusion we have two variables per face
  uint nConstraints = 0;
  uint nEntries = 0;

  Math1D::Vector<uint> pair_var_base(edge_pairs.size()+1);
  Math1D::Vector<uint> pair_con_base(edge_pairs.size());

  uint nQuadCliques = 0;
  uint nTripleCliques = 0;
  uint nDoubleCliques = 0;
  for (uint p=0; p < edge_pairs.size(); p++) {

    pair_var_base[p] = nVars;
    pair_con_base[p] = nConstraints;

    uint first_edge = edge_pairs[p].first_edge_idx_;
    uint second_edge = edge_pairs[p].second_edge_idx_;    

    std::set<uint> adjacent_regions;

    const std::vector<uint>& adj1 = mesh.adjacent_faces(first_edge);
    for (std::vector<uint>::const_iterator it = adj1.begin(); it != adj1.end(); it++)
      adjacent_regions.insert(*it);

    const std::vector<uint>& adj2 = mesh.adjacent_faces(second_edge);
    for (std::vector<uint>::const_iterator it = adj2.begin(); it != adj2.end(); it++)
      adjacent_regions.insert(*it);

    if (adjacent_regions.size() == 4) {
      nVars += 16;
      nConstraints += 4*2;
      nEntries += 16*4 + 4*2;

      nQuadCliques++;
    }
    else if (adjacent_regions.size() == 3) {
      nVars += 8;
      nConstraints += 3*2;
      nEntries += 8*3 + 3*2;
      
      nTripleCliques++;
    }
    else if (adjacent_regions.size() == 2) {
      nVars += 4; 
      nConstraints += 2*2; 
      nEntries += 4*2 + 2*2;
      nDoubleCliques++;
    }
    
    //NOTE: even a count of 1 is possible, at the image corners. 
    // The cost can then be directly added to the cost of the face, so there is no need for extra variables
  }

  pair_var_base[edge_pairs.size()] = nVars;

  Math1D::Vector<double> cost(nVars,0.0);
  
  for (uint i=0; i<mesh.nFaces(); ++i) {
    cost[2*i] = calculate_data_term(i, mesh, data_term);
    
    if (options.fix_regions_) {
      TODO("implement fixing of regions");
    }
  }

  //SparseMatrixDescription<double> lp_descr(nEntries, nConstraints, nVars); //for the current comp. of lower bounds
  SparseMatrixDescription<char> lp_descr(nEntries, nConstraints, nVars);

  /*** joint calculation of cost and constraint system ***/
  for (uint p=0; p < edge_pairs.size(); p++) {

    const uint first_edge = edge_pairs[p].first_edge_idx_;
    const uint second_edge = edge_pairs[p].second_edge_idx_;    

    const std::vector<uint>& adj1 = mesh.adjacent_faces(first_edge);
    const std::vector<uint>& adj2 = mesh.adjacent_faces(second_edge);

    double mul1 = (mesh.edge(first_edge).to_idx_ == edge_pairs[p].common_point_idx_) ? 1.0 : -1.0;
    double mul2 = (mesh.edge(second_edge).from_idx_ == edge_pairs[p].common_point_idx_) ? 1.0 : -1.0;

    double curvature_weight = gamma * curv_weight(mesh,edge_pairs[p],2.0,bruckstein);
    if (adj1.size() == 1 && adj2.size() == 1)
      curvature_weight = 0.0; //don't penalize the image corners

    double boundary_cost = curvature_weight;

    //no length for the image border
    if (adj1.size() > 1)
      boundary_cost += 0.5*lambda*mesh.edge_length(first_edge); 
    if (adj2.size() > 1)
      boundary_cost += 0.5*lambda*mesh.edge_length(second_edge);

    std::set<uint> adjacent_regions_set;
    std::map<uint,double> match1;
    std::map<uint,double> match2;

    for (std::vector<uint>::const_iterator it = adj1.begin(); it != adj1.end(); it++) {
      adjacent_regions_set.insert(*it);
      match1[*it] = mesh.match(*it,first_edge);
    }

    for (std::vector<uint>::const_iterator it = adj2.begin(); it != adj2.end(); it++) {
      adjacent_regions_set.insert(*it);
      match2[*it] = mesh.match(*it,second_edge);
    }

    const uint cur_var_base = pair_var_base[p];
    const uint cur_con_base = pair_con_base[p];
  
    Math1D::Vector<uint> adjacent_regions(adjacent_regions_set.size());
    uint k=0;
    for (std::set<uint>::iterator it = adjacent_regions_set.begin(); it != adjacent_regions_set.end(); it++) {
      adjacent_regions[k] = *it;
      k++;
    }

    if (adjacent_regions.size() == 4) {

      for (uint f=0; f < 4; f++) {
	
	lp_descr.add_entry(cur_con_base + 2*f, 2*adjacent_regions[f], -1);
	lp_descr.add_entry(cur_con_base + 2*f+1, 2*adjacent_regions[f]+1, -1);
      }
      
      for (uint v_addon = 0; v_addon < 16; v_addon++) {

	double sum_e1 = 0.0;
	double sum_e2 = 0.0;

	for (uint f=0; f < 4; f++) {

	  uint mask = 1 << f;
	  assert(mask != 0);

	  //uint offs = (v_addon & mask) >> f;
	  uint offs = 1 - ((v_addon & mask) >> f);
	  lp_descr.add_entry(cur_con_base + 2*f+offs, cur_var_base + v_addon, 1);

	  sum_e1 += (1 - offs) * match1[adjacent_regions[f]];
	  sum_e2 += (1 - offs) * match2[adjacent_regions[f]];	  
	}

	sum_e1 *= mul1;
	sum_e2 *= mul2;

	if (fabs(sum_e1) > 0.5 && fabs(sum_e2) > 0.5 && sum_e1*sum_e2 > 0.0) {
	  cost[cur_var_base + v_addon] = boundary_cost;
	}
      }
    }
    else if (adjacent_regions.size() == 3) {
      
      for (uint f=0; f < 3; f++) {
	
	lp_descr.add_entry(cur_con_base + 2*f, 2*adjacent_regions[f], -1);
	lp_descr.add_entry(cur_con_base + 2*f+1, 2*adjacent_regions[f]+1, -1);
      }
      
      for (uint v_addon = 0; v_addon < 8; v_addon++) {

	double sum_e1 = 0.0;
	double sum_e2 = 0.0;

	for (uint f=0; f < 3; f++) {
	
	  uint mask = 1 << f;
	  assert(mask != 0);

	  //uint offs = (v_addon & mask) >> f;
	  uint offs = 1 - ((v_addon & mask) >> f);
	  lp_descr.add_entry(cur_con_base + 2*f+offs, cur_var_base + v_addon, 1);
	  
	  sum_e1 += (1 - offs) * match1[adjacent_regions[f]];
	  sum_e2 += (1 - offs) * match2[adjacent_regions[f]];
	}

	sum_e1 *= mul1;
	sum_e2 *= mul2;

	if (fabs(sum_e1) > 0.5 && fabs(sum_e2) > 0.5 && sum_e1*sum_e2 > 0.0) {

	  //at the image boundary there can be triple cliques that do not share a region
	  if (adj1.size() > 1 && adj2.size() > 1) {
	    
	    assert(sum_e1*sum_e2 > 0.0);
	  }
	  
	  cost[cur_var_base + v_addon] = boundary_cost;
	}
      }
    }
    else if (adjacent_regions.size() == 2) {
      
      //here we use the construction via the absolute differences

      for (uint f=0; f < 2; f++) {
	
	lp_descr.add_entry(cur_con_base + 2*f, 2*adjacent_regions[f], -1);
	lp_descr.add_entry(cur_con_base + 2*f+1, 2*adjacent_regions[f]+1, -1);
      }
      
      for (uint v_addon = 0; v_addon < 4; v_addon++) {

	double sum_e1 = 0.0;
	double sum_e2 = 0.0;

	for (uint f=0; f < 2; f++) {
	
	  uint mask = 1 << f;
	  assert(mask != 0);

	  //uint offs = (v_addon & mask) >> f;
	  uint offs = 1 - ((v_addon & mask) >> f);
	  lp_descr.add_entry(cur_con_base + 2*f+offs, cur_var_base + v_addon, 1);
	  
	  sum_e1 += (1-offs) * match1[adjacent_regions[f]];
	  sum_e2 += (1-offs) * match2[adjacent_regions[f]];
	}

	sum_e1 *= mul1;
	sum_e2 *= mul2;

	if (fabs(sum_e1) > 0.5 && fabs(sum_e2) > 0.5 && sum_e1*sum_e2 > 0.0) {

	  //at the image boundary there can be double cliques that do not share a region
	  //if (adj1.size() > 1 && adj2.size() > 1) {
	  //  assert(sum_e1*sum_e2 > 0.0);
	  //}
	  
	  cost[cur_var_base + v_addon] = boundary_cost;
	}
      }
    }
  }

  Math1D::Vector<uint> row_start;
  lp_descr.sort_by_row(row_start,true);

  //we only need the column indices, so to save memory we copy them, then release the sparse matrix
  Math1D::Vector<uint> col_idx(lp_descr.nEntries());
  for (uint i=0; i < col_idx.size(); i++)
    col_idx[i] = lp_descr.col_indices()[i];
  
  lp_descr.reset(0);

  double lower_bound = 0.0;

  std::cerr.precision(12);
  
  /**** problem constructed, now start the minsum-diffusion ****/
  for (uint iter = 1; iter <= 12000; iter++) {
    
    std::cerr << "******** iteration " << iter << std::endl;
    
    //a) compute current lower bound
    lower_bound = 0.0;
    for (uint f=0; f < mesh.nFaces(); f++) {
      
      lower_bound += std::min(cost[2*f],cost[2*f+1]);
    }
    std::cerr << "unaries: " << lower_bound << std::endl;
    for (uint p=0; p < edge_pairs.size(); p++) {
	
      const uint cur_var_base = pair_var_base[p];
      const uint cur_var_end = pair_var_base[p+1];
	
      double min_val = MAX_DOUBLE;
      for (uint v=cur_var_base; v < cur_var_end; v++) {
	double val = cost[v];
	if (val < min_val)
	  min_val = val;
      }
	
      lower_bound += min_val;
    }    
      
    std::cerr << "lower bound: " << lower_bound << "   (" << (lower_bound+energy_offset) << ")" << std::endl;
      
    //b) update the cost parameters
    for (uint row = 0; row < row_start.size()-1; row++) {
      
      //std::cerr << "row: " << row << std::endl;
      
      double a = MAX_DOUBLE;
      double b = MAX_DOUBLE;
      
      for (uint k = row_start[row]; k < row_start[row+1];  k++) {
	
	uint var_idx = col_idx[k]; 
	
	if (var_idx < 2*mesh.nFaces()) {
	  
	  assert(b == MAX_DOUBLE);
	  b = cost[var_idx];
	}
	else if (cost[var_idx] < a)
	  a = cost[var_idx];
      }
      
      double shift = 0.5*(b-a);
	
      for (uint k = row_start[row]; k < row_start[row+1];  k++) {
	  
	uint var_idx = col_idx[k]; 
	  
	if (var_idx < 2*mesh.nFaces()) {
	  cost[var_idx] -= shift;
	}
	else 
	  cost[var_idx] += shift;
      }
    }    
  }

  /**** analyzation stage ****/
  uint nNegligible = 0;
  for (uint p=0; p < edge_pairs.size(); p++) {
    
    const uint first_edge = edge_pairs[p].first_edge_idx_;
    const uint second_edge = edge_pairs[p].second_edge_idx_;
    
    std::set<uint> adjacent_regions_set;
      
    const std::vector<uint>& adj1 = mesh.adjacent_faces(first_edge);
    for (std::vector<uint>::const_iterator it = adj1.begin(); it != adj1.end(); it++) {
      adjacent_regions_set.insert(*it);
    }
      
    const std::vector<uint>& adj2 = mesh.adjacent_faces(second_edge);
    for (std::vector<uint>::const_iterator it = adj2.begin(); it != adj2.end(); it++) {
      adjacent_regions_set.insert(*it);
    }
      
    Math1D::Vector<uint> adjacent_regions(adjacent_regions_set.size());
    uint k=0;
    for (std::set<uint>::iterator it = adjacent_regions_set.begin(); it != adjacent_regions_set.end(); it++) {
      adjacent_regions[k] = *it;
      k++;
    }
      
    if (adjacent_regions.size() > 1) {
	
      uint nCurParams = 1 << adjacent_regions.size();
	
      const uint cur_var_base = pair_var_base[p];
      
      double sum = 0.0;
      for (uint v = cur_var_base; v < cur_var_base + nCurParams; v++)
	sum += fabs(cost[v]);
      
      if (sum < 1e-3) 
	nNegligible++;
    }
  }
    
  std::cerr << nNegligible << " cliques are negligible after reparameterization" << std::endl;

  /*** extract solution ***/
  Math1D::Vector<double> face_solution(mesh.nFaces(),0.0);
  Math1D::Vector<uint> face_uint_solution(mesh.nFaces(),0);
  
  for (uint f=0; f < mesh.nFaces(); f++) {
    
    if (cost[2*f] < cost[2*f+1]) {
      face_solution[f] = 1.0;
      face_uint_solution[f] = 1;
    }
  }

  double disc_energy = clique_curv_energy(mesh, edge_pairs, data_term, options, face_uint_solution)
    + energy_offset;
  
  std::cerr << "solution found by stand-alone routine: " << disc_energy << std::endl;

  //compute cost of solution
  //currently based on calling an lp-solver
  // {

  //   Math1D::Vector<double> var_lb(nVars,0.0);
  //   Math1D::Vector<double> var_ub(nVars,1.0);

  //   for (uint f=0; f < mesh.nFaces(); f++) {

  //     if (face_solution[f] == 1) {
  // 	var_lb[2*f] = 1.0;
  // 	var_ub[2*f+1] = 0.0;
  //     }
  //     else {
  // 	var_ub[2*f] = 0.0;
  // 	var_lb[2*f+1] = 1.0;
  //     }
  //   }


  //   CoinPackedMatrix coinMatrix(false,(int*) lp_descr.row_indices(),(int*) lp_descr.col_indices(),
  // 				lp_descr.value(),lp_descr.nEntries());

  //   Math1D::Vector<double> rhs(nConstraints,0.0); //only needed to compute the cost of the thresholded solution
    
  //   ClpSimplex lpSolver;
  //   lpSolver.loadProblem(coinMatrix, var_lb.direct_access(), var_ub.direct_access(),   
  // 			 cost.direct_access(), rhs.direct_access(), rhs.direct_access());

  //   coinMatrix.cleanMatrix();
    
  //   for (uint v=0; v < nVars; v++) {
  //     lpSolver.setInteger(v);
  //   }
  
  //   lpSolver.dual();

  //   uint nNonInt = 0;
  //   for (uint v=0; v < nVars; v++) {
  //     double val = lpSolver.primalColumnSolution()[v];

  //     if (val > 0.01 && val < 0.99)
  // 	nNonInt++;
  //   }


  //   double int_energy = lpSolver.getObjValue() + energy_offset;

  //   std::cerr.precision(10);
  //   std::cerr << "energy of thresholded solution: " << int_energy 
  // 	      << "(" << nNonInt << " non-integral variables)" << std::endl;

  // }
  
  uint out_factor = options.output_factor_;
  mesh.enlarge(out_factor,out_factor);
    
  segmentation.resize(xDim*out_factor,yDim*out_factor);

  Math1D::Vector<int> labels(mesh.nFaces());
  
  Math2D::Matrix<double> output(xDim*out_factor,yDim*out_factor,0);
  for (uint i=0;i<mesh.nFaces();++i) {
    add_grid_output(i,face_solution[i],mesh,output);	
  }
  for (uint y=0; y < yDim*out_factor; y++) {
    for (uint x=0; x < xDim*out_factor; x++) {
      segmentation(x,y) = uint(output(x,y)*255.0);
    }
  }
  
  if (mesh.nFaces() > 20000) {
    //statusFailed();
  }
  else {
    mesh.draw_labels_with_pairs("mesh_lp.svg",face_solution.direct_access(),edge_pairs,xDim,yDim);
  }
  
  return lower_bound;
}

