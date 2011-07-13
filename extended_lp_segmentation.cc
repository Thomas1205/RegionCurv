/*** written by Thomas Schoenemann as an employee of Lund University, Sweden, 2010 - 2011 ***/

#include "extended_lp_segmentation.hh"
#include "mesh2D.hh"
#include "sparse_matrix_description.hh"
#include "timing.hh"
#include "curvature.hh"


#include "segmentation_common.hh"
#include "stl_out.hh"

#include <coin/ClpSimplex.hpp>
#include <coin/ClpPlusMinusOneMatrix.hpp>
#include <coin/ClpFactorization.hpp>

#include "conv_lp_solving.hh"

#define USE_PM_ONE_MATRIX

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

    if (curv_weight(mesh,edge_pairs[p],2.0,bruckstein) < curv_threshold)
      continue;

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
  for (uint v=mesh.nFaces(); v < mesh.nFaces()+2*mesh.nEdges(); v++) {
    if (mesh.adjacent_faces((v-mesh.nFaces())/2).size() == 2)
      cost[v] = lambda;
  }
  
  SparseMatrixDescription<double> lp_descr(nEntries, nConstraints, nVars);

  for (uint e = 0; e < mesh.nEdges(); e++) {

    for (uint f=0; f < mesh.adjacent_faces(e).size(); f++) {

      uint var = mesh.adjacent_faces(e)[f];
      lp_descr.add_entry(e, var, mesh.match(var,e));
    }

    lp_descr.add_entry(e,mesh.nFaces()+2*e,1.0);
    lp_descr.add_entry(e,mesh.nFaces()+2*e+1,-1.0);
  }
  
  /*** joint calculation of cost and constraint system ***/
  for (uint p=0; p < edge_pairs.size(); p++) {

    const uint first_edge = edge_pairs[p].first_edge_idx_;
    const uint second_edge = edge_pairs[p].second_edge_idx_;    

    double mul1 = (mesh.edge(first_edge).to_idx_ == edge_pairs[p].common_point_idx_) ? 1.0 : -1.0;
    double mul2 = (mesh.edge(second_edge).from_idx_ == edge_pairs[p].common_point_idx_) ? 1.0 : -1.0;

    double curvature_weight = gamma * curv_weight(mesh,edge_pairs[p],2.0,bruckstein);

    if ((curvature_weight / gamma) < curv_threshold)
      continue; //length-based cost are handled via two-cliques

    double boundary_cost = curvature_weight;

    std::set<uint> adjacent_regions_set;
    std::map<uint,double> match1;
    std::map<uint,double> match2;

    const std::vector<uint>& adj1 = mesh.adjacent_faces(first_edge);
    for (std::vector<uint>::const_iterator it = adj1.begin(); it != adj1.end(); it++) {
      adjacent_regions_set.insert(*it);
      match1[*it] = mesh.match(*it,first_edge);
    }

    const std::vector<uint>& adj2 = mesh.adjacent_faces(second_edge);
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

	  uint offs = (v_addon & mask) >> f;
	  lp_descr.add_entry(cur_con_base + 2*f+offs, cur_var_base + v_addon, 1.0);

	  sum_e1 += offs * match1[adjacent_regions[f]];
	  sum_e2 += offs * match2[adjacent_regions[f]];	  
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

	  uint offs = (v_addon & mask) >> f;
	  lp_descr.add_entry(cur_con_base + 2*f+offs, cur_var_base + v_addon, 1.0);
	  
	  sum_e1 += offs * match1[adjacent_regions[f]];
	  sum_e2 += offs * match2[adjacent_regions[f]];
	}

	sum_e1 *= mul1;
	sum_e2 *= mul2;

	if (fabs(sum_e1) > 0.5 && fabs(sum_e2) > 0.5) {

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

      //here we use the construction via the absolute differences
      lp_descr.add_entry(cur_con_base, adjacent_regions[0], 1.0);
      lp_descr.add_entry(cur_con_base, adjacent_regions[1], -1.0);
      
      lp_descr.add_entry(cur_con_base, cur_var_base, 1.0);
      lp_descr.add_entry(cur_con_base, cur_var_base+1, -1.0);

      cost[cur_var_base] = boundary_cost;
      cost[cur_var_base+1] = boundary_cost;
    }
    else {
      assert(adjacent_regions.size() == 1);
      cost[adjacent_regions[0]] += boundary_cost;
    }
  }

  Math1D::Vector<uint> row_start(nConstraints+1);
  lp_descr.sort_by_row(row_start,true);

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

  uint nNonInt = 0;
  for (uint v=0; v < mesh.nFaces(); v++) {
    
    double val = lp_solution[v];

    if (val > 0.01 && val < 0.99)
      nNonInt++;
    
    if (val < 0.5) 
      var_ub[v] = 0.0;
    else
      var_lb[v] = 1.0;
  }
  
  double energy = lp_energy;
  
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


double clique_lp_segment_curvreg_maxsum_diffusion(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset, 
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

  SparseMatrixDescription<char> lp_descr(nEntries, nConstraints, nVars);
  
  /*** joint calculation of cost and constraint system ***/
  for (uint p=0; p < edge_pairs.size(); p++) {

    const uint first_edge = edge_pairs[p].first_edge_idx_;
    const uint second_edge = edge_pairs[p].second_edge_idx_;    

    double mul1 = (mesh.edge(first_edge).to_idx_ == edge_pairs[p].common_point_idx_) ? 1.0 : -1.0;
    double mul2 = (mesh.edge(second_edge).from_idx_ == edge_pairs[p].common_point_idx_) ? 1.0 : -1.0;

    double curvature_weight = gamma * curv_weight(mesh,edge_pairs[p],2.0,bruckstein);

    double boundary_cost = 0.5*lambda*(mesh.edge_length(first_edge) + mesh.edge_length(second_edge)) + curvature_weight;

    std::set<uint> adjacent_regions_set;
    std::map<uint,double> match1;
    std::map<uint,double> match2;

    const std::vector<uint>& adj1 = mesh.adjacent_faces(first_edge);
    for (std::vector<uint>::const_iterator it = adj1.begin(); it != adj1.end(); it++) {
      adjacent_regions_set.insert(*it);
      match1[*it] = mesh.match(*it,first_edge);
    }

    const std::vector<uint>& adj2 = mesh.adjacent_faces(second_edge);
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

	  uint offs = (v_addon & mask) >> f;
	  lp_descr.add_entry(cur_con_base + 2*f+offs, cur_var_base + v_addon, 1);

	  sum_e1 += offs * match1[adjacent_regions[f]];
	  sum_e2 += offs * match2[adjacent_regions[f]];	  
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

	  uint offs = (v_addon & mask) >> f;
	  lp_descr.add_entry(cur_con_base + 2*f+offs, cur_var_base + v_addon, 1);
	  
	  sum_e1 += offs * match1[adjacent_regions[f]];
	  sum_e2 += offs * match2[adjacent_regions[f]];
	}

	sum_e1 *= mul1;
	sum_e2 *= mul2;

	if (fabs(sum_e1) > 0.5 && fabs(sum_e2) > 0.5) {

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

	  uint offs = (v_addon & mask) >> f;
	  lp_descr.add_entry(cur_con_base + 2*f+offs, cur_var_base + v_addon, 1);
	  
	  sum_e1 += offs * match1[adjacent_regions[f]];
	  sum_e2 += offs * match2[adjacent_regions[f]];
	}

	sum_e1 *= mul1;
	sum_e2 *= mul2;

	if (fabs(sum_e1) > 0.5 && fabs(sum_e2) > 0.5) {

	  //at the image boundary there can be triple cliques that do not share a region
	  if (adj1.size() > 1 && adj2.size() > 1) {
	    
	    assert(sum_e1*sum_e2 > 0.0);
	  }
	  
	  cost[cur_var_base + v_addon] = boundary_cost;
	}
      }
    }
    else {
      assert(adjacent_regions.size() == 1);
      cost[2*adjacent_regions[0]] += boundary_cost;
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

  for (uint f=0; f < mesh.nFaces(); f++) {

    if (cost[2*f] < cost[2*f+1])
      face_solution[f] = 1.0;
  }

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


