/*** first versions written by Thomas Schoenemann as an employee of Lund University, Sweden, 2010 - 2011 ***/
/*** refined by Thomas Schoenemann as an employee of the University of Pisa, Italy, 2011 ****/

#include "extended_lp_segmentation.hh"
#include "mesh2D.hh"
#include "sparse_matrix_description.hh"
#include "line_drawing.hh"
#include "timing.hh"
#include "curvature.hh"


#include "segmentation_common.hh"
#include "stl_out.hh"

#include <coin/ClpSimplex.hpp>
#include <coin/ClpPlusMinusOneMatrix.hpp>
#include <coin/ClpFactorization.hpp>
#include <coin/OsiClpSolverInterface.hpp>
#include <coin/CbcModel.hpp>

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

double multi_point_energy(const Mesh2D& mesh, uint point, const LPSegOptions& options, 
			  const Storage1D<uint>& integral_solution, uint nRegions,
			  const std::vector<Mesh2DEdgePair>& edge_pairs, 
			  const NamedStorage1D<std::vector<uint> >& point_pair,
			  const NamedStorage1D<std::vector<uint> >& point_edge) {


  double energy = 0.0;

  //std::cerr << "**** point_energy(" << point << ") ****" << std::endl;
  //std::cerr << "nRegions: " << nRegions << std::endl;

  double lambda = 0.5*options.lambda_;
  double gamma = 0.5*options.gamma_;
  bool bruckstein = options.bruckstein_;
  bool crossings_allowed = !options.prevent_crossings_;
  
  Storage1D<Math1D::Vector<double> > drop(nRegions);
  for (uint i=0; i < nRegions; i++)
    drop[i].resize(point_edge[point].size(),0.0);
  
  Math1D::Vector<double> abs_drop_sum(nRegions,0.0);
  
  std::map<uint,uint> drop_idx;

  uint nEdges = point_edge[point].size();

  for (uint e=0; e < nEdges; e++) {

    uint edge = point_edge[point][e];
    drop_idx[edge] = e;

    //std::cerr << "edge " << edge << std::endl;

    const std::vector<uint>& adjacent_faces = mesh.adjacent_faces(edge);

    Math1D::Vector<double> cur_drop(nRegions,0.0);

    for (uint i=0; i < adjacent_faces.size(); i++) {
      uint face = adjacent_faces[i];

      uint cur_label = integral_solution[face];

      //std::cerr << "face " << face << " ---> label " << cur_label << std::endl;

      cur_drop[cur_label] += mesh.match(face,edge);
    }

    for (uint r=0; r < nRegions; r++) {
      drop[r][e] = cur_drop[r];
      
      abs_drop_sum[r] += fabs(cur_drop[r]);
    }
  }

  //std::cerr << "abs_drop_sum: " << abs_drop_sum << std::endl;

  //can now solve for each region independently
  for (uint r=0; r < nRegions; r++) {

    if (abs_drop_sum[r] > 0.0) {

      //std::cerr << "pair present" << std::endl;

      int ads = int (abs_drop_sum[r]+0.5);
      assert((ads % 2) == 0);

      if (ads == 2) {
	//simple case
	
	for (uint j=0; j < point_pair[point].size(); j++) {
	  
	  const Mesh2DEdgePair& cur_edge_pair = edge_pairs[point_pair[point][j]];

	  uint first = cur_edge_pair.first_edge_idx_;
	  uint second = cur_edge_pair.second_edge_idx_;
	  
	  if (drop[r][drop_idx[first]] != 0.0 && drop[r][drop_idx[second]] != 0.0) {

	    uint nFirstAdjacent = uint( mesh.adjacent_faces(first).size() );
	    uint nSecondAdjacent = uint( mesh.adjacent_faces(second).size() );
	    
	    double weight = 0.0;
	    if (nFirstAdjacent > 1)
	      weight += 0.5*lambda*mesh.edge_length(first);
	    if (nSecondAdjacent > 1)
	      weight += 0.5*lambda*mesh.edge_length(second);
	    //do not penalize the image corners for their curvature
	    if (nFirstAdjacent > 1 || nSecondAdjacent > 1)
	      weight += gamma * curv_weight(mesh,cur_edge_pair,2.0,bruckstein);
	    
	    //since the cost are identical for both orientations, we don't care about orientation here
	    energy += weight;
	  }
	}
      }
      else {
	//use Integer Linear Programming
	
	uint nVars = 2* point_pair[point].size();
	
	//std::cerr << "ILP, r=" << r << std::endl;
	//std::cerr << "nVars: " << nVars << std::endl;
	//std::cerr << "drop_sum: " << abs_drop_sum[r] << std::endl;
	//std::cerr << "drop: " << drop[r] << std::endl;
	
	Math1D::Vector<double> var_lb(nVars,0.0);
	Math1D::Vector<double> var_ub(nVars,1.0);
	Math1D::Vector<double> lp_cost(nVars,0.0);
	
	//set cost and remove pairs that cannot be active
	for (uint j=0; j < point_pair[point].size(); j++) {
	  
	  const Mesh2DEdgePair& cur_edge_pair = edge_pairs[point_pair[point][j]];

	  uint first = cur_edge_pair.first_edge_idx_;
	  uint second = cur_edge_pair.second_edge_idx_;
	  
	  uint nFirstAdjacent = uint( mesh.adjacent_faces(first).size() );
	  uint nSecondAdjacent = uint( mesh.adjacent_faces(second).size() );
	  
	  double weight = 0.0;
	  if (nFirstAdjacent > 1)
	    weight += 0.5*lambda*mesh.edge_length(first);
	  if (nSecondAdjacent > 1)
	    weight += 0.5*lambda*mesh.edge_length(second);
	  //do not penalize the image corners for their curvature
	  if (nFirstAdjacent > 1 || nSecondAdjacent > 1)
	    weight += gamma * curv_weight(mesh,cur_edge_pair,2.0,bruckstein);
	  
	  lp_cost[2*j]   = weight;
	  lp_cost[2*j+1] = weight;
	  
	  if (drop[r][drop_idx[first]] == 0.0 || drop[r][drop_idx[second]] == 0.0) {
	    var_ub[2*j] = 0.0;
	    var_ub[2*j+1] = 0.0;
	  }
	}

	std::vector<std::pair<uint,uint> > conflicts;

	//find conflicting pairs
	for (uint p1=0; p1 < point_pair[point].size(); p1++) {
	  for (uint p2=p1+1; p2 < point_pair[point].size(); p2++) {
	    
	    uint pair1 = point_pair[point][p1];
	    uint pair2 = point_pair[point][p2];
	    
	    if (line_pairs_with_meeting_point_cross(mesh, edge_pairs[pair1], edge_pairs[pair2]) ) {
	      conflicts.push_back(std::make_pair(p1,p2));
	    }
	  }
	}
	
	if (crossings_allowed && ads != 4)
	  conflicts.clear();

	//std::cerr << conflicts.size() << " conflicts" << std::endl;
	
	uint nConstraints = 2*nEdges + conflicts.size();
	uint nMatrixEntries = 4*point_pair.size() + 4*conflicts.size();
	
	SparseMatrixDescription<double> lp_descr(nMatrixEntries, nConstraints, nVars);

	Math1D::Vector<double> rhs_lower(nConstraints,0.0);
	Math1D::Vector<double> rhs_upper(nConstraints,0.0);
	
	for (uint e=0; e < nEdges; e++) {
	  rhs_lower[e] = -drop[r][e];
	  rhs_upper[e] = -drop[r][e];
	  rhs_upper[nEdges+e] = 1.0;
	}
	
	/*** code constraint system ***/
	for (uint j=0; j < point_pair[point].size(); j++) {
	  
	  const Mesh2DEdgePair& cur_edge_pair = edge_pairs[point_pair[point][j]];
	  
	  uint first_edge = cur_edge_pair.first_edge_idx_;
	  uint second_edge = cur_edge_pair.second_edge_idx_;
	  
	  uint middle_point = cur_edge_pair.common_point_idx_;
	  
	  lp_descr.add_entry(nEdges+drop_idx[first_edge],2*j,1);
	  lp_descr.add_entry(nEdges+drop_idx[first_edge],2*j+1,1);
	  
	  lp_descr.add_entry(nEdges+drop_idx[second_edge],2*j,1);
	  lp_descr.add_entry(nEdges+drop_idx[second_edge],2*j+1,1);
	  
	  if (mesh.edge(first_edge).to_idx_ == middle_point) {
	    lp_descr.add_entry(drop_idx[first_edge],2*j,1);
	    lp_descr.add_entry(drop_idx[first_edge],2*j+1,-1);
	  }
	  else {
	    lp_descr.add_entry(drop_idx[first_edge],2*j,-1);
	    lp_descr.add_entry(drop_idx[first_edge],2*j+1,1);
	  }
	  
	  if (mesh.edge(second_edge).to_idx_ == middle_point) {
	    lp_descr.add_entry(drop_idx[second_edge],2*j,-1);
	    lp_descr.add_entry(drop_idx[second_edge],2*j+1,1);
	  }
	  else {
	    lp_descr.add_entry(drop_idx[second_edge],2*j,1);
	    lp_descr.add_entry(drop_idx[second_edge],2*j+1,-1);
	  }
	}
	
	for (uint k=0; k < conflicts.size(); k++) {

	  uint row = 2*point_edge[point].size() + k;
	  rhs_upper[row] = 1.0;
	  
	  uint first = conflicts[k].first;
	  uint second = conflicts[k].second;
	  
	  lp_descr.add_entry(row, 2*first, 1.0);
	  lp_descr.add_entry(row, 2*first+1, 1.0);
	  lp_descr.add_entry(row, 2*second, 1.0);
	  lp_descr.add_entry(row, 2*second+1, 1.0);
	}

	OsiClpSolverInterface lpSolver;
	lpSolver.messageHandler()->setLogLevel(5);
	
	CoinPackedMatrix coinMatrix(false,(int*) lp_descr.row_indices(),(int*) lp_descr.col_indices(),
				    lp_descr.value(),lp_descr.nEntries());
	
	lpSolver.loadProblem (coinMatrix, var_lb.direct_access(), var_ub.direct_access(),   
			      lp_cost.direct_access(), rhs_lower.direct_access(), rhs_upper.direct_access());
	
	for (int v=0; v < lpSolver.getNumCols(); v++) {
	  lpSolver.setInteger(v);
	}
	
	//ILP
	CbcModel cbc_model(lpSolver);
	cbc_model.setLogLevel(0);
	//cbc_model.solver()->messageHandler()->setLogLevel(5);

	cbc_model.branchAndBound();

	energy += cbc_model.getObjValue();
	
	if (!cbc_model.isProvenOptimal()) {
	  
	  std::cerr << "ERROR: the optimal solution could not be found. Exiting..." << std::endl;
	  exit(1);
	}
	if (cbc_model.isProvenInfeasible()) {
	  
	  std::cerr << "ERROR: problem marked as infeasible. Exiting..." << std::endl;
	  exit(1);
	}
	if (cbc_model.isAbandoned()) {
	  
	  std::cerr << "ERROR: problem was abandoned. Exiting..." << std::endl;
	  exit(1);
	}
      }
    }
  }

  //std::cerr << "leaving" << std::endl;

  return energy;
}


double multi_curv_energy(const Mesh2D& mesh, const Storage1D<uint>& integral_solution, 
			 const Math1D::Vector<double>& cost, const LPSegOptions& options, uint nRegions) {

  double energy = 0.0;

  for (uint f=0; f < mesh.nFaces(); f++) {
    uint label = integral_solution[f];
    energy += cost[f*nRegions+label];
  }
    
  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  //std::cerr << "unary terms: " << energy << std::endl;
  
  NamedStorage1D<std::vector<uint> > point_pair(mesh.nPoints(),MAKENAME(point_pair));
  NamedStorage1D<std::vector<uint> > point_edge(mesh.nEdges(),MAKENAME(point_edge));

  for (uint e=0; e < edge_pairs.size(); e++) {

    uint point = edge_pairs[e].common_point_idx_;

    point_pair[point].push_back(e);
    uint first = edge_pairs[e].first_edge_idx_;
    uint second = edge_pairs[e].second_edge_idx_;

    if (std::find(point_edge[point].begin(), point_edge[point].end(), first) == point_edge[point].end())
      point_edge[point].push_back(first);
    if (std::find(point_edge[point].begin(), point_edge[point].end(), second) == point_edge[point].end())
      point_edge[point].push_back(second);
  }

  for (uint p=0; p < mesh.nPoints(); p++) {

    energy += multi_point_energy(mesh, p, options, integral_solution, nRegions,
				 edge_pairs, point_pair, point_edge);
  }

  return energy;
}

double multi_curv_icm(const Mesh2D& mesh, const Storage1D<uint>& integral_solution, 
		      const Math1D::Vector<double>& cost, const LPSegOptions& options, uint nRegions) {

  double energy = 0.0;

  for (uint f=0; f < mesh.nFaces(); f++) {
    uint label = integral_solution[f];
    energy += cost[f*nRegions+label];
  }
    
  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);
  
  //std::cerr << "unary terms: " << energy << std::endl;
  //std::cerr << "nRegions: " << nRegions << std::endl;
  
  NamedStorage1D<std::vector<uint> > point_pair(mesh.nPoints(),MAKENAME(point_pair));
  NamedStorage1D<std::vector<uint> > point_edge(mesh.nEdges(),MAKENAME(point_edge));
  
  for (uint e=0; e < edge_pairs.size(); e++) {
    
    uint point = edge_pairs[e].common_point_idx_;
    
    point_pair[point].push_back(e);
    uint first = edge_pairs[e].first_edge_idx_;
    uint second = edge_pairs[e].second_edge_idx_;
    
    if (std::find(point_edge[point].begin(), point_edge[point].end(), first) == point_edge[point].end())
      point_edge[point].push_back(first);
    if (std::find(point_edge[point].begin(), point_edge[point].end(), second) == point_edge[point].end())
      point_edge[point].push_back(second);
  }
  
  Math1D::Vector<double> cur_point_energy(mesh.nPoints());
  
  for (uint p=0; p < mesh.nPoints(); p++) {
  
    cur_point_energy[p] = multi_point_energy(mesh, p, options, integral_solution, nRegions,
					     edge_pairs, point_pair, point_edge);
    
    energy += cur_point_energy[p];
  }
  
  
  std::cerr << "ICM. Initial energy: " << energy << std::endl;
  
  bool changes = true;
  
  for (uint iter=1; changes && iter <= 15; iter++) {
    
    changes = false;
    uint nChanges = 0;
    
    std::cerr << "ICM iteration " << iter << std::endl;
    
    for (uint f=0; f < mesh.nFaces(); f++) {
      
      for (uint hyp_label = 0; hyp_label < nRegions; hyp_label++) {
	
	uint cur_label = integral_solution[f];
	if (cur_label == hyp_label)
	  continue;
	
	double cur_energy = 0.0;
	cur_energy += cost[f*nRegions+cur_label];
	
	std::vector<uint> point_indices;
	mesh.get_polygon_points(f, point_indices);
	
	for (uint k=0; k < point_indices.size(); k++) 
	  cur_energy += cur_point_energy[point_indices[k]];
	
	//temporarily modify the solution
	integral_solution[f] = hyp_label;
	
	Math1D::Vector<double> hyp_point_cost(point_indices.size());
	
	double hyp_energy = 0.0;
	hyp_energy += cost[f*nRegions+hyp_label];
	
	for (uint p=0; p < point_indices.size(); p++) {
	  hyp_point_cost[p] = multi_point_energy(mesh, p, options, integral_solution, nRegions,
						 edge_pairs, point_pair, point_edge);
	  hyp_energy += hyp_point_cost[p];
	}
	
	if (cur_energy <= hyp_energy) {
	  integral_solution[f] = cur_label;
	}
	else {
	  //update point cost
	  for (uint k=0; k < point_indices.size(); k++) {
	    uint idx = point_indices[k];
	    cur_point_energy[idx] = hyp_point_cost[k];
	  }
	  changes = true;
	  nChanges++;
	  
	  energy -= (cur_energy - hyp_energy);
	}
      }
    }

    std::cerr << "energy " << energy << " (" << nChanges << " changes)"<< std::endl;
  }

  return energy;
}

double multi_curv_icm(const Math3D::Tensor<float>& data_term, const LPSegOptions& options,
		      Math2D::Matrix<uint>& segmentation) {

  uint nRegions = data_term.zDim();

  int neighborhood = options.neighborhood_;
  
  assert(neighborhood <= 16); 

  uint xDim = uint( data_term.xDim() );
  uint yDim = uint( data_term.yDim() );

  Mesh2D mesh;  
  create_mesh(options, data_term, 0, mesh);

  Math1D::Vector<double> region_cost(nRegions * mesh.nFaces());
  Math1D::Vector<uint> face_label(mesh.nFaces(),0);

  for (uint f=0; f < mesh.nFaces(); f++) {

    double best = 1e300;
    uint arg_best = 0;

    for (uint l=0; l < nRegions; l++) {

      double hyp = calculate_data_term(f, l, mesh, data_term);

      region_cost[f*nRegions+l] = hyp;

      if (hyp < best) {
	
	best = hyp;
	arg_best = l;
      }
    }
    
    face_label[f] = arg_best;
  }

  //std::cerr << "calling main for " << data_term.zDim() << " regions" << std::endl;

  double energy = multi_curv_icm(mesh, face_label, region_cost, options, data_term.zDim());

  const uint out_factor = options.output_factor_;

  segmentation.resize(xDim*out_factor,yDim*out_factor);

  mesh.enlarge(out_factor,out_factor);
  
  Math2D::Matrix<double> output(xDim*out_factor,yDim*out_factor,0);
  for (uint i=0; i < mesh.nFaces(); ++i) {
    add_grid_output(i,face_label[i],mesh,output);	
  }
  for (uint y=0; y < yDim*out_factor; y++) {
    for (uint x=0; x < xDim*out_factor; x++) {
      segmentation(x,y) = uint(output(x,y)*(255 / (nRegions-1)));
    }
  }

  return energy;
}

double lp_segment_pottscurvreg(const Math3D::Tensor<float>& data_term, const LPSegOptions& options, Math2D::Matrix<uint>& segmentation) {

  std::cerr << "CURVATURE POTTS MODEL" << std::endl;

  double lambda = options.lambda_;
  double gamma = options.gamma_;
  int neighborhood = options.neighborhood_;
  bool enforce_consistent_boundaries = options.enforce_consistent_boundaries_;
  bool enforce_regionedge = options.enforce_regionedge_;
  bool bruckstein = options.bruckstein_;
  bool reduce_edge_pairs = options.reduce_edge_pairs_;
  
  if (options.light_constraints_) {
    std::cerr << "WARNING: light constraints currently not implemented" << std::endl;
  }

  uint light_factor = 2;

  assert(neighborhood <= 16); 

  uint xDim = uint( data_term.xDim() );
  uint yDim = uint( data_term.yDim() );
  uint nRegions = uint( data_term.zDim() );

  Mesh2D mesh;  
  create_mesh(options, data_term, 0, mesh);

  Storage1D<PixelFaceRelation> shares;
  Math1D::Vector<uint> share_start;
  compute_pixel_shares(mesh, xDim, yDim, shares, share_start);

  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  size_t nRemoved = 0;
  if (reduce_edge_pairs) {
    nRemoved = filter_edge_pairs(mesh, edge_pairs); 
    std::cerr << "removed " << nRemoved << " edge pairs." << std::endl;
  }


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

  //uint nStandardEntries = lp_descr.nEntries();

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

  //threshold, compute energy and post-process
  Math1D::Vector<uint> face_labeling(mesh.nFaces(),0);
  for (uint i=0; i < mesh.nFaces();++i) {
    
    uint arg_max = MAX_UINT;
    double max_val = 0.0;
    
    for (uint r=0; r < nRegions; r++) {
      
      double val = lp_solution[i*nRegions + r];
      if (val > max_val) {

	max_val = val;
	arg_max = r;
      }
    }
    
    face_labeling[i] = arg_max;
  }

  double thresh_energy = multi_curv_energy(mesh, face_labeling, cost, options, nRegions);
  std::cerr << "energy of thresholded solution: " << thresh_energy << std::endl;

  double icm_energy = multi_curv_icm(mesh, face_labeling, cost, options, nRegions);
  std::cerr << "energy of ICM solution: " << icm_energy << std::endl;

  /**** generate output ****/

  Math1D::Vector<int> labels(mesh.nFaces(),-1);
  
  uint out_factor = options.output_factor_;

  segmentation.resize(xDim*out_factor,yDim*out_factor);

  Math3D::Tensor<double> output(xDim*out_factor,yDim*out_factor,nRegions,0.0);
  mesh.enlarge(out_factor,out_factor);

  //re-compute pixel shares for the now larger mesh
  compute_pixel_shares(mesh, out_factor*xDim, out_factor*yDim, shares, share_start);

  uint seg_factor = 255 / (nRegions-1);

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

double factor_energy(const Mesh2DEdgePair& edge_pair, const Mesh2D& mesh, const LPSegOptions& options,
		     const Math1D::Vector<uint>& face_label, uint nLabels) {


  double lambda = options.lambda_;
  double gamma = options.gamma_;
  bool bruckstein = options.bruckstein_;

  double energy = 0.0;

  const uint first_edge = edge_pair.first_edge_idx_;
  const uint second_edge = edge_pair.second_edge_idx_;    
    
  double mul1 = (mesh.edge(first_edge).to_idx_ == edge_pair.common_point_idx_) ? 1.0 : -1.0;
  double mul2 = (mesh.edge(second_edge).from_idx_ == edge_pair.common_point_idx_) ? 1.0 : -1.0;
    
  const std::vector<uint>& adj1 = mesh.adjacent_faces(first_edge);
  const std::vector<uint>& adj2 = mesh.adjacent_faces(second_edge);

  double curvature_weight = gamma * curv_weight(mesh,edge_pair,2.0,bruckstein);
  
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
    return 0.0;
  
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
    
    Storage1D< Math3D::Tensor<float> > cost(nLabels);
    for (uint i=0; i < nLabels; i++)
      cost[i].resize(nLabels,nLabels,nLabels,0.0);

    if (nLabels == 2) {
      cost[0](1,0,1) = pair_cost;
      cost[1](0,1,0) = pair_cost;
    }
    else {
      for (uint l=0; l < nLabels; l++) {

	for (uint ll=0; ll < nLabels; ll++) {
	  for (uint lll=0; lll < nLabels; lll++) {
	    
	    if (ll != l && lll != l) {
	      cost[l](ll,l,lll) += 0.5 * pair_cost;
	    } 
	  }
	}
      }
    }

    
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
  }
  else if (adjacent_regions.size() == 3) {

    Math3D::Tensor<float> cost(nLabels,nLabels,nLabels,0.0);
      
    if (adj1.size() == 1) {
      //boundary edge pair
	
      uint v1 = adj1[0];
      double m1 = match1[v1]*mul1;

      uint v2 = adj2[0];
      uint v3 = adj2[1];

      if (match2[v2]*mul2 != m1)
	std::swap(v2,v3);

      if (nLabels == 2) 
	cost(1,1,0) = pair_cost;
      else {

	for (uint l=0; l < nLabels; l++) {
	  for (uint ll=0; ll < nLabels; ll++) {
	    if (ll != l)
	      cost(l,l,ll) = pair_cost;
	  }
	}
      }

      energy += cost(face_label[v1],face_label[v2],face_label[v3]);
    }
    else if (adj2.size() == 1) {

      //boundary edge pair
	
      uint v1 = adj2[0];
      double m2 = match2[v1]*mul2;
	
      uint v2 = adj1[0];
      uint v3 = adj1[1];

      if (match1[v2]*mul1 != m2)
	std::swap(v2,v3);

      if (nLabels == 2) 
	cost(1,1,0) = pair_cost;
      else {
	for (uint l=0; l < nLabels; l++) {
	  for (uint ll=0; ll < nLabels; ll++) {
	    if (ll != l)
	      cost(ll,l,l) = pair_cost;
	  }
	}
      }

      energy += cost(face_label[v1],face_label[v2],face_label[v3]);
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

      if (nLabels == 2) {
	cost(0,1,1) = pair_cost;
	cost(1,0,0) = pair_cost;
      }
      else {
	for (uint l=0; l < nLabels; l++) {
	  for (uint ll=0; ll < nLabels; ll++) {
	    if (ll != l) {
	      cost(ll,l,l) = pair_cost;
	      cost(l,ll,ll) = pair_cost;
	    }
	  }
	}
      }
	
      assert(shared < MAX_UINT);
      assert(v1 < MAX_UINT);
      assert(v2 < MAX_UINT);	  
      	
      energy += cost(face_label[shared],face_label[v1],face_label[v2]);
    }
  }
  else if (adjacent_regions.size() == 2) {

    Math2D::Matrix<float> cost(nLabels,nLabels,0.0);

    if (adj1.size() == 1 && adj2.size() == 1) {

      double val1 = match1[adj1[0]] * mul1;
      double val2 = match2[adj2[0]] * mul2;
	
      if (val1 == val2) { //both regions are on the same side of the line pair
	if (nLabels == 2)
	  cost(1,1) = pair_cost;
      }
      else {
	std::cerr << "should this happen???" << std::endl;
	assert(false);
	cost(0,1) = pair_cost;
	cost(1,0) = pair_cost;
      }
    }
    else {

      //there should only be one shared region
      assert(adj1.size() == 1 || adj2.size() == 1);

      if (nLabels == 2) {
	if (find(adj1.begin(),adj1.end(),adjacent_regions[0]) != adj1.end()
	    && find(adj2.begin(),adj2.end(),adjacent_regions[0]) != adj2.end()) {
	  cost(1,0) = pair_cost;
	}
	else
	  cost(0,1) = pair_cost;
      }
    }

    energy += cost(face_label[adjacent_regions[0]],face_label[adjacent_regions[1]]);
  }

  return energy;
}

double factor_curv_energy(Mesh2D& mesh, const std::vector<Mesh2DEdgePair>& edge_pairs,
			  const Math2D::Matrix<float>& data_term, const LPSegOptions& options,
			  const Math1D::Vector<uint>& face_label) {

  double energy = 0.0;

  for (uint k=0; k < face_label.size(); k++) {

    if (face_label[k] == 1)
      energy += calculate_data_term(k, mesh, data_term);
  }

  // process factors
  for (uint p=0; p < edge_pairs.size(); p++) {
    
    energy += factor_energy(edge_pairs[p],mesh,options,face_label,2);
  }

  return energy;
}

double factor_curv_energy(Mesh2D& mesh, const std::vector<Mesh2DEdgePair>& edge_pairs,
			  const Math3D::Tensor<float>& data_term, const LPSegOptions& options,
			  const Math1D::Vector<uint>& face_label) {

  double energy = 0.0;

  for (uint k=0; k < face_label.size(); k++) {

    energy += calculate_data_term(k, face_label[k], mesh, data_term);
  }

  // process factors
  for (uint p=0; p < edge_pairs.size(); p++) {
    
    energy += factor_energy(edge_pairs[p],mesh,options,face_label,data_term.zDim());
  }

  return energy;
}

double factor_curv_icm(Mesh2D& mesh, const std::vector<Mesh2DEdgePair>& edge_pairs,
		       const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset,
		       Math1D::Vector<uint>& face_label) {

  double energy = 0.0;
  
  for (uint f=0; f < mesh.nFaces(); f++) {

    if (face_label[f] == 1)
      energy += calculate_data_term(f, mesh, data_term);
  }

  NamedStorage1D<std::vector<uint> > point_pair(mesh.nPoints(),MAKENAME(point_pair));

  for (uint e=0; e < edge_pairs.size(); e++) {

    uint point = edge_pairs[e].common_point_idx_;

    point_pair[point].push_back(e);
  }

  //1.) compute initial energy
  Math1D::Vector<double> cur_factor_energy(edge_pairs.size(),0.0);

  for (uint p=0; p < edge_pairs.size(); p++) {

    cur_factor_energy[p] = factor_energy(edge_pairs[p],mesh,options,face_label,2);
    energy += cur_factor_energy[p];
  }

  std::cerr << "initial energy: " << energy << std::endl;

  for (uint iter = 1; iter <= 5; iter++) {

    double last_energy = energy;

    std::cerr << "**** ICM iteration " << iter << " ******" << std::endl;

    for (uint f = 0; f < mesh.nFaces(); f++) {

      double cur_contrib = 0.0;

      uint org_label = face_label[f];

      if (org_label == 1)
	cur_contrib += calculate_data_term(f, mesh, data_term);

      //temporarily change label to compute the cost
      face_label[f] = 1 - org_label;
      double hyp_contrib = 0.0;
      if (face_label[f] == 1) 
	hyp_contrib += calculate_data_term(f, mesh, data_term);

      std::vector<uint> point_indices;
      mesh.get_polygon_points(f, point_indices);
      
      for (uint k=0; k < point_indices.size(); k++) {
	for (uint l=0; l < point_pair[k].size(); l++) {

	  uint pair = point_pair[k][l];
	  cur_contrib += cur_factor_energy[pair];
	  hyp_contrib += factor_energy(edge_pairs[pair],mesh,options,face_label,2);
	}
      }

      if (hyp_contrib >= cur_contrib)
	face_label[f] = org_label;
      else {
	//maybe better to store these values above
	for (uint k=0; k < point_indices.size(); k++) {
	  for (uint l=0; l < point_pair[k].size(); l++) {
	    
	    uint pair = point_pair[k][l];
	    cur_factor_energy[pair] = factor_energy(edge_pairs[pair],mesh,options,face_label,2);
	  }
	}

	energy += hyp_contrib - cur_contrib;
      }

    }

    std::cerr << "energy: " << energy << " (" << (energy + energy_offset) << ")" << std::endl;
    if (last_energy == energy)
      break;
  }

  return energy;
}


double factor_curv_icm(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset,
		       Math2D::Matrix<uint>& segmentation) {

  //Get the options
  int neighborhood = options.neighborhood_;
  
  assert(neighborhood <= 16); 

  uint xDim = uint( data_term.xDim() );
  uint yDim = uint( data_term.yDim() );

  Mesh2D mesh;  
  create_mesh(options, data_term, 0, mesh);

  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  Math1D::Vector<uint> face_label(mesh.nFaces(),0);

  for (uint f = 0; f < mesh.nFaces(); f++) {

    double cur_dt = calculate_data_term(f, mesh, data_term);

    if (cur_dt < 0)
      face_label[f] = 1;
  }

  double energy = factor_curv_icm(mesh, edge_pairs, data_term, options, energy_offset, face_label);

  const uint out_factor = options.output_factor_;

  segmentation.resize(xDim*out_factor,yDim*out_factor);

  mesh.enlarge(out_factor,out_factor);
  
  Math2D::Matrix<double> output(xDim*out_factor,yDim*out_factor,0);
  for (uint i=0; i < mesh.nFaces(); ++i) {
    add_grid_output(i,face_label[i],mesh,output);	
  }
  for (uint y=0; y < yDim*out_factor; y++) {
    for (uint x=0; x < xDim*out_factor; x++) {
      segmentation(x,y) = uint(output(x,y)*255.0);
    }
  }

  return energy;
}


double factor_curv_icm(Mesh2D& mesh, const std::vector<Mesh2DEdgePair>& edge_pairs,
		       const Math3D::Tensor<float>& data_term, const LPSegOptions& options,
		       Math1D::Vector<uint>& face_label) {

  uint nLabels = data_term.zDim();

  double energy = 0.0;

  for (uint k=0; k < face_label.size(); k++) {

    energy += calculate_data_term(k, face_label[k], mesh, data_term);
  }

  NamedStorage1D<std::vector<uint> > point_pair(mesh.nPoints(),MAKENAME(point_pair));

  for (uint e=0; e < edge_pairs.size(); e++) {

    uint point = edge_pairs[e].common_point_idx_;

    point_pair[point].push_back(e);
  }

  //1.) compute initial energy
  Math1D::Vector<double> cur_factor_energy(edge_pairs.size(),0.0);

  for (uint p=0; p < edge_pairs.size(); p++) {

    cur_factor_energy[p] = factor_energy(edge_pairs[p],mesh,options,face_label,nLabels);
    energy += cur_factor_energy[p];
  }

  std::cerr << "initial energy: " << energy << std::endl;

  for (uint iter = 1; iter <= 5; iter++) {

    double last_energy = energy;

    std::cerr << "**** ICM iteration " << iter << " ******" << std::endl;

    for (uint f = 0; f < mesh.nFaces(); f++) {

      for (uint l=0; l < nLabels; l++) {
	uint org_label = face_label[f];
	
	if (l == org_label)
	  continue;

	double cur_contrib = 0.0;
	
	cur_contrib += calculate_data_term(f, org_label, mesh, data_term);
	
	//temporarily change label to compute the cost
	face_label[f] = l;
	double hyp_contrib = 0.0;
	hyp_contrib += calculate_data_term(f, l, mesh, data_term);
	
	std::vector<uint> point_indices;
	mesh.get_polygon_points(f, point_indices);
	
	for (uint k=0; k < point_indices.size(); k++) {
	  for (uint l=0; l < point_pair[k].size(); l++) {
	    
	    uint pair = point_pair[k][l];
	    cur_contrib += cur_factor_energy[pair];
	    hyp_contrib += factor_energy(edge_pairs[pair],mesh,options,face_label,nLabels);
	  }
	}

	if (hyp_contrib >= cur_contrib)
	  face_label[f] = org_label;
	else {
	  //maybe better to store these values above
	  for (uint k=0; k < point_indices.size(); k++) {
	    for (uint l=0; l < point_pair[k].size(); l++) {
	      
	      uint pair = point_pair[k][l];
	      cur_factor_energy[pair] = factor_energy(edge_pairs[pair],mesh,options,face_label,nLabels);
	    }
	  }

	  energy += hyp_contrib - cur_contrib;
	}
      }

    }

    std::cerr << "energy: " << energy << std::endl;
    if (last_energy == energy)
      break;
  }

  return energy;
}

double factor_curv_icm(const Math3D::Tensor<float>& data_term, const LPSegOptions& options,
		       Math2D::Matrix<uint>& segmentation) {

  //Get the options
  int neighborhood = options.neighborhood_;
  
  assert(neighborhood <= 16); 

  uint xDim = uint( data_term.xDim() );
  uint yDim = uint( data_term.yDim() );

  Mesh2D mesh;  
  create_mesh(options, data_term, 0, mesh);

  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  Math1D::Vector<uint> face_label(mesh.nFaces(),0);

  for (uint f=0; f < mesh.nFaces(); f++) {
    
    double best = 1e300;
    uint best_label = 0;

    for (uint l=0; l < data_term.zDim(); l++) {
      double hyp = calculate_data_term(f, l, mesh, data_term);

      if (hyp < best) {
	best = hyp;
	best_label = l;
      }
    }
    
    face_label[f] = best_label;
  }

  double energy = factor_curv_icm(mesh, edge_pairs, data_term, options, face_label);

  const uint out_factor = options.output_factor_;

  segmentation.resize(xDim*out_factor,yDim*out_factor);

  mesh.enlarge(out_factor,out_factor);
  
  Math2D::Matrix<double> output(xDim*out_factor,yDim*out_factor,0);
  for (uint i=0; i < mesh.nFaces(); ++i) {
    add_grid_output(i,face_label[i],mesh,output);	
  }
  for (uint y=0; y < yDim*out_factor; y++) {
    for (uint x=0; x < xDim*out_factor; x++) {
      segmentation(x,y) = uint(output(x,y)*(255  / (data_term.zDim() - 1)));
    }
  }

  return energy;
}


double factor_lp_segment_curvreg(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset, 
				 Math2D::Matrix<uint>& segmentation, const Math2D::Matrix<int>* fixed_labels) {
  
  //Get the options
  double lambda = options.lambda_;
  double gamma = options.gamma_;
  int neighborhood = options.neighborhood_;
  bool bruckstein = options.bruckstein_;
  bool reduce_edge_pairs = options.reduce_edge_pairs_;

  assert(neighborhood <= 16); 

  uint xDim = uint( data_term.xDim() );
  uint yDim = uint( data_term.yDim() );

  Mesh2D mesh;  
  create_mesh(options, data_term, fixed_labels, mesh);

  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  size_t nRemoved = 0;
  if (reduce_edge_pairs) {
    nRemoved = filter_edge_pairs(mesh, edge_pairs); 
    std::cerr << "removed " << nRemoved << " edge pairs." << std::endl;
  }

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

  if (lambda > 0.0) {
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

      //note: if there is no direction change, we will never get here (see continue-statement above)

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

      //note: if there is no direction change, we will never get here (see continue-statement above)
      
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

  lp_descr.reset(0);
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

  double disc_energy = factor_curv_energy(mesh, edge_pairs, data_term, options, face_uint_solution)
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

double factor_lp_segment_pottscurvreg(const Math3D::Tensor<float>& data_term, const LPSegOptions& options, 
				      Math2D::Matrix<uint>& segmentation) {


  //Get the options
  double lambda = options.lambda_;
  double gamma = options.gamma_;
  int neighborhood = options.neighborhood_;
  bool bruckstein = options.bruckstein_;
  bool reduce_edge_pairs = options.reduce_edge_pairs_;

  lambda *= 0.5;
  gamma *= 0.5;

  assert(neighborhood <= 16); 

  uint xDim = uint( data_term.xDim() );
  uint yDim = uint( data_term.yDim() );
  uint nLabels = data_term.zDim();

  Mesh2D mesh;  
  create_mesh(options, data_term, 0, mesh);

  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  size_t nRemoved = 0;
  if (reduce_edge_pairs) {
    nRemoved = filter_edge_pairs(mesh, edge_pairs); 
    std::cerr << "removed " << nRemoved << " edge pairs." << std::endl;
  }

  std::cerr << edge_pairs.size() << " edge pairs." << std::endl;

  uint nVars = nLabels*mesh.nFaces() + nLabels*2*mesh.nEdges();  //region vars and length cliques
  uint nConstraints = nLabels*mesh.nEdges();
  nConstraints += mesh.nFaces(); // simplex constraints
  uint nEntries = nLabels*4*mesh.nEdges();

  Math1D::Vector<uint> pair_var_base(edge_pairs.size());
  Math1D::Vector<uint> pair_con_base(edge_pairs.size());

  uint abs_var_offs = nLabels*mesh.nFaces();

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
      uint nNewVars = nLabels*nLabels*nLabels*nLabels;
      nVars += nNewVars;
      nConstraints += 4*nLabels;
      nEntries += nNewVars*4 + 4*nLabels;

      nQuadCliques++;
    }
    else if (adjacent_regions.size() == 3) {
      uint nNewVars = nLabels*nLabels*nLabels;
      nVars += nNewVars;
      nConstraints += 3*nLabels;
      nEntries += nNewVars*3 + 3*nLabels;
      
      nTripleCliques++;
    }
    else if (adjacent_regions.size() == 2) {
      uint nNewVars = nLabels*nLabels;
      nVars += nNewVars; 
      nConstraints += 2*nLabels; 
      nEntries += nNewVars*2 + 2*nLabels;
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
    for (uint r = 0; r < nLabels; r++) 
      cost[i*nLabels+r] = calculate_data_term(i, r, mesh, data_term);
  }


  for (uint v=nLabels*mesh.nFaces(); v < nLabels*(mesh.nFaces()+2*mesh.nEdges()); v++) {
    if (mesh.adjacent_faces((v-mesh.nFaces())/(2*nLabels)).size() == 2)
      cost[v] = lambda;
  }
  
  SparseMatrixDescription<double> lp_descr(nEntries, nConstraints, nVars);
  

  //add simplex constraints
  for (uint f=0; f < mesh.nFaces(); f++) {
    
    rhs[f] = 1.0;

    for (uint l=0; l < nLabels; l++)
      lp_descr.add_entry(f,f*nLabels+l,1.0);
  }

  //handle length regularity
  if (lambda > 0.0) {
    for (uint e = 0; e < mesh.nEdges(); e++) {
      
      if (mesh.adjacent_faces(e).size() == 2) { 
	for (uint l=0; l < nLabels; l++) {
	  for (uint f=0; f < mesh.adjacent_faces(e).size(); f++) {
	  
	    uint var = mesh.adjacent_faces(e)[f];
	    lp_descr.add_entry(mesh.nFaces()+e, nLabels*var+l, mesh.match(var,e));
	  }
	
	  lp_descr.add_entry(mesh.nFaces()+e,abs_var_offs+2*(nLabels*e+l),1.0);
	  lp_descr.add_entry(mesh.nFaces()+e,abs_var_offs+2*(nLabels*e+l)+1,-1.0);
	}
      }
      else {
	//don't penalize the image border
      }
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

      uint nValues = nLabels*nLabels*nLabels*nLabels;

      //note: if there is no direction change, we will never get here (see continue-statement above)

      for (uint f=0; f < 4; f++) {

	for (uint l=0; l < nLabels; l++) {
	  lp_descr.add_entry(cur_con_base + nLabels*f+l, nLabels*adjacent_regions[f]+l, -1.0);
	}
      }
      
      for (uint v_addon = 0; v_addon < nValues; v_addon++) {

	Math1D::Vector<double> sum_e1(nLabels,0.0);
	Math1D::Vector<double> sum_e2(nLabels,0.0);

	uint inter_value = v_addon;

	for (uint f=0; f < 4; f++) {

	  uint l = inter_value % nLabels;

	  inter_value /= nLabels;

	  lp_descr.add_entry(cur_con_base + nLabels*f+l, cur_var_base + v_addon, 1.0);

	  sum_e1[l] += match1[adjacent_regions[f]];
	  sum_e2[l] += match2[adjacent_regions[f]];	  
	}

	sum_e1 *= mul1;
	sum_e2 *= mul2;

	cost[cur_var_base + v_addon] = 0.0;

	for (uint l=0; l < nLabels; l++) {
	  if (fabs(sum_e1[l]) > 0.5 && fabs(sum_e2[l]) > 0.5 && sum_e1[l]*sum_e2[l] > 0.0) {
	    cost[cur_var_base + v_addon] += boundary_cost;
	  }
	}
      }
    }
    else if (adjacent_regions.size() == 3) {

      uint nValues = nLabels*nLabels*nLabels;

      //note: if there is no direction change, we will never get here (see continue-statement above)

      for (uint f=0; f < 3; f++) {

	for (uint l=0; l < nLabels; l++) {
	  lp_descr.add_entry(cur_con_base + nLabels*f+l, nLabels*adjacent_regions[f]+l, -1.0);
	}
      }
      
      for (uint v_addon = 0; v_addon < nValues; v_addon++) {

	Math1D::Vector<double> sum_e1(nLabels,0.0);
	Math1D::Vector<double> sum_e2(nLabels,0.0);

	uint inter_value = v_addon;

	for (uint f=0; f < 3; f++) {

	  uint l = inter_value % nLabels;

	  inter_value /= nLabels;

	  lp_descr.add_entry(cur_con_base + nLabels*f+l, cur_var_base + v_addon, 1.0);

	  sum_e1[l] += match1[adjacent_regions[f]];
	  sum_e2[l] += match2[adjacent_regions[f]];	  
	}

	sum_e1 *= mul1;
	sum_e2 *= mul2;

	cost[cur_var_base + v_addon] = 0.0;

	for (uint l=0; l < nLabels; l++) {
	  if (fabs(sum_e1[l]) > 0.5 && fabs(sum_e2[l]) > 0.5 && sum_e1[l]*sum_e2[l] > 0.0) {
	    cost[cur_var_base + v_addon] += boundary_cost;
	  }
	}
      }
    }
    else {

      uint nValues = nLabels*nLabels;

      //note: if there is no direction change, we will never get here (see continue-statement above)

      for (uint f=0; f < 2; f++) {

	for (uint l=0; l < nLabels; l++) {
	  lp_descr.add_entry(cur_con_base + nLabels*f+l, nLabels*adjacent_regions[f]+l, -1.0);
	}
      }
      
      for (uint v_addon = 0; v_addon < nValues; v_addon++) {

	Math1D::Vector<double> sum_e1(nLabels,0.0);
	Math1D::Vector<double> sum_e2(nLabels,0.0);

	uint inter_value = v_addon;

	for (uint f=0; f < 2; f++) {

	  uint l = inter_value % nLabels;

	  inter_value /= nLabels;

	  lp_descr.add_entry(cur_con_base + nLabels*f+l, cur_var_base + v_addon, 1.0);

	  sum_e1[l] += match1[adjacent_regions[f]];
	  sum_e2[l] += match2[adjacent_regions[f]];	  
	}

	sum_e1 *= mul1;
	sum_e2 *= mul2;

	cost[cur_var_base + v_addon] = 0.0;

	for (uint l=0; l < nLabels; l++) {
	  if (fabs(sum_e1[l]) > 0.5 && fabs(sum_e2[l]) > 0.5 && sum_e1[l]*sum_e2[l] > 0.0) {
	    cost[cur_var_base + v_addon] += boundary_cost;
	  }
	}
      }
    }
  }

  Math1D::Vector<uint> row_start(nConstraints+1);
  lp_descr.sort_by_row(row_start,true);

  CoinPackedMatrix coinMatrix(false,(int*) lp_descr.row_indices(),(int*) lp_descr.col_indices(),
			      lp_descr.value(),lp_descr.nEntries());

  ClpSimplex lpSolver;
  lpSolver.loadProblem(coinMatrix, var_lb.direct_access(), var_ub.direct_access(),   
		       cost.direct_access(), rhs.direct_access(), rhs.direct_access());

  lp_descr.reset(0);
  coinMatrix.cleanMatrix();

  for (uint v=0; v < nVars; v++) {
    lpSolver.setInteger(v);
  }
  
  lpSolver.dual();

  const double* lp_solution = lpSolver.getColSolution();

  double lp_energy = lpSolver.getObjValue();

  std::cerr.precision(10);
  std::cerr << "original relaxed energy: " << (lp_energy) << std::endl;

  Math1D::Vector<uint> face_uint_solution(mesh.nFaces(),0);
  
  //extract solution
  Storage1D<PixelFaceRelation> shares;
  Math1D::Vector<uint> share_start;

  Math1D::Vector<int> labels(mesh.nFaces(),-1);
  
  uint out_factor = options.output_factor_;

  segmentation.resize(xDim*out_factor,yDim*out_factor);

  Math3D::Tensor<double> output(xDim*out_factor,yDim*out_factor,nLabels,0.0);
  mesh.enlarge(out_factor,out_factor);

  //re-compute pixel shares for the now larger mesh
  compute_pixel_shares(mesh, out_factor*xDim, out_factor*yDim, shares, share_start);

  uint seg_factor = 255 / (nLabels-1);

  for (uint y=0; y < yDim*out_factor; y++) {
    for (uint x=0; x < xDim*out_factor; x++) {

      Math1D::Vector<double> votes(nLabels);

      for (uint r=0; r < nLabels; r++) {
	double sum = 0.0;
	for (uint k= share_start[y*(xDim*out_factor)+x]; k < share_start[y*(xDim*out_factor)+x+1]; k++) {
	  uint face = shares[k].face_idx_;
	  double cur_share = shares[k].share_;
	  sum += mesh.convex_area(face) * cur_share * lp_solution[face*nLabels + r];
	}

	votes[r] = sum;
      }

      double max_val = 0.0;
      uint arg_max = MAX_UINT;

      for (uint r=0; r < nLabels; r++) {

	double val = votes[r];
	if (val > max_val) {
	  max_val = val;
	  arg_max = r;
	}
      }

      segmentation(x,y) = arg_max * seg_factor;
    }
  }

  return lp_energy; 
}

double factor_lp_segment_pottscurvreg_layered(const Math3D::Tensor<float>& data_term, const LPSegOptions& options, 
					      Math2D::Matrix<uint>& segmentation) {


  //Get the options
  double lambda = options.lambda_;
  double gamma = options.gamma_;
  int neighborhood = options.neighborhood_;
  bool bruckstein = options.bruckstein_;
  bool reduce_edge_pairs = options.reduce_edge_pairs_;

  lambda *= 0.5;
  gamma *= 0.5;

  assert(neighborhood <= 16); 

  uint xDim = uint( data_term.xDim() );
  uint yDim = uint( data_term.yDim() );
  uint nLabels = data_term.zDim();

  Mesh2D mesh;  
  create_mesh(options, data_term, 0, mesh);

  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  size_t nRemoved = 0;
  if (reduce_edge_pairs) {
    nRemoved = filter_edge_pairs(mesh, edge_pairs); 
    std::cerr << "removed " << nRemoved << " edge pairs." << std::endl;
  }


  std::cerr << edge_pairs.size() << " edge pairs." << std::endl;

  uint nVars = nLabels*mesh.nFaces() + nLabels*2*mesh.nEdges();  //region vars and length cliques
  uint nConstraints = nLabels*mesh.nEdges();
  nConstraints += mesh.nFaces(); // simplex constraints
  uint nEntries = nLabels*4*mesh.nEdges();

  Math1D::Vector<uint> pair_var_base(edge_pairs.size());
  Math1D::Vector<uint> pair_con_base(edge_pairs.size());

  uint abs_var_offs = nLabels*mesh.nFaces();

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
      nVars += nLabels*16;
      nConstraints += nLabels*4*2;
      nEntries += nLabels*(16*4 + 4*2);

      nQuadCliques++;
    }
    else if (adjacent_regions.size() == 3) {
      nVars += nLabels*8;
      nConstraints += nLabels*3*2;
      nEntries += nLabels*(8*3 + 3*2);
      
      nTripleCliques++;
    }
    else if (adjacent_regions.size() == 2) {
      nVars += nLabels*2; //here we use the construction via absolutes
      nConstraints += nLabels; 
      nEntries += nLabels*4;
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
    for (uint r = 0; r < nLabels; r++) 
      cost[i*nLabels+r] = calculate_data_term(i, r, mesh, data_term);
  }

  for (uint v=nLabels*mesh.nFaces(); v < nLabels*(mesh.nFaces()+2*mesh.nEdges()); v++) {
    if (mesh.adjacent_faces((v-mesh.nFaces())/(2*nLabels)).size() == 2)
      cost[v] = lambda;
  }
  
  SparseMatrixDescription<double> lp_descr(nEntries, nConstraints, nVars);  

  //add simplex constraints
  for (uint f=0; f < mesh.nFaces(); f++) {
    
    rhs[f] = 1.0;

    for (uint l=0; l < nLabels; l++)
      lp_descr.add_entry(f,f*nLabels+l,1.0);
  }

  //handle length regularity
  if (lambda > 0.0) {
    for (uint e = 0; e < mesh.nEdges(); e++) {
      
      if (mesh.adjacent_faces(e).size() == 2) { 
	for (uint l=0; l < nLabels; l++) {
	  for (uint f=0; f < mesh.adjacent_faces(e).size(); f++) {
	  
	    uint var = mesh.adjacent_faces(e)[f];
	    lp_descr.add_entry(mesh.nFaces()+e, nLabels*var+l, mesh.match(var,e));
	  }
	
	  lp_descr.add_entry(mesh.nFaces()+e,abs_var_offs+2*(nLabels*e+l),1.0);
	  lp_descr.add_entry(mesh.nFaces()+e,abs_var_offs+2*(nLabels*e+l)+1,-1.0);
	}
      }
      else {
	//don't penalize the image border
      }
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

      //note: if there is no direction change, we will never get here (see continue-statement above)

      for (uint f=0; f < 4; f++) {

	for (uint l=0; l < nLabels; l++) {
	  lp_descr.add_entry(cur_con_base + 8*l + 2*f, nLabels*adjacent_regions[f] + l, -1.0);
	  lp_descr.add_entry(cur_con_base + 8*l + 2*f+1, nLabels*adjacent_regions[f] + l, 1.0);
	  rhs[cur_con_base + 8*l + 2*f+1] = 1.0;
	}
      }
      
      for (uint v_addon = 0; v_addon < 16; v_addon++) {

	double sum_e1 = 0.0;
	double sum_e2 = 0.0;

	for (uint f=0; f < 4; f++) {

	  uint mask = 1 << f;
	  assert(mask != 0);

	  uint offs = 1 - ((v_addon & mask) >> f);
	  for (uint l=0; l < nLabels; l++) 
	    lp_descr.add_entry(cur_con_base + 8*l + 2*f+offs, cur_var_base + 16*l + v_addon, 1.0);

	  sum_e1 += (1-offs) * match1[adjacent_regions[f]];
	  sum_e2 += (1-offs) * match2[adjacent_regions[f]];	  
	}

	sum_e1 *= mul1;
	sum_e2 *= mul2;

	if (fabs(sum_e1) > 0.5 && fabs(sum_e2) > 0.5 && sum_e1*sum_e2 > 0.0) {
	  for (uint l=0; l < nLabels; l++) 
	    cost[cur_var_base + 16*l + v_addon] = boundary_cost;
	}
      }
    }
    else if (adjacent_regions.size() == 3) {

      //note: if there is no direction change, we will never get here (see continue-statement above)
      
      for (uint f=0; f < 3; f++) {
	
	for (uint l=0; l < nLabels; l++) {
	  lp_descr.add_entry(cur_con_base + 6*l + 2*f, nLabels*adjacent_regions[f]+l, -1.0);
	  lp_descr.add_entry(cur_con_base + 6*l + 2*f+1, nLabels*adjacent_regions[f]+l, 1.0);
	  rhs[cur_con_base + 6*l + 2*f+1] = 1.0;
	}
      }
      
      for (uint v_addon = 0; v_addon < 8; v_addon++) {

	double sum_e1 = 0.0;
	double sum_e2 = 0.0;

	for (uint f=0; f < 3; f++) {
	
	  uint mask = 1 << f;
	  assert(mask != 0);

	  uint offs = 1 - ((v_addon & mask) >> f);
	  for (uint l=0; l < nLabels; l++) 
	    lp_descr.add_entry(cur_con_base + 6*l + 2*f+offs, cur_var_base + 8*l + v_addon, 1.0);
	  
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
	  
	  for (uint l=0; l < nLabels; l++) 
	    cost[cur_var_base + 8*l + v_addon] = boundary_cost;
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
	for (uint l=0; l < nLabels; l++)  {
	  lp_descr.add_entry(cur_con_base+l, nLabels*adjacent_regions[0]+l, 1.0);
	  lp_descr.add_entry(cur_con_base+l, nLabels*adjacent_regions[1]+l, -1.0);
	
	  lp_descr.add_entry(cur_con_base+l, cur_var_base + 2*l, 1.0);
	  lp_descr.add_entry(cur_con_base+l, cur_var_base + 2*l +1, -1.0);
	  
	  cost[cur_var_base] = c01;
	  cost[cur_var_base+1] = c10;
	}
      }
    }
  }
  
  Math1D::Vector<uint> row_start(nConstraints+1);
  lp_descr.sort_by_row(row_start,true);

  CoinPackedMatrix coinMatrix(false,(int*) lp_descr.row_indices(),(int*) lp_descr.col_indices(),
			      lp_descr.value(),lp_descr.nEntries());

  ClpSimplex lpSolver;
  lpSolver.loadProblem(coinMatrix, var_lb.direct_access(), var_ub.direct_access(),   
		       cost.direct_access(), rhs.direct_access(), rhs.direct_access());

  lp_descr.reset(0);
  coinMatrix.cleanMatrix();

  for (uint v=0; v < nVars; v++) {
    lpSolver.setInteger(v);
  }
  
  lpSolver.dual();

  const double* lp_solution = lpSolver.getColSolution();

  double lp_energy = lpSolver.getObjValue();

  std::cerr.precision(10);
  std::cerr << "original relaxed energy: " << (lp_energy) << std::endl;

  Math1D::Vector<uint> face_uint_solution(mesh.nFaces(),0);
  
  //extract solution
  Storage1D<PixelFaceRelation> shares;
  Math1D::Vector<uint> share_start;

  Math1D::Vector<int> labels(mesh.nFaces(),-1);
  
  uint out_factor = options.output_factor_;

  segmentation.resize(xDim*out_factor,yDim*out_factor);

  Math3D::Tensor<double> output(xDim*out_factor,yDim*out_factor,nLabels,0.0);
  mesh.enlarge(out_factor,out_factor);

  //re-compute pixel shares for the now larger mesh
  compute_pixel_shares(mesh, out_factor*xDim, out_factor*yDim, shares, share_start);

  uint seg_factor = 255 / (nLabels-1);

  for (uint y=0; y < yDim*out_factor; y++) {
    for (uint x=0; x < xDim*out_factor; x++) {

      Math1D::Vector<double> votes(nLabels);

      for (uint r=0; r < nLabels; r++) {
	double sum = 0.0;
	for (uint k= share_start[y*(xDim*out_factor)+x]; k < share_start[y*(xDim*out_factor)+x+1]; k++) {
	  uint face = shares[k].face_idx_;
	  double cur_share = shares[k].share_;
	  sum += mesh.convex_area(face) * cur_share * lp_solution[face*nLabels + r];
	}

	votes[r] = sum;
      }

      double max_val = 0.0;
      uint arg_max = MAX_UINT;

      for (uint r=0; r < nLabels; r++) {

	double val = votes[r];
	if (val > max_val) {
	  max_val = val;
	  arg_max = r;
	}
      }

      segmentation(x,y) = arg_max * seg_factor;
    }
  }

  return lp_energy;
}



double factor_lp_segment_curvreg_minsum_diffusion(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset, 
						  Math2D::Matrix<uint>& segmentation, const Math2D::Matrix<int>* fixed_labels) {

  //NOTE: fixed_labels is currently ignored

  //Get the options
  double lambda = options.lambda_;
  double gamma = options.gamma_;
  int neighborhood = options.neighborhood_;
  bool bruckstein = options.bruckstein_;
  bool reduce_edge_pairs = options.reduce_edge_pairs_;
  
  assert(neighborhood <= 16); 

  uint xDim = uint( data_term.xDim() );
  uint yDim = uint( data_term.yDim() );

  Mesh2D mesh;  
  create_mesh(options, data_term, fixed_labels, mesh);

  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  size_t nRemoved = 0;
  if (reduce_edge_pairs) {
    nRemoved = filter_edge_pairs(mesh, edge_pairs); 
    std::cerr << "removed " << nRemoved << " edge pairs." << std::endl;
  }

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

  double disc_energy = factor_curv_energy(mesh, edge_pairs, data_term, options, face_uint_solution)
    + energy_offset;
  
  std::cerr << "solution found by stand-alone routine: " << disc_energy << std::endl;

  double icm_energy = factor_curv_icm(mesh, edge_pairs, data_term, options, energy_offset, face_uint_solution)
    + energy_offset;

  std::cerr << "energy of ICM solution: " << icm_energy << std::endl;

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

  return lower_bound;
}

double factor_lp_segment_curvreg_minsum_diffusion_memsave(const Math3D::Tensor<float>& data_term, 
							  const LPSegOptions& options, double energy_offset, 
							  Math2D::Matrix<uint>& segmentation, 
							  const Math2D::Matrix<int>* fixed_labels) {
  //NOTE: fixed_labels is currently ignored

  uint nLabels = data_term.zDim();

  //Get the options
  double lambda = options.lambda_;
  double gamma = options.gamma_;
  int neighborhood = options.neighborhood_;
  bool bruckstein = options.bruckstein_;
  bool reduce_edge_pairs = options.reduce_edge_pairs_;

  assert(neighborhood <= 16); 

  uint xDim = uint( data_term.xDim() );
  uint yDim = uint( data_term.yDim() );

  Mesh2D mesh;  
  create_mesh(options, data_term, fixed_labels, mesh);

  Math2D::Matrix<double> var_cost(mesh.nFaces(),nLabels,0.0);
  for (uint f=0; f < mesh.nFaces(); f++) {
    
    for (uint l=0; l < nLabels; l++)
      var_cost(f,l) = calculate_data_term(f, l, mesh, data_term);
  }
  
  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  size_t nRemoved = 0;
  if (reduce_edge_pairs) {
    nRemoved = filter_edge_pairs(mesh, edge_pairs); 
    std::cerr << "removed " << nRemoved << " edge pairs." << std::endl;
  }

  std::cerr << edge_pairs.size() << " edge pairs." << std::endl;

  Math1D::Vector<float> factor_weight(edge_pairs.size());
  Math1D::Vector<uint> repar_start(edge_pairs.size()+1,MAX_UINT);
  Math1D::Vector<uchar> factor_type(edge_pairs.size(),255);

  uint repar_size = 0;

  for (uint p=0; p < edge_pairs.size(); p++) {

    const uint first_edge = edge_pairs[p].first_edge_idx_;
    const uint second_edge = edge_pairs[p].second_edge_idx_;    

    const std::vector<uint>& adj1 = mesh.adjacent_faces(first_edge);
    const std::vector<uint>& adj2 = mesh.adjacent_faces(second_edge);

    double curvature_weight = gamma * curv_weight(mesh,edge_pairs[p],2.0,bruckstein);
    if (adj1.size() == 1 && adj2.size() == 1)
      curvature_weight = 0.0; //don't penalize the image corners

    double boundary_cost = curvature_weight;

    //no length for the image border
    if (adj1.size() > 1)
      boundary_cost += 0.5*lambda*mesh.edge_length(first_edge); 
    if (adj2.size() > 1)
      boundary_cost += 0.5*lambda*mesh.edge_length(second_edge);

    factor_weight[p] = boundary_cost;

    if (fabs(factor_weight[p]) < 0.001)
      continue;

    std::set<uint> adjacent_regions_set;

    for (std::vector<uint>::const_iterator it = adj1.begin(); it != adj1.end(); it++) {
      adjacent_regions_set.insert(*it);
    }

    for (std::vector<uint>::const_iterator it = adj2.begin(); it != adj2.end(); it++) {
      adjacent_regions_set.insert(*it);
    }

    repar_size += adjacent_regions_set.size();
  }

  Math1D::Vector<uint> factor_face(repar_size);
  Math1D::Vector<double> repar(nLabels*repar_size,0.0);

  uint cur_pos = 0;

  Storage1D< Math3D::Tensor<float> > cost0(nLabels);
  for (uint i=0; i < nLabels; i++)
    cost0[i].resize(nLabels,nLabels,nLabels,0.0);
  
  if (nLabels == 2) {
    cost0[0](1,0,1) = 1.0;
    cost0[1](0,1,0) = 1.0;
  }
  else {
    for (uint l=0; l < nLabels; l++) {

      for (uint ll=0; ll < nLabels; ll++) {
	for (uint lll=0; lll < nLabels; lll++) {

	  if (ll != l && lll != l) {
	    cost0[l](ll,l,lll) += 0.5;
	  }
	}
      }
    }
  }

  Math3D::Tensor<float> cost1(nLabels,nLabels,nLabels,0.0);
  if (nLabels == 2) {
    cost1(1,1,0) = 1.0;
  }
  else {
    for (uint l=0; l < nLabels; l++) {
      for (uint ll=0; ll < nLabels; ll++) {
	if (ll != l)
	  cost1(l,l,ll) = 1.0;
      }
    }
  }

  Math3D::Tensor<float> cost2(nLabels,nLabels,nLabels,0.0);
  if (nLabels == 2) {
    cost2(0,1,1) = 1.0;
    cost2(1,0,0) = 1.0;
  }
  else {
    for (uint l=0; l < nLabels; l++) {
      for (uint ll=0; ll < nLabels; ll++) {
	if (ll != l) {
	  cost2(ll,l,l) = 1.0;
	  cost2(l,ll,ll) = 1.0;
	}
      }
    }
  }

  Math2D::Matrix<float> cost3(nLabels,nLabels,0.0);
  if (nLabels == 2)
    cost3(1,1) = 1.0;
  else {
    ; //should there be cost here for multi-region?
  }

  Math2D::Matrix<float> cost4(nLabels,nLabels,0.0);
  if (nLabels == 2) {
    cost4(0,1) = 1.0;
    cost4(1,0) = 1.0;
  }
  else {
    for (uint l=0; l < nLabels; l++) {
      for (uint ll=0; ll < nLabels; ll++) {
	if (ll != l) {
	  cost4(l,ll) = 1.0;
	}
      }
    }
  }

  Math2D::Matrix<float> cost5(nLabels,nLabels,0.0);
  if (nLabels == 2)
    cost5(1,0) = 1.0;
  else {
    ; //should there be cost here for multi-region?
  }

  double bound = 1e300;

  for (uint p=0; p < edge_pairs.size(); p++) {

    repar_start[p] = cur_pos;

    if (fabs(factor_weight[p]) < 0.001)
      continue;

    //std::cerr << "p: " << p << std::endl;
    
    const uint first_edge = edge_pairs[p].first_edge_idx_;
    const uint second_edge = edge_pairs[p].second_edge_idx_;    
    
    double mul1 = (mesh.edge(first_edge).to_idx_ == edge_pairs[p].common_point_idx_) ? 1.0 : -1.0;
    double mul2 = (mesh.edge(second_edge).from_idx_ == edge_pairs[p].common_point_idx_) ? 1.0 : -1.0;
    
    const std::vector<uint>& adj1 = mesh.adjacent_faces(first_edge);
    const std::vector<uint>& adj2 = mesh.adjacent_faces(second_edge);
    
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
    
    //now really add the factors
    if (adjacent_regions.size() == 4) {
      
      assert(adj1.size() == 2);
      assert(adj2.size() == 2);

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

      factor_face[cur_pos + 0] = v1;
      factor_face[cur_pos + 1] = v2;
      factor_face[cur_pos + 2] = v3;
      factor_face[cur_pos + 3] = v4;

      factor_type[p] = 0;
      cur_pos += 4;
    }
    else if (adjacent_regions.size() == 3) {

      if (adj1.size() == 1) {
	//boundary edge pair
	
	//std::cerr << "A" << std::endl;
	
	uint v1 = adj1[0];
	double m1 = match1[v1]*mul1;

	uint v2 = adj2[0];
	uint v3 = adj2[1];

	if (match2[v2]*mul2 != m1)
	  std::swap(v2,v3);

	factor_face[cur_pos + 0] = v1;
	factor_face[cur_pos + 1] = v2;
	factor_face[cur_pos + 2] = v3;
	
	factor_type[p] = 1;
      }
      else if (adj2.size() == 1) {

	//boundary edge pair
	
	uint v1 = adj2[0];
	double m2 = match2[v1]*mul2;
	
	uint v2 = adj1[0];
	uint v3 = adj1[1];

	if (match1[v2]*mul1 != m2)
	    std::swap(v2,v3);

	//this is symmetric to the above case

	factor_face[cur_pos + 0] = v1;
	factor_face[cur_pos + 1] = v2;
	factor_face[cur_pos + 2] = v3;
	
	factor_type[p] = 1;
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
	
	assert(shared < MAX_UINT);
	assert(v1 < MAX_UINT);
	assert(v2 < MAX_UINT);	  

	factor_face[cur_pos + 0] = shared;
	factor_face[cur_pos + 1] = v1;
	factor_face[cur_pos + 2] = v2;

	factor_type[p] = 2;	
      }

      cur_pos += 3;
    }
    else if (adjacent_regions.size() == 2) {

      uint v1 = adjacent_regions[0];
      uint v2 = adjacent_regions[0];
      

      if (adj1.size() == 1 && adj2.size() == 1) {

	double val1 = match1[adj1[0]] * mul1;
	double val2 = match2[adj2[0]] * mul2;
	
	if (val1 == val2) { //both regions are on the same side of the line pair
	  factor_type[p] = 3;
	}
	else {
	  std::cerr << "should this happen???" << std::endl;
	  factor_type[p] = 4;
	}
      }
      else {

	//std::cerr << "should this happen?" << std::endl;
	//std::cerr << "adj1: " << adj1 << std::endl;
	//std::cerr << "adj2: " << adj2 << std::endl;

	//there should only be one shared region
	assert(adj1.size() == 1 || adj2.size() == 1);

	factor_type[p] = 5;

	if (find(adj1.begin(),adj1.end(),adjacent_regions[0]) != adj1.end()
	    && find(adj2.begin(),adj2.end(),adjacent_regions[0]) != adj2.end()) {
	  //order is already correct
	}
	else
	  std::swap(v1,v2);
      }

      factor_face[cur_pos + 0] = v1;
      factor_face[cur_pos + 1] = v2;

      cur_pos += 2;
    }
  }

  Math1D::Vector<uint> face_uint_solution(mesh.nFaces(),0);
  
  for (uint iter=1; iter <= 100; iter++) {

    std::cerr << "**** iter: " << iter << " *****" << std::endl;

    //a) compute current cost
    double cur_bound = 0.0;
    for (uint f=0; f < mesh.nFaces(); f++) {
      double best = 1e300;
      for (uint l=0; l < nLabels; l++) {

	double hyp = var_cost(f,l);
	
	if (hyp < best) {
	  best = hyp;
	  face_uint_solution[f] = l;
	}
      }

      cur_bound += best;
    }

    for (uint p=0; p < edge_pairs.size(); p++) {

      uint start_pos = repar_start[p];

      uchar cur_type = factor_type[p];

      double* cur_repar = repar.direct_access() + nLabels*start_pos;

      double cur_weight = factor_weight[p];

      double best = 1e300;

      if (cur_type == 0) {
	//4-clique

	for (uint l1=0; l1 < nLabels; l1++) {
	  for (uint l2=0; l2 < nLabels; l2++) {
	    for (uint l3=0; l3 < nLabels; l3++) {
	      for (uint l4=0; l4 < nLabels; l4++) {

		double hyp= cur_weight * cost0[l1](l2,l3,l4)
		  + cur_repar[l1] + cur_repar[nLabels + l2] + cur_repar[2*nLabels + l3] + cur_repar[3*nLabels + l4];
		
		if (hyp < best)
		  best = hyp;
	      }
	    }
	  }
	} 

      }
      else if (cur_type <= 2) {
	//3-clique
	const Math3D::Tensor<float>& cost = (cur_type == 1) ? cost1 : cost2;


	for (uint l1=0; l1 < nLabels; l1++) {
	  for (uint l2=0; l2 < nLabels; l2++) {
	    for (uint l3=0; l3 < nLabels; l3++) {

	      double hyp= cur_weight * cost(l1,l2,l3)
		+ cur_repar[l1] + cur_repar[nLabels + l2] + cur_repar[2*nLabels + l3];

	      if (hyp < best)
		best = hyp;
	    }
	  }
	} 
	
      }
      else if (cur_type <= 5) {
	//2-clique
	const Math2D::Matrix<float>& cost = (cur_type == 3) ? cost3 : ((cur_type == 4) ? cost4 : cost5);


	for (uint l1=0; l1 < nLabels; l1++) {
	  for (uint l2=0; l2 < nLabels; l2++) {

	    double hyp= cur_weight * cost(l1,l2)
	      + cur_repar[l1] + cur_repar[nLabels + l2];

	    if (hyp < best)
	      best = hyp;
	  }
	} 
      }
      else {
	//std::cerr << "WARNING: unknown type " << cur_type << std::endl;
	best = 0.0;
      }

      cur_bound += best;
    }

    std::cerr << "bound: " << cur_bound << std::endl;
    bound = cur_bound;

    // b) reparamterize

    for (uint p=0; p < edge_pairs.size(); p++) {

      uint start_pos = repar_start[p];

      uchar cur_type = factor_type[p];

      double* cur_repar = repar.direct_access() + nLabels*start_pos;

      double cur_weight = factor_weight[p];

      if (cur_type == 0) {
	//4-clique

	uint v1 = factor_face[start_pos];
	uint v2 = factor_face[start_pos+1];
	uint v3 = factor_face[start_pos+2];
	uint v4 = factor_face[start_pos+3];

	//repar v1
	for (uint l1=0; l1 < nLabels; l1++) {

	  double b = var_cost(v1,l1);
	  double a = 1e300;
	    
	  for (uint l2=0; l2 < nLabels; l2++) {
	    for (uint l3=0; l3 < nLabels; l3++) {
	      for (uint l4=0; l4 < nLabels; l4++) {
	  
		double hyp = cur_weight * cost0[l1](l2,l3,l4)
		  + cur_repar[nLabels + l2] + cur_repar[2*nLabels + l3] + cur_repar[3*nLabels + l4];

		if (hyp < a)
		  a = hyp;
	      }
	    }
	  }

	  a += cur_repar[l1];

	  double shift = 0.5*(b-a);
	  var_cost(v1,l1) -= shift;
	  cur_repar[l1] += shift;
	}

	//repar v2
	for (uint l2=0; l2 < nLabels; l2++) {

	  double b = var_cost(v2,l2);
	  double a = 1e300;

	  for (uint l1=0; l1 < nLabels; l1++) {
	    for (uint l3=0; l3 < nLabels; l3++) {
	      for (uint l4=0; l4 < nLabels; l4++) {
	  
		double hyp = cur_weight * cost0[l1](l2,l3,l4)
		  + cur_repar[l1] + cur_repar[2*nLabels + l3] + cur_repar[3*nLabels + l4];

		if (hyp < a)
		  a = hyp;
	      }
	    }
	  }

	  a += cur_repar[nLabels + l2];

	  double shift = 0.5*(b-a);
	  var_cost(v2,l2) -= shift;
	  cur_repar[nLabels + l2] += shift;
	}

	//repar v3
	for (uint l3=0; l3 < nLabels; l3++) {

	  double b = var_cost(v3,l3);
	  double a = 1e300;

	  for (uint l1=0; l1 < nLabels; l1++) {
	    for (uint l2=0; l2 < nLabels; l2++) {
	      for (uint l4=0; l4 < nLabels; l4++) {
	  
		double hyp = cur_weight * cost0[l1](l2,l3,l4)
		  + cur_repar[l1] + cur_repar[nLabels + l2] + cur_repar[3*nLabels + l4];

		if (hyp < a)
		  a = hyp;
	      }
	    }
	  }

	  a += cur_repar[2*nLabels + l3];
	  
	  double shift = 0.5*(b-a);
	  var_cost(v3,l3) -= shift;
	  cur_repar[2*nLabels + l3] += shift;
	}

	//repar v4
	for (uint l4=0; l4 < nLabels; l4++) {

	  double b = var_cost(v4,l4);
	  double a = 1e300;

	  for (uint l1=0; l1 < nLabels; l1++) {
	    for (uint l2=0; l2 < nLabels; l2++) {
	      for (uint l3=0; l3 < nLabels; l3++) {
	  
		double hyp = cur_weight * cost0[l1](l2,l3,l4)
		  + cur_repar[l1] + cur_repar[nLabels + l2] + cur_repar[2*nLabels + l3];

		if (hyp < a)
		  a = hyp;
	      }
	    }
	  }

	  a += cur_repar[3*nLabels + l4];

	  double shift = 0.5*(b-a);
	  var_cost(v4,l4) -= shift;
	  cur_repar[3*nLabels + l4] += shift;
	}
      }
      else if (cur_type <= 2) {
	//3-clique
	const Math3D::Tensor<float>& cost = (cur_type == 1) ? cost1 : cost2;

	uint v1 = factor_face[start_pos];
	uint v2 = factor_face[start_pos+1];
	uint v3 = factor_face[start_pos+2];
	
	//repar v1
	for (uint l1=0; l1 < nLabels; l1++) {

	  double b = var_cost(v1,l1);
	  double a = 1e300;
	    
	  for (uint l2=0; l2 < nLabels; l2++) {
	    for (uint l3=0; l3 < nLabels; l3++) {
	  
	      double hyp = cur_weight * cost(l1,l2,l3)
	        + cur_repar[nLabels + l2] + cur_repar[2*nLabels + l3];
	      
	      if (hyp < a)
		a = hyp;
	    }
	  }

	  a += cur_repar[l1];

	  double shift = 0.5*(b-a);
	  var_cost(v1,l1) -= shift;
	  cur_repar[l1] += shift;
	}

	//repar v2
	for (uint l2=0; l2 < nLabels; l2++) {

	  double b = var_cost(v2,l2);
	  double a = 1e300;
	    
	  for (uint l1=0; l1 < nLabels; l1++) {
	    for (uint l3=0; l3 < nLabels; l3++) {
	  
	      double hyp = cur_weight * cost(l1,l2,l3)
		+ cur_repar[l1] + cur_repar[2*nLabels + l3];
	      
	      if (hyp < a)
		a = hyp;
	    }
	  }

	  a += cur_repar[nLabels + l2];

	  double shift = 0.5*(b-a);
	  var_cost(v2,l2) -= shift;
	  cur_repar[nLabels + l2] += shift;
	}

	//repar v3
	for (uint l3=0; l3 < nLabels; l3++) {

	  double b = var_cost(v3,l3);
	  double a = 1e300;
	    
	  for (uint l1=0; l1 < nLabels; l1++) {
	    for (uint l2=0; l2 < nLabels; l2++) {
	  
	      double hyp = cur_weight * cost(l1,l2,l3)
		+ cur_repar[l1] + cur_repar[nLabels + l2];
	      
	      if (hyp < a)
		a = hyp;
	    }
	  }

	  a += cur_repar[2*nLabels + l3];

	  double shift = 0.5*(b-a);
	  var_cost(v3,l3) -= shift;
	  cur_repar[2*nLabels + l3] += shift;
	}
      }
      else if (cur_type <= 5) {
	//2-clique
	const Math2D::Matrix<float>& cost = (cur_type == 3) ? cost3 : ((cur_type == 4) ? cost4 : cost5);

	uint v1 = factor_face[start_pos];
	uint v2 = factor_face[start_pos+1];

	//repar v1
	for (uint l1=0; l1 < nLabels; l1++) {

	  double b = var_cost(v1,l1);
	  double a = 1e300;
	    
	  for (uint l2=0; l2 < nLabels; l2++) {
	  
	    double hyp = cur_weight * cost(l1,l2)
	      + cur_repar[nLabels + l2];
	      
	    if (hyp < a)
	      a = hyp;
	  }

	  a += cur_repar[l1];

	  double shift = 0.5*(b-a);
	  var_cost(v1,l1) -= shift;
	  cur_repar[l1] += shift;
	}


	//repar v2
	for (uint l2=0; l2 < nLabels; l2++) {

	  double b = var_cost(v2,l2);
	  double a = 1e300;
	    
	  for (uint l1=0; l1 < nLabels; l1++) {
	  
	    double hyp = cur_weight * cost(l1,l2)
	      + cur_repar[l1];
	      
	    if (hyp < a)
	      a = hyp;
	  }

	  a += cur_repar[nLabels + l2];

	  double shift = 0.5*(b-a);
	  var_cost(v2,l2) -= shift;
	  cur_repar[nLabels + l2] += shift;
	}
      }
      else {
	//std::cerr << "WARNING: unknown type " << cur_type << std::endl;
      }

    }
  }

  std::cerr << "offset: " << energy_offset << std::endl;

  double disc_energy = factor_curv_energy(mesh, edge_pairs, data_term, options, face_uint_solution)
    + energy_offset;
  
  std::cerr << "solution found by stand-alone routine: " << disc_energy << std::endl;

  double icm_energy = factor_curv_icm(mesh, edge_pairs, data_term, options, face_uint_solution)
    + energy_offset;

  std::cerr << "energy of ICM solution: " << icm_energy << std::endl;

  //extract solution
  Storage1D<PixelFaceRelation> shares;
  Math1D::Vector<uint> share_start;

  Math1D::Vector<int> labels(mesh.nFaces(),-1);
  
  uint out_factor = options.output_factor_;

  segmentation.resize(xDim*out_factor,yDim*out_factor);

  Math3D::Tensor<double> output(xDim*out_factor,yDim*out_factor,nLabels,0.0);
  mesh.enlarge(out_factor,out_factor);

  //re-compute pixel shares for the now larger mesh
  compute_pixel_shares(mesh, out_factor*xDim, out_factor*yDim, shares, share_start);

  uint seg_factor = 255 / (nLabels-1);

  for (uint y=0; y < yDim*out_factor; y++) {
    for (uint x=0; x < xDim*out_factor; x++) {

      Math1D::Vector<double> votes(nLabels,0.0);

      for (uint k= share_start[y*(xDim*out_factor)+x]; k < share_start[y*(xDim*out_factor)+x+1]; k++) {
	uint face = shares[k].face_idx_;
	double cur_share = shares[k].share_;
	votes[face_uint_solution[face]] += mesh.convex_area(face) * cur_share;
      }

      double max_val = 0.0;
      uint arg_max = MAX_UINT;

      for (uint r=0; r < nLabels; r++) {

	double val = votes[r];
	if (val > max_val) {
	  max_val = val;
	  arg_max = r;
	}
      }

      segmentation(x,y) = arg_max * seg_factor;
    }
  }

  return bound;
}
