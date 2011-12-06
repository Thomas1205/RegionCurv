/*** First version written by Thomas Schoenemann as a private person without employment, October 2009 ***/
/*** continued by Thomas Schoenemann as an employee of Lund University, Sweden, January 2010 ***/
/*** extended by Thomas Schoenemann and Petter Strandmark as employees of Lund University, Sweden, September 2010-2011 ***/

#include "lp_segmentation.hh"
#include "mesh2D.hh"
#include "sparse_matrix_description.hh"
#include "timing.hh"
#include "curvature.hh"

#include "segmentation_common.hh"
#include "stl_out.hh"

#include <coin/ClpSimplex.hpp>
#include <coin/ClpPlusMinusOneMatrix.hpp>
#include <coin/CbcModel.hpp>
#include <coin/OsiClpSolverInterface.hpp>

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

//solves a segmentation problem with length regularity via an LP
//@param lambda: the weight for the length term
double lp_segment_lenreg(const Math2D::Matrix<float>& data_term, const LPSegOptions& options,
  double energy_offset, Math2D::Matrix<uint>& segmentation, const Math2D::Matrix<int>* fixed_labels) {

    double lambda = options.lambda_;
    int neighborhood = options.neighborhood_;

    uint xDim = uint( data_term.xDim() );
    uint yDim = uint( data_term.yDim() );

    Mesh2D mesh;

    uint nAreasPerPixel = 1;
    if (neighborhood == 8)
      nAreasPerPixel = 4;
    if (neighborhood == 16)
      nAreasPerPixel = 32;

    generate_mesh(xDim,yDim,neighborhood,mesh);

    Storage1D<PixelFaceRelation> shares;
    Math1D::Vector<uint> share_start;
    compute_pixel_shares(mesh, xDim, yDim, shares, share_start);

    uint nVars = mesh.nFaces() + 2*mesh.nEdges();
    uint nConstraints = mesh.nEdges();

    Math1D::NamedVector<double> cost(nVars,0.0,MAKENAME(cost));
    for (uint y=0; y < yDim; y++) {
      for (uint x=0; x < xDim; x++) {

        uint base = (y*xDim+x)*nAreasPerPixel;
        double cur_data = data_term(x,y); 

        for (uint i=0; i < nAreasPerPixel; i++) {
          double area = mesh.convex_area(base+i);
          cost[base+i] = area * cur_data;
        }
      }
    }

    uint edge_offset = mesh.nFaces();
    for (uint i=0; i < mesh.nEdges(); i++) {
      double length_weight = lambda * mesh.edge_length(i);

      if (mesh.adjacent_faces(i).size() >= 2) {
        cost[edge_offset + 2*i]   = length_weight;
        cost[edge_offset + 2*i+1] = length_weight;
      }
    }

    Math1D::NamedVector<double> rhs(nConstraints,0.0,MAKENAME(rhs));

    Math1D::NamedVector<double> var_lb(nVars,0.0,MAKENAME(var_lb));
    Math1D::NamedVector<double> var_ub(nVars,1.0,MAKENAME(var_ub));

    if (fixed_labels != 0) {

      bool has_warned = false;

      for (uint y=0; y < yDim; y++) {
        for (uint x=0; x < xDim; x++) {

          int fixed = (*fixed_labels)(x,y);
          if (fixed == 0) {
            for (uint v=share_start[y*xDim+x]; v < share_start[y*xDim+x+1]; v++) {
              uint f = shares[v].face_idx_;
              float area_share = std::min(1.0f, shares[v].share_);
              if (area_share >= 0.95)
                var_ub[f] = 0.0;
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

    std::cerr << "coding matrix" << std::endl;

    uint nEntries = 4*mesh.nEdges(); 
    SparseMatrixDescription<double> lp_descr(nEntries, nConstraints, nVars);

    for (uint j=0; j < mesh.nEdges(); j++) {
      lp_descr.add_entry(j,edge_offset+2*j,1);
      lp_descr.add_entry(j,edge_offset+2*j+1,-1);

      for (std::vector<uint>::const_iterator it = mesh.adjacent_faces(j).begin();
        it != mesh.adjacent_faces(j).end(); it++) {

          lp_descr.add_entry(j,*it,mesh.match(*it,j));
      }
    }

    std::cerr << nVars << " variables, " << nConstraints << " constraints" << std::endl;

    std::cerr << "converting matrix" << std::endl;

    std::clock_t tStartCLP, tEndCLP;  

    CoinPackedMatrix coinMatrix(false,(int*) lp_descr.row_indices(),(int*) lp_descr.col_indices(),
      lp_descr.value(),lp_descr.nEntries());

    ClpSimplex lpSolver;
    lpSolver.loadProblem (coinMatrix, var_lb.direct_access(), var_ub.direct_access(),   
      cost.direct_access(), rhs.direct_access(), rhs.direct_access());

    coinMatrix.cleanMatrix();

    tStartCLP = std::clock();

    int error = lpSolver.dual();
    //lpSolver.initialSolve();
    //int error = 1 - lpSolver.isProvenOptimal();

    tEndCLP = std::clock();

    std::cerr << "CLP-time: " << diff_seconds(tEndCLP,tStartCLP) << " seconds. " << std::endl;

    if (error)
      std::cerr << "WARNING: solving the LP-relaxation failed!!!" << std::endl;

    const double* lp_solution = lpSolver.primalColumnSolution();

#ifdef HAS_VNK_GRAPH 
    Graph<double,double,double> graph(mesh.nFaces(), 2*mesh.nEdges());

    graph.add_node(mesh.nFaces());

    for (uint f=0; f < mesh.nFaces(); f++) {

      graph.add_tweights(f, cost[f], 0.0);

      if (var_lb[f] == 1.0)
        graph.add_tweights(f, 0.0, 1e50);
      if (var_ub[f] == 0.0) 
        graph.add_tweights(f, 1e50,0.0);
    }
    for (uint e=0; e < mesh.nEdges(); e++) {

      const std::vector<uint>& adjacent = mesh.adjacent_faces(e);

      if (adjacent.size() == 1) {

        graph.add_tweights(adjacent[0],lambda*mesh.edge_length(e),0.0);
      }
      else {
        assert(adjacent.size() == 2);

        graph.add_edge(adjacent[0],adjacent[1],lambda*mesh.edge_length(e),lambda*mesh.edge_length(e));
      }
    }

    double gc_energy = graph.maxflow();

    std::cerr << "graph cut energy: " << gc_energy << std::endl;

    Math1D::Vector<double> gc_solution(nVars,0.0);

    for (uint f=0; f < mesh.nFaces(); f++) {
      int seg = graph.what_segment(f);
      gc_solution[f] = (seg == Graph<double,double,double>::SOURCE) ? 0 : 1;
    }

    lp_solution = gc_solution.direct_access();
#endif


    double energy = energy_offset;
    for (uint i=0; i < nVars; i++)
      energy += cost[i] * lp_solution[i];

    std::cerr << "energy: " << energy << std::endl;

    uint out_factor = options.output_factor_;

    segmentation.resize(xDim*out_factor,yDim*out_factor);

    mesh.enlarge(out_factor,out_factor);

    //re-compute pixel shares for the now larger mesh
    compute_pixel_shares(mesh, out_factor*xDim, out_factor*yDim, shares, share_start);

    for (uint y=0; y < yDim*out_factor; y++) {
      for (uint x=0; x < xDim*out_factor; x++) {

        double sum = 0.0;
        for (uint k= share_start[y*(xDim*out_factor)+x]; k < share_start[y*(xDim*out_factor)+x+1]; k++) {
          uint face = shares[k].face_idx_;
          sum += mesh.convex_area(face) * shares[k].share_ * lp_solution[face];
        }
        segmentation(x,y) = uint(sum*255.0);
      }
    }

    return energy;
}

/***********************************************************************************************************************/

double point_energy(uint point, const uint* solution, const Mesh2D& mesh, 
		    const std::vector<Mesh2DEdgePair>& edge_pairs, 
		    const NamedStorage1D<std::vector<uint> >& point_pair,
		    const NamedStorage1D<std::vector<uint> >& point_edge, 
		    double lambda, double gamma, bool bruckstein, bool crossings_allowed = false) {

  //NOTE: if crossings are ALLOWED, this routine only APPROXIMATES the energy

  double energy = 0.0;

  Math1D::Vector<double> drop(point_edge[point].size(),0.0);

  double abs_drop_sum = 0.0;

  std::map<uint,uint> drop_idx;

  for (uint e=0; e < drop.size(); e++) {

    uint edge = point_edge[point][e];
    drop_idx[edge] = e;

    const std::vector<uint>& adjacent_faces = mesh.adjacent_faces(edge);

    double cur_drop = 0.0;

    for (uint i=0; i < adjacent_faces.size(); i++) {
      uint face = adjacent_faces[i];

      uint cur_label = solution[face];

      if (cur_label != 0) 
        cur_drop += mesh.match(face,edge);
    }

    drop[e] = cur_drop;

    abs_drop_sum += fabs(cur_drop);
  }

  if (abs_drop_sum > 0.0) {

    //std::cerr << "pair present" << std::endl;

    int ads = int (abs_drop_sum+0.5);
    assert((ads % 2) == 0);

    if (ads == 2) {
      //simple case

      for (uint j=0; j < point_pair[point].size(); j++) {

        const Mesh2DEdgePair& cur_edge_pair = edge_pairs[point_pair[point][j]];

        uint first = cur_edge_pair.first_edge_idx_;
        uint second = cur_edge_pair.second_edge_idx_;

        if (drop[drop_idx[first]] != 0.0 && drop[drop_idx[second]] != 0.0) {

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

      //std::cerr << "nVars: " << nVars << std::endl;

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

        if (drop[drop_idx[first]] == 0.0 || drop[drop_idx[second]] == 0.0) {
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

      uint nConstraints = 2*point_edge[point].size() + conflicts.size();
      uint nMatrixEntries = 4*point_pair.size() + 4*conflicts.size();

      SparseMatrixDescription<double> lp_descr(nMatrixEntries, nConstraints, nVars);

      Math1D::Vector<double> rhs_lower(nConstraints,0.0);
      Math1D::Vector<double> rhs_upper(nConstraints,0.0);

      for (uint e=0; e < drop.size(); e++) {
        rhs_lower[e] = -drop[e];
        rhs_upper[e] = -drop[e];
        rhs_upper[drop.size()+e] = 1.0;
      }

      /*** code constraint system ***/
      for (uint j=0; j < point_pair[point].size(); j++) {

        const Mesh2DEdgePair& cur_edge_pair = edge_pairs[point_pair[point][j]];

        uint first_edge = cur_edge_pair.first_edge_idx_;
        uint second_edge = cur_edge_pair.second_edge_idx_;

        uint middle_point = cur_edge_pair.common_point_idx_;

        lp_descr.add_entry(drop.size()+drop_idx[first_edge],2*j,1);
        lp_descr.add_entry(drop.size()+drop_idx[first_edge],2*j+1,1);

        lp_descr.add_entry(drop.size()+drop_idx[second_edge],2*j,1);
        lp_descr.add_entry(drop.size()+drop_idx[second_edge],2*j+1,1);

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

  return energy;

}

//currently assuming that crossing prevention is active
double curv_energy(const uint* solution, const double* region_cost, const Mesh2D& mesh,
		   double lambda, double gamma, bool bruckstein, bool crossings_allowed = false) { 

  double energy = 0.0;
  
  for (uint f=0; f < mesh.nFaces(); f++)
    energy += region_cost[f] * double(solution[f]);

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

    energy += point_energy(p, solution, mesh, edge_pairs, point_pair, point_edge, 
			   lambda, gamma, bruckstein, crossings_allowed);
  }

  //std::cerr << "final energy: " << energy << std::endl;
  
  return energy;
}

double curv_icm(uint* solution, const double* region_cost, const Mesh2D& mesh,
		double lambda, double gamma, bool bruckstein, bool crossings_allowed = false) { 

  double energy = 0.0;

  uint edge_pair_var_offs = mesh.nFaces();

  for (uint f=0; f < mesh.nFaces(); f++)
    energy += region_cost[f] * double(solution[f]);

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

  Math1D::Vector<double> cur_point_energy(mesh.nPoints());

  for (uint p=0; p < mesh.nPoints(); p++) {

    cur_point_energy[p] = point_energy(p, solution, mesh, edge_pairs, point_pair, point_edge, 
      lambda, gamma, bruckstein, crossings_allowed);

    energy += cur_point_energy[p];
  }


  bool changes = true;

  for (uint iter=1; changes && iter <= 15; iter++) {

    changes = false;
    uint nChanges = 0;

    std::cerr << "ICM iteration " << iter << std::endl;

    for (uint f=0; f < mesh.nFaces(); f++) {

      uint cur_label = solution[f];

      double cur_energy = 0.0;
      if (cur_label == 1)
        cur_energy += region_cost[f];

      std::vector<uint> point_indices;
      mesh.get_polygon_points(f, point_indices);

      for (uint k=0; k < point_indices.size(); k++) 
        cur_energy += cur_point_energy[point_indices[k]];

      uint hyp_label = 1 - cur_label;

      //temporarily modify the solution
      solution[f] = hyp_label;

      Math1D::Vector<double> hyp_point_cost(point_indices.size());

      double hyp_energy = 0.0;
      if (hyp_label == 1)
        hyp_energy += region_cost[f];

      for (uint k=0; k < point_indices.size(); k++) {
        hyp_point_cost[k] = point_energy(point_indices[k], solution, mesh, edge_pairs, point_pair, point_edge, 
          lambda, gamma, bruckstein, crossings_allowed);
        hyp_energy += hyp_point_cost[k];
      }

      if (cur_energy <= hyp_energy) {
        solution[f] = cur_label;
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

    std::cerr << "energy " << energy << "(" << nChanges << " changes)"<< std::endl;
  }

  return energy;
}


/***********************************************************************************************************************/

//solves a segmentation problem with length and curvature regularity via an LP
double lp_segment_curvreg(const Math2D::Matrix<float>& data_term, const LPSegOptions& options, double energy_offset, 
			  Math2D::Matrix<uint>& segmentation, const Math2D::Matrix<int>* fixed_labels) {

  //Get the options
  double lambda = options.lambda_;
  double gamma = options.gamma_;
  int neighborhood = options.neighborhood_;
  bool enforce_consistent_boundaries = options.enforce_consistent_boundaries_;
  bool enforce_consistent_points = options.enforce_consistent_points_;
  bool enforce_regionedge = options.enforce_regionedge_;
  bool prevent_crossings = options.prevent_crossings_;
  if (options.convex_prior_) {
    prevent_crossings = false;
  }
  std::vector<double> convex_sol;

  bool light_constraints = options.light_constraints_;
  bool bruckstein = options.bruckstein_;
  bool reduce_edge_pairs = options.reduce_edge_pairs_;

  std::string solver = options.solver_;

  std::cerr.precision(10);

  uint light_factor = (light_constraints) ? 1 : 2;

  std::cerr << "light constraints: " << light_constraints << std::endl;

  assert(neighborhood <= 16); 

  uint xDim = uint( data_term.xDim() );
  uint yDim = uint( data_term.yDim() );

  Mesh2D mesh;  
  if (options.gridtype_ == options.Square) {
    if (options.adaptive_mesh_n_ < 0) {
      double xfac = double(xDim) / double(options.griddim_xDim_);
      double yfac = double(yDim) / double(options.griddim_yDim_);
      generate_mesh( options.griddim_xDim_, options.griddim_yDim_, neighborhood, mesh, false, fixed_labels);
      mesh.enlarge(xfac,yfac);
    }
    else {
      //Adaptive mesh
      generate_adaptive_mesh(data_term, mesh, neighborhood, options.adaptive_mesh_n_);
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
      generate_adaptive_hexagonal_mesh(data_term, mesh, neighborhood, options.adaptive_mesh_n_);
    }
  }

  std::vector<Mesh2DEdgePair> edge_pairs;
  Petter::statusTry("Generating edge pairs...");
  mesh.generate_edge_pair_list(edge_pairs);
  size_t nRemoved = 0;
  if (reduce_edge_pairs) {
    nRemoved = filter_edge_pairs(mesh, edge_pairs); 
  }

  Petter::statusOK();

  if (reduce_edge_pairs) {
    std::cerr << "removed " << nRemoved << " edge pairs." << std::endl;
  }

  std::cerr << edge_pairs.size() << " edge pairs." << std::endl;

  uint nVars = uint( mesh.nFaces() + 2*edge_pairs.size() );
  uint nConstraints = 3*mesh.nEdges();

  const uint consistency_con_offs = nConstraints;
  if (enforce_consistent_boundaries) {
    nConstraints += light_factor*mesh.nEdges();
  }
  const uint point_consistency_con_offs = nConstraints;
  if (enforce_consistent_points) {
    nConstraints += mesh.nPoints();
  }
  const uint regionedge_constraints_offs = nConstraints;
  if (enforce_regionedge) {
    nConstraints += 2*light_factor*mesh.nEdges();
  }

  Math1D::NamedVector<double> var_lb(nVars,0.0,MAKENAME(var_lb));
  Math1D::NamedVector<double> var_ub(nVars,1.0,MAKENAME(var_ub));
  Math1D::NamedVector<double> cost(nVars,0.0,MAKENAME(cost));

  Petter::statusTry("Creating data terms...");
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
  Petter::statusOK();

  uint edge_pair_var_offs = mesh.nFaces();

  Petter::statusTry("Line pair costs...");
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

    cost[edge_pair_var_offs+2*j] = weight;
    cost[edge_pair_var_offs+2*j+1] = weight;

    /*** check if (at the image border) one of the edge pairs is impossible. if so, set its upper bound to 0 ***/
    uint edge = first;
    if (mesh.adjacent_faces(edge).size() == 1) {
      int match = mesh.match(mesh.adjacent_faces(edge)[0],edge);

      uint y = edge_pair_var_offs+2*j;

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

      uint y = edge_pair_var_offs+2*j;

      if (edge_pairs[j].common_point_idx_ == mesh.edge(edge).from_idx_)
	      match *= -1;

      if (match == -1)
	      var_ub[y] = 0.0;
      else if (match == 1)
	      var_ub[y+1] = 0.0;
    }
  }
  Petter::statusOK();

  Math1D::NamedVector<double> rhs_lower(nConstraints,0.0,MAKENAME(rhs_lower));
  Math1D::NamedVector<double> rhs_upper(nConstraints,0.0,MAKENAME(rhs_upper));

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
	          if (area_share >= 0.95)
	            var_ub[f] = 0.0;
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

  uint nEntries = uint( 2*mesh.nEdges() + 6*edge_pairs.size() );  
  if (enforce_consistent_boundaries) {
    //we also allocate space for the slack variables used with the convex solver
    nEntries += uint( light_factor*2*edge_pairs.size() + light_factor*mesh.nEdges()); 
  }
  if (enforce_consistent_points) {
    nEntries += uint( 2*edge_pairs.size() );
  }
  if (enforce_regionedge) {
    nEntries += uint( 8*light_factor*edge_pairs.size() //Note: not exact number
		      + 2*light_factor*mesh.nEdges() ); //we also allocate space for the slack variables used with the convex solver
  }

  SparseMatrixDescription<double> lp_descr(nEntries, nConstraints, nVars);

  /**** a) code surface continuation constraints *****/

  Petter::statusTry("Coding surface cont. constr....");
  for (uint j=0; j < mesh.nEdges(); j++) {

    const std::vector<uint>& adjacent_faces = mesh.adjacent_faces(j);
    for (std::vector<uint>::const_iterator it = adjacent_faces.begin();
    it != adjacent_faces.end(); it++) {

      lp_descr.add_entry(j,*it,mesh.match(*it,j));
    }
  }

  for (uint j=0; j < edge_pairs.size(); j++) {

    uint first_edge = edge_pairs[j].first_edge_idx_;
    uint second_edge = edge_pairs[j].second_edge_idx_;

    uint middle_point = edge_pairs[j].common_point_idx_;

    if (mesh.edge(first_edge).to_idx_ == middle_point) {
      lp_descr.add_entry(first_edge,edge_pair_var_offs+2*j,1);
    }
    else {
      lp_descr.add_entry(first_edge,edge_pair_var_offs+2*j,-1);
    }

    if (mesh.edge(second_edge).to_idx_ == middle_point) {
      lp_descr.add_entry(second_edge,edge_pair_var_offs+2*j+1,1);
    }
    else {
      lp_descr.add_entry(second_edge,edge_pair_var_offs+2*j+1,-1);
    }
  }
  Petter::statusOK();

  /**** b) code boundary continuation constraints *****/
  Petter::statusTry("Coding boundary cont. constr....");
  uint boundary_con_offset = mesh.nEdges();

  for (uint j=0; j < edge_pairs.size(); j++) {

    uint first_edge = edge_pairs[j].first_edge_idx_;
    uint second_edge = edge_pairs[j].second_edge_idx_;

    uint middle_point = edge_pairs[j].common_point_idx_;

    if (mesh.edge(first_edge).to_idx_ == middle_point) {      
      lp_descr.add_entry(boundary_con_offset + 2*first_edge, edge_pair_var_offs+2*j, 1);
      lp_descr.add_entry(boundary_con_offset + 2*first_edge+1, edge_pair_var_offs+2*j+1, -1);
    }
    else {
      lp_descr.add_entry(boundary_con_offset + 2*first_edge+1, edge_pair_var_offs+2*j, 1);
      lp_descr.add_entry(boundary_con_offset + 2*first_edge, edge_pair_var_offs+2*j+1, -1);
    }

    if (mesh.edge(second_edge).from_idx_ == middle_point) {
      lp_descr.add_entry(boundary_con_offset + 2*second_edge, edge_pair_var_offs+2*j, -1);
      lp_descr.add_entry(boundary_con_offset + 2*second_edge+1, edge_pair_var_offs+2*j+1, 1);
    }
    else {
      lp_descr.add_entry(boundary_con_offset + 2*second_edge+1, edge_pair_var_offs+2*j, -1);
      lp_descr.add_entry(boundary_con_offset + 2*second_edge, edge_pair_var_offs+2*j+1, 1);
    }
  }
  Petter::statusOK();
  uint nStandardEntries = lp_descr.nEntries();

  if (enforce_consistent_boundaries) {

    Petter::statusTry("Coding consistent boundaries...");

    //constraints in words: for each oriented edge, the pairs that start with this oriented edge 
    // and the pairs that end in the oppositely oriented edge may not sum to more than 1.0
    // (i.e. they are mutually exclusive)

    for (uint c=consistency_con_offs; c < consistency_con_offs + light_factor*mesh.nEdges(); c += light_factor) {

      rhs_upper[c] = 1.0;
      if (!light_constraints)
        rhs_upper[c+1] = 1.0;
    }

    for (uint j=0; j < edge_pairs.size(); j++) {

      uint middle_point = edge_pairs[j].common_point_idx_;

      uint first_edge = edge_pairs[j].first_edge_idx_;
      uint second_edge = edge_pairs[j].second_edge_idx_;


      if (mesh.edge(first_edge).to_idx_ == middle_point) { 
        lp_descr.add_entry(consistency_con_offs + light_factor*first_edge, edge_pair_var_offs+2*j, 1);
        lp_descr.add_entry(consistency_con_offs + light_factor*first_edge, edge_pair_var_offs+2*j+1, 1);
      }
      else {
        if (!light_constraints) {
          lp_descr.add_entry(consistency_con_offs + light_factor*first_edge+1, edge_pair_var_offs+2*j, 1);
          lp_descr.add_entry(consistency_con_offs + light_factor*first_edge+1, edge_pair_var_offs+2*j+1, 1);
        }
      }

      if (mesh.edge(second_edge).to_idx_ == middle_point) {
        lp_descr.add_entry(consistency_con_offs + light_factor*second_edge, edge_pair_var_offs+2*j, 1);
        lp_descr.add_entry(consistency_con_offs + light_factor*second_edge, edge_pair_var_offs+2*j+1, 1);
      }
      else {
        if (!light_constraints) {
          lp_descr.add_entry(consistency_con_offs + light_factor*second_edge+1, edge_pair_var_offs+2*j, 1);
          lp_descr.add_entry(consistency_con_offs + light_factor*second_edge+1, edge_pair_var_offs+2*j+1, 1);
        }
      }
    }

    Petter::statusOK();
  }


  if (enforce_consistent_points) {
    Petter::statusTry("Coding consistent points...");

    for (uint c=point_consistency_con_offs; c < nConstraints; c++) {
      rhs_upper[c] = 1.0;
    }

    for (uint j=0; j < edge_pairs.size(); j++) {
      uint middle_point = edge_pairs[j].common_point_idx_;

      lp_descr.add_entry(point_consistency_con_offs + middle_point, edge_pair_var_offs+2*j, 1);
      lp_descr.add_entry(point_consistency_con_offs + middle_point, edge_pair_var_offs+2*j+1, 1);
    }
    Petter::statusOK();
  }

  if (enforce_regionedge) {
    Petter::statusTry("Coding boundary/edge pair constr....");
    uint rowoff = regionedge_constraints_offs;

    for (uint edge=0; edge < mesh.nEdges(); edge++) {

      //Get the two adjacent faces
      if (mesh.adjacent_faces(edge).size() != 2)
      {
        //One of the edges is at the border of the image
        continue;
      }

      uint x1 = mesh.adjacent_faces(edge)[0];
      uint x2 = mesh.adjacent_faces(edge)[1];

      lp_descr.add_entry(rowoff+2*light_factor*edge  , x1, 1);
      lp_descr.add_entry(rowoff+2*light_factor*edge  , x2, 1);
      lp_descr.add_entry(rowoff+2*light_factor*edge+1, x1, -1);
      lp_descr.add_entry(rowoff+2*light_factor*edge+1, x2, -1);
      if (!light_constraints) {
        lp_descr.add_entry(rowoff+2*light_factor*edge+2, x1, 1);
        lp_descr.add_entry(rowoff+2*light_factor*edge+2, x2, 1);
        lp_descr.add_entry(rowoff+2*light_factor*edge+3, x1, -1);
        lp_descr.add_entry(rowoff+2*light_factor*edge+3, x2, -1);
      }

      //NOTE (important for the construction with slacks: the binding constraints are in both cases the upper bounds)
      rhs_lower[rowoff+2*light_factor*edge]   = 0;
      rhs_upper[rowoff+2*light_factor*edge]   = 2;
      rhs_lower[rowoff+2*light_factor*edge+1] = -2;
      rhs_upper[rowoff+2*light_factor*edge+1] = 0;

      if (!light_constraints) {
        rhs_lower[rowoff+2*light_factor*edge+2] = 0;
        rhs_upper[rowoff+2*light_factor*edge+2] = 2;
        rhs_lower[rowoff+2*light_factor*edge+3] = -2;
        rhs_upper[rowoff+2*light_factor*edge+3] = 0;
      }
    }

    for (uint j=0; j < edge_pairs.size(); j++) {
      uint first = edge_pairs[j].first_edge_idx_;
      uint second = edge_pairs[j].second_edge_idx_;	

      uint y = mesh.nFaces() + 2*j;

      uint edge = first;
      if (mesh.adjacent_faces(edge).size() == 2) {
        lp_descr.add_entry(rowoff+2*light_factor*edge   ,y  , 1);
        lp_descr.add_entry(rowoff+2*light_factor*edge+1 ,y  , 1);
        if (!light_constraints) {
          lp_descr.add_entry(rowoff+2*light_factor*edge+2 ,y+1, 1);
          lp_descr.add_entry(rowoff+2*light_factor*edge+3 ,y+1, 1);
        }
      }

      edge = second;
      if (mesh.adjacent_faces(edge).size() == 2) {
        lp_descr.add_entry(rowoff+2*light_factor*edge   ,y+1, 1);
        lp_descr.add_entry(rowoff+2*light_factor*edge+1 ,y+1, 1);
        if (!light_constraints) {
          lp_descr.add_entry(rowoff+2*light_factor*edge+2 ,y  , 1);
          lp_descr.add_entry(rowoff+2*light_factor*edge+3 ,y  , 1);
        }
      }
    }

    Petter::statusOK();
  }



  if (options.convex_prior_) {
    Petter::statusTry("Adding convexity constraints...");

    for (uint j=0; j < edge_pairs.size(); j++) {
      uint first = edge_pairs[j].first_edge_idx_;
      uint second = edge_pairs[j].second_edge_idx_;

      uint middle_point = edge_pairs[j].common_point_idx_;
      uint first_edge = edge_pairs[j].first_edge_idx_;
      uint second_edge = edge_pairs[j].second_edge_idx_;

      //Vector of the first edge 
      double from1_x = mesh.point(mesh.edge(first).from_idx_).x_;
      double from1_y = mesh.point(mesh.edge(first).from_idx_).y_;
      double to1_x = mesh.point(mesh.edge(first).to_idx_).x_;
      double to1_y = mesh.point(mesh.edge(first).to_idx_).y_;
      double dx1 = to1_x - from1_x;
      double dy1 = to1_y - from1_y;
      

      //Vector of the second edge 
      double from2_x = mesh.point(mesh.edge(second).from_idx_).x_;
      double from2_y = mesh.point(mesh.edge(second).from_idx_).y_;
      double to2_x = mesh.point(mesh.edge(second).to_idx_).x_;
      double to2_y = mesh.point(mesh.edge(second).to_idx_).y_;
      double dx2 = to2_x - from2_x;
      double dy2 = to2_y - from2_y;

      bool ok;
      bool bothok=false;

      // First edge determines direction
      if (mesh.edge(first).from_idx_ == edge_pairs[j].common_point_idx_) {
        //CASE A: First edge away from the center

        if (mesh.edge(second).from_idx_ == edge_pairs[j].common_point_idx_) {
          dx2 = -dx2;
          dy2 = -dy2;
        }

        //    <-----o<-----

        double det = dx2*dy1 - dx1*dy2;
        if (det >= 1e-6) {
          ok = false;
        }
        else if (det <= -1e-6) {
          ok = true;
        }
        else {
          bothok = true;
        }

      }
      else {
        //CASE B: First edge towards the center
        
        if (mesh.edge(second).from_idx_ == edge_pairs[j].common_point_idx_) {
          dx2 = -dx2;
          dy2 = -dy2;
        }

        //    ----->o----->

        double det = dx1*dy2 - dx2*dy1;
        if (det >= 1e-6) {
          ok = false;
        }
        else if (det <= -1e-6) {
          ok = true;
        }
        else {
          bothok = true;
        }

      }


      if (!bothok) {
        if (ok) {
          // The other pair is non-convex and not allowed
          var_lb[edge_pair_var_offs+2*j+1] = 0;
          var_ub[edge_pair_var_offs+2*j+1] = 0;
        }
        else {
          // The first pair is non-convex and not allowed
          var_lb[edge_pair_var_offs+2*j] = 0;
          var_ub[edge_pair_var_offs+2*j] = 0;
        }
      }


    }


    /*for (int i=0; i<nVars; ++i) { 
      if (cost[i] < 0) {
        var_lb[i] = 1;
        var_ub[i] = 1;
      }
      else {
        var_lb[i] = 0;
        var_ub[i] = 0;
      }
    }*/

    Petter::statusOK();
  }




  Math1D::Vector<uint> row_start(nConstraints+1);
  lp_descr.sort_by_row(row_start);

  bool solver_known = false;
  const double* lp_solution = 0;

  int error = 0;

#ifdef HAS_GUROBI
  Math1D::Vector<double> GUROBI_solution;

  GRBenv   *grb_env   = NULL;
  GRBmodel *grb_model = NULL;

  if (solver == "gurobi") {

    solver_known = true;

    /* Create environment */

    error = GRBloadenv(&grb_env,NULL);
    GRBsetintparam(grb_env, GRB_INT_PAR_METHOD, GRB_METHOD_BARRIER);
    GRBsetdblparam(grb_env, "BarConvTol", 1e-10);
    GRBsetintparam(grb_env, "Crossover", 0);
    //GRBsetintparam(grb_env, "CrossoverBasis", 1);
    GRBsetintparam(grb_env, "Presolve", 1);
    GRBsetintparam(grb_env, "PrePasses", 2);

    assert (!error && grb_env != NULL);

    /* Create an empty model */
    error = GRBnewmodel(grb_env, &grb_model, "curv-lp", 0, NULL, NULL, NULL, NULL, NULL);
    assert(!error);

    Storage1D<char> vtype(nVars,GRB_CONTINUOUS);

    error = GRBaddvars(grb_model,nVars,0,NULL,NULL,NULL,cost.direct_access(),var_lb.direct_access(),
		       var_ub.direct_access(),vtype.direct_access(),NULL);
    assert(!error);

    error = GRBupdatemodel(grb_model);
    assert(!error);

    for (uint c=0; c < nConstraints; c++) {

      if (rhs_lower[c] == rhs_upper[c])
	error = GRBaddconstr(grb_model, row_start[c+1]-row_start[c], ((int*) lp_descr.col_indices()) + row_start[c], 
			     lp_descr.value() + row_start[c], GRB_EQUAL, rhs_upper[c], NULL);
      else
	GRBaddrangeconstr(grb_model, row_start[c+1]-row_start[c], ((int*) lp_descr.col_indices()) + row_start[c], 
			  lp_descr.value() + row_start[c], rhs_lower[c], rhs_upper[c], NULL);      

      assert(!error);
    }

    //deallocate memory where possible
    vtype.resize(0);
    row_start.resize(0);
    lp_descr.reset(0);

    /* Optimize model */
    error = GRBoptimize(grb_model);
    assert(!error);

    //extract solution
    GUROBI_solution.resize(nVars);

    for (uint v=0; v < nVars; v++)
      GRBgetdblattrelement(grb_model,"X",v, GUROBI_solution.direct_access()+v);

    lp_solution = GUROBI_solution.direct_access();
  }    
#endif

#ifdef HAS_CPLEX
  Math1D::Vector<double> CPLEX_solution;

  CPXENVptr     cp_env = NULL;
  CPXLPptr      cp_lp = NULL;

  if (solver == "cplex") {

    solver_known = true;

    int status = 0;

    /* Initialize the CPLEX environment */

    cp_env = CPXopenCPLEX (&status);
    CPXsetintparam(cp_env, CPX_PARAM_BARCROSSALG, -1);
    //CPXsetintparam(cp_env, CPX_PARAM_PREIND, CPX_OFF);
    //CPXsetintparam(cp_env, CPX_PARAM_PREPASS, 0);

    /* If an error occurs, the status value indicates the reason for
       failure.  A call to CPXgeterrorstring will produce the text of
       the error message.  Note that CPXopenCPLEX produces no output,
       so the only way to see the cause of the error is to use
       CPXgeterrorstring.  For other CPLEX routines, the errors will
       be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON.  */

    if ( cp_env == NULL ) {
      char  errmsg[1024];
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (cp_env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      exit(1);
    }

    /* Turn on output to the screen */

    status = CPXsetintparam (cp_env, CPX_PARAM_SCRIND, CPX_ON);
    if ( status ) {
      fprintf (stderr,"Failure to turn on screen indicator, error %d.\n", status);
      exit(1);
    }

    //set problem data 
    cp_lp = CPXcreateprob (cp_env, &status, "curv-lp");

    /* A returned pointer of NULL may mean that not enough memory
       was available or there was some other problem.  In the case of
       failure, an error message will have been written to the error
       channel from inside CPLEX.  In this example, the setting of
       the parameter CPX_PARAM_SCRIND causes the error message to
       appear on stdout.  */

    if ( cp_lp == NULL ) {
      fprintf (stderr, "Failed to create LP.\n");
      exit(1);
    }

    /* Now copy the problem data into the lp */

    char* row_sense = new char[nConstraints];
    for (uint c=0; c < nConstraints; c++) {

      if (rhs_lower[c] == rhs_upper[c]) {
	row_sense[c] = 'E';
      }
      else {
	row_sense[c] = 'R';
      }
    }

    int* row_count = new int[nConstraints];
    for (uint c=0; c < nConstraints; c++)
      row_count[c] = row_start[c+1] - row_start[c];

    status = CPXnewcols (cp_env, cp_lp, nVars, cost.direct_access(), var_lb.direct_access(), 
			 var_ub.direct_access(), NULL, NULL);

    std::cerr << "copying" << std::endl;

    if ( status )  
      exit(1);

    std::cerr << "adding rows" << std::endl;

    CPXaddrows(cp_env, cp_lp, 0, nConstraints, lp_descr.nEntries(), rhs_lower.direct_access(), row_sense, 
	       (int*) row_start.direct_access(), (int*) lp_descr.col_indices(), lp_descr.value(),
	       NULL, NULL);

    uint count = 0;
    Math1D::Vector<int> idx(nConstraints);
    Math1D::Vector<double> range(nConstraints);
    for (int c=0; c < (int) nConstraints; c++) {

      if (row_sense[c] == 'R') {

	idx[count] = c;
	range[count] = rhs_upper[c] - rhs_lower[c];
	count++;
      }
    }
    CPXchgrngval (cp_env, cp_lp, count, idx.direct_access(), range.direct_access());

    delete[] row_sense;
    delete[] row_count;

    std::cerr << "calling optimize" << std::endl;

    //status = CPXlpopt (cp_env, cp_lp);
    status = CPXbaropt (cp_env, cp_lp);

    if ( status ) {
      fprintf (stderr, "Failed to optimize MIP.\n");
      exit(1);
    }

    CPLEX_solution.resize_dirty(nVars);

    CPXsolution (cp_env, cp_lp, NULL, NULL, CPLEX_solution.direct_access(), NULL, NULL, NULL);

    lp_solution = CPLEX_solution.direct_access();
  }
#endif

#ifdef HAS_XPRESS

  Math1D::Vector<double> XPRESS_solution;
  
  XPRSprob xp_prob = NULL;

  if (solver == "xpress") {

    solver_known = true;

    int nReturn;
    char banner[256];

    nReturn=XPRSinit("/opt/xpressmp/");

    if (nReturn != 0) {

      char msg[512];
      XPRSgetlicerrmsg(msg,512);

      std::cerr << "error message: " << msg << std::endl;
    }

    assert(nReturn == 0);

    XPRSgetbanner(banner); printf("banner: %s \n",banner);

    nReturn=XPRScreateprob(&xp_prob);
    XPRSsetintcontrol(xp_prob,XPRS_CROSSOVER,0);

    row_start.resize(0);

    Math1D::Vector<uint> col_start(nVars+1);
    lp_descr.sort_by_column(col_start);

    char* row_sense = new char[nConstraints];
    double* row_range = new double[nConstraints];

    for (uint c=0; c < nConstraints; c++) {

      if (rhs_lower[c] == rhs_upper[c]) {
	row_sense[c] = 'E';
      }
      else {
	row_sense[c] = 'R';
	row_range[c] = rhs_upper[c] - rhs_lower[c];
      }
    }

    nReturn = XPRSloadlp(xp_prob, "curv-ilp", nVars, nConstraints, row_sense,
			 rhs_upper.direct_access(), row_range, cost.direct_access(), 
			 (int*) col_start.direct_access(), NULL, (int*) lp_descr.row_indices(), lp_descr.value(),
			 var_lb.direct_access(), var_ub.direct_access());

    delete[] row_sense;
    delete[] row_range;

    XPRSlpoptimize(xp_prob,"b");

    XPRESS_solution.resize_dirty(nVars);

    XPRSgetlpsol(xp_prob, XPRESS_solution.direct_access(), 0, 0, 0); 

    lp_solution = XPRESS_solution.direct_access();
  }
#endif

  ClpSimplex lpSolver;
  double solverTime = -1;

  if (!solver_known)
    solver = "clp";

  if (solver == "clp") {
    std::clock_t tStartCLP,tEndCLP;

#ifndef USE_PM_ONE_MATRIX
    CoinPackedMatrix coinMatrix(false,(int*) lp_descr.row_indices(),(int*) lp_descr.col_indices(),
				lp_descr.value(),lp_descr.nEntries());

    lpSolver.loadProblem (coinMatrix, var_lb.direct_access(), var_ub.direct_access(),   
			  cost.direct_access(), rhs_lower.direct_access(), rhs_upper.direct_access());

    coinMatrix.cleanMatrix();
#else

    uint nPos = 0;
    uint nNeg = 0;

    for (uint k=0; k < lp_descr.nEntries(); k++) {

      double entry = lp_descr.value()[k];

      if (entry == -1.0)
        nNeg++;
      else if (entry == 1.0)
        nPos++;
      else {
        INTERNAL_ERROR << " cannot create a +1/-1 matrix from this matrix. exiting..." << std::endl;
        exit(1);
      }
    }

    Math1D::NamedVector<int> pmone_idx(nPos+nNeg,MAKENAME(pmone_idx));

    Math1D::NamedVector<int> row_pos_start(row_start.size(),MAKENAME(row_pos_start));
    Math1D::NamedVector<int> row_neg_start(row_start.size(),MAKENAME(row_neg_start));

    uint cur_pos = 0;

    for (uint r=0; r < row_start.size()-1; r++) {

      row_pos_start[r] = cur_pos;

      for (uint k=row_start[r]; k < row_start[r+1]; k++) {

        if (lp_descr.value()[k] == 1.0) {
          pmone_idx[cur_pos] = lp_descr.col_indices()[k];
          cur_pos++;
        }
      }

      row_neg_start[r] = cur_pos;
      for (uint k=row_start[r]; k < row_start[r+1]; k++) {

        if (lp_descr.value()[k] == -1.0) {
          pmone_idx[cur_pos] = lp_descr.col_indices()[k];
          cur_pos++;
        }
      }    
    }

    row_pos_start[row_start.size()-1] = cur_pos;
    row_neg_start[row_start.size()-1] = cur_pos;

    ClpPlusMinusOneMatrix pmone_matrix(row_start.size()-1,nVars,false,pmone_idx.direct_access(),
                                       row_pos_start.direct_access(),row_neg_start.direct_access());

    pmone_idx.resize(0);
    row_pos_start.resize(0);
    row_neg_start.resize(0);

    std::cerr << "matrix created" << std::endl;

    lpSolver.loadProblem (pmone_matrix, var_lb.direct_access(), var_ub.direct_access(),   
			  cost.direct_access(), rhs_lower.direct_access(), rhs_upper.direct_access());
#endif

    //outcomment this when you are debugging
    lp_descr.reset(0);


    //lpSolver.writeMps("curv.mps");

    tStartCLP = std::clock();
    

    if (options.convex_prior_) {

#ifdef ALLOW_CONVEX_PRIORS
      //
      // We need an integer solver for this task.
      // The LP solution will essentially be worthless. Unfortunately, this
      // means that the convex constraints are of limited usefulness in practice.
      //

      OsiClpSolverInterface solver1;
      CoinPackedMatrix coinMatrix(false,(int*) lp_descr.row_indices(),(int*) lp_descr.col_indices(),
				  lp_descr.value(),lp_descr.nEntries());
      solver1.loadProblem(coinMatrix, var_lb.direct_access(), var_ub.direct_access(),   
			    cost.direct_access(), rhs_lower.direct_access(), rhs_upper.direct_access());
      for (int i=0;i<nVars;++i) {
        solver1.setInteger(i);
      }
      solver1.setObjSense(1.0);

      // Save the MPS file for other solvers to attack
      //solver1.writeMps("convex","mps",1.0);

      CbcModel model(solver1);
      model.branchAndBound();
      bool is_optimal = model.isProvenOptimal();
      
      if (is_optimal) {
        std::cerr << "Optimal solution found!" << std::endl;
      }

      convex_sol.resize(nVars);
      for (int i=0; i < nVars; ++i) {
        convex_sol[i] = model.getColSolution()[i];
      }
      lp_solution = &convex_sol[0];
#else
      std::cerr << "Define ALLOW_CONVEX_PRIORS in order to compile with an integer solver." << std::endl;
      std::exit(1);
#endif
    }
    else {
      lpSolver.dual();
      lp_solution = lpSolver.primalColumnSolution();

      error = 1 - lpSolver.isProvenOptimal();

      if (error != 0)
        std::cerr << "!!!!!!!!!!!!!!LP-solver failed!!!!!!!!!!!!!!!!!!!" << std::endl;
    }

    //ClpSolve solve_options;
    //solve_options.setSolveType(ClpSolve::useDual);
    //solve_options.setSolveType(ClpSolve::useBarrier);
    //solve_options.setPresolveType(ClpSolve::presolveNumber,5);
    //lpSolver.initialSolve(solve_options);

    tEndCLP = std::clock();

    if (error != 0)
      std::cerr << "!!!!!!!!!!!!!!LP-solver failed!!!!!!!!!!!!!!!!!!!" << std::endl;

    std::cerr << "CLP-time: " << diff_seconds(tEndCLP,tStartCLP) << " seconds. " << std::endl;
    solverTime = diff_seconds(tEndCLP,tStartCLP);
    if (mesh.nFaces() <= 20000) {
      Petter::statusTry("Saving SVG...");
      mesh.draw_labels_with_pairs(options.base_filename_ + ".lp.svg",lp_solution,edge_pairs,xDim,yDim);
      Petter::statusOK();
    }
  }
  
  //list of edge pairs that have the respective point as the common point
  Storage1D<std::vector<uint> > point_pairs(mesh.nPoints());
  
  std::cerr << "checking for crossing line pairs" << std::endl;
  
  for (uint j=0; j < edge_pairs.size(); j++) {
    
    uint common_point = edge_pairs[j].common_point_idx_;
    point_pairs[common_point].push_back(j);
  }
  
  {
    uint nNonInt = 0;
    for (uint i=0; i < mesh.nFaces(); i++) {
      
      double val = lp_solution[i];
      if (val > 0.01 && val < 0.99) {
        nNonInt++;
      }
    }

    std::cerr << "before adding violated constraints: " << nNonInt << " variables are fractional" << std::endl;
  }

  uint nIter = 0;
  while (prevent_crossings) {

    nIter++;

    std::cerr << "##### constraint generation iter #" << nIter << std::endl;

    uint nConstraintsAdded = 0;
    for (uint p=0; p < mesh.nPoints(); p++) {

      std::set<uint> active_pair;

      double sum = 0.0;
      for (std::vector<uint>::iterator it = point_pairs[p].begin(); it != point_pairs[p].end(); it++) {
	double contrib = lp_solution[edge_pair_var_offs +  2*(*it)] + lp_solution[edge_pair_var_offs +  2*(*it) + 1];
	if (fabs(contrib) >= 0.0025) {
	  sum += contrib;
	  active_pair.insert(*it);
	}
      }
      
      if (sum >= 1.0025 && active_pair.size() > 1) {
	std::cerr << "point #" << p << " might be a crossing point. " 
		  << active_pair.size() << " active pairs. " << std::endl;

	// this handles all pairwise and triple-constellations exactly,
	// as well as SOME of the higher ones. Note that we don't need to add
	// low order interactions if they are subsumed in a high order one
	
	uint nPrevAdded = nConstraintsAdded;
	
	for (uint k1=0; k1 < point_pairs[p].size()-1; k1++) {
	  
	  uint pair1 = point_pairs[p][k1];
	  
	  if (active_pair.find(pair1) == active_pair.end())
	    continue;

	  for (uint k2=k1+1; k2 < point_pairs[p].size(); k2++) {
	    
	    uint pair2 = point_pairs[p][k2];
	    
	    if (active_pair.find(pair2) == active_pair.end())
	      continue;

	    if (!line_pairs_with_meeting_point_cross(mesh, edge_pairs[pair1], edge_pairs[pair2]))
	      continue;
	    
	    double sum = lp_solution[edge_pair_var_offs +  2*pair1] + lp_solution[edge_pair_var_offs +  2*pair1 + 1]
	      + lp_solution[edge_pair_var_offs +  2*pair2] + lp_solution[edge_pair_var_offs +  2*pair2 + 1];
	    
	    //std::cerr << "checking pair" << std::endl;
	    
	    std::vector<uint> base(2);
	    base[0] = pair1;
	    base[1] = pair2;
	    
	    std::vector<uint> addons;
	    for (uint k3=k2+1; k3 < point_pairs[p].size(); k3++) {
	      
	      uint pair3 = point_pairs[p][k3];
	      
	      if (active_pair.find(pair3) != active_pair.end()
		  && line_pairs_with_meeting_point_cross(mesh, edge_pairs[pair1], edge_pairs[pair3])
		  && line_pairs_with_meeting_point_cross(mesh, edge_pairs[pair2], edge_pairs[pair3]))
		addons.push_back(point_pairs[p][k3]);
	    }

	    std::vector<std::vector<uint> > constraint_list;
	    if (addons.empty()) {
	      if (sum >= 1.005)
		constraint_list.push_back(base);
	    }
	    else {
#if 1
	      //handle higher interactions:
	      //NOTE: currently some of the generated constraints can be redundant to others
	      for (uint k=0; k < addons.size(); k++) {

		std::vector<uint> new_list = base;
		new_list.push_back(addons[k]);

		double cur_sum = sum + lp_solution[edge_pair_var_offs +  2*addons[k]] 
		  + lp_solution[edge_pair_var_offs +  2*addons[k] + 1];

		for (uint l=k+1; l < addons.size(); l++) {

		  uint p1 = addons[l];

		  bool compatible = true;
		  for (uint j=2; j < new_list.size(); j++) {

		    uint p2 = new_list[j];
		    
		    if (!line_pairs_with_meeting_point_cross(mesh, edge_pairs[p1], edge_pairs[p2])) {
		      compatible = false;
		      break;
		    }
		  }

		  if (compatible) {
		    new_list.push_back(addons[l]);
		    cur_sum += lp_solution[edge_pair_var_offs +  2*addons[l]] 
		      + lp_solution[edge_pair_var_offs +  2*addons[l] + 1];
		  }
		}

		if (cur_sum >= 1.005) {
		  std::cerr << "sum: " << cur_sum << std::endl;
		  constraint_list.push_back(new_list);
		}
	      }
#else
	      //fallback to pairwise:
	      if (sum >= 1.005)
		constraint_list.push_back(base);
#endif
	    }
	    
	    for (uint i=0; i < constraint_list.size(); i++) {
	      
	      std::vector<uint> pairs = constraint_list[i];

	      //add the constraint
	      int* cols = new int[2*pairs.size()];
	      double* coeffs = new double[2*pairs.size()];
	      for (uint k=0; k < 2*pairs.size(); k++)
		coeffs[k] = 1.0;
	      
	      for (uint p=0; p < pairs.size(); p++) {
		cols[2*p] = edge_pair_var_offs +  2*pairs[p];
		cols[2*p+1] = edge_pair_var_offs +  2*pairs[p]+1;
	      }
	      
	      //note: adding constraints separately is VERY inefficient
	      if (solver == "clp")
		lpSolver.addRow(2*pairs.size(), cols, coeffs, 0.0, 1.0);
#ifdef HAS_GUROBI
	      if (solver == "gurobi") {
		GRBaddconstr(grb_model,2*pairs.size(),cols,coeffs,'L',1.0,NULL);
	      }
#endif
#ifdef HAS_XPRESS
	      if (solver == "xpress") {
		
		double new_rhs[1] = {1.0};
		double new_range[1] = {0.0};
		int new_start[2] = {0,2*pairs.size()};
		XPRSaddrows(xp_prob, 1, 2*pairs.size(), "L", new_rhs,new_range,new_start,cols,coeffs);
	      }
#endif
#ifdef HAS_CPLEX
	      if (solver == "cplex") {
		
		double new_rhs[1] = {1.0};
		int new_start[2] = {0,2*pairs.size()};
		
		CPXaddrows(cp_env, cp_lp, 0, 1, 2*pairs.size(), new_rhs, "L", new_start, cols, coeffs,  NULL, NULL);
	      }
#endif		
	      nConstraintsAdded++;

	      delete cols;
	      delete coeffs;
	    }
	  }
	}
	
	if (nPrevAdded != nConstraintsAdded)
	  std::cerr << "added " << (nConstraintsAdded - nPrevAdded) << " constraints" << std::endl;
      }
    }

    std::cerr << "added " << nConstraintsAdded << " constraints" << std::endl;

    if (nConstraintsAdded > 0) {

      if (solver == "clp") {
        lpSolver.dual();
        lp_solution = lpSolver.primalColumnSolution();
      }
#ifdef HAS_GUROBI
      if (solver == "gurobi") {

        int error = GRBupdatemodel(grb_model);
        error = GRBoptimize(grb_model);
        assert(!error);

        for (uint v=0; v < nVars; v++)
          GRBgetdblattrelement(grb_model,"X",v, GUROBI_solution.direct_access()+v);

        lp_solution = GUROBI_solution.direct_access();
      }
#endif
#ifdef HAS_XPRESS
      if (solver == "xpress") {

        XPRSlpoptimize(xp_prob,"b");      
        XPRSgetlpsol(xp_prob, XPRESS_solution.direct_access(), 0, 0, 0); 

        lp_solution = XPRESS_solution.direct_access();
      }
#endif
#ifdef HAS_CPLEX
      if (solver == "cplex") {

        int cpx_status = CPXbaropt(cp_env,cp_lp);
        //int cpx_status = CPXlpopt(cp_env,cp_lp);

        if ( cpx_status ) {
          fprintf (stderr, "Failed to optimize MIP.\n");
          exit(1);
        }

        CPXsolution (cp_env, cp_lp, NULL, NULL, CPLEX_solution.direct_access(), NULL, NULL, NULL);

        lp_solution = CPLEX_solution.direct_access();
      }
#endif
    }
    else
      break;
  }      
  

  double lp_energy = 0;
  for (uint i=0; i < nVars; i++) {
    lp_energy += cost[i] * lp_solution[i];    
  }
  double energy = energy_offset + lp_energy;
  
  std::cerr << "lp energy: " << lp_energy << std::endl; 
  std::cerr << "original relaxation energy: " << energy << std::endl;

  Math1D::Vector<uint> labeling(mesh.nFaces(),0);

  uint nNonInt = 0;
  for (uint i=0; i < mesh.nFaces(); i++) {

    double val = lp_solution[i];
    if (val > 0.01 && val < 0.99) {
      nNonInt++;
    }

    if (val < 0.5)
      var_ub[i] = 0.0;
    else {
      var_lb[i] = 1.0;
      labeling[i] = 1;
    }
  }

  uint nNonIntAuxVars = 0;
  for (uint i= mesh.nFaces(); i < nVars; i++) {
    double val = lp_solution[i];
    if (val > 0.01 && val < 0.99) {
      nNonIntAuxVars++;
    }
  }

  std::cerr << nNonInt << " region variables are fractional" << std::endl;
  std::cerr << nNonIntAuxVars << " auxiliary (non-region) variables are fractional" << std::endl;
  
  double integral_energy = curv_energy(labeling.direct_access(), cost.direct_access(), mesh,
				       lambda, gamma, bruckstein, !prevent_crossings)
    + energy_offset;

  std::cerr << "integral energy according to independent routine: " << integral_energy << std::endl;

  double icm_energy = curv_icm(labeling.direct_access(), cost.direct_access(), mesh,
			       lambda, gamma, bruckstein, !prevent_crossings) + energy_offset;

  std::cerr << "energy after ICM: " << icm_energy << std::endl;


  uint nOpposingLinePairs = 0;

  for (uint v=mesh.nFaces(); v < nVars; v+=2) {

    double pair_sum = lp_solution[v] + lp_solution[v+1]; 

    //if (pair_sum > 1.02) {
    if (lp_solution[v] > 0.05 && lp_solution[v+1] > 0.05) {
      nOpposingLinePairs++;
    }
  }

  std::cerr << nOpposingLinePairs << " opposing line pairs" << std::endl;
  

  Math1D::Vector<double> frac_solution(mesh.nFaces());
  for (uint k=0; k < mesh.nFaces(); k++) {
    frac_solution[k] = lp_solution[k];
  }

  if (mesh.nFaces() <= 20000) {
    Petter::statusTry("Saving SVG...");
    mesh.draw_labels_with_pairs(options.base_filename_ + ".lp_simple.svg",lp_solution,edge_pairs,xDim,yDim);
    Petter::statusOK();
  }

  if (nNonInt > 0 && !options.convex_prior_) { // Thresholding does not work with convex prior, because 
                                               // the problem will most likely be infeasible

    std::cerr << nNonInt << " non-integral region variables. Computing thresholded solution" << std::endl;

    //fix the region variables (according to the above thresholding) and re-run the solver
    // to get the boundary variables

    if (enforce_consistent_boundaries) {

      //in this mode we want to make sure that no fractional boundary variables remain.
      std::cerr << "disable pairs that do not agree with the fixed segmentation" << std::endl;

      for (uint p=0; p < edge_pairs.size(); p++) {

        uint edge1 = edge_pairs[p].first_edge_idx_;
        uint edge2 = edge_pairs[p].second_edge_idx_;

        double drop1 = 0.0;
        for (uint i=0; i < mesh.adjacent_faces(edge1).size(); i++) {
          uint face = mesh.adjacent_faces(edge1)[i];
          drop1 += var_ub[face] * mesh.match(face,edge1);      
        }

        double drop2 = 0.0;
        for (uint i=0; i < mesh.adjacent_faces(edge2).size(); i++) {
          uint face = mesh.adjacent_faces(edge2)[i];
          drop2 += var_ub[face] * mesh.match(face,edge2);      
        }

        if (drop1 == 0.0 || drop2 == 0.0) {
          var_ub[edge_pair_var_offs + 2*p] = 0.0;
          var_ub[edge_pair_var_offs + 2*p+1] = 0.0;
        }
        else {

          if (mesh.edge(edge1).from_idx_ == edge_pairs[p].common_point_idx_) {

            drop1 *= -1.0;
            drop2 *= -1.0;
          }

          if (drop1 < 0.0)
            var_ub[edge_pair_var_offs + 2*p+1] = 0.0;
          else
            var_ub[edge_pair_var_offs + 2*p] = 0.0;
        }
      }
    }

    std::cerr << "tightening of bounds done" << std::endl;

    if (solver == "clp") {
      for (uint i=0; i < (uint) lpSolver.getNumCols(); i++) {

        if (! (var_lb[i] <= var_ub[i]))
          std::cerr << i << ", lb: " << var_lb[i] << ", ub: " << var_ub[i] << std::endl;

        assert(var_lb[i] <= var_ub[i]);
        lpSolver.setColumnBounds(i,var_lb[i],var_ub[i]);
      }
    }
#ifdef HAS_GUROBI
    if (solver == "gurobi") {	
      for (uint v=0; v < mesh.nFaces(); v++) {
        GRBsetdblattrelement(grb_model,"LB",v,var_lb[v]);
        GRBsetdblattrelement(grb_model,"UB",v,var_ub[v]);
      }
    }
#endif
#ifdef HAS_XPRESS
    if (solver == "xpress") {
      Math1D::Vector<int> indices(mesh.nFaces());
      for (uint v=0; v < mesh.nFaces(); v++)
        indices[v] = v;

      Storage1D<char> indicator(mesh.nFaces(),'L');

      XPRSchgbounds(xp_prob, mesh.nFaces(), indices.direct_access(), indicator.direct_access(), var_lb.direct_access());
      indicator.set_constant('U');
      XPRSchgbounds(xp_prob, mesh.nFaces(), indices.direct_access(), indicator.direct_access(), var_ub.direct_access());
    }
#endif
#ifdef HAS_CPLEX
    if (solver == "cplex") {

      Math1D::Vector<int> indices(mesh.nFaces());
      for (uint v=0; v < mesh.nFaces(); v++)
        indices[v] = v;

      Storage1D<char> indicator(mesh.nFaces(),'L');
      CPXchgbds(cp_env, cp_lp, mesh.nFaces(), indices.direct_access(), indicator.direct_access(), var_lb.direct_access());
      indicator.set_constant('U');
      CPXchgbds(cp_env, cp_lp, mesh.nFaces(), indices.direct_access(), indicator.direct_access(), var_ub.direct_access());
    }
#endif

    if (solver == "clp") {
      lpSolver.dual();
      lp_solution = lpSolver.primalColumnSolution();

      error = 1 - lpSolver.isProvenOptimal();
      assert(!error || enforce_consistent_points);
    }

#ifdef HAS_GUROBI
    if (solver == "gurobi") {

      int error = GRBupdatemodel(grb_model);
      error = GRBoptimize(grb_model);
      assert(!error);

      for (uint v=0; v < nVars; v++)
        GRBgetdblattrelement(grb_model,"X",v, GUROBI_solution.direct_access()+v);

      lp_solution = GUROBI_solution.direct_access();
    }
#endif
#ifdef HAS_XPRESS
    if (solver == "xpress") {

      XPRSlpoptimize(xp_prob,"b");      
      XPRSgetlpsol(xp_prob, XPRESS_solution.direct_access(), 0, 0, 0); 

      lp_solution = XPRESS_solution.direct_access();
    }
#endif
#ifdef HAS_CPLEX
    if (solver == "cplex") {

      int cpx_status = CPXbaropt(cp_env,cp_lp);
      //int cpx_status = CPXlpopt(cp_env,cp_lp);

      if ( cpx_status ) {
        fprintf (stderr, "Failed to optimize MIP.\n");
        exit(1);
      }

      CPXsolution (cp_env, cp_lp, NULL, NULL, CPLEX_solution.direct_access(), NULL, NULL, NULL);

      lp_solution = CPLEX_solution.direct_access();
    }
#endif

    uint nIter = 0;
    while (prevent_crossings) {

      nIter++;

      std::cerr << "##### constraint generation (for thresholded solution) iter #" << nIter << std::endl;

      uint nConstraintsAdded = 0;
      for (uint p=0; p < mesh.nPoints(); p++) {

	
	std::set<uint> active_pair;
	
	double sum = 0.0;
	for (std::vector<uint>::iterator it = point_pairs[p].begin(); it != point_pairs[p].end(); it++) {
	  double contrib = lp_solution[edge_pair_var_offs +  2*(*it)] + lp_solution[edge_pair_var_offs +  2*(*it) + 1];
	  if (fabs(contrib) >= 0.0025) {
	    sum += contrib;
	    active_pair.insert(*it);
	  }
	}
	
	if (sum >= 1.0025 && active_pair.size() > 1) {
	  std::cerr << "point #" << p << " might be a crossing point. " 
		    << active_pair.size() << " active pairs. " << std::endl;
	  
	  // this handles all pairwise and triple-constellations exactly,
	  // as well as SOME of the higher ones. Note that we don't need to add
	  // low order interactions if they are subsumed in a high order one
	  
	  uint nPrevAdded = nConstraintsAdded;
	  
	  for (uint k1=0; k1 < point_pairs[p].size()-1; k1++) {
	    
	    uint pair1 = point_pairs[p][k1];
	    
	    if (active_pair.find(pair1) == active_pair.end())
	      continue;
	    
	    for (uint k2=k1+1; k2 < point_pairs[p].size(); k2++) {
	    
	      uint pair2 = point_pairs[p][k2];
	      
	      if (active_pair.find(pair2) == active_pair.end())
		continue;
	      
	      if (!line_pairs_with_meeting_point_cross(mesh, edge_pairs[pair1], edge_pairs[pair2]))
		continue;
	      
	      double sum = lp_solution[edge_pair_var_offs +  2*pair1] + lp_solution[edge_pair_var_offs +  2*pair1 + 1]
		+ lp_solution[edge_pair_var_offs +  2*pair2] + lp_solution[edge_pair_var_offs +  2*pair2 + 1];
	      
	      //std::cerr << "checking pair" << std::endl;
	      
	      std::vector<uint> base(2);
	      base[0] = pair1;
	      base[1] = pair2;
	      
	      std::vector<uint> addons;
	      for (uint k3=k2+1; k3 < point_pairs[p].size(); k3++) {
		
		uint pair3 = point_pairs[p][k3];
		
		if (active_pair.find(pair3) != active_pair.end()
		    && line_pairs_with_meeting_point_cross(mesh, edge_pairs[pair1], edge_pairs[pair3])
		    && line_pairs_with_meeting_point_cross(mesh, edge_pairs[pair2], edge_pairs[pair3]))
		  addons.push_back(point_pairs[p][k3]);
	      }
	      
	      std::vector<std::vector<uint> > constraint_list;
	      if (addons.empty()) {
		if (sum >= 1.005)
		  constraint_list.push_back(base);
	      }
	      else {
#if 1
		//handle higher interactions:
		//NOTE: currently some of the generated constraints can be redundant to others
		for (uint k=0; k < addons.size(); k++) {
		  
		  std::vector<uint> new_list = base;
		  new_list.push_back(addons[k]);
		  
		  double cur_sum = sum + lp_solution[edge_pair_var_offs +  2*addons[k]] 
		    + lp_solution[edge_pair_var_offs +  2*addons[k] + 1];
		  
		  for (uint l=k+1; l < addons.size(); l++) {
		    
		    uint p1 = addons[l];
		    
		    bool compatible = true;
		    for (uint j=2; j < new_list.size(); j++) {
		      
		      uint p2 = new_list[j];
		      
		      if (!line_pairs_with_meeting_point_cross(mesh, edge_pairs[p1], edge_pairs[p2])) {
			compatible = false;
			break;
		      }
		    }

		    if (compatible) {
		      new_list.push_back(addons[l]);
		      cur_sum += lp_solution[edge_pair_var_offs +  2*addons[l]] 
			+ lp_solution[edge_pair_var_offs +  2*addons[l] + 1];
		    }
		  }
		  
		  if (cur_sum >= 1.005) {
		    std::cerr << "sum: " << cur_sum << std::endl;
		    constraint_list.push_back(new_list);
		  }
		}
#else
		//fallback to pairwise:
		if (sum >= 1.005)
		  constraint_list.push_back(base);
#endif
	      }
	    
	      for (uint i=0; i < constraint_list.size(); i++) {
	      
		std::vector<uint> pairs = constraint_list[i];
		
		//add the constraint
		int* cols = new int[2*pairs.size()];
		double* coeffs = new double[2*pairs.size()];
		for (uint k=0; k < 2*pairs.size(); k++)
		  coeffs[k] = 1.0;
		
		for (uint p=0; p < pairs.size(); p++) {
		  cols[2*p] = edge_pair_var_offs +  2*pairs[p];
		  cols[2*p+1] = edge_pair_var_offs +  2*pairs[p]+1;
		}
		
		//note: adding constraints separately is VERY inefficient
		if (solver == "clp")
		  lpSolver.addRow(2*pairs.size(), cols, coeffs, 0.0, 1.0);
#ifdef HAS_GUROBI
		if (solver == "gurobi") {
		  GRBaddconstr(grb_model,2*pairs.size(),cols,coeffs,'L',1.0,NULL);
		}
#endif
#ifdef HAS_XPRESS
		if (solver == "xpress") {
		
		  double new_rhs[1] = {1.0};
		  double new_range[1] = {0.0};
		  int new_start[2] = {0,2*pairs.size()};
		  XPRSaddrows(xp_prob, 1, 2*pairs.size(), "L", new_rhs,new_range,new_start,cols,coeffs);
		}
#endif
#ifdef HAS_CPLEX
		if (solver == "cplex") {
		
		  double new_rhs[1] = {1.0};
		  int new_start[2] = {0,2*pairs.size()};
		  
		  CPXaddrows(cp_env, cp_lp, 0, 1, 2*pairs.size(), new_rhs, "L", new_start, cols, coeffs,  NULL, NULL);
		}
#endif		
		nConstraintsAdded++;
		
		delete cols;
		delete coeffs;
	      }
	    }
	  }
	}
      }

      std::cerr << "added " << nConstraintsAdded << " constraints" << std::endl;

      if (nConstraintsAdded > 0) {

        if (solver == "clp") {
          lpSolver.dual();
          lp_solution = lpSolver.primalColumnSolution();
        }

#ifdef HAS_GUROBI
        if (solver == "gurobi") {

          int error = GRBupdatemodel(grb_model);
          error = GRBoptimize(grb_model);
          assert(!error);

          for (uint v=0; v < nVars; v++)
            GRBgetdblattrelement(grb_model,"X",v, GUROBI_solution.direct_access()+v);

          lp_solution = GUROBI_solution.direct_access();
        }
#endif
#ifdef HAS_XPRESS
        if (solver == "xpress") {

          XPRSlpoptimize(xp_prob,"b");      
          XPRSgetlpsol(xp_prob, XPRESS_solution.direct_access(), 0, 0, 0); 

          lp_solution = XPRESS_solution.direct_access();
        }
#endif
#ifdef HAS_CPLEX
        if (solver == "cplex") {

          int cpx_status = CPXbaropt(cp_env,cp_lp);
          //int cpx_status = CPXlpopt(cp_env,cp_lp);

          if ( cpx_status ) {
            fprintf (stderr, "Failed to optimize MIP.\n");
            exit(1);
          }

          CPXsolution (cp_env, cp_lp, NULL, NULL, CPLEX_solution.direct_access(), NULL, NULL, NULL);

          lp_solution = CPLEX_solution.direct_access();
        }
#endif
      }
      else
        break;
    }

    double thresh_energy = energy_offset;
    uint nNonIntThresh = 0;
    for (uint i=0; i < nVars; i++) {
      double val = lp_solution[i]; 
      if (val > 0.01 && val < 0.99)
        nNonIntThresh++;
      thresh_energy += cost[i] * val;
    }

    std::cerr << "energy of thresholded solution: " << thresh_energy 
      << "  (" << nNonIntThresh << " non-integral variables)" << std::endl;
  }

  int nOpposingLinePairs_after = 0;

  for (uint v=mesh.nFaces(); v < nVars; v+=2) {

    double pair_sum = lp_solution[v] + lp_solution[v+1]; 

    //if (pair_sum > 1.02) {
    if (lp_solution[v] > 0.05 && lp_solution[v+1] > 0.05) {
      nOpposingLinePairs_after++;
    }
  }

  std::cerr << nOpposingLinePairs_after << " opposing line pairs" << std::endl;

  if (mesh.nFaces() <= 20000) {
    Petter::statusTry("Saving SVG...");
    mesh.draw_labels_with_pairs(options.base_filename_ + ".final.svg",lp_solution,edge_pairs,xDim,yDim);
    Petter::statusOK();
  }

  Petter::statusTry("Building output...");

  Math1D::Vector<int> labels(mesh.nFaces());

  uint out_factor = options.output_factor_;

  segmentation.resize(xDim*out_factor,yDim*out_factor);
  Math2D::Matrix<double> frac_seg(xDim*out_factor,yDim*out_factor);

  mesh.enlarge(out_factor,out_factor);

  Storage1D<PixelFaceRelation> shares;
  Math1D::Vector<uint> share_start;

  //re-compute pixel shares for the now larger mesh
  compute_pixel_shares(mesh, out_factor*xDim, out_factor*yDim, shares, share_start);

  for (uint i=0;i<mesh.nFaces();++i) {
    if (lp_solution[i] > 0.99) {
      labels[i] = 1;
    }
    else if (lp_solution[i] < 0.01) {
      labels[i] = 0;
    }
    else {
      labels[i] = -1;
    }
  }

  for (uint y=0; y < yDim*out_factor; y++) {
    for (uint x=0; x < xDim*out_factor; x++) {

      double sum = 0.0;
      double frac_sum = 0.0;
      for (uint k= share_start[y*(xDim*out_factor)+x]; k < share_start[y*(xDim*out_factor)+x+1]; k++) {
        uint face = shares[k].face_idx_;
        sum += mesh.convex_area(face) * shares[k].share_ * lp_solution[face];
        frac_sum += mesh.convex_area(face) * shares[k].share_ * frac_solution[face];
      }
      double seg = int(sum*255.0 + 0.5);
      if (seg > 255) seg = 255;
      segmentation(x,y) = seg;
      frac_seg(x,y) = frac_sum * 255.0;
    }
  }

  frac_seg.savePGM(options.base_filename_ + ".frac.pgm",255);

  Petter::statusOK();

  

#ifdef HAS_GUROBI
  if (grb_model != NULL)
    GRBfreemodel(grb_model);
  if (grb_env != NULL)
    GRBfreeenv(grb_env);
#endif
#ifdef HAS_CPLEX
  if (cp_lp != NULL)
    CPXfreeprob (cp_env, &cp_lp);
  if (cp_env != NULL)
    CPXcloseCPLEX (&cp_env);
#endif
#ifdef HAS_XPRESS
  if (xp_prob != NULL) {
    XPRSdestroyprob(xp_prob);
    XPRSfree();
  }
#endif



  //
  // Open the log file (append mode)
  //
  // <lambda> <gamma> <Clp time (s)> <fractional (%)> <fractional pairs (#)> <pairs overlapping (#)>
  //
  std::string logfile_name = options.base_filename_ + ".lplog";
  std::ofstream logfile(logfile_name.c_str(), std::ios::app);
  logfile << options.lambda_ << " "     // 0 
          << options.gamma_  << " "     // 1
          << solverTime      << " "     // 2
		  << 100.0*double(nNonInt + nNonIntAuxVars)/double(nVars) << " " // 3
          << nNonIntAuxVars  << ' '     // 4
          << nOpposingLinePairs << ' '  // 5
		  << mesh.nFaces()   << ' '     // 6
		  << edge_pairs.size()  << ' '  // 7
		  << lp_energy       << ' '     // 8
		  ;


  //Close logfile
  logfile << std::endl;

  return energy;
}
