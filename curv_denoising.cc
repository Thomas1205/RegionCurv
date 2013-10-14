/******** written by Thomas Schoenemann as an employee of Lund University, Sweden, October 2010 ****/

#include "curv_denoising.hh"
#include "mesh2D.hh"
#include "sparse_matrix_description.hh"
#include <coin/ClpFactorization.hpp>
#include <coin/OsiClpSolverInterface.hpp>


#include "tensor.hh"
#include "timing.hh"

//#define USE_GUROBI
//#define USE_CPLEX
//#define USE_XPRESS

#ifdef HAS_GUROBI
#include "gurobi_c++.h"
#endif

#ifdef HAS_CPLEX
#include <ilcplex/cplex.h>
#endif

#ifdef HAS_XPRESS
#include "xprs.h" 
#endif

void curv_denoise(const Math3D::Tensor<float>& image, const LPSegOptions& options,
                  Math3D::Tensor<double>& denoised_image, uint nBins) {

  double lambda = options.lambda_;
  double gamma = options.gamma_;
  int neighborhood = options.neighborhood_;
  bool enforce_consistent_boundaries = options.enforce_consistent_boundaries_;
  bool enforce_regionedge = options.enforce_regionedge_;
  bool bruckstein = options.bruckstein_;
  bool light_constraints = options.light_constraints_;

  std::string solver = options.solver_;
  bool solver_known = false;
  
  assert(neighborhood <= 16); 

  double im_const = image.min();

  uint xDim = image.xDim();
  uint yDim = image.yDim();
  uint zDim = image.zDim();

  if (zDim > 1) {
    std::cerr << "WARNING: currently only the first channel is denoised" << std::endl;
  }

  Math2D::Matrix<float> gray_image(xDim,yDim);
  for (uint y=0; y < yDim; y++)
    for (uint x=0; x < xDim; x++)
      gray_image(x,y) = image(x,y,0);

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
      generate_adaptive_mesh(gray_image, mesh, neighborhood, options.adaptive_mesh_n_);
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
      generate_adaptive_mesh(gray_image, mesh, neighborhood, options.adaptive_mesh_n_);
    }
  }

  uint light_factor = (light_constraints) ? 1 : 2;

  Storage1D<PixelFaceRelation> shares;
  Math1D::Vector<uint> share_start;
  compute_pixel_shares(mesh, xDim, yDim, shares, share_start);

  Math1D::NamedVector<double> face_intensity(mesh.nFaces(),0.0,MAKENAME(face_intensity));
  Math1D::NamedVector<double> face_weight(mesh.nFaces(),0.0,MAKENAME(face_intensity));
  for (uint i=0; i < xDim*yDim; i++) {
    
    double cur_intensity = image.direct_access(i) - im_const;
    for (uint k=share_start[i]; k < share_start[i+1]; k++) {
      
      uint face_idx = shares[k].face_idx_;
      double fac = shares[k].share_ * mesh.convex_area(face_idx);

      face_intensity[face_idx] += fac * cur_intensity;
      face_weight[face_idx] += fac;
    }
  }

  for (uint f=0; f < mesh.nFaces(); f++)
    face_intensity[f] /= face_weight[f];

  std::vector<Mesh2DEdgePair> edge_pairs;
  mesh.generate_edge_pair_list(edge_pairs);

  std::cerr << edge_pairs.size() << " edge pairs." << std::endl;

  double bin_size = (image.max() - image.min()) / nBins;

  uint nVars = (nBins+2)*mesh.nFaces() + 2*nBins*edge_pairs.size();
  uint nConstraints = 3*nBins*mesh.nEdges() //surface and boundary
    + mesh.nFaces() //for the absolutes
    + (nBins-1)*mesh.nFaces(); //level constraints

  const uint abs_var_offs = nBins*mesh.nFaces();
  const uint edge_pair_var_offs = (nBins+2)*mesh.nFaces();

  const uint abs_con_offs = 3*nBins*mesh.nEdges();
  const uint level_con_offs = abs_con_offs + mesh.nFaces();

  const uint consistency_con_offs = nConstraints;
  if (enforce_consistent_boundaries) {
    nConstraints += nBins*light_constraints*mesh.nEdges();
  }
  const uint regionedge_constraints_offs = nConstraints;
  if (enforce_regionedge) {
    nConstraints += 2*nBins*light_constraints*mesh.nEdges();
  }

  Math1D::NamedVector<double> var_lb(nVars,0.0,MAKENAME(var_lb));
  Math1D::NamedVector<double> var_ub(nVars,(image.max()-im_const)/nBins,MAKENAME(var_ub));
  Math1D::NamedVector<double> cost(nVars,0.0,MAKENAME(cost));

  //set cost for the absolute variables
  for (uint f=0; f < mesh.nFaces(); f++) {

    uint v = abs_var_offs + 2*f;

    cost[v] = mesh.convex_area(f);
    var_ub[v] = image.max() - im_const; //TODO: set bounds based on corresponding image value
    cost[v+1] = mesh.convex_area(f);
    var_ub[v+1] = image.max() - im_const; //TODO: set bounds based on corresponding image value    
  }

  uint edge_var_shift = 2*edge_pairs.size();
  uint edge_con_shift = mesh.nEdges();
  uint var_shift = mesh.nFaces();

  for (uint j=0; j < edge_pairs.size(); j++) {
    
    uint first = edge_pairs[j].first_edge_idx_;
    uint nFirstAdjacent = mesh.adjacent_faces(first).size();
    
    uint second = edge_pairs[j].second_edge_idx_;
    uint nSecondAdjacent = mesh.adjacent_faces(second).size();

    double weight = 0.0;
    if (nFirstAdjacent > 1)
      weight += 0.5*lambda*mesh.edge_length(first);
    if (nSecondAdjacent > 1)
      weight += 0.5*lambda*mesh.edge_length(second);
    
    //do not penalize the image corners for their curvature
    if (nFirstAdjacent > 1 || nSecondAdjacent > 1)
      weight += gamma * curv_weight(mesh,edge_pairs[j],2.0,bruckstein);

    //double weight = 0.5*lambda*(mesh.edge_length(first) + mesh.edge_length(second))
    //+ gamma * curv_weight(mesh,edge_pairs[j],2.0,bruckstein);

    for (uint b=0; b < nBins; b++) {
      cost[edge_pair_var_offs+b*edge_var_shift+2*j] = weight;
      cost[edge_pair_var_offs+b*edge_var_shift+2*j+1] = weight;
    }

    /*** check if (at the image border) one of the edge pairs is impossible. if so, set its upper bound to 0 ***/
    uint edge = first;
    if (mesh.adjacent_faces(edge).size() == 1) {
      int match = mesh.match(mesh.adjacent_faces(edge)[0],edge);

      uint y = edge_pair_var_offs+2*j;
      
      if (edge_pairs[j].common_point_idx_ == mesh.edge(edge).from_idx_)
        match *= -1;

      if (match == -1) {
        for (uint b=0; b < nBins; b++)
          var_ub[y+1+b*edge_var_shift] = 0.0;
      }
      else if (match == 1) {
        for (uint b=0; b < nBins; b++)
          var_ub[y+b*edge_var_shift] = 0.0;
      }
    }
    
    edge = second;
    if (mesh.adjacent_faces(edge).size() == 1) {
      int match = mesh.match(mesh.adjacent_faces(edge)[0],edge);

      uint y = edge_pair_var_offs+2*j;
      
      if (edge_pairs[j].common_point_idx_ == mesh.edge(edge).from_idx_)
        match *= -1;
      
      if (match == -1) {
        for (uint b=0; b < nBins; b++)
          var_ub[y+b*edge_var_shift] = 0.0;
      }
      else if (match == 1) {
        for (uint b=0; b < nBins; b++)
          var_ub[y+1+b*edge_var_shift] = 0.0;
      }
    }
  }

  Math1D::NamedVector<double> rhs_lower(nConstraints,0.0,MAKENAME(rhs_lower));
  Math1D::NamedVector<double> rhs_upper(nConstraints,0.0,MAKENAME(rhs_upper));

  for (uint f=0; f < mesh.nFaces(); f++) {
    uint cur_con_offs = abs_con_offs + f;
    rhs_lower[cur_con_offs] = face_intensity[f];
    rhs_upper[cur_con_offs] = face_intensity[f];
  }

  uint nEntries = 2*nBins*mesh.nEdges() + (nBins+2)*mesh.nFaces() + 8*nBins*edge_pairs.size(); 
  if (enforce_consistent_boundaries) {
    nEntries += nBins*(2*light_constraints*edge_pairs.size() + light_constraints*mesh.nEdges());
    //we also allocate space for the slack variables used with the convex solver
  }
  if (enforce_regionedge) {
    nEntries += nBins*light_constraints*(8*edge_pairs.size() //TODO: not exact number
                                         + 2*mesh.nEdges()); //we also allocate space for the slack variables used with the convex solver
  }

  SparseMatrixDescription<double> lp_descr(nEntries, nConstraints, nVars);
  
  /**** a) code surface continuation constraints *****/

  for (uint j=0; j < mesh.nEdges(); j++) {

    const std::vector<uint>& adjacent_faces = mesh.adjacent_faces(j);
    for (std::vector<uint>::const_iterator it = adjacent_faces.begin();
         it != adjacent_faces.end(); it++) {

      for (uint b=0; b < nBins; b++) {	
        lp_descr.add_entry(j+b*edge_con_shift,*it+b*var_shift,mesh.match(*it,j));
      }
    }
  }

  std::cerr << "A" << std::endl;

  for (uint j=0; j < edge_pairs.size(); j++) {

    uint first_edge = edge_pairs[j].first_edge_idx_;
    uint second_edge = edge_pairs[j].second_edge_idx_;

    uint middle_point = edge_pairs[j].common_point_idx_;

    if (mesh.edge(first_edge).to_idx_ == middle_point) {
      for (uint b=0; b < nBins; b++) {	
        lp_descr.add_entry(first_edge+b*edge_con_shift,edge_pair_var_offs+2*j+b*edge_var_shift,1);
      }
    }
    else {
      for (uint b=0; b < nBins; b++) {	
        lp_descr.add_entry(first_edge+b*edge_con_shift,edge_pair_var_offs+2*j+b*edge_var_shift,-1);
      }
    }

    if (mesh.edge(second_edge).to_idx_ == middle_point) {
      for (uint b=0; b < nBins; b++) {
        lp_descr.add_entry(second_edge+b*edge_con_shift,edge_pair_var_offs+2*j+1+b*edge_var_shift,1);
      }
    }
    else {
      for (uint b=0; b < nBins; b++) {	
        lp_descr.add_entry(second_edge+b*edge_con_shift,edge_pair_var_offs+2*j+1+b*edge_var_shift,-1);
      }
    }
  }

  std::cerr << "B" << std::endl;

  /***** b) constraints to model the absolutes  ******/
  for (uint f=0; f < mesh.nFaces(); f++) {

    uint row = abs_con_offs + f;
    for (uint b=0; b < nBins; b++) {
      lp_descr.add_entry(row, f + b*var_shift, 1.0);
    }
    lp_descr.add_entry(row, abs_var_offs + 2*f, 1.0);
    lp_descr.add_entry(row, abs_var_offs + 2*f+1, -1.0);
  }

  std::cerr << "C" << std::endl;

  /**** c) code boundary continuation constraints *****/
  uint boundary_con_offset = nBins*mesh.nEdges();

  for (uint j=0; j < edge_pairs.size(); j++) {
    
    uint first_edge = edge_pairs[j].first_edge_idx_;
    uint second_edge = edge_pairs[j].second_edge_idx_;
    
    uint middle_point = edge_pairs[j].common_point_idx_;
    
    if (mesh.edge(first_edge).to_idx_ == middle_point) {      
      for (uint b=0; b < nBins; b++) {	
        lp_descr.add_entry(boundary_con_offset + 2*first_edge + b*2*edge_con_shift, 
                           edge_pair_var_offs+2*j + b*edge_var_shift, 1);
        lp_descr.add_entry(boundary_con_offset + 2*first_edge+1 + b*2*edge_con_shift, 
                           edge_pair_var_offs+2*j+1 + b*edge_var_shift, -1);
      }
    }
    else {
      for (uint b=0; b < nBins; b++) {	
        lp_descr.add_entry(boundary_con_offset + 2*first_edge+1 + b*2*edge_con_shift,
                           edge_pair_var_offs+2*j + b*edge_var_shift, 1);
        lp_descr.add_entry(boundary_con_offset + 2*first_edge + b*2*edge_con_shift, 
                           edge_pair_var_offs+2*j+1 + b*edge_var_shift, -1);
      }
    }
    
    if (mesh.edge(second_edge).from_idx_ == middle_point) {
      for (uint b=0; b < nBins; b++) {	
        lp_descr.add_entry(boundary_con_offset + 2*second_edge + b*2*edge_con_shift, 
                           edge_pair_var_offs+2*j + b*edge_var_shift, -1);
        lp_descr.add_entry(boundary_con_offset + 2*second_edge+1 + b*2*edge_con_shift, 
                           edge_pair_var_offs+2*j+1 + b*edge_var_shift, 1);
      }
    }
    else {
      for (uint b=0; b < nBins; b++) {	
        lp_descr.add_entry(boundary_con_offset + 2*second_edge+1 + b*2*edge_con_shift, 
                           edge_pair_var_offs+2*j + b*edge_var_shift, -1);
        lp_descr.add_entry(boundary_con_offset + 2*second_edge + b*2*edge_con_shift, 
                           edge_pair_var_offs+2*j+1 + b*edge_var_shift, 1);
      }
    }
  }

  std::cerr << "D" << std::endl;

  //level constraints
  if (nBins > 1) {

    for (uint b=0; b < nBins-1; b++) {
      for (uint f=0; f < mesh.nFaces(); f++) {
        uint row = level_con_offs + f + b*mesh.nFaces();
        rhs_lower[row] = -1e300;
        lp_descr.add_entry(row,f+b*mesh.nFaces(),-1.0);
        lp_descr.add_entry(row,f+(b+1)*mesh.nFaces(),1.0);
      }
    }
  }

  uint nStandardEntries = lp_descr.nEntries();

  if (enforce_consistent_boundaries) {

    //constraints in words: for each oriented edge, the pairs that start with this oriented edge 
    // and the pairs that end in the oppositely oriented edge may not sum to more than 1.0
    // (i.e. they are mutually exclusive)

    double range = bin_size;

    for (uint c=consistency_con_offs; c < regionedge_constraints_offs; c+=light_factor) { 
      rhs_upper[c] = range;
      if (!light_constraints)
        rhs_upper[c+1] = range;
    }

    for (uint b=0; b < nBins; b++) {

      uint cur_con_offs = consistency_con_offs + light_factor*b*mesh.nEdges();

      for (uint j=0; j < edge_pairs.size(); j++) {
	
        uint middle_point = edge_pairs[j].common_point_idx_;
	
        uint first_edge = edge_pairs[j].first_edge_idx_;
        uint second_edge = edge_pairs[j].second_edge_idx_;
	
        if (mesh.edge(first_edge).to_idx_ == middle_point) {      
          lp_descr.add_entry(cur_con_offs + light_factor*first_edge, edge_pair_var_offs+2*j + b*edge_var_shift, 1);
          lp_descr.add_entry(cur_con_offs + light_factor*first_edge, edge_pair_var_offs+2*j+1 + b*edge_var_shift, 1);
        }
        else if (!light_constraints) {
          lp_descr.add_entry(cur_con_offs + light_factor*first_edge+1, edge_pair_var_offs+2*j + b*edge_var_shift, 1);
          lp_descr.add_entry(cur_con_offs + light_factor*first_edge+1, edge_pair_var_offs+2*j+1 + b*edge_var_shift, 1);
        }
	
        if (mesh.edge(second_edge).to_idx_ == middle_point) {
          lp_descr.add_entry(cur_con_offs + light_factor*second_edge, edge_pair_var_offs+2*j + b*edge_var_shift, 1);
          lp_descr.add_entry(cur_con_offs + light_factor*second_edge, edge_pair_var_offs+2*j+1 + b*edge_var_shift, 1);
        }
        else if (!light_constraints) {
          lp_descr.add_entry(cur_con_offs + light_factor*second_edge+1, edge_pair_var_offs+2*j + b*edge_var_shift, 1);
          lp_descr.add_entry(cur_con_offs + light_factor*second_edge+1, edge_pair_var_offs+2*j+1 + b*edge_var_shift, 1);
        }
      }
    }
  }
  
  if (enforce_regionedge) {
    Petter::statusTry("Coding boundary/edge pair constr....");

    for (uint b=0; b < nBins; b++) {

      uint rowoff = regionedge_constraints_offs +2*b*light_factor*mesh.nEdges();
      
      for (uint edge=0; edge < mesh.nEdges(); edge++) {
	
        //Get the two adjacent faces
        if (mesh.adjacent_faces(edge).size() != 2)
          {
            //One of the edges is at the border of the image
	  
            continue;
          }
      
        uint x1 = mesh.adjacent_faces(edge)[0] + b*var_shift;
        uint x2 = mesh.adjacent_faces(edge)[1] + b*var_shift;
	
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
	
        //NOTE (important for the construction with slacks): the binding constraints are in both cases the upper bounds
        rhs_lower[rowoff+2*light_factor*edge]   = 0;
        rhs_upper[rowoff+2*light_factor*edge]   = 2*bin_size;
        rhs_lower[rowoff+2*light_factor*edge+1] = -2*bin_size;
        rhs_upper[rowoff+2*light_factor*edge+1] = 0;
        if (!light_constraints) {
          rhs_lower[rowoff+2*light_factor*edge+2] = 0;
          rhs_upper[rowoff+2*light_factor*edge+2] = 2*bin_size;
          rhs_lower[rowoff+2*light_factor*edge+3] = -2*bin_size;
          rhs_upper[rowoff+2*light_factor*edge+3] = 0;
        }
      }
      
      for (uint j=0; j < edge_pairs.size(); j++) {
        uint first = edge_pairs[j].first_edge_idx_;
        uint second = edge_pairs[j].second_edge_idx_;	
	
        uint y = edge_pair_var_offs + 2*j + b*edge_var_shift;
	
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
    }
    
    Petter::statusOK();
  }

  std::cerr << "sorting" << std::endl;
  
  Math1D::NamedVector<uint> row_start(nConstraints+1,MAKENAME(row_start));
  lp_descr.sort_by_row(row_start);

  std::cerr << "sorted" << std::endl;

  //const double* lp_solution = 0;

  int error;

  Math1D::Vector<double> lp_solution(nVars);

#ifdef HAS_GUROBI

  if (solver == "gurobi") {

    solver_known = true;

    GRBenv   *grb_env   = NULL;
    GRBmodel *grb_model = NULL;

    /* Create environment */
    
    error = GRBloadenv(&grb_env,NULL);
    //GRBsetintparam(grb_env, GRB_INT_PAR_METHOD, GRB_LPMETHOD_DUAL);
    GRBsetintparam(grb_env, GRB_INT_PAR_METHOD, GRB_METHOD_BARRIER);
    //GRBsetdblparam(grb_env, "BarConvTol", 1e-10);
    GRBsetintparam(grb_env, "Crossover", 0);
    //GRBsetintparam(grb_env, "CrossoverBasis", 1);
    GRBsetintparam(grb_env, "Presolve", 1);
    GRBsetintparam(grb_env, "PrePasses", 2);

    assert (!error && grb_env != NULL);

    /* Create an empty model */
    error = GRBnewmodel(grb_env, &grb_model, "curv-denoise-lp", 0, NULL, NULL, NULL, NULL, NULL);
    assert(!error);
    
    NamedStorage1D<char> vtype(nVars,GRB_CONTINUOUS,MAKENAME(vtype));
  
    error = GRBaddvars(grb_model,nVars,0,NULL,NULL,NULL,cost.direct_access(),var_lb.direct_access(),
                       var_ub.direct_access(),vtype.direct_access(),NULL);
    assert(!error);
    
    error = GRBupdatemodel(grb_model);
    assert(!error);
    
    std::cerr << "F" << std::endl;
    
    for (uint c=0; c < row_start.size()-1; c++) {
      
      if (rhs_lower[c] == rhs_upper[c])
        error = GRBaddconstr(grb_model, row_start[c+1]-row_start[c], ((int*) lp_descr.col_indices()) + row_start[c], 
                             lp_descr.value() + row_start[c], GRB_EQUAL, rhs_upper[c], NULL);
      else
        GRBaddrangeconstr(grb_model, row_start[c+1]-row_start[c], ((int*) lp_descr.col_indices()) + row_start[c], 
                          lp_descr.value() + row_start[c], rhs_lower[c], rhs_upper[c], NULL);      
      
      assert(!error);
    }
    
    std::cerr << "starting Gurobi optimization" << std::endl;
    
    /* Optimize model */
    error = GRBoptimize(grb_model);
    assert(!error);
    
    for (uint v=0; v < nVars; v++) {
      GRBgetdblattrelement(grb_model,"X",v, lp_solution.direct_access()+v);
    }
    
    GRBfreemodel(grb_model);
    GRBfreeenv(grb_env);
  }
#endif
#ifdef HAS_CPLEX

  if (solver == "cplex") {

    solver_known = true;

    CPXENVptr     env = NULL;
    CPXLPptr      lp = NULL;
    int status = 0;

    /* Initialize the CPLEX environment */
    
    env = CPXopenCPLEX (&status);
    //CPXsetintparam(env, CPX_PARAM_STARTALG, CPX_ALG_BARRIER);
    //CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);
    //CPXsetintparam(env, CPX_PARAM_PREIND, CPX_OFF);
    //CPXsetintparam(env, CPX_PARAM_PREPASS, 0);
    CPXsetintparam(env, CPX_PARAM_BARCROSSALG, -1);
    
    /* If an error occurs, the status value indicates the reason for
       failure.  A call to CPXgeterrorstring will produce the text of
       the error message.  Note that CPXopenCPLEX produces no output,
       so the only way to see the cause of the error is to use
       CPXgeterrorstring.  For other CPLEX routines, the errors will
       be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON.  */
    
    if ( env == NULL ) {
      char  errmsg[1024];
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      exit(1);
    }
    
    /* Turn on output to the screen */
  
    status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
    if ( status ) {
      fprintf (stderr,
               "Failure to turn on screen indicator, error %d.\n", status);
      exit(1);
    }

    //necessary when using own cut generator (or heuristic??) with CPLEX
    status = CPXsetintparam (env, CPX_PARAM_PREIND, CPX_OFF);
    
    //set problem data
    
    lp = CPXcreateprob (env, &status, "curv-lp");

    /* A returned pointer of NULL may mean that not enough memory
       was available or there was some other problem.  In the case of
       failure, an error message will have been written to the error
       channel from inside CPLEX.  In this example, the setting of
       the parameter CPX_PARAM_SCRIND causes the error message to
       appear on stdout.  */
  
    if ( lp == NULL ) {
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
    for (uint c=0; c < row_start.size()-1; c++)
      row_count[c] = row_start[c+1] - row_start[c];
    
    status = CPXnewcols (env, lp, nVars, cost.direct_access(), var_lb.direct_access(), 
                         var_ub.direct_access(), NULL, NULL);
    
    std::cerr << "copying" << std::endl;
    // status = CPXcopylp (env, lp, 0, nConstraints, CPX_MIN, cost.direct_access(), 
    //  		      rhs_upper.direct_access(), row_sense, 
    //  		      (int*) row_start.direct_access(), row_count, (int*) lp_descr.col_indices(), lp_descr.value(),
    //  		      var_lb.direct_access(), var_ub.direct_access(), NULL);
    
    if ( status )  
      exit(1);
    
    std::cerr << "adding rows" << std::endl;
    
    CPXaddrows(env, lp, 0, row_start.size()-1, lp_descr.nEntries(), rhs_lower.direct_access(), row_sense, 
               (int*) row_start.direct_access(), (int*) lp_descr.col_indices(), lp_descr.value(),
               NULL, NULL);
    
    std::cerr << "done adding rows" << std::endl;
    
    uint count = 0;
    Math1D::Vector<int> idx(nConstraints);
    Math1D::Vector<double> range(nConstraints);
    for (int c=0; c < (int) row_start.size()-1; c++) {
      
      if (row_sense[c] == 'R') {
	
        idx[count] = c;
        range[count] = rhs_upper[c] - rhs_lower[c];
        count++;
      }
    }
    CPXchgrngval (env, lp, count, idx.direct_access(), range.direct_access());
    
    delete[] row_sense;
    delete[] row_count;
    
    std::cerr << "calling optimize" << std::endl;
  
    //status = CPXlpopt (env, lp);
    status = CPXbaropt (env, lp);
    
    if ( status ) {
      fprintf (stderr, "Failed to optimize MIP.\n");
      exit(1);
    }


    CPXsolution (env, lp, NULL, NULL, lp_solution.direct_access(), NULL, NULL, NULL);
  
    CPXfreeprob (env, &lp);
    CPXcloseCPLEX (&env);
  }
#endif

#ifdef HAS_XPRESS


  if (solver == "xpress") {

    solver_known = true;

    int nReturn;
    XPRSprob xp_prob;
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
  
    SparseMatrixDescription<double> lp_copy(lp_descr);

    Math1D::Vector<uint> col_start(nVars+1);
    lp_copy.sort_by_column(col_start);
    
    for (uint k=0; k < lp_copy.nEntries(); k++) {
      assert(lp_copy.row_indices()[k] < nConstraints);
    }
    
    char* row_sense = new char[nConstraints];
    double* row_range = new double[nConstraints];
    
    for (uint c=0; c < row_start.size()-1; c++) {
      
      if (rhs_lower[c] == rhs_upper[c]) {
        row_sense[c] = 'E';
      }
      else {
        row_sense[c] = 'R';
        row_range[c] = rhs_upper[c] - rhs_lower[c];
      }
    }
    
    XPRSsetintcontrol(xp_prob,XPRS_CROSSOVER,0);

    nReturn = XPRSloadlp(xp_prob, "curvdenoise-lp", col_start.size()-1, row_start.size()-1, row_sense,
                         rhs_upper.direct_access(), row_range, cost.direct_access(), 
                         (int*) col_start.direct_access(), NULL, (int*) lp_copy.row_indices(), lp_copy.value(),
                         var_lb.direct_access(), var_ub.direct_access());
 
    delete[] row_sense;
    delete[] row_range;

    XPRSlpoptimize(xp_prob,"b");  //barrier
    //XPRSlpoptimize(xp_prob,"d");  //dual simplex
    
    XPRSgetlpsol(xp_prob, lp_solution.direct_access(), 0, 0, 0); 
    
    nReturn=XPRSdestroyprob(xp_prob);
    nReturn=XPRSfree();
  }
#endif

  if (!solver_known && solver != "clp") {

    std::cerr << "WARNING: solver\"" << solver << "\" unknown. Taking Clp instead" << std::endl;
    solver = "clp";
  }
    
  if ( solver == "clp") {

    CoinPackedMatrix coinMatrix(false,(int*) lp_descr.row_indices(),(int*) lp_descr.col_indices(),
                                lp_descr.value(),lp_descr.nEntries());
    

    OsiClpSolverInterface lpSolver;
    lpSolver.loadProblem(coinMatrix, var_lb.direct_access(), var_ub.direct_access(),   
                         cost.direct_access(), rhs_lower.direct_access(), rhs_upper.direct_access());
    
    coinMatrix.cleanMatrix();

    std::clock_t tStartCLP,tEndCLP;
    
    tStartCLP = std::clock();
    
    //lpSolver.resolve();
    
    ClpSolve solve_options;
    solve_options.setSolveType(ClpSolve::useDual);
    //solve_options.setSolveType(ClpSolve::useBarrier);
    solve_options.setPresolveType(ClpSolve::presolveNumber,5);
    lpSolver.setSolveOptions(solve_options);
    lpSolver.initialSolve();
    
    tEndCLP = std::clock();
    
    error = 1 - lpSolver.isProvenOptimal();
    
    if (error != 0)
      std::cerr << "!!!!!!!!!!!!!!LP-solver failed!!!!!!!!!!!!!!!!!!!" << std::endl;
    
    std::cerr << "CLP-time: " << diff_seconds(tEndCLP,tStartCLP) << " seconds. " << std::endl;

    uint nSolverVars = lpSolver.getNumCols();

    for (uint v=0; v < nSolverVars; v++) {

      lp_solution[v] = lpSolver.getColSolution()[v];
    }

    assert(nSolverVars == nVars);
  }
  
  denoised_image.resize(xDim,yDim,zDim);
  denoised_image.set_constant(im_const);

  for (uint y=0; y < yDim; y++) {
    for (uint x=0; x < xDim; x++) {

      uint i = y*xDim+x;
    
      for (uint k=share_start[i]; k < share_start[i+1]; k++) {
	
        uint face_idx = shares[k].face_idx_;
        double fac = shares[k].share_;

        double lp_val = 0.0;
        for (uint b=0; b < nBins; b++) {
          lp_val += lp_solution[b*mesh.nFaces() + face_idx];
        }
	
        for (uint z=0; z < zDim; z++)
          denoised_image(x,y,z) += lp_val * fac * mesh.convex_area(face_idx);
      }
    }
  }
}
