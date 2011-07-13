/**** written by Thomas Schoenemann as an employee of Lund University, Sweden, 2010 *****/

#include "lp_segmenter.hh"
#include "conv_lp_solving.hh"

LpSegmenter::LpSegmenter(const Math2D::Matrix<float>& image, const  LPSegOptions& options, uint nRegions, bool shared_boundary_vars) 
  : image_(image), options_(options), nRegions_(nRegions), shared_boundary_vars_(shared_boundary_vars), lp_descr_(0) {
  
  xDim_ = uint( image_.xDim() );
  yDim_ = uint( image_.yDim() );

  solver_ = options.solver_;

  bool solver_known = false;

  if (solver_ == "convex-solver" || solver_ == "own-conv") {
    solver_ = "convex-solver";
    solver_known = true;
  }
  if (solver_ == "clp")
    solver_known = true;

#ifdef HAS_GUROBI
  if (solver_ == "gurobi")
    solver_known = true;
  grb_env_ = 0;
  grb_model_ = 0;
#endif
#ifdef HAS_XPRESS
  if (solver_ == "xpress")
    solver_known = true;
  xp_prob_ = 0;
#endif

  if (!solver_known) {
    std::cerr << "WARNING: solver unknown. Taking Clp instead." << std::endl;
    solver_ = "clp";
  }

  if (shared_boundary_vars_) {

    std::cerr << "ERROR: there is currently no correct implementation of shared boundaries. It is not clear if one exists." << std::endl;
    std::cerr << "exiting...." << std::endl;
    exit(1);
  }

  /**** initialize data term ****/

  data_term_.resize(xDim_,yDim_,nRegions_);
  mean_.resize(nRegions_);

  Math1D::Vector<float> intensity(image.size());
  for (uint i=0; i < image.size(); i++)
    intensity[i] = image.direct_access(i);
  std::sort(intensity.direct_access(),intensity.direct_access()+image.size());

  //spread percentiles equally
  double share = 1.0 / nRegions;
  for (uint i=0; i < nRegions_; i++) {
    mean_[i] = intensity[std::size_t(image.size() * (i+0.5)*share)];
  }
  
  for (uint y=0; y < yDim_; y++) {
    for (uint x=0; x < xDim_; x++) {
      float cur = image_(x,y);
      for (uint r=0; r < nRegions_; r++) {
	float data = float((cur-mean_[r])*(cur-mean_[r]));
	data_term_(x,y,r) = data;
      }
    }
  }

  uint out_factor = options.output_factor_;
  segmentation_.resize(xDim_*out_factor,yDim_*out_factor);

  double lambda = options.lambda_;
  double gamma = options.gamma_;
  int neighborhood = options.neighborhood_;
  bool enforce_consistent_boundaries = options.enforce_consistent_boundaries_;
  bool enforce_regionedge = options.enforce_regionedge_;
  bool bruckstein = options.bruckstein_;
  bool light_constraints = options_.light_constraints_;

  uint light_factor = (light_constraints) ? 1 : 2;

  if (options.gridtype_ == options.Square) {
    if (options.adaptive_mesh_n_ < 0) {
      double xfac = double(xDim_) / double(options.griddim_xDim_);
      double yfac = double(yDim_) / double(options.griddim_yDim_);
      generate_mesh( options_.griddim_xDim_, options.griddim_yDim_, neighborhood, mesh_, false, 0);
      mesh_.enlarge(xfac,yfac);
    }
    else {
      //Adaptive mesh
      TODO("combination of curvature potts model and adaptive meshes");
      //generate_adaptive_mesh(data_term, mesh, neighborhood, options.adaptive_mesh_n_);
    }
  }
  else {

    if (options.adaptive_mesh_n_ < 0) {
      double xfac = double(xDim_) / double(options.griddim_xDim_);
      double yfac = double(yDim_) / double(options.griddim_yDim_);
      
      generate_hexagonal_mesh( xDim_, yDim_, 0.5*(xfac+yfac), neighborhood, mesh_); //TODO: proper handling of the factor
    }
    else {
      //Adaptive mesh
      TODO("combination of curvature potts model and adaptive hexagonal meshes");
      //generate_adaptive_hexagonal_mesh(data_term, mesh, neighborhood, options.adaptive_mesh_n_);
    }
  }

  compute_pixel_shares(mesh_, xDim_, yDim_, shares_, share_start_);
  mesh_.generate_edge_pair_list(edge_pairs_);

  std::cerr << edge_pairs_.size() << " edge pairs." << std::endl;

  nVars_ = nRegions_*mesh_.nFaces();
  if (!shared_boundary_vars_)
    nVars_ += uint( 2*nRegions_*edge_pairs_.size() );
  else
    nVars_ += uint( 2*edge_pairs_.size() ); 

  nSlacks_ = 0;

  nConstraints_ = mesh_.nFaces() //simplex constraints  
    + nRegions_*mesh_.nEdges(); //surface continuation constraints
  nConstraints_ += nRegions_*2*mesh_.nEdges(); //boundary continuation constraints
  
  uint surface_con_offs = mesh_.nFaces();
  uint boundary_con_offset = surface_con_offs + nRegions_*mesh_.nEdges();

  const uint consistency_con_offs = nConstraints_;
  if (enforce_consistent_boundaries) {
    nConstraints_ += nRegions_*light_factor*mesh_.nEdges();
  }
  const uint regionedge_constraints_offs = nConstraints_;
  if (enforce_regionedge) {
    nConstraints_ += nRegions_*2*light_factor*mesh_.nEdges();
  }

  var_lb_.resize(nVars_,0.0);
  var_ub_.resize(nVars_,1.0);
  cost_.resize(nVars_,0.0);
  
  for (uint i=0; i<mesh_.nFaces(); ++i) {
    for (uint r = 0; r < nRegions_; r++) 
      cost_[i*nRegions+r] = calculate_data_term(i, r, mesh_, data_term_);
  }

  uint edge_pair_var_offs = nRegions_*mesh_.nFaces();

  for (uint j=0; j < edge_pairs_.size(); j++) {
    
    uint first = edge_pairs_[j].first_edge_idx_;
    uint nFirstAdjacent = uint( mesh_.adjacent_faces(first).size() );
    
    uint second = edge_pairs_[j].second_edge_idx_;
    uint nSecondAdjacent = uint( mesh_.adjacent_faces(second).size() );

    double weight = 0.0;
    if (nFirstAdjacent > 1)
      weight += 0.5*lambda*mesh_.edge_length(first);
    if (nSecondAdjacent > 1)
      weight += 0.5*lambda*mesh_.edge_length(second);
    
    //do not penalize the image corners for their curvature
    if (nFirstAdjacent > 1 || nSecondAdjacent > 1)
      weight += gamma * curv_weight(mesh_,edge_pairs_[j],2.0,bruckstein);
    
    uint nRelevantRegions = nRegions_; //might be different for shared vars

    //0.5 as each region border is counted twice
    //TODO: this is not the correct factor at the image border (but in the end it is just a constant offset)
    double boundary_fac = 0.5;

    for (uint r = 0; r < nRelevantRegions; r++) {

      cost_[edge_pair_var_offs + 2*r*edge_pairs_.size() + 2*j]   = boundary_fac*weight;
      cost_[edge_pair_var_offs + 2*r*edge_pairs_.size() + 2*j+1] = boundary_fac*weight;
      
      /*** check if (at the image border) one of the edge pairs is impossible. if so, set its upper bound to 0 ***/
      uint edge = first;
      if (mesh_.adjacent_faces(edge).size() == 1) {
	int match = mesh_.match(mesh_.adjacent_faces(edge)[0],edge);
	
	uint y = uint( edge_pair_var_offs+2*r*edge_pairs_.size()+2*j );
	
	if (edge_pairs_[j].common_point_idx_ == mesh_.edge(edge).from_idx_)
	  match *= -1;
	
	if (match == -1)
	  var_ub_[y+1] = 0.0;
	else if (match == 1)
	  var_ub_[y] = 0.0;
      }
    
      edge = second;
      if (mesh_.adjacent_faces(edge).size() == 1) {
	int match = mesh_.match(mesh_.adjacent_faces(edge)[0],edge);
	
	uint y = uint( edge_pair_var_offs + 2*r*edge_pairs_.size() + 2*j );
	
	if (edge_pairs_[j].common_point_idx_ == mesh_.edge(edge).from_idx_)
	  match *= -1;
	
	if (match == -1)
	  var_ub_[y] = 0.0;
	else if (match == 1)
	  var_ub_[y+1] = 0.0;
      }
    }
  }

  rhs_lower_.resize(nConstraints_,0.0);
  rhs_upper_.resize(nConstraints_,0.0);

  //prepare the simplex constraints
  for (uint c=0; c < mesh_.nFaces(); c++) {
    rhs_lower_[c] = 1.0;
    rhs_upper_[c] = 1.0;
  }

  std::cerr << "coding matrix" << std::endl;

  nEntries_ = uint( nRegions_*(mesh_.nFaces() + 2*mesh_.nEdges() + 2*edge_pairs_.size()) + 4*edge_pairs_.size() );  
  nEntries_ += uint( (nRegions_ - 1)*4*edge_pairs_.size() );

  if (enforce_consistent_boundaries) {
    nEntries_ += uint( 2*light_factor*nRegions_*edge_pairs_.size() 
		       + light_factor*nRegions_*mesh_.nEdges()); //we also allocate space for the slack variables used with the convex solver
  }
  if (enforce_regionedge) {
    nEntries_ += uint( nRegions_*(8*light_factor*edge_pairs_.size() 
			  + 2*light_factor*mesh_.nEdges()) ); //we also allocate space for the slack variables used with the convex solver
  }

  SparseMatrixDescription<double> lp_descr(nEntries_, nConstraints_, nVars_);

  /*** a) region simplex constraints ****/
  for (uint i=0; i < mesh_.nFaces(); i++) {

    for (uint r=0; r < nRegions_; r++)
      lp_descr.add_entry(i, i*nRegions_ + r, 1.0);
  }

  uint nSimplexConEntries = mesh_.nFaces()*nRegions_;
  assert(nSimplexConEntries == lp_descr.nEntries());

  /*** b) surface continuation constraints ****/
  for (uint r=0; r < nRegions_; r++) {
    uint cur_con_offs = surface_con_offs+r*mesh_.nEdges();
    
    for (uint j=0; j < mesh_.nEdges(); j++) {
      
      const std::vector<uint>& adjacent_faces = mesh_.adjacent_faces(j);
      for (std::vector<uint>::const_iterator it = adjacent_faces.begin();
	   it != adjacent_faces.end(); it++) {
	
	lp_descr.add_entry(cur_con_offs + j, nRegions_*(*it)+r, mesh_.match(*it,j));
      }
    }
    
    uint cur_edge_var_offs = edge_pair_var_offs;
    cur_edge_var_offs += uint( 2*r*edge_pairs_.size() );

    for (uint j=0; j < edge_pairs_.size(); j++) {
      
      uint first_edge = edge_pairs_[j].first_edge_idx_;
      uint second_edge = edge_pairs_[j].second_edge_idx_;

      double first_surf_bound_coeff  = 1.0;
      double second_surf_bound_coeff = 1.0;
      
      uint middle_point = edge_pairs_[j].common_point_idx_;
      
      if (mesh_.edge(first_edge).to_idx_ == middle_point) {
	lp_descr.add_entry(cur_con_offs+first_edge,cur_edge_var_offs+2*j, first_surf_bound_coeff);
      }
      else {
	lp_descr.add_entry(cur_con_offs+first_edge,cur_edge_var_offs+2*j, -first_surf_bound_coeff);
      }
      
      if (mesh_.edge(second_edge).to_idx_ == middle_point) {
	lp_descr.add_entry(cur_con_offs+second_edge,cur_edge_var_offs+2*j+1, second_surf_bound_coeff);
      }
      else {
	lp_descr.add_entry(cur_con_offs+second_edge,cur_edge_var_offs+2*j+1, -second_surf_bound_coeff);
      }
    }
  }

  /*** c) boundary continuation constraints ****/
  for (uint j=0; j < edge_pairs_.size(); j++) {
    
    uint first_edge = edge_pairs_[j].first_edge_idx_;
    uint second_edge = edge_pairs_[j].second_edge_idx_;
    
    uint middle_point = edge_pairs_[j].common_point_idx_;

    uint nRelevantRegions = nRegions_; //might be different for shared boundary vars
    
    for (uint r=0; r < nRelevantRegions; r++) {

      uint cur_edge_offs = uint( edge_pair_var_offs + 2*r*edge_pairs_.size() );

      if (mesh_.edge(first_edge).to_idx_ == middle_point) {      
	lp_descr.add_entry(boundary_con_offset + 2*(nRelevantRegions*first_edge+r), cur_edge_offs+2*j, 1);
	lp_descr.add_entry(boundary_con_offset + 2*(nRelevantRegions*first_edge+r)+1, cur_edge_offs+2*j+1, -1);
      }
      else {
	lp_descr.add_entry(boundary_con_offset + 2*(nRelevantRegions*first_edge+r)+1, cur_edge_offs+2*j, 1);
	lp_descr.add_entry(boundary_con_offset + 2*(nRelevantRegions*first_edge+r), cur_edge_offs+2*j+1, -1);
      }
      
      if (mesh_.edge(second_edge).from_idx_ == middle_point) {
	lp_descr.add_entry(boundary_con_offset + 2*(nRelevantRegions*second_edge+r), cur_edge_offs+2*j, -1);
	lp_descr.add_entry(boundary_con_offset + 2*(nRelevantRegions*second_edge+r)+1, cur_edge_offs+2*j+1, 1);
      }
      else {
	lp_descr.add_entry(boundary_con_offset + 2*(nRelevantRegions*second_edge+r)+1, cur_edge_offs+2*j, -1);
	lp_descr.add_entry(boundary_con_offset + 2*(nRelevantRegions*second_edge+r), cur_edge_offs+2*j+1, 1);
      }
    }
  }

  uint nStandardEntries = lp_descr.nEntries();

  if (enforce_consistent_boundaries) {

    assert(!shared_boundary_vars_);

    //constraints in words: for each oriented edge, the pairs that start with this oriented edge 
    // and the pairs that end in the oppositely oriented edge may not sum to more than 1.0
    // (i.e. they are mutually exclusive)

    for (uint c=consistency_con_offs; c < consistency_con_offs + light_factor*nRegions_*mesh_.nEdges(); c+= light_factor) {
      rhs_upper_[c] = 1.0;
      if (!light_constraints)
	rhs_upper_[c+1] = 1.0;
    }

    for (uint j=0; j < edge_pairs_.size(); j++) {

      uint middle_point = edge_pairs_[j].common_point_idx_;
      
      uint first_edge = edge_pairs_[j].first_edge_idx_;
      uint second_edge = edge_pairs_[j].second_edge_idx_;

      for (uint r = 0; r < nRegions_; r++) {

	uint cur_con_offs = consistency_con_offs + light_factor*r*mesh_.nEdges();
	uint cur_edge_offs = uint( edge_pair_var_offs + 2*r*edge_pairs_.size() );

	if (mesh_.edge(first_edge).to_idx_ == middle_point) {      
	  lp_descr.add_entry(cur_con_offs + light_factor*first_edge, cur_edge_offs+2*j, 1);
	  lp_descr.add_entry(cur_con_offs + light_factor*first_edge, cur_edge_offs+2*j+1, 1);
	}
	else if (!light_constraints) {
	  lp_descr.add_entry(cur_con_offs + light_factor*first_edge+1, cur_edge_offs+2*j, 1);
	  lp_descr.add_entry(cur_con_offs + light_factor*first_edge+1, cur_edge_offs+2*j+1, 1);
	}
	
	if (mesh_.edge(second_edge).to_idx_ == middle_point) {
	  lp_descr.add_entry(cur_con_offs + light_factor*second_edge, cur_edge_offs+2*j, 1);
	  lp_descr.add_entry(cur_con_offs + light_factor*second_edge, cur_edge_offs+2*j+1, 1);
	}
	else if (!light_constraints) {
	  lp_descr.add_entry(cur_con_offs + light_factor*second_edge+1, cur_edge_offs+2*j, 1);
	  lp_descr.add_entry(cur_con_offs + light_factor*second_edge+1, cur_edge_offs+2*j+1, 1);
	}
      }
    }
  }
  if (enforce_regionedge) {

    uint rowoff = regionedge_constraints_offs;

    for (uint edge=0; edge < mesh_.nEdges(); edge++) {
      
      //Get the two adjacent faces
      if (mesh_.adjacent_faces(edge).size() != 2) {
	  //One of the edges is at the border of the image
      
	continue;
      }
      
      uint x1 = mesh_.adjacent_faces(edge)[0];
      uint x2 = mesh_.adjacent_faces(edge)[1];

      for (uint r=0; r < nRegions_; r++) {

	uint cur_row = rowoff + 2*light_factor*r*mesh_.nEdges() + 2*light_factor*edge; 
	
	lp_descr.add_entry(cur_row , nRegions_*x1+r, 1);
	lp_descr.add_entry(cur_row , nRegions_*x2+r, 1);
	lp_descr.add_entry(cur_row+1 , nRegions_*x1+r, -1);
	lp_descr.add_entry(cur_row+1, nRegions_*x2+r, -1);
	if (!light_constraints) {
	  lp_descr.add_entry(cur_row+2, nRegions_*x1+r, 1);
	  lp_descr.add_entry(cur_row+2, nRegions_*x2+r, 1);
	  lp_descr.add_entry(cur_row+3, nRegions_*x1+r, -1);
	  lp_descr.add_entry(cur_row+3, nRegions_*x2+r, -1);
	}

	//NOTE (important for the construction with slacks: the binding constraints are in both cases the upper bounds)
	rhs_lower_[cur_row]   = 0;
	rhs_upper_[cur_row]   = 2;
	rhs_lower_[cur_row+1] = -2;
	rhs_upper_[cur_row+1] = 0;
	
	if (!light_constraints) {
	  rhs_lower_[cur_row+2] = 0;
	  rhs_upper_[cur_row+2] = 2;
	  rhs_lower_[cur_row+3] = -2;
	  rhs_upper_[cur_row+3] = 0;
	}
      }
    }

    for (uint r = 0; r < nRegions_; r++) {

      uint cur_edge_offs = uint( edge_pair_var_offs + 2*r*edge_pairs_.size() );

      for (uint j=0; j < edge_pairs_.size(); j++) {
	uint first = edge_pairs_[j].first_edge_idx_;
	uint second = edge_pairs_[j].second_edge_idx_;	
	
	uint y = cur_edge_offs + 2*j;
	
	uint edge = first;
	if (mesh_.adjacent_faces(edge).size() == 2) {
	  
	  uint cur_row = rowoff + 2*light_factor*r*mesh_.nEdges() + 2*light_factor*edge; 

	  lp_descr.add_entry(cur_row   ,y  , 1);
	  lp_descr.add_entry(cur_row+1 ,y  , 1);
	  if (!light_constraints) {
	    lp_descr.add_entry(cur_row+2 ,y+1, 1);
	    lp_descr.add_entry(cur_row+3 ,y+1, 1);
	  }
	}
      
	edge = second;
	if (mesh_.adjacent_faces(edge).size() == 2) {
	
	  uint cur_row = rowoff + 2*light_factor*r*mesh_.nEdges() + 2*light_factor*edge; 

	  lp_descr.add_entry(cur_row   ,y+1, 1);
	  lp_descr.add_entry(cur_row+1 ,y+1, 1);
	  if (!light_constraints) {
	    lp_descr.add_entry(cur_row+2 ,y  , 1);
	    lp_descr.add_entry(cur_row+3 ,y  , 1);
	  }
	}
      }
    }
  }

  std::cerr << "sorting" << std::endl;

  Math1D::Vector<uint> row_start(nConstraints_+1);
  lp_descr.sort_by_row(row_start);

#ifdef HAS_GUROBI
  int error;
  grb_env_ = 0;
  grb_model_ = 0;

  if (solver_ == "gurobi") {

    /* Create environment */
    
    error = GRBloadenv(&grb_env_,NULL);
    //GRBsetintparam(grb_env_, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
    GRBsetintparam(grb_env_, GRB_INT_PAR_METHOD, GRB_METHOD_BARRIER);
    //GRBsetdblparam(grb_env_, "BarConvTol", 1e-10);
    GRBsetintparam(grb_env_, "Crossover", 0);
    //GRBsetintparam(grb_env_, "CrossoverBasis", 1);
    GRBsetintparam(grb_env_, "Presolve", 0);
    GRBsetintparam(grb_env_, "PrePasses", 1);
    
    assert (!error && grb_env_ != NULL);

    /* Create an empty model */

    error = GRBnewmodel(grb_env_, &grb_model_, "potts-curv-lp", 0, NULL, NULL, NULL, NULL, NULL);
    assert(!error);

    Storage1D<char> vtype(nVars_,GRB_CONTINUOUS);
  
    error = GRBaddvars(grb_model_,nVars_,0,NULL,NULL,NULL,cost_.direct_access(),var_lb_.direct_access(),
		       var_ub_.direct_access(),vtype.direct_access(),NULL);
    assert(!error);

    error = GRBupdatemodel(grb_model_);
    assert(!error);
    
    for (uint c=0; c < nConstraints_; c++) {
      
      if (rhs_lower_[c] == rhs_upper_[c])
	error = GRBaddconstr(grb_model_, row_start[c+1]-row_start[c], ((int*) lp_descr.col_indices()) + row_start[c], 
			     lp_descr.value() + row_start[c], GRB_EQUAL, rhs_upper_[c], NULL);
      else
	GRBaddrangeconstr(grb_model_, row_start[c+1]-row_start[c], ((int*) lp_descr.col_indices()) + row_start[c], 
			  lp_descr.value() + row_start[c], rhs_lower_[c], rhs_upper_[c], NULL);      
      
      assert(!error);
    }    
  }
#endif
#ifdef HAS_XPRESS
  if (solver_ == "xpress") {

    char banner[256];
    
    int nReturn=XPRSinit("/opt/xpressmp/");
    
    if (nReturn != 0) {
      
      char msg[512];
      XPRSgetlicerrmsg(msg,512);
      
      std::cerr << "error message: " << msg << std::endl;
    }
    
    assert(nReturn == 0);
    
    XPRSgetbanner(banner); printf("banner: %s \n",banner);
    
    nReturn=XPRScreateprob(&xp_prob_);

    Math1D::Vector<uint> col_start(nVars_+1);
    lp_descr.sort_by_column(col_start);
    
    char* row_sense = new char[nConstraints_];
    double* row_range = new double[nConstraints_];

    for (uint c=0; c < nConstraints_; c++) {
      
      if (rhs_lower_[c] == rhs_upper_[c]) {
	row_sense[c] = 'E';
      }
      else {
	row_sense[c] = 'R';
	row_range[c] = rhs_upper_[c] - rhs_lower_[c];
      }
    }
    
    XPRSsetintcontrol(xp_prob_,XPRS_CROSSOVER,0);

    nReturn = XPRSloadlp(xp_prob_, "mcurvseg-ilp", nVars_, nConstraints_, row_sense,
			 rhs_upper_.direct_access(), row_range, cost_.direct_access(), 
			 (int*) col_start.direct_access(), NULL, (int*) lp_descr.row_indices(), lp_descr.value(),
			 var_lb_.direct_access(), var_ub_.direct_access());
    
    delete[] row_sense;
    delete[] row_range;
  }
#endif
  
  if (solver_ == "clp") {
    coinMatrix_ = CoinPackedMatrix (false,(int*) lp_descr.row_indices(),(int*) lp_descr.col_indices(),
				    lp_descr.value(),lp_descr.nEntries());

    lpSolver_.loadProblem (coinMatrix_, var_lb_.direct_access(), var_ub_.direct_access(),   
			   cost_.direct_access(), rhs_lower_.direct_access(), rhs_upper_.direct_access());
    //std::cerr << "setting precision" << std::endl;
    //lpSolver_.messageHandler()->setPrecision(12);
    //std::cerr << "precision set" << std::endl;

    std::cerr << "#estimated vars: " << nVars_ << std::endl;
    std::cerr << "solver knows " << lpSolver_.getNumCols() << " variables." << std::endl;
  }
  if (solver_ == "convex-solver") {

    simplex_start_.resize(mesh_.nFaces()+1,0);
    for (uint i=0; i <= mesh_.nFaces(); i++)
      simplex_start_[i] = i*nRegions_;      
    
    lp_descr_.reset(lp_descr.nReservedEntries()-nSimplexConEntries, nConstraints_-mesh_.nFaces(), nVars_+nSlacks_);

    for (uint n=nSimplexConEntries; n < lp_descr.nEntries(); n++) {
      
      lp_descr_.add_entry(lp_descr.row_indices()[n]-mesh_.nFaces(), lp_descr.col_indices()[n], (char) lp_descr.value()[n]);
    }

    if (enforce_consistent_boundaries) {
      
      uint next_slack_var = nVars_ + nSlacks_;
      nSlacks_ += 2*nRegions_*mesh_.nEdges();

      lp_descr_.increase_nColumns(2*nRegions_*mesh_.nEdges());

      var_lb_.resize(nVars_ + nSlacks_, 0.0);
      var_ub_.resize(nVars_ + nSlacks_, 1.0);
      cost_.resize(nVars_ + nSlacks_, 0.0);

      for (uint c=consistency_con_offs; c < consistency_con_offs + light_factor*nRegions_*mesh_.nEdges(); c++) {
	lp_descr_.add_entry(c - mesh_.nFaces(), next_slack_var, 1);
	next_slack_var++;
      }
    }

    //make sure this matrix is sorted (important for the cuda-based solver)
    Math1D::Vector<uint> row_start(nConstraints_-mesh_.nFaces()-1,MAX_UINT);
    lp_descr_.sort_by_row(row_start);
  }

  output_.resize(xDim_*out_factor,yDim_*out_factor,nRegions_,0.0);
  mesh_.enlarge(out_factor,out_factor);
}

LpSegmenter::~LpSegmenter() {

#ifdef HAS_GUROBI
  if (grb_model_ != 0)
    GRBfreemodel(grb_model_);
  if (grb_env_ != 0)
    GRBfreeenv(grb_env_);
#endif
#ifdef HAS_XPRESS
  if (xp_prob_ != 0) {
    int nReturn=XPRSdestroyprob(xp_prob_);
    nReturn=XPRSfree();
  }
#endif
}

const Math2D::Matrix<uint>& LpSegmenter::segmentation() {

  return segmentation_;
}


double LpSegmenter::segment(uint nIter) {

  int error = 0;
  double energy = 0.0;

  Math1D::Vector<double> solver_solution;

  if (solver_ != "clp")
    solver_solution.resize(nVars_+nSlacks_,0.0);

  for (uint iter=1; iter <= nIter; iter++) {

    std::cerr << "******************** iteration " << iter << std::endl;

    std::cerr << "mean values: " << mean_ << std::endl;

    const double* lp_solution = 0;

#ifdef HAS_GUROBI
    if (solver_ == "gurobi") {
    
      //to be sure we reset the variable bounds of the region variables
      for (uint v=0; v < nVars_; v++) {
	GRBsetdblattrelement(grb_model_,"LB",v, var_lb_.direct_access(v));
	GRBsetdblattrelement(grb_model_,"UB",v, var_ub_.direct_access(v));
      }

      error = GRBupdatemodel(grb_model_);
      error = GRBoptimize(grb_model_);
      assert(!error);

      solver_solution.resize(nVars_);
      for (uint v=0; v < nVars_; v++) {
	GRBgetdblattrelement(grb_model_,"X",v, solver_solution.direct_access()+v);
      }
      lp_solution = solver_solution.direct_access();
    }
#endif
#ifdef HAS_XPRESS
    if (solver_ == "xpress") {
      int nReturn;

      Math1D::Vector<int> var_idx(nVars_,0);
      for (uint v=0; v < nVars_; v++) {
	var_idx[v] = v;
      }
      Math1D::Vector<char> btype(nVars_,'L');
      XPRSchgbounds(xp_prob_,nVars_, var_idx.direct_access(), btype.direct_access(), var_lb_.direct_access() );
      btype.set_constant('U');
      XPRSchgbounds(xp_prob_,nVars_, var_idx.direct_access(), btype.direct_access(), var_ub_.direct_access() );

      nReturn = XPRSlpoptimize(xp_prob_,"b");  //barrier

      XPRSgetlpsol(xp_prob_, solver_solution.direct_access(), 0, 0, 0); 
      lp_solution = solver_solution.direct_access();      
    }
#endif

    if (solver_ == "clp") {
      //to be sure we reset the variables bounds of the region variables
      lpSolver_.setColLower(var_lb_.direct_access());
      lpSolver_.setColUpper(var_ub_.direct_access());
      
      std::clock_t tStartCLP,tEndCLP;
      tStartCLP = std::clock();
      
      lpSolver_.resolve();
      
      tEndCLP = std::clock();
      
      lp_solution = lpSolver_.getColSolution();

      //ClpSolve solve_options;
      //solve_options.setSolveType(ClpSolve::useDual);
      //solve_options.setSolveType(ClpSolve::useBarrier);
      //solve_options.setPresolveType(ClpSolve::presolveNumber,5);
      //lpSolver.initialSolve(solve_options);
      
      error = 1 - lpSolver_.isProvenOptimal();

      tEndCLP = std::clock();
    
      if (error != 0)
	std::cerr << "!!!!!!!!!!!!!!LP-solver failed!!!!!!!!!!!!!!!!!!!" << std::endl;
      
      std::cerr << "CLP-time: " << diff_seconds(tEndCLP,tStartCLP) << " seconds. " << std::endl;
    }    
    if (solver_ == "convex-solver") {

#if 1
      if (iter == 1) {

	for (uint i=0; i < mesh_.nFaces(); ++i) {
	  double best_cost = 1e300;
	  uint best_region = MAX_UINT;
	  for (uint r = 0; r < nRegions_; r++) {
	    if (cost_[i*nRegions_+r] < best_cost) {
	      best_cost = cost_[i*nRegions_+r];
	      best_region = r;
	    }
	  }
	  solver_solution[i*nRegions_+best_region] = 1.0;
	}

	uint edge_pair_var_offs = nRegions_*mesh_.nFaces();

	for (uint r = 0; r < nRegions_; r++) {

	  uint cur_edge_offs = uint( edge_pair_var_offs + 2*r*edge_pairs_.size() );

	  for (uint p=0; p < edge_pairs_.size(); p++) {
	    
	    uint edge1 = edge_pairs_[p].first_edge_idx_;
	    uint edge2 = edge_pairs_[p].second_edge_idx_;
	    
	    double drop1 = 0.0;
	    for (uint i=0; i < mesh_.adjacent_faces(edge1).size(); i++) {
	      uint face = mesh_.adjacent_faces(edge1)[i];
	      drop1 += solver_solution[face*nRegions_+r] * mesh_.match(face,edge1);      
	    }
	    
	    double drop2 = 0.0;
	    for (uint i=0; i < mesh_.adjacent_faces(edge2).size(); i++) {
	      uint face = mesh_.adjacent_faces(edge2)[i];
	      drop2 += solver_solution[face*nRegions_+r] * mesh_.match(face,edge2);      
	    }
	    
	    if (fabs(drop1) >= 0.1 && fabs(drop2) > 0.1) {
	      
	      if (mesh_.edge(edge1).from_idx_ == edge_pairs_[p].common_point_idx_) {
		
		drop1 *= -1.0;
		drop2 *= -1.0;
	      }
	    
	      if (drop1 > 0.0)
		solver_solution[cur_edge_offs + 2*p] = fabs(drop1);
	      else
		solver_solution[cur_edge_offs + 2*p+1] = fabs(drop1);
	    }
	  }
	}

	for (uint v=nVars_; v < nVars_+nSlacks_; v++) {
	  solver_solution[v] = 1.0;
	}
      }
      
#endif

      double stepsize_coeff = 0.5;
      if (nSlacks_ > 0)
	stepsize_coeff = 0.2;

      eq_and_simplex_constrained_lp_solve_auglagr_nesterov(nVars_+nSlacks_, nConstraints_-mesh_.nFaces(), cost_.direct_access(), 
							   var_lb_.direct_access(), var_ub_.direct_access(),
							   lp_descr_, rhs_upper_.direct_access()+mesh_.nFaces(), 
							   mesh_.nFaces(),simplex_start_.direct_access(),solver_solution.direct_access(),
							   options_.gamma_*50.0, stepsize_coeff, 2500, 30, 1.1);

      lp_solution = solver_solution.direct_access();
    }

    uint nFrac = 0;
    energy = 0.0;
    for (uint i=0; i < nVars_; i++) {
      double val = lp_solution[i];
      energy += cost_[i] * val;
      
      if (val > 0.01 && val < 0.99)
	nFrac++;
    }
  
    std::cerr << nFrac << "/" << nVars_ << " variables are fractional" << std::endl;

    //threshold and re-run the solver.
    if (iter == nIter) {

      //a) threshold
      for (uint i=0; i < mesh_.nFaces();++i) {

	uint arg_max = MAX_UINT;
	double max_val = 0.0;

	for (uint r=0; r < nRegions_; r++) {

	  double val = lp_solution[i*nRegions_ + r];
	  if (val > max_val) {

	    max_val = val;
	    arg_max = r;
	  }
	  var_ub_[i*nRegions_ + r] = 0.0;
	}
	
	var_ub_[i*nRegions_ + arg_max] = 1.0;
	var_lb_[i*nRegions_ + arg_max] = 1.0;
      }
      
      if (options_.enforce_consistent_boundaries_) {

	uint edge_pair_var_offs = nRegions_*mesh_.nFaces();

	for (uint p=0; p < edge_pairs_.size(); p++) {

	  uint edge1 = edge_pairs_[p].first_edge_idx_;
	  uint edge2 = edge_pairs_[p].second_edge_idx_;
	  
	  for (uint r=0; r < nRegions_; r++) {
	    
	    double drop1 = 0.0;
	    for (uint i=0; i < mesh_.adjacent_faces(edge1).size(); i++) {
	      uint face = mesh_.adjacent_faces(edge1)[i];
	      drop1 += var_ub_[face*nRegions_+r] * mesh_.match(face,edge1);      
	    }
	    
	    double drop2 = 0.0;
	    for (uint i=0; i < mesh_.adjacent_faces(edge2).size(); i++) {
	      uint face = mesh_.adjacent_faces(edge2)[i];
	      drop2 += var_ub_[face*nRegions_+r] * mesh_.match(face,edge2);      
	    }
	    
	    if (drop1 == 0.0 || drop2 == 0.0) {
	      var_ub_[edge_pair_var_offs + 2*r*edge_pairs_.size() + 2*p] = 0.0;
	      var_ub_[edge_pair_var_offs + 2*r*edge_pairs_.size() + 2*p+1] = 0.0;
	    }
	    else {
	      
	      if (mesh_.edge(edge1).from_idx_ == edge_pairs_[p].common_point_idx_) {
		
		drop1 *= -1.0;
		drop2 *= -1.0;
	      }
	      
	      if (drop1 < 0.0)
		var_ub_[edge_pair_var_offs + 2*r*edge_pairs_.size() + 2*p+1] = 0.0;
	      else
		var_ub_[edge_pair_var_offs + 2*r*edge_pairs_.size() + 2*p] = 0.0;
	    }
	  }
	}
      }

      std::cerr << "----------------------------computing energy of thresholded solution" << std::endl;

      if (solver_ == "clp") {

	//reset the variables bounds of the region variables
	lpSolver_.setColLower(var_lb_.direct_access());
	lpSolver_.setColUpper(var_ub_.direct_access());
	
	std::clock_t tStartCLP,tEndCLP;
	tStartCLP = std::clock();
	
	lpSolver_.resolve();
	
	tEndCLP = std::clock();
	
	lp_solution = lpSolver_.getColSolution();
	
	error = 1 - lpSolver_.isProvenOptimal();
	
	tEndCLP = std::clock();
	
	if (error != 0)
	  std::cerr << "!!!!!!!!!!!!!!LP-solver failed!!!!!!!!!!!!!!!!!!!" << std::endl;
	
	std::cerr << "CLP-time for thresholded solution: " << diff_seconds(tEndCLP,tStartCLP) << " seconds. " << std::endl;
      }    
#ifdef HAS_GUROBI
      if (solver_ == "gurobi") {
    
	//to be sure we reset the variables bounds of the region variables
	for (uint v=0; v < nVars_; v++) {
	  GRBsetdblattrelement(grb_model_,"LB",v, var_lb_.direct_access(v));
	  GRBsetdblattrelement(grb_model_,"UB",v, var_ub_.direct_access(v));
	}
	
	error = GRBresetmodel(grb_model_);
	error = GRBupdatemodel(grb_model_);
	error = GRBoptimize(grb_model_);
	assert(!error);
	
	solver_solution.resize(nVars_);
	for (uint v=0; v < nVars_; v++) {
	  GRBgetdblattrelement(grb_model_,"X",v, solver_solution.direct_access()+v);
	}
	lp_solution = solver_solution.direct_access();
      }
#endif

#ifdef HAS_XPRESS
      if (solver_ == "xpress") {

	Math1D::Vector<int> var_idx(nVars_,0);
	for (uint v=0; v < nVars_; v++) {
	  var_idx[v] = v;
	}
	Math1D::Vector<char> btype(nVars_,'L');
	XPRSchgbounds(xp_prob_,nVars_, var_idx.direct_access(), btype.direct_access(), var_lb_.direct_access() );
	btype.set_constant('U');
	XPRSchgbounds(xp_prob_,nVars_, var_idx.direct_access(), btype.direct_access(), var_ub_.direct_access() );
	
	int nReturn = XPRSlpoptimize(xp_prob_,"b");  //barrier

	XPRSgetlpsol(xp_prob_, solver_solution.direct_access(), 0, 0, 0); 
	lp_solution = solver_solution.direct_access();      
      }
#endif

      uint nFrac = 0;
      energy = 0.0;
      for (uint i=0; i < nVars_; i++) {
	double val = lp_solution[i];
	energy += cost_[i] * val;
	
	if (val > 0.01 && val < 0.99)
	  nFrac++;
      }
      
      std::cerr << "check: now " <<  nFrac << "/" << nVars_ << " variables are fractional" << std::endl;


      /*** reset variable bounds ***/
      //NOTE: this will not work if seed nodes are given
      for (uint v=0; v < nVars_; v++) {
	var_lb_[v] = 0.0;
	var_ub_[v] = 1.0;
      }
    }

    output_.set_constant(0.0);

    /*** extract solution ****/
    for (uint i=0; i < mesh_.nFaces();++i) {
      
      add_grid_output_mreg(i, lp_solution+i*nRegions_, nRegions_, mesh_, output_);
    }

    uint seg_factor = 255 / nRegions_;
    uint out_factor = options_.output_factor_;

    for (uint y=0; y < yDim_*out_factor; y++) {
      for (uint x=0; x < xDim_*out_factor; x++) {

	double max_val = 0.0;
	uint arg_max = MAX_UINT;
	
	for (uint r=0; r < nRegions_; r++) {
	  
	  double val = output_(x,y,r);
	  if (val > max_val) {
	    max_val = val;
	    arg_max = r;
	  }
	}
	
	segmentation_(x,y) = arg_max * seg_factor;
      }
    }

    mesh_.enlarge(1.0/out_factor,1.0/out_factor);
    
    /*** reestimate means and data terms ****/
    mean_.set_constant(0.0);
    Math1D::Vector<double> mean_support(nRegions_,0.0);

    for (uint y=0; y < yDim_; y++) {
      for (uint x=0; x < xDim_; x++) {

	double cur_intensity = image_(x,y);
	
	for (uint v=share_start_[y*xDim_+x]; v < share_start_[y*xDim_+x+1]; v++) {
	  uint f = shares_[v].face_idx_;
	  float area = std::min(1.0f, shares_[v].share_) * float(mesh_.convex_area(f));

	  for (uint r=0; r < nRegions_; r++) {

	    double val = lp_solution[f*nRegions_+r];
	    
	    mean_[r] += val*area*cur_intensity;
	    mean_support[r] += val*area;
	  }
	}
      }
    }

    for (uint r=0; r < nRegions_; r++) {

      if (mean_support[r] > 1e-50)
	mean_[r] /= mean_support[r];
    }

    std::cerr << "new mean values: " << mean_ << std::endl;
    
    //reset data term
    for (uint y=0; y < yDim_; y++) {
      for (uint x=0; x < xDim_; x++) {
	float cur = image_(x,y);
	for (uint r=0; r < nRegions_; r++) {
	  float data = float( (cur-mean_[r])*(cur-mean_[r]) );
	  data_term_(x,y,r) = data;
	  assert(data >= 0.0);
	}
      }
    }

    //reset cost function
    for (uint i=0; i < mesh_.nFaces(); ++i) {
      for (uint r = 0; r < nRegions_; r++) 
	cost_[i*nRegions_+r] = calculate_data_term(i, r, mesh_, data_term_);
    }

    mesh_.enlarge(out_factor,out_factor);

    for (uint v=0; v < nVars_; v++) {
      assert(cost_[v] >= 0.0);
    }

#ifdef HAS_GUROBI
    if (solver_ == "gurobi") {
      for (uint v=0; v < nVars_; v++) {
	GRBsetdblattrelement(grb_model_,"Obj",v, cost_.direct_access(v));
      }
    }
#endif
#ifdef HAS_XPRESS
    if (solver_ == "xpress") {

      Math1D::Vector<int> var_idx(nVars_,0);
      for (uint v=0; v < nVars_; v++) {
	var_idx[v] = v;
      }

      XPRSchgobj(xp_prob_, nVars_, var_idx.direct_access(), cost_.direct_access() );      
    }
#endif

    if (solver_ == "clp")
      lpSolver_.setObjective(cost_.direct_access());
  }

  std::cerr << "leaving segment" << std::endl;

  return energy;
}
