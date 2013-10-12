/******* written by Thomas Schoenemann. Started as a private person without employment, July 2011, *****/
/******* continued at the University of Pisa, Italy, October - December 2011 ***/
/******* continued at the University of Düsseldorf, Germany, October - December 2011 ***/

#include "factorDualOpt.hh"
#include <map>

DualVariableNode::DualVariableNode(const Math1D::Vector<float>& cost) : cost_(cost) {
  dual_var_.resize(cost.size(),0);
}

void DualVariableNode::add_cost(const Math1D::Vector<float>& add_cost) {
  cost_ += add_cost;
}
  
void DualVariableNode::add_factor(DualFactorNode* node) {

  uint nPrevFactors = neighboring_factor_.size();
  neighboring_factor_.resize(nPrevFactors+1,node);
  dual_var_.resize(nLabels(),nPrevFactors+1,0.0); 
}

const double* DualVariableNode::get_dual_var_start() const {

  return dual_var_.direct_access();
}

double* DualVariableNode::get_dual_var_start() {

  return dual_var_.direct_access();
}

  
double* DualVariableNode::get_dual_vars(const DualFactorNode* node) {
  
  double* ptr = dual_var_.direct_access();

  const uint label_size = dual_var_.xDim();
  const uint nFactors = neighboring_factor_.size();

  bool found = false;
  for (uint i=0; i < nFactors; i++) {
    if (neighboring_factor_[i] == node) {
      found = true;
      break;
    }
    else
      ptr += label_size;
  }

  if (!found) {
    INTERNAL_ERROR << "node not found" << std::endl;
    exit(1);
  }

  return ptr;
}

const double* DualVariableNode::get_dual_vars(const DualFactorNode* node) const {
  
  const double* ptr = dual_var_.direct_access();

  const uint label_size = dual_var_.xDim();
  const uint nFactors = neighboring_factor_.size();

  bool found = false;
  for (uint i=0; i < nFactors; i++) {
    if (neighboring_factor_[i] == node) {
      found = true;
      break;
    }
    else
      ptr += label_size;
  }

  if (!found) {
    INTERNAL_ERROR << "node not found" << std::endl;
    exit(1);
  }

  return ptr;
}


  
uint DualVariableNode::nLabels() const {
  return dual_var_.xDim();
}

void DualVariableNode::init_dual_vars() {
  dual_var_.set_constant(0.0);
}

const Storage1D<DualFactorNode*>& DualVariableNode::neighboring_factor() const {
  return neighboring_factor_;
}

  
void DualVariableNode::compute_message(const DualFactorNode* node, Math1D::Vector<double>& msg) const {

  const uint label_size = cost_.size();
  const uint nFactors = neighboring_factor_.size();

  msg.resize_dirty(label_size);
  
  for (uint l=0; l < label_size; l++)
    msg[l] = cost_[l];

  for (uint k=0; k < nFactors; k++) {

    if (neighboring_factor_[k] != node) {

      for (uint l=0; l < label_size; l++)
        msg[l] += dual_var_(l,k);
    }
  }
}

double DualVariableNode::cost(uint label) const {
  return cost_[label];
}

const Math1D::Vector<float>& DualVariableNode::cost() const {
  return cost_;
}

double DualVariableNode::dual_value(uint& arg_min) const {

  Math1D::Vector<double> msg;
  this->compute_message(0,msg);

  double min_val = 1e300;
  arg_min = MAX_UINT;

  for (uint k=0; k < msg.size(); k++) {
    if (msg[k] < min_val) {
      min_val = msg[k];
      arg_min = k;
    }
  }

  return min_val;
}

/**********************************/

DualFactorNode::DualFactorNode(const Storage1D<DualVariableNode*>& participating_vars)
  : participating_var_(participating_vars) {

  for (uint i=0; i < participating_var_.size(); i++)
    participating_var_[i]->add_factor(this);
}

/*virtual*/ DualFactorNode::~DualFactorNode() {}
  
//factors can modify the current labeling, e.g. to satisfy some constraints
/*virtual*/ bool DualFactorNode::process_labeling(Math1D::Vector<uint>& /*labels*/) const {
  return false;
}

const Storage1D<DualVariableNode*>& DualFactorNode::participating_nodes() const {
  return participating_var_;
}

/*virtual*/ double DualFactorNode::dual_value() const {

  Math1D::Vector<uint> labels;
  return compute_minimizer(labels);
}

/**********************************/

BinaryDualFactorNodeBase::BinaryDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars) :
  DualFactorNode(participating_vars) {

  assert(participating_vars.size() == 2);

  if (participating_vars.size() != 2 ){
    INTERNAL_ERROR << "attempt to instantiate a binary factor with " << participating_vars.size() << " variables. Exiting." << std::endl;
    exit(1);
  }
}

void BinaryDualFactorNodeBase::update_duals(const Math2D::Matrix<float>& cost, DualBCAMode mode) {


  const uint nLabels1 = cost.xDim();
  const uint nLabels2 = cost.yDim();

  //std::cerr << "bin" << std::endl;

  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());

  double* dp0 = participating_var_[0]->get_dual_vars(this);
  double* dp1 = participating_var_[1]->get_dual_vars(this);

  Math1D::Vector<double> msg0(nLabels1);
  Math1D::Vector<double> msg1(nLabels2);

  if (mode == DUAL_BCA_MODE_MPLP) {

    for (uint l1=0; l1 < nLabels1; l1++)
      dp0[l1] = 1e300;
    for (uint l2=0; l2 < nLabels2; l2++)
      dp1[l2] = 1e300;

    participating_var_[0]->compute_message(this,msg0);
    participating_var_[1]->compute_message(this,msg1);

    for (uint l1=0; l1 < nLabels1; l1++) {

      const double inter_val = msg0[l1]; 

      for (uint l2=0; l2 < nLabels2; l2++) {
	
	const double hyp = cost(l1,l2) + inter_val + msg1[l2]; 
	  
	if (hyp < dp0[l1])
	  dp0[l1] = hyp;
	if (hyp < dp1[l2])
	  dp1[l2] = hyp;
      }
    }

    for (uint l1=0; l1 < nLabels1; l1++)
      dp0[l1] = 0.5 * dp0[l1] - msg0[l1];
    for (uint l2=0; l2 < nLabels2; l2++)
      dp1[l2] = 0.5 * dp1[l2] - msg1[l2];

  }
  else {
    //MSD-Mode

    //var 1
    {
      participating_var_[0]->compute_message(this, msg0);

      for (uint l1=0; l1 < nLabels1; l1++) {

	double min_cost = 1e300;

	for (uint l2=0; l2 < nLabels2; l2++) {
	  
	  double hyp = cost(l1,l2) - dp1[l2];
	    
	  if (hyp < min_cost)
	    min_cost = hyp;
	}

	dp0[l1] = 0.5 * (min_cost - msg0[l1]);
      }
    }
    
    //var 2 
    {
      participating_var_[1]->compute_message(this, msg1);

      for (uint l2=0; l2 < nLabels2; l2++) {
	
	double min_cost = 1e300;

	for (uint l1=0; l1 < nLabels1; l1++) {
	  
	  double hyp = cost(l1,l2) - dp0[l1];
	    
          if (hyp < min_cost)
            min_cost = hyp;
	}

	dp1[l2] = 0.5 * (min_cost - msg1[l2]);
      }
    }
  }
}
  
double BinaryDualFactorNodeBase::dual_value(const Math2D::Matrix<float>& cost) const {

  NamedStorage1D<const double*> dual_ptr(2, MAKENAME(dual_ptr));

  uint nLabels1 = cost.xDim();
  uint nLabels2 = cost.yDim();
  
  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());

  for (uint v=0; v < 2; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  double min_cost = 1e300;
  
  for (uint l1=0; l1 < nLabels1; l1++) {

    double dual1 = dual_ptr[0][l1];

    for (uint l2=0; l2 < nLabels2; l2++) {
      double cur_cost = cost(l1,l2) - dual1 - dual_ptr[1][l2];
      if (cur_cost < min_cost)
        min_cost = cur_cost;
    }
  }

  return min_cost;
}

double BinaryDualFactorNodeBase::compute_minimizer(const Math2D::Matrix<float>& cost, Math1D::Vector<uint>& min_labels) const {

  min_labels.resize_dirty(2);

  NamedStorage1D<const double*> dual_ptr(2, MAKENAME(dual_ptr));

  uint nLabels1 = cost.xDim();
  uint nLabels2 = cost.yDim();
  
  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());

  for (uint v=0; v < 2; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  double min_cost = 1e300;
  
  for (uint l1=0; l1 < nLabels1; l1++) {
    for (uint l2=0; l2 < nLabels2; l2++) {
      double cur_cost = cost(l1,l2) - dual_ptr[0][l1] - dual_ptr[1][l2];
      if (cur_cost < min_cost) {
        min_cost = cur_cost;
        min_labels[0] = l1;
        min_labels[1] = l2;
      }
    }
  }

  return min_cost;
}

/**********************************/

PottsDualFactorNode::PottsDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, float lambda) :
  DualFactorNode(participating_vars), lambda_(lambda) {

  if (participating_vars.size() != 2 ){
    INTERNAL_ERROR << "attempt to instantiate a (standard) Potts factor with " 
                   << participating_vars.size() << " variables. Exiting." << std::endl;
    exit(1);
  }
}

/*virtual*/ PottsDualFactorNode::~PottsDualFactorNode() {}

/*virtual*/ double PottsDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  return (labels[0] == labels[1]) ? 0.0 : lambda_;
}

/*virtual*/ void PottsDualFactorNode::update_duals(DualBCAMode mode) {

  NamedStorage1D<Math1D::Vector<double> > msg(2, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(2, MAKENAME(dual_ptr));

  const uint nLabels1 = participating_var_[0]->nLabels();
  const uint nLabels2 = participating_var_[1]->nLabels();

  for (uint v=0; v < 2; v++) {
    if (mode == DUAL_BCA_MODE_MPLP)
      participating_var_[v]->compute_message(this, msg[v]);
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  double* dp0 = dual_ptr[0];
  double* dp1 = dual_ptr[1];

  if (mode == DUAL_BCA_MODE_MPLP) {

    Math1D::Vector<double> msg_min(2); //note: already includes lambda

    for (uint k=0; k < 2; k++) {
      
      const Math1D::Vector<double>& cur_msg = msg[k];
      
      const uint cur_size = participating_var_[k]->nLabels();
      
      double cur_min = 1e300;
      for (uint l=0; l < cur_size; l++) {
        if (cur_msg[l] < cur_min)
          cur_min = cur_msg[l];
      }
      
      msg_min[k] = cur_min + lambda_;
    }      

    //NOTE: this code assumes that nLabels1 == nLabels2

    //Message 1
    for (uint l=0; l < nLabels1; l++) {
	
      dp0[l] = 0.5 * (std::min(msg_min[1], msg[1][l]) - msg[0][l]);
    }
    
    //Message 2
    for (uint k=0; k < nLabels2; k++) {

      dp1[k] = 0.5 * (std::min(msg_min[0], msg[0][k]) - msg[1][k]);
    }  
  }
  else {
    // MSD-mode
    //NOTE: this code assumes that nLabels1 == nLabels2

    //variable 1
    participating_var_[0]->compute_message(this, msg[0]);

    double min2 = 1e300;
    for (uint l2=0; l2 < nLabels2; l2++) {
      if (-dp1[l2] < min2)
        min2 = -dp1[l2];
    }

    for (uint l1=0; l1 < nLabels1; l1++) {
      dp0[l1] = 0.5 * ( std::min(dp1[l1],min2+lambda_) - msg[0][l1]);
    }

    //variable 2
    participating_var_[1]->compute_message(this, msg[1]);

    double min1 = 1e300;
    for (uint l1=0; l1 < nLabels1; l1++) {
      if (-dp0[l1] < min1)
        min1 = -dp0[l1];
    }
    
    for (uint l2=0; l2 < nLabels2; l2++) {
      dp1[l2] = 0.5 * ( std::min(-dp0[l2],min1+lambda_) - msg[1][l2]);
    }
  }
}

/*virtual*/ double PottsDualFactorNode::dual_value() const {

  NamedStorage1D<const double*> dual_ptr(2, MAKENAME(dual_ptr));

  const uint nLabels1 = participating_var_[0]->nLabels();
  const uint nLabels2 = participating_var_[1]->nLabels();

  for (uint v=0; v < 2; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }
  
  Math1D::Vector<double> neg_dual_min(2);

  for (uint k=0; k < 2; k++) {

    double cur_min = 1e300;
    for (uint l=0; l < participating_var_[k]->nLabels(); l++) {
      double hyp = -dual_ptr[k][l];

      if (hyp < cur_min)
        cur_min = hyp;
    }
    
    neg_dual_min[k] = cur_min;
  }

  double min_cost = neg_dual_min.sum() + lambda_;

  for (uint l=0; l < std::min(nLabels1,nLabels2); l++) {

    double hyp =  -dual_ptr[0][l] - dual_ptr[1][l];

    if (hyp < min_cost)
      min_cost = hyp;
  }
  
  return min_cost;
}

/*virtual*/ double PottsDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {

  min_labels.resize(2);

  NamedStorage1D<const double*> dual_ptr(2, MAKENAME(dual_ptr));

  const uint nLabels1 = participating_var_[0]->nLabels();
  const uint nLabels2 = participating_var_[1]->nLabels();

  for (uint v=0; v < 2; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }
  
  Math1D::Vector<double> neg_dual_min(2);

  for (uint k=0; k < 2; k++) {

    double cur_min = 1e300;
    for (uint l=0; l < participating_var_[k]->nLabels(); l++) {
      double hyp = -dual_ptr[k][l];

      if (hyp < cur_min) {
        cur_min = hyp;
        min_labels[k] = l;
      }
    }
    
    neg_dual_min[k] = cur_min;
  }

  double min_cost = neg_dual_min.sum() + lambda_;

  for (uint l=0; l < std::min(nLabels1,nLabels2); l++) {

    double hyp =  -dual_ptr[0][l] - dual_ptr[1][l];

    if (hyp < min_cost) {
      min_cost = hyp;
      
      for (uint k=0; k < 2; k++)
        min_labels[k] = l;
    }
  }
  
  return min_cost;
}


/**********************************/

BinaryDualFactorNode::BinaryDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                                           const Math2D::Matrix<float>& cost) :
  BinaryDualFactorNodeBase(participating_vars), cost_(cost) {

  if (cost_.xDim() < participating_var_[0]->nLabels() || cost_.yDim() < participating_var_[1]->nLabels()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
  }
}

/*virtual*/ BinaryDualFactorNode::~BinaryDualFactorNode() {}

/*virtual*/ void BinaryDualFactorNode::update_duals(DualBCAMode mode) {
  BinaryDualFactorNodeBase::update_duals(cost_, mode);
}

/*virtual*/ double BinaryDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  return cost_(labels[0],labels[1]);
}

/*virtual*/ double BinaryDualFactorNode::dual_value() const {
  return BinaryDualFactorNodeBase::dual_value(cost_);
}

/*virtual*/ double BinaryDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {
  return BinaryDualFactorNodeBase::compute_minimizer(cost_,min_labels);
}

/**********************************/

BinaryDualRefFactorNode::BinaryDualRefFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                                                 const Math2D::Matrix<float>& cost) :
  BinaryDualFactorNodeBase(participating_vars), cost_(cost) {
}

/*virtual*/ BinaryDualRefFactorNode::~BinaryDualRefFactorNode() {}

/*virtual*/ void BinaryDualRefFactorNode::update_duals(DualBCAMode mode) {

  assert(cost_.xDim() >= participating_var_[0]->nLabels());
  assert(cost_.yDim() >= participating_var_[1]->nLabels());
  BinaryDualFactorNodeBase::update_duals(cost_, mode);
}

/*virtual*/ double BinaryDualRefFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  return cost_(labels[0],labels[1]);
}

/*virtual*/ double BinaryDualRefFactorNode::dual_value() const {
  return BinaryDualFactorNodeBase::dual_value(cost_);
}

/*virtual*/ double BinaryDualRefFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {
  return BinaryDualFactorNodeBase::compute_minimizer(cost_,min_labels);
}


/**********************************/

GenericDualFactorNode::GenericDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                                             const VarDimStorage<float>& cost) :
  DualFactorNode(participating_vars), cost_(cost)
{
  if (cost.nDims() != participating_vars.size()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
    exit(1);
  }

  for (uint v = 0; v < participating_vars.size(); v++) {
    if (cost.dim(v) < participating_vars[v]->nLabels()) {
      INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
      exit(1);
    }
  }
}

/*virtual*/ GenericDualFactorNode::~GenericDualFactorNode() {}

/*virtual*/ void GenericDualFactorNode::update_duals(DualBCAMode mode) {

  assert(mode == DUAL_BCA_MODE_MPLP);

  const uint nVars = participating_var_.size();

  NamedStorage1D<Math1D::Vector<double> > msg(nVars, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  Math1D::NamedVector<uint> nLabels(nVars, MAKENAME(nLabels));

  for (uint v=0; v < nVars; v++) {
    participating_var_[v]->compute_message(this, msg[v]);
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
    nLabels[v] = participating_var_[v]->nLabels();

    for (uint l=0; l < nLabels[v]; l++) 
      dual_ptr[v][l] = 1e300; 
  }

  Math1D::Vector<size_t> labeling(nVars,0);

  while (true) {
    
    double cost = cost_(labeling);
    for (uint v=0; v < nVars; v++) {
      cost += msg[v][labeling[v]];
    }

    for (uint v=0; v < nVars; v++) {

      if (cost < dual_ptr[v][labeling[v]])
        dual_ptr[v][labeling[v]] = cost; 
    }

    //increase labeling
    uint l;
    for (l=0; l < nVars; l++) {

      labeling[l] = (labeling[l] + 1) % nLabels[l];
      if (labeling[l] != 0)
        break;
    }

    if (l == nVars) //all zero after increase => cycle completed
      break;
  }

  for (uint v=0; v < nVars; v++) {
    for (uint l=0; l < nLabels[v]; l++) {
      dual_ptr[v][l] /= double(nVars);
      dual_ptr[v][l] -= msg[v][l];
    }
  }
}

/*virtual*/ double GenericDualFactorNode::dual_value() const {

  uint nVars = participating_var_.size();

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  Math1D::NamedVector<uint> nLabels(nVars, MAKENAME(nLabels));

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
    nLabels[v] = participating_var_[v]->nLabels();
  }

  Math1D::Vector<size_t> labeling(nVars,0);

  double min_val = 1e300;

  while (true) {
    
    double cost = cost_(labeling);
    for (uint v=0; v < nVars; v++) {
      cost -= dual_ptr[v][labeling[v]];
    }

    if (cost < min_val)
      min_val = cost;

    //increase labeling
    uint l;
    for (l=0; l < nVars; l++) {

      labeling[l] = (labeling[l] + 1) % nLabels[l];
      if (labeling[l] != 0)
        break;
    }

    if (l == nVars) //all zero after increase => cycle completed
      break;
  }

  return min_val;
}

/*virtual*/ double GenericDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {

  const uint nVars = participating_var_.size();

  min_labels.resize_dirty(nVars);

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  Math1D::NamedVector<uint> nLabels(nVars, MAKENAME(nLabels));

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
    nLabels[v] = participating_var_[v]->nLabels();
  }

  Math1D::Vector<size_t> labeling(nVars,0);

  double min_val = 1e300;

  while (true) {
    
    double cost = cost_(labeling);
    for (uint v=0; v < nVars; v++) {
      cost -= dual_ptr[v][labeling[v]];
    }

    if (cost < min_val) {
      min_val = cost;
      for (uint v=0; v < nVars; v++)
        min_labels[v] = labeling[v];
    }

    //increase labeling
    uint l;
    for (l=0; l < nVars; l++) {

      labeling[l] = (labeling[l] + 1) % nLabels[l];
      if (labeling[l] != 0)
        break;
    }

    if (l == nVars) //all zero after increase => cycle completed
      break;
  }

  return min_val;
}

/*virtual*/ double GenericDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  Math1D::Vector<size_t> size_t_labels(labels.size());

  for (uint k=0; k < labels.size(); k++)
    size_t_labels[k] = labels[k];

  return cost_(size_t_labels);
}

/**********************************/

TernaryDualFactorNodeBase::TernaryDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars)
  : DualFactorNode(participating_vars) {

  if (participating_vars.size() != 3){
    INTERNAL_ERROR << "attempt to instantiate a binary factor with " << participating_vars.size() << " variables. Exiting." << std::endl;
    exit(1);
  }
}

void TernaryDualFactorNodeBase::update_duals(const Math3D::Tensor<float>& cost, DualBCAMode mode) {

  const uint nLabels1 = cost.xDim();
  const uint nLabels2 = cost.yDim();
  const uint nLabels3 = cost.zDim();

  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());
  assert(nLabels3 == participating_var_[2]->nLabels());

  double* dp0 = participating_var_[0]->get_dual_vars(this);
  double* dp1 = participating_var_[1]->get_dual_vars(this);
  double* dp2 = participating_var_[2]->get_dual_vars(this);

  Math1D::Vector<double> msg0(nLabels1);
  Math1D::Vector<double> msg1(nLabels2);
  Math1D::Vector<double> msg2(nLabels3);


  if (mode == DUAL_BCA_MODE_MPLP) {

    for (uint l1=0; l1 < nLabels1; l1++)
      dp0[l1] = 1e300;
    for (uint l2=0; l2 < nLabels2; l2++)
      dp1[l2] = 1e300;
    for (uint l3=0; l3 < nLabels3; l3++)
      dp2[l3] = 1e300;

    participating_var_[0]->compute_message(this,msg0);
    participating_var_[1]->compute_message(this,msg1);
    participating_var_[2]->compute_message(this,msg2);
    

    for (uint l1=0; l1 < nLabels1; l1++) {

      const double inter_val = msg0[l1];

      for (uint l2=0; l2 < nLabels2; l2++) {
	
        const double part_sum = inter_val + msg1[l2];
	
        for (uint l3=0; l3 < nLabels3; l3++) {
	  
          const double hyp = cost(l1,l2,l3) + part_sum + msg2[l3];
	  
          if (hyp < dp0[l1])
            dp0[l1] = hyp;
          if (hyp < dp1[l2])
            dp1[l2] = hyp;
          if (hyp < dp2[l3])
            dp2[l3] = hyp;
        }
      }
    }

    for (uint l1=0; l1 < nLabels1; l1++)
      dp0[l1] = (dp0[l1] / 3.0) - msg0[l1];
    for (uint l2=0; l2 < nLabels2; l2++)
      dp1[l2] = (dp1[l2] / 3.0) - msg1[l2];
    for (uint l3=0; l3 < nLabels3; l3++)
      dp2[l3] = (dp2[l3] / 3.0) - msg2[l3];
  }
  else {
    //MSD-Mode

    //var 1
    {

      participating_var_[0]->compute_message(this, msg0);

      for (uint l1=0; l1 < nLabels1; l1++) {

	double min_cost = 1e300;
	
	for (uint l2=0; l2 < nLabels2; l2++) {
	  
	  const double inter_val = dp1[l2];
	  
	  for (uint l3=0; l3 < nLabels3; l3++) {
	    
	    double hyp = cost(l1,l2,l3) - inter_val - dp2[l3];
	    
	    if (hyp < min_cost)
	      min_cost = hyp;
	  }
	}
      
	dp0[l1] = 0.5 * (min_cost - msg0[l1]);
      }
    }
    
    //var 2 
    {
      participating_var_[1]->compute_message(this, msg1);
	
      for (uint l2=0; l2 < nLabels2; l2++) {
	
	double min_cost = 1e300;
	
	for (uint l1=0; l1 < nLabels1; l1++) {
	  
	  const double inter_val = dp0[l1];
	  
	  for (uint l3=0; l3 < nLabels3; l3++) {
	    
	    double hyp = cost(l1,l2,l3) - inter_val - dp2[l3];
	    
          if (hyp < min_cost)
            min_cost = hyp;
	  }
	}

	dp1[l2] = 0.5 * (min_cost - msg1[l2]);
      }
    }

    //var 3
    {
      participating_var_[2]->compute_message(this, msg2);

      for (uint l3=0; l3 < nLabels3; l3++) {
	
	double min_cost = 1e300;
	
	for (uint l1=0; l1 < nLabels1; l1++) {
	  
	  const double inter_val = dp0[l1];
	  
	  for (uint l2=0; l2 < nLabels2; l2++) {
	    
	    double hyp = cost(l1,l2,l3) - inter_val - dp1[l2];
	    
	    if (hyp < min_cost)
	      min_cost = hyp;
	  }
	}
	
	dp2[l3] = 0.5 * (min_cost - msg2[l3]);
      }
    }
  }

}

double TernaryDualFactorNodeBase::dual_value(const Math3D::Tensor<float>& cost) const {

  NamedStorage1D<const double*> dual_ptr(3, MAKENAME(dual_ptr));

  for (uint v=0; v < 3; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  const uint nLabels1 = cost.xDim();
  const uint nLabels2 = cost.yDim();
  const uint nLabels3 = cost.zDim();

  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());
  assert(nLabels3 == participating_var_[2]->nLabels());

  double min_val = 1e300;

  for (uint l1=0; l1 < nLabels1; l1++) {

    const double inter_val = -dual_ptr[0][l1];
    
    for (uint l2=0; l2 < nLabels2; l2++) {
      
      const double part_sum = inter_val - dual_ptr[1][l2];

      for (uint l3=0; l3 < nLabels3; l3++) {

        const double hyp = cost(l1,l2,l3) + part_sum - dual_ptr[2][l3];

        if (hyp < min_val)
          min_val = hyp;
      }
    }
  }  

  return min_val;
}

double TernaryDualFactorNodeBase::compute_minimizer(const Math3D::Tensor<float>& cost, Math1D::Vector<uint>& min_labels) const {

  min_labels.resize_dirty(3);

  NamedStorage1D<const double*> dual_ptr(3, MAKENAME(dual_ptr));

  for (uint v=0; v < 3; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  const uint nLabels1 = cost.xDim();
  const uint nLabels2 = cost.yDim();
  const uint nLabels3 = cost.zDim();

  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());
  assert(nLabels3 == participating_var_[2]->nLabels());

  double min_val = 1e300;

  for (uint l1=0; l1 < nLabels1; l1++) {

    for (uint l2=0; l2 < nLabels2; l2++) {
      
      for (uint l3=0; l3 < nLabels3; l3++) {

        const double hyp = cost(l1,l2,l3) - dual_ptr[0][l1] - dual_ptr[1][l2] - dual_ptr[2][l3];

        if (hyp < min_val) {
          min_val = hyp;
          min_labels[0] = l1;
          min_labels[1] = l2;
          min_labels[2] = l3;
        }
      }
    }
  }  

  return min_val;
}

/**********************************/

TernaryDualFactorNode::TernaryDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                                             const Math3D::Tensor<float>& cost) :
  TernaryDualFactorNodeBase(participating_vars), cost_(cost) {

  if (cost_.xDim() < participating_var_[0]->nLabels() || cost_.yDim() < participating_var_[1]->nLabels()
      || cost_.zDim() < participating_var_[2]->nLabels()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
  }
}

/*virtual*/ TernaryDualFactorNode::~TernaryDualFactorNode() {}

/*virtual*/ void TernaryDualFactorNode::update_duals(DualBCAMode mode) {
  TernaryDualFactorNodeBase::update_duals(cost_, mode);
}

/*virtual*/ double TernaryDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  return cost_(labels[0],labels[1],labels[2]);
}

/*virtual*/ double TernaryDualFactorNode::dual_value() const {
  return TernaryDualFactorNodeBase::dual_value(cost_);
}

/*virtual*/ double TernaryDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {
  return TernaryDualFactorNodeBase::compute_minimizer(cost_,min_labels);
}

/**********************************/

TernaryDualRefFactorNode::TernaryDualRefFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                                                   const Math3D::Tensor<float>& cost) :
  TernaryDualFactorNodeBase(participating_vars), cost_(cost) {}

/*virtual*/ TernaryDualRefFactorNode::~TernaryDualRefFactorNode() {}

/*virtual*/ void TernaryDualRefFactorNode::update_duals(DualBCAMode mode) {

  assert(cost_.xDim() >= participating_var_[0]->nLabels());
  assert(cost_.yDim() >= participating_var_[1]->nLabels());
  assert(cost_.zDim() >= participating_var_[2]->nLabels());

  TernaryDualFactorNodeBase::update_duals(cost_, mode);
}

/*virtual*/ double TernaryDualRefFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  return cost_(labels[0],labels[1],labels[2]);
}

/*virtual*/ double TernaryDualRefFactorNode::dual_value() const {
  return TernaryDualFactorNodeBase::dual_value(cost_);
}

/*virtual*/ double TernaryDualRefFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {
  return TernaryDualFactorNodeBase::compute_minimizer(cost_,min_labels);
}

/**********************************/

SecondDiffDualFactorNode::SecondDiffDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, float lambda) :
  DualFactorNode(participating_vars), lambda_(lambda) {

  if (participating_vars.size() != 3){
    INTERNAL_ERROR << "attempt to instantiate a second difference factor with " 
                   << participating_vars.size() << " variables. Exiting." << std::endl;
    exit(1);
  }
}  

/*virtual*/ SecondDiffDualFactorNode::~SecondDiffDualFactorNode() {}

/*virtual*/
void SecondDiffDualFactorNode::update_duals(DualBCAMode mode) {

  const uint nLabels1 = participating_var_[0]->nLabels();
  const uint nLabels2 = participating_var_[1]->nLabels();
  const uint nLabels3 = participating_var_[2]->nLabels();

  NamedStorage1D<Math1D::Vector<double> > msg(3, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(3, MAKENAME(dual_ptr));

  for (uint v=0; v < 3; v++) {

    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);

    if (mode == DUAL_BCA_MODE_MPLP) {
      participating_var_[v]->compute_message(this, msg[v]);

      for (uint l=0; l < participating_var_[v]->nLabels() ; l++) {
        dual_ptr[v][l] = 1e300;
      }
    }
    else
      msg[v].resize(participating_var_[v]->nLabels());
  }

  if (mode == DUAL_BCA_MODE_MPLP) {

    Math1D::Vector<double> min_msg(3);
    for (uint v=0; v < 3; v++) {
      min_msg[v] = msg[v].min();
    }

    for (int l1=0; l1 < int(nLabels1); l1++) {

      double best = min_msg[1] + min_msg[2] + 3*lambda_;

      for (int l2 = std::max(0,l1-1); l2 <= std::min<int>(nLabels2-1,l1+1); l2++) {

        const double m2 = msg[1][l2];

        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {
	
          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = m2+ msg[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = m2 + msg[2][l3] + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      dual_ptr[0][l1] = best + msg[0][l1];
    }

    for (int l2=0; l2 < int(nLabels2); l2++) {

      double best = min_msg[0] + min_msg[2] + 3*lambda_;

      for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {

        const double m1 = msg[0][l1];

        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = m1 + msg[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = m1 + msg[2][l3] + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      dual_ptr[1][l2] = best + msg[1][l2];
    }

    for (int l3=0; l3 < int(nLabels3); l3++) {

      double best = min_msg[0] + min_msg[1] + 3*lambda_;

      for (int l2 = std::max(0,l3-1); l2 <= std::min<int>(nLabels2-1,l3+1); l2++) {

         const double m2 = msg[1][l2];

        for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = msg[0][l1] + m2;
          }
          else if (abs(so_diff) <= 1) {
            hyp = msg[0][l1] + m2 + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      dual_ptr[2][l3] = best + msg[2][l3];
    }

    for (uint k=0; k < 3; k++) {
      for (uint l=0; l < participating_var_[k]->nLabels() ; l++) {
        dual_ptr[k][l] *= 1.0 / 3.0;
        dual_ptr[k][l] -= msg[k][l];
      }
    }
  }
  else {
    //MSD-Mode

    //var 1
    participating_var_[0]->compute_message(this, msg[0]);

    double base_cost = 3*lambda_;

    double mincost2 = 1e300;
    for (uint l2=0; l2 < nLabels2; l2++) {

      if (-dual_ptr[1][l2] < mincost2)
        mincost2 = -dual_ptr[1][l2];
    }

    base_cost += mincost2;

    double mincost3 = 1e300;
    for (uint l3=0; l3 < nLabels3; l3++) {

      if (-dual_ptr[2][l3] < mincost3)
        mincost3 = -dual_ptr[2][l3];
    }

    base_cost += mincost3;

    for (int l1=0; l1 < int(nLabels1); l1++) {

      double min_cost = base_cost;

      for (int l2 = std::max(0,l1-1); l2 <= std::min<int>(nLabels2-1,l1+1); l2++) {

        const double m2 = -dual_ptr[1][l2];

        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {
	
          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = m2 - dual_ptr[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = m2 - dual_ptr[2][l3] + lambda_;
          }

          if (hyp < min_cost)
            min_cost = hyp;
        }
      }      
      
      dual_ptr[0][l1] = 0.5 * (min_cost - msg[0][l1]);
    }

    //var 2
    participating_var_[1]->compute_message(this, msg[1]);

    base_cost = 3*lambda_ + mincost3;

    double mincost1 = 1e300;
    for (uint l1=0; l1 < nLabels1; l1++) {

      if (-dual_ptr[0][l1] < mincost1)
        mincost1 = -dual_ptr[0][l1];
    }

    base_cost += mincost1;
    
    for (int l2=0; l2 < int(nLabels2); l2++) {

      double min_cost = base_cost;

      for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {

        const double m1 = -dual_ptr[0][l1];
        
        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = m1 - dual_ptr[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = m1 - dual_ptr[2][l3] + lambda_;
          }

          if (hyp < min_cost)
            min_cost = hyp;
        }
      }

      dual_ptr[1][l2] = 0.5 * (min_cost - msg[1][l2]);      
    }    
    
    //var 3
    participating_var_[2]->compute_message(this, msg[2]);

    base_cost = 3*lambda_ + mincost1;

    mincost2 = 1e300;

    for (uint l2=0; l2 < nLabels2; l2++) {

      if (-dual_ptr[1][l2] < mincost2)
        mincost2 = -dual_ptr[1][l2];
    }

    base_cost += mincost2;

    for (int l3=0; l3 < int(nLabels3); l3++) {

      double min_cost = base_cost;

      for (int l2 = std::max(0,l3-1); l2 <= std::min<int>(nLabels2-1,l3+1); l2++) {

        const double m2 = - dual_ptr[1][l2];

        for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = l3 - 2*l2 + l1;
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -dual_ptr[0][l1] + m2;
          }
          else if (abs(so_diff) <= 1) {
            hyp = -dual_ptr[0][l1] + m2 + lambda_;
          }

          if (hyp < min_cost)
            min_cost = hyp;
        }
      }

      dual_ptr[2][l3] = 0.5 * (min_cost - msg[2][l3]);
    }
  }

}

/*virtual*/ 
double SecondDiffDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  int diff1 = labels[1] - labels[0];
  int diff2 = labels[2] - labels[1];
  
  int so_diff = diff2 - diff1;
  
  if (abs(diff1) <= 1 && abs(diff2) <= 1 && so_diff == 0)
    return 0.0; //no cost
  else if (abs(diff1) <= 1 && abs(diff2) <= 1 && abs(so_diff) == 1)
    return lambda_;

  return 3*lambda_;
}

/*virtual*/ 
double SecondDiffDualFactorNode::dual_value() const {

  NamedStorage1D<const double*> dual_ptr(3, MAKENAME(dual_ptr));

  for (uint v=0; v < 3; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  double best = 3*lambda_;
  for (uint v=0; v < 3; v++) {
    double cur_min = 1e300;

    for (uint l=0; l < participating_var_[v]->nLabels(); l++) {
      if (-dual_ptr[v][l] < cur_min)
        cur_min = -dual_ptr[v][l];
    }

    best += cur_min;
  }

  int nLabels1 = participating_var_[0]->nLabels();
  int nLabels2 = participating_var_[1]->nLabels();
  int nLabels3 = participating_var_[2]->nLabels();

  for (int l1=0; l1 < nLabels1; l1++) {

    double best_inter = 1e300;

    for (int l2=std::max(0,l1-1); l2 < std::min(l1+1,nLabels2-1); l2++) {

      const double dual2 = -dual_ptr[1][l2];

      for (int l3=std::max(0,l2-1); l3 < std::min(l2+1,nLabels3-1); l3++) {
      
        assert(abs(l2-l1) <= 1);
        assert(abs(l3-l2) <= 1);

        int so_diff = l3 - 2*l2 + l1;
        
        double hyp = 1e300;
	
        if (so_diff == 0) {
          hyp = dual2 - dual_ptr[2][l3];
        }
        else if (abs(so_diff) <= 1) {
          hyp = dual2 - dual_ptr[2][l3] + lambda_;
        }

        if (hyp < best_inter)
          best_inter = hyp;
      }
    }

    if (best_inter - dual_ptr[0][l1] < best)
      best = best_inter - dual_ptr[0][l1];
  }

  return best;
}

/*virtual*/
double SecondDiffDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {

  min_labels.resize_dirty(3);

  NamedStorage1D<const double*> dual_ptr(3, MAKENAME(dual_ptr));

  for (uint v=0; v < 3; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  double best = 3*lambda_;
  for (uint v=0; v < 3; v++) {
    double cur_min = 1e300;

    for (uint l=0; l < participating_var_[v]->nLabels(); l++) {
      if (-dual_ptr[v][l] < cur_min) {
        cur_min = -dual_ptr[v][l];
        min_labels[v] = l;
      }
    }

    best += cur_min;
  }

  int nLabels1 = participating_var_[0]->nLabels();
  int nLabels2 = participating_var_[1]->nLabels();
  int nLabels3 = participating_var_[2]->nLabels();

  for (int l1=0; l1 < nLabels1; l1++) {

    double best_inter = 1e300;

    int opt_l2 = MAX_UINT;
    int opt_l3 = MAX_UINT;


    for (int l2=std::max(0,l1-1); l2 < std::min(l1+1,nLabels2-1); l2++) {
      for (int l3=std::max(0,l2-1); l3 < std::min(l2+1,nLabels3-1); l3++) {
      
        assert(abs(l2-l1) <= 1);
        assert(abs(l3-l2) <= 1);

        int so_diff = l3 - 2*l2 + l1;
	
        double hyp = 1e300;
	
        if (so_diff == 0) {
          hyp = -dual_ptr[1][l2] - dual_ptr[2][l3];
        }
        else if (abs(so_diff) <= 1) {
          hyp = -dual_ptr[1][l2] - dual_ptr[2][l3] + lambda_;
        }

        if (hyp < best_inter) {
          best_inter = hyp;
          opt_l2 = l2;
          opt_l3 = l3;
        }
      }
    }

    if (best_inter - dual_ptr[0][l1] < best) {
      best = best_inter - dual_ptr[0][l1];
      min_labels[0] = l1;
      min_labels[1] = opt_l2;
      min_labels[2] = opt_l3;
    }
  }

  return best;
}


/**********************************/

FourthOrderDualFactorNodeBase::FourthOrderDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars) :
  DualFactorNode(participating_vars) {

  if (participating_vars.size() != 4){
    INTERNAL_ERROR << "attempt to instantiate a fourth order factor with " 
                   << participating_vars.size() << " variables. Exiting." << std::endl;
    exit(1);
  }
}

void FourthOrderDualFactorNodeBase::update_duals(const Storage1D<Math3D::Tensor<float> >& cost, DualBCAMode mode) {

  const uint nLabels1 = cost.size();
  const uint nLabels2 = cost[0].xDim();
  const uint nLabels3 = cost[0].yDim();
  const uint nLabels4 = cost[0].zDim();

  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());
  assert(nLabels3 == participating_var_[2]->nLabels());
  assert(nLabels4 == participating_var_[3]->nLabels());

  double* dp0 = participating_var_[0]->get_dual_vars(this);
  double* dp1 = participating_var_[1]->get_dual_vars(this);
  double* dp2 = participating_var_[2]->get_dual_vars(this);
  double* dp3 = participating_var_[3]->get_dual_vars(this);

  Math1D::Vector<double> msg0(nLabels1);
  Math1D::Vector<double> msg1(nLabels2);
  Math1D::Vector<double> msg2(nLabels3);
  Math1D::Vector<double> msg3(nLabels4);


  if (mode == DUAL_BCA_MODE_MPLP) {

    for (uint l1=0; l1 < nLabels1; l1++)
      dp0[l1] = 1e300;
    for (uint l2=0; l2 < nLabels2; l2++)
      dp1[l2] = 1e300;
    for (uint l3=0; l3 < nLabels3; l3++)
      dp2[l3] = 1e300;
    for (uint l4=0; l4 < nLabels4; l4++)
      dp3[l4] = 1e300;


    participating_var_[0]->compute_message(this,msg0);
    participating_var_[1]->compute_message(this,msg1);
    participating_var_[2]->compute_message(this,msg2);
    participating_var_[3]->compute_message(this,msg3);


    for (uint l1=0; l1 < nLabels1; l1++) {

      double* cur_dp0 = dp0 + l1;

      const Math3D::Tensor<float>& cur_cost = cost[l1];

      const double inter1 = msg0[l1];

      for (uint l2=0; l2 < nLabels2; l2++) {
	
        double* cur_dp1 = dp1 + l2;
	
        const double part_sum1 = inter1 + msg1[l2];
	
        for (uint l3=0; l3 < nLabels3; l3++) {
	  
          const double part_sum2 = part_sum1 + msg2[l3];
	  
          for (uint l4=0; l4 < nLabels4; l4++) {
	    
            const double hyp = cur_cost(l2,l3,l4) + part_sum2 + msg3[l4];
	    
            if (hyp < (*cur_dp0))
              *cur_dp0 = hyp;
            if (hyp < (*cur_dp1))
              *cur_dp1 = hyp;
            if (hyp < dp2[l3])
              dp2[l3] = hyp;
            if (hyp < dp3[l4])
              dp3[l4] = hyp;
          }
        }
      }
    }

    for (uint l1=0; l1 < nLabels1; l1++)
      dp0[l1] = 0.25 * dp0[l1] - msg0[l1];
    for (uint l2=0; l2 < nLabels2; l2++)
      dp1[l2] = 0.25 * dp1[l2] - msg1[l2];
    for (uint l3=0; l3 < nLabels3; l3++)
      dp2[l3] = 0.25 * dp2[l3] - msg2[l3];
    for (uint l4=0; l4 < nLabels4; l4++)
      dp3[l4] = 0.25 * dp3[l4] - msg3[l4];
  }
  else {
    
    //var 1
    participating_var_[0]->compute_message(this, msg0);

    for (uint l1=0; l1 < nLabels1; l1++) {

      double min_cost = 1e300;

      const Math3D::Tensor<float>& cur_cost = cost[l1];

      for (uint l2=0; l2 < nLabels2; l2++) {
	
        const double inter1 = dp1[l2];

        for (uint l3=0; l3 < nLabels3; l3++) {

          const double inter2 = inter1 + dp2[l3];

          for (uint l4=0; l4 < nLabels4; l4++) {
	    
            double hyp = cur_cost(l2,l3,l4) - inter2 - dp3[l4];
	  
            if (hyp < min_cost)
              min_cost = hyp;
          }
        }
      }
      
      dp0[l1] = 0.5 * (min_cost - msg0[l1]);
    }

    //var 2
    participating_var_[1]->compute_message(this, msg1);

    for (uint l2=0; l2 < nLabels2; l2++) {

      double min_cost = 1e300;

      for (uint l1=0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost[l1];

        const double inter1 = dp0[l1];

        for (uint l3=0; l3 < nLabels3; l3++) {

          const double inter2 = inter1 + dp2[l3];

          for (uint l4=0; l4 < nLabels4; l4++) {

            double hyp = cur_cost(l2,l3,l4) - inter2 - dp3[l4];
	  
            if (hyp < min_cost)
              min_cost = hyp;
          }
        }
      }

      dp1[l2] = 0.5 * (min_cost - msg1[l2]);
    }

    //var 3
    participating_var_[2]->compute_message(this, msg2);
    
    for (uint l3=0; l3 < nLabels3; l3++) {

      double min_cost = 1e300;

      for (uint l1=0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost[l1];

        const double inter1 = dp0[l1];

        for (uint l2=0; l2 < nLabels2; l2++) {

          const double inter2 = inter1 + dp1[l2];

          for (uint l4=0; l4 < nLabels4; l4++) {

            double hyp = cur_cost(l2,l3,l4) - inter2 - dp3[l4];
	  
            if (hyp < min_cost)
              min_cost = hyp;
          }
        }
      }

      dp2[l3] = 0.5 * (min_cost - msg2[l3]);
    }

    //var 4
    participating_var_[3]->compute_message(this, msg3);
    
    for (uint l4=0; l4 < nLabels4; l4++) {

      double min_cost = 1e300;

      for (uint l1=0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost[l1];

        const double inter1 = dp0[l1];

        for (uint l2=0; l2 < nLabels2; l2++) {

          const double inter2 = inter1 + dp1[l2];

          for (uint l3=0; l3 < nLabels3; l3++) {

            double hyp = cur_cost(l2,l3,l4) - inter2 - dp2[l3];
	  
            if (hyp < min_cost)
              min_cost = hyp;
          }
        }
      }

      dp3[l4] = 0.5 * (min_cost - msg3[l4]);
    }
  }

}

double FourthOrderDualFactorNodeBase::dual_value(const Storage1D<Math3D::Tensor<float> >& cost) const {

  NamedStorage1D<const double*> dual_ptr(4, MAKENAME(dual_ptr));

  for (uint v=0; v < 4; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  const uint nLabels1 = cost.size();
  const uint nLabels2 = cost[0].xDim();
  const uint nLabels3 = cost[0].yDim();
  const uint nLabels4 = cost[0].zDim();

  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());
  assert(nLabels3 == participating_var_[2]->nLabels());
  assert(nLabels4 == participating_var_[3]->nLabels());

  double min_val = 1e300;

  for (uint l1=0; l1 < nLabels1; l1++) {

    const Math3D::Tensor<float>& cur_cost = cost[l1];

    const double dual1 = - dual_ptr[0][l1];

    for (uint l2=0; l2 < nLabels2; l2++) {

      const double dual2 = dual1 - dual_ptr[1][l2];
      
      for (uint l3=0; l3 < nLabels3; l3++) {
	
        const double dual3 = dual2 - dual_ptr[2][l3];

        for (uint l4=0; l4 < nLabels4; l4++) {
	  
          const double hyp = cur_cost(l2,l3,l4) + dual3 - dual_ptr[3][l4];
	  
          if (hyp < min_val)
            min_val = hyp;
        }
      }
    }
  }  

  return min_val;
}

double FourthOrderDualFactorNodeBase::compute_minimizer(const Storage1D<Math3D::Tensor<float> >& cost, 
                                                        Math1D::Vector<uint>& min_labels) const {

  min_labels.resize_dirty(4);

  NamedStorage1D<const double*> dual_ptr(4, MAKENAME(dual_ptr));

  for (uint v=0; v < 4; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  const uint nLabels1 = cost.size();
  const uint nLabels2 = cost[0].xDim();
  const uint nLabels3 = cost[0].yDim();
  const uint nLabels4 = cost[0].zDim();

  assert(nLabels1 == participating_var_[0]->nLabels());
  assert(nLabels2 == participating_var_[1]->nLabels());
  assert(nLabels3 == participating_var_[2]->nLabels());
  assert(nLabels4 == participating_var_[3]->nLabels());

  double min_val = 1e300;

  for (uint l1=0; l1 < nLabels1; l1++) {

    const Math3D::Tensor<float>& cur_cost = cost[l1];

    for (uint l2=0; l2 < nLabels2; l2++) {
      
      for (uint l3=0; l3 < nLabels3; l3++) {
	
        for (uint l4=0; l4 < nLabels4; l4++) {
	  
          const double hyp = cur_cost(l2,l3,l4) - dual_ptr[0][l1] - dual_ptr[1][l2] - dual_ptr[2][l3] - dual_ptr[3][l4];
	  
          if (hyp < min_val) {
            min_val = hyp;

            min_labels[0] = l1;
            min_labels[1] = l2;
            min_labels[2] = l3;
            min_labels[3] = l4;
          }
        }
      }
    }
  }  

  return min_val;
}


/**********************************/

FourthOrderDualFactorNode::FourthOrderDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                                                     const Storage1D<Math3D::Tensor<float> >& cost) :
  FourthOrderDualFactorNodeBase(participating_vars), cost_(cost) {

  if (cost_[0].size() < participating_var_[0]->nLabels() || cost_[0].xDim() < participating_var_[1]->nLabels()
      || cost_[0].yDim() < participating_var_[2]->nLabels() || cost_[0].zDim() < participating_var_[3]->nLabels()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
  }
}

/*virtual*/ FourthOrderDualFactorNode::~FourthOrderDualFactorNode() {}

/*virtual*/ void FourthOrderDualFactorNode::update_duals(DualBCAMode mode) {
  FourthOrderDualFactorNodeBase::update_duals(cost_, mode);
}

/*virtual*/ double FourthOrderDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  return cost_[labels[0]](labels[1],labels[2],labels[3]);
}

/*virtual*/ double FourthOrderDualFactorNode::dual_value() const {
  return FourthOrderDualFactorNodeBase::dual_value(cost_);
}

/*virtual*/ double FourthOrderDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {
  return FourthOrderDualFactorNodeBase::compute_minimizer(cost_,min_labels);
}


/**********************************/

FourthOrderDualRefFactorNode::FourthOrderDualRefFactorNode(const Storage1D<DualVariableNode*>& participating_vars, 
                                                           const Storage1D<Math3D::Tensor<float> >& cost) :
  FourthOrderDualFactorNodeBase(participating_vars), cost_(cost) {
}

/*virtual*/ FourthOrderDualRefFactorNode::~FourthOrderDualRefFactorNode() {}

/*virtual*/ void FourthOrderDualRefFactorNode::update_duals(DualBCAMode mode) {

  assert(cost_[0].size() >= participating_var_[0]->nLabels());
  assert(cost_[0].xDim() >= participating_var_[1]->nLabels());
  assert(cost_[0].yDim() >= participating_var_[2]->nLabels()); 
  assert(cost_[0].zDim() >= participating_var_[3]->nLabels());
  
  FourthOrderDualFactorNodeBase::update_duals(cost_, mode);
}

/*virtual*/ double FourthOrderDualRefFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  return cost_[labels[0]](labels[1],labels[2],labels[3]);
}

/*virtual*/ double FourthOrderDualRefFactorNode::dual_value() const {
  return FourthOrderDualFactorNodeBase::dual_value(cost_);
}

/*virtual*/ double FourthOrderDualRefFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {
  return FourthOrderDualFactorNodeBase::compute_minimizer(cost_,min_labels);
}


/**********************************/

OneOfNDualFactorNode::OneOfNDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars) :
  DualFactorNode(participating_vars) {

  for (uint v=0; v < participating_vars.size(); v++) {
    if (participating_vars[v]->nLabels() != 2) {
      INTERNAL_ERROR << "instantiation of a 1-of-N factor with non-binary variables. Exiting." << std::endl;
      exit(1);
    }
  }
}

/*virtual*/ OneOfNDualFactorNode::~OneOfNDualFactorNode() {}

/*virtual*/ double OneOfNDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  return (labels.sum() == 1) ? 0.0 : 1e75;
}

/*virtual*/ void OneOfNDualFactorNode::update_duals(DualBCAMode mode) {

  const uint nVars = participating_var_.size();

  NamedStorage1D<Math1D::Vector<double> > msg(nVars, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    if (mode == DUAL_BCA_MODE_MPLP)
      participating_var_[v]->compute_message(this, msg[v]);
    else
      msg[v].resize(2);
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  if (mode == DUAL_BCA_MODE_MPLP) {

    Math1D::Vector<double> rel_msg(nVars);

    double offs = 0.0;
    
    for (uint k=0; k < nVars; k++) {
      
      rel_msg[k] = msg[k][1] - msg[k][0];
      offs += msg[k][0];
    }
    
    double best = 1e15;
    uint arg_best = MAX_UINT;
    double second_best = 1e15;
    
    for (uint k=0; k < nVars; k++) {
      
      if (rel_msg[k] < best) {
        
        second_best = best;
        
        best = rel_msg[k];
        arg_best = k;
      }
      else if (rel_msg[k] < second_best) {
      
        second_best = rel_msg[k];
      }
    }
    
    best += offs;
    best /= nVars;
    second_best += offs;
    second_best /= nVars;
    
    for (uint k=0; k < nVars; k++) {
      
      //msg0
      if (k == arg_best)
        dual_ptr[k][0] = second_best - msg[k][0];
      else
        dual_ptr[k][0] = best - msg[k][0];
      
      //msg 1
      dual_ptr[k][1] = ((offs + rel_msg[k]) / nVars) - msg[k][1];
    }
  }
  else {

    assert(nVars >= 2);

    double offs = 0.0;

    double best = 1e300;
    uint arg_best = MAX_UINT;
    double second_best = 1e300;
    uint arg_second_best = MAX_UINT;

    for (uint v=0; v < nVars; v++) {
      double cur_diff = dual_ptr[v][0]-dual_ptr[v][1];
      if (cur_diff < best) {
        second_best = best;
        arg_second_best = arg_best;
        best = cur_diff;
        arg_best = v;
      }
      else if (cur_diff < second_best) {
        second_best = cur_diff;
        arg_second_best = v;
      }

      offs -= dual_ptr[v][0];
    }

    for (uint v=0; v < nVars; v++) {

      participating_var_[v]->compute_message(this, msg[v]);

      offs += dual_ptr[v][0]; //the current dual vars don't enter!

      double prev_diff = dual_ptr[v][0] - dual_ptr[v][1];
      dual_ptr[v][1] = 0.5 * (offs - msg[v][1]);
      if (v == arg_best) 
        dual_ptr[v][0] = 0.5 * (offs + second_best - msg[v][0]);
      else
        dual_ptr[v][0] = 0.5 * (offs + best - msg[v][0]);

      offs -= dual_ptr[v][0];

      double new_diff = dual_ptr[v][0] - dual_ptr[v][1];
      
      if (new_diff <= prev_diff) {
        if (v == arg_best)
          best = new_diff;
        else if (v == arg_second_best) {
          if (new_diff < best) {
            std::swap(best,second_best);
            std::swap(arg_best,arg_second_best);
          }
          else
            second_best = new_diff;
        }
      }
      else {
        if (v == arg_best && new_diff < second_best)
          best = new_diff;
        else if (v == arg_best || v == arg_second_best) {
          //recompute
          best = 1e300;
          second_best = 1e300;
          for (uint vv=0; vv < nVars; vv++) {
            double cur_diff = dual_ptr[vv][0] - dual_ptr[vv][1];
            if (cur_diff < best) {
              second_best = best;
              arg_second_best = arg_best;
              best = cur_diff;
              arg_best = vv; 
            }
            else if (cur_diff < second_best) {
              second_best = cur_diff;
              arg_second_best = vv;
            }
          }
        }
      }
    }    
  }
}

/*virtual*/ double OneOfNDualFactorNode::dual_value() const {

  uint nVars = participating_var_.size();

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  double dual_zero_sum = 0.0;

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
    dual_zero_sum -= dual_ptr[v][0];
  }

  double min_val = 1e300;

  for (uint v=0; v < nVars; v++) {

    double hyp = dual_zero_sum + dual_ptr[v][0] - dual_ptr[v][1];

    if (hyp < min_val)
      min_val = hyp;
  }
  
  return min_val;
}

/*virtual*/ double OneOfNDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {

  const uint nVars = participating_var_.size();

  min_labels.resize_dirty(nVars);
  min_labels.set_constant(0);

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  double dual_zero_sum = 0.0;

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
    dual_zero_sum -= dual_ptr[v][0];
  }

  double min_val = 1e300;
  uint arg_min = MAX_UINT;

  for (uint v=0; v < nVars; v++) {

    double hyp = dual_zero_sum + dual_ptr[v][0] - dual_ptr[v][1];

    if (hyp < min_val) {
      min_val = hyp;
      arg_min = v;
    }
  }
  
  min_labels[arg_min] = 1;
  return min_val;
}

/**********************************/

CardinalityDualFactorNodeBase::CardinalityDualFactorNodeBase(const Storage1D<DualVariableNode*>& participating_vars) 
  : DualFactorNode(participating_vars) {

  for (uint v=0; v < participating_vars.size(); v++) {
    if (participating_vars[v]->nLabels() != 2) {
      INTERNAL_ERROR << "instantiation of a cardinality factor with non-binary variables. Exiting." << std::endl;
      exit(1);
    }
  }
}

void CardinalityDualFactorNodeBase::update_duals(DualBCAMode mode, const Math1D::Vector<float>& card_cost) {

  const uint nVars = participating_var_.size();

  NamedStorage1D<Math1D::Vector<double> > msg(nVars, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    if (mode == DUAL_BCA_MODE_MPLP)     
      participating_var_[v]->compute_message(this, msg[v]);
    else
      msg[v].resize(2);
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }
  
  if (mode == DUAL_BCA_MODE_MPLP) {

    const double inv_nVars = 1.0 / nVars;

    //TEMPORARY
    // for (uint idx=0; idx < nVars; idx++) {

    //   Math1D::NamedVector<double> rel_param(nVars-1,MAKENAME(rel_param));
    
    //   double offs = 0.0;
    
    //   uint next = 0;
    //   for (uint k=0; k < nVars; k++) {
    
    //     if (k != idx) {
    // 	rel_param[next] = (msg[k][1]) - (msg[k][0]);
    // 	offs += msg[k][0];
    
    // 	next++;
    //     }
    //   }
    
    //   std::sort(rel_param.direct_access(), rel_param.direct_access() + nVars-1);
  
    //   double cum_sum = 0.0;
    
    //   Math1D::Vector<double> message(2,1e300);
    
    //   for (uint c=0; c < nVars; c++) {
    
    //     double hyp0 = cum_sum + card_cost[c];
    //     if (hyp0 < message[0])
    // 	message[0] = hyp0;
    
    //     double hyp1 = cum_sum + card_cost[c+1];
    //     if (hyp1 < message[1])
    // 	message[1] = hyp1;
    
    //     if (c+1 < nVars) 
    // 	cum_sum += rel_param[c];
    //   }
    
    //   for (uint l=0; l < 2; l++)
    //     message[l] += offs + msg[idx][l];
    
    //   for (uint l=0; l < 2; l++) {
    
    //     dual_ptr[idx][l] = (message[l] / double(nVars)) - msg[idx][l]; 
    //   }
    // }
    //END_TEMPORARY


    Math1D::Vector<std::pair<double,uint> > rel_msg(nVars);

    double offs = 0.0;
    
    for (uint k=0; k < nVars; k++) {
      
      rel_msg[k].first = msg[k][1] - msg[k][0];
      rel_msg[k].second = k;
      offs += msg[k][0];
    }
    
    std::sort(rel_msg.direct_access(), rel_msg.direct_access() + nVars);
    
    Math1D::Vector<uint> order(nVars);
    for (uint k=0; k < nVars; k++) {
      order[rel_msg[k].second] = k + 1;
    }
    
    Math1D::Vector<double> cum_sum(nVars+1);
    cum_sum[0] = rel_msg[0].first;
    for (uint k=1; k < nVars; k++) {
      assert(rel_msg[k].first >= rel_msg[k-1].first);
      cum_sum[k] = cum_sum[k-1] + rel_msg[k].first;
    }

    Math1D::Vector<double> cum_best(nVars + 1);
    cum_best[0] = card_cost[0]; 
    
    for (uint k=1; k <= nVars; k++) {
      
      double hyp = card_cost[k] + cum_sum[k-1];
      cum_best[k] = std::min(hyp, cum_best[k-1]);
    }
    
    
    Math1D::Vector<double> cum_best_m1(nVars+1);
    cum_best_m1[0] = 1e300;
    cum_best_m1[1] = card_cost[1];
    
    for (uint k=2; k <= nVars; k++) {
      cum_best_m1[k] = std::min(cum_best_m1[k-1], card_cost[k] + cum_sum[k-2]  ); 
    }
    
    Math1D::Vector<double> rev_cum_best(nVars + 1);
    rev_cum_best[nVars] = card_cost[nVars] + cum_sum[nVars-1];
    for (int k=nVars-1; k >= 1; k--) {
      
      double hyp = card_cost[k] + cum_sum[k-1];
      rev_cum_best[k] = std::min(hyp, rev_cum_best[k+1]);
    }
    rev_cum_best[0] = 1e300;
    
    Math1D::Vector<double> rev_cum_best_p1(nVars+1);
    rev_cum_best_p1[nVars] = 1e300;
    for (int k=nVars-1; k >= 0; k--) {
      rev_cum_best_p1[k] = std::min(rev_cum_best_p1[k+1], card_cost[k] + cum_sum[k]);
    }
    
    for (uint k=0; k < nVars; k++) {

      double* cur_dp = dual_ptr[k];
      const Math1D::Vector<double>& cur_msg = msg[k];
      
      const uint cur_order = order[k];
      
      const double cur_rel_msg = cur_msg[1] - cur_msg[0];
      
      const double val0 = std::min(cum_best[cur_order-1], rev_cum_best_p1[cur_order] - cur_rel_msg);
      cur_dp[0] = val0 + offs;

      const double val1 = std::min(rev_cum_best[cur_order], cum_best_m1[cur_order-1] + cur_rel_msg );
      cur_dp[1] = val1 + offs;
      
      //DEBUG
#if 0
      //a) check msg0
      double best0 = card_cost[0];
      uint arg_best0 = 0;
      for (uint l=1; l < nVars; l++) {
        
        double hyp = card_cost_[l];
        if (cur_order > l)
          hyp += cum_sum[l-1];
        else {
          hyp += cum_sum[l] - cur_rel_msg;
        }
        
        if (hyp < best0) {
          best0 = hyp;
          arg_best0 = l;
        }
      }
      if ( ! (fabs(best0-val0) < 1e-2)) {
        std::cerr << "best0: " << best0 << ", msg0: " << val0 << ", arg best: "  << arg_best0 
                  << ", order: " << cur_order << ", factor size: " << nVars << std::endl;
        
        std::cerr << "cum_best value: " << cum_best[cur_order-1] << std::endl;
        
        std::cerr << "p1 value: " << (rev_cum_best_p1[cur_order] - cur_rel_msg) << std::endl;
        
        std::cerr << "card_cost[arg0]: " << card_cost_[arg_best0] << std::endl;
        std::cerr << "cum_sum[arg0]: " << cum_sum[arg_best0] << std::endl;
        std::cerr << "rel msg: " << cur_rel_msg << std::endl;
      }
      assert(fabs(best0-val0) < 1e-2);
      
      //b) check msg1
      double best1 = 1e300;
      uint arg_best1 = 0;
      for (uint l=1; l <= nVars; l++) {
        
        double hyp = card_cost[l];
        if (cur_order <= l)
          hyp += cum_sum[l-1] - cur_rel_msg;
        else {
          if (l >= 2)
            hyp += cum_sum[l-2];
        }

        if (hyp < best1) {
          best1 = hyp;
          arg_best1 = l;
        }
      }
      
      if ( ! (fabs(best1-val1) < 1e-2)) {

        std::cerr << "best1: " << best1 << ", msg1: " << val1 << ", arg_best1: " << arg_best1 
                  << ", factor size: " << nVars << std::endl;
        std::cerr << "rev_cum_best: " << rev_cum_best[cur_order] << std::endl;
        std::cerr << "cum_best_m1 value: " << cum_best_m1[cur_order-1] << std::endl;
        
        
        std::cerr << "cur_order: " << cur_order << std::endl;
        std::cerr << "cur_rel_msg: " << cur_rel_msg << std::endl;
      }
      
      assert(fabs(best1-val1) < 1e-2);
#endif
      //END_DEBUG

      for (uint i=0; i < 2; i++) {
        cur_dp[i] *= inv_nVars;
        cur_dp[i] -= cur_msg[i];
      }
    }
  }
  else { //MSD

    assert(nVars >= 2);

    Storage1D<std::pair<double,uint> > values(nVars);

    double offs = 0.0;

    for (uint v=0; v < nVars; v++) {
      const double* cur_dp = dual_ptr[v];
      values[v] = std::make_pair(cur_dp[0]-cur_dp[1],v);
      offs -= cur_dp[0];
    }

    //std::cerr << "offs: " << offs << std::endl;

    std::sort(values.direct_access(),values.direct_access()+nVars);

    Storage1D<uint> order(nVars);
    for (uint v=0; v < nVars; v++) {
      order[values[v].second] = v;
    }

    for (uint v=0; v < nVars; v++) {

      Math1D::Vector<double>& cur_msg = msg[v];
      double* cur_dp = dual_ptr[v];      

      participating_var_[v]->compute_message(this, cur_msg);

      offs += cur_dp[0]; //the current dual vars don't enter!

      uint cur_order = order[v];

      double best_zero = card_cost[0];
      double cum_sum = 0.0;
      for (uint c=1; c < nVars; c++) {

        if (cur_order < c)
          cum_sum += values[c].first;
        else
          cum_sum += values[c-1].first;

        double hyp = cum_sum + card_cost[c];
        if (hyp < best_zero)
          best_zero  = hyp;
      }


      cum_sum = 0.0; // the current var doesn't enter!! 
      double best_one = card_cost[1];
      for (uint c=2; c <= nVars; c++) {
        if (cur_order < c-1)
          cum_sum += values[c-1].first;
        else
          cum_sum += values[c-2].first;
        
        double hyp = cum_sum + card_cost[c];
        if (hyp < best_one)
          best_one  = hyp;
      }

      cur_dp[1] = 0.5 * (best_one + offs - cur_msg[1]);
      cur_dp[0] = 0.5 * (best_zero + offs - cur_msg[0]);
        
      offs -= cur_dp[0];

      //update values and order
      values[cur_order].first = cur_dp[0]-cur_dp[1];
      double val = values[cur_order].first;

#if 0
      while(cur_order > 0 && values[cur_order-1].first > val) {
        uint other_var = values[cur_order-1].second;
        //assert(order[other_var] == cur_order-1);
        std::swap(values[cur_order-1],values[cur_order]);
        order[other_var]++;
        cur_order--;
      }
      while (cur_order+1 < nVars && values[cur_order+1].first < val) {
        uint other_var = values[cur_order+1].second;
        //assert(order[other_var] == cur_order+1);
        std::swap(values[cur_order+1],values[cur_order]);
        order[other_var]--;
        cur_order++;
      }
#else
      uint k=cur_order;
      while (k>0 && values[k-1].first > val)
        k--;
      if (k != cur_order) {
        for (uint kk=cur_order; kk > k; kk--) {
          values[kk] = values[kk-1];
          order[values[kk].second] = kk;
        }
      }
      else {
        while (k+1 < nVars && values[k+1].first < val)
          k++;
        for (uint kk=cur_order; kk < k; kk++) {
          values[kk] = values[kk+1];
          order[values[kk].second] = kk;
        }
      }
      values[k].first = val;
      values[k].second = v;
#endif
      //order[v] = cur_order; //not that we will ever need this again....
    }
  }

}

double CardinalityDualFactorNodeBase::dual_value(const Math1D::Vector<float>& card_cost) const {

  const uint nVars = participating_var_.size();

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  Math1D::Vector<double> rel_msg(nVars);

  double offs = 0.0;

  for (uint k=0; k < nVars; k++) {

    const double val0 = - dual_ptr[k][0];
    const double val1 = - dual_ptr[k][1];

    rel_msg[k] = val1 - val0;
    offs += val0;
  }

  std::sort(rel_msg.direct_access(), rel_msg.direct_access() + nVars);

  double min_val = card_cost[0];

  double cum_sum = 0.0;
  for (uint c=1; c <= nVars; c++) {
    cum_sum += rel_msg[c-1];

    const double hyp = cum_sum + card_cost[c];

    if (hyp < min_val) {
      min_val = hyp;
    }
  }

  return min_val + offs;
}

double CardinalityDualFactorNodeBase::compute_minimizer(Math1D::Vector<uint>& min_labels, 
                                                        const Math1D::Vector<float>& card_cost) const {


  const uint nVars = participating_var_.size();

  min_labels.resize_dirty(nVars);

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  Storage1D<std::pair<double,uint> > rel_msg(nVars);

  double offs = 0.0;

  for (uint k=0; k < nVars; k++) {

    const double val0 = - dual_ptr[k][0];
    const double val1 = - dual_ptr[k][1];

    rel_msg[k] = std::make_pair(val1 - val0,k);
    offs += val0;
  }

  std::sort(rel_msg.direct_access(), rel_msg.direct_access() + nVars);

  double min_val = card_cost[0];
  uint arg_min = 0;
  
  double cum_sum = 0.0;
  for (uint c=1; c <= nVars; c++) {
    cum_sum += rel_msg[c-1].first;
    
    const double hyp = cum_sum + card_cost[c];

    if (hyp < min_val) {
      min_val = hyp;
      arg_min = c;
    }
  }

  min_labels.set_constant(0);
  for (uint k=0; k < arg_min; k++)
    min_labels[rel_msg[k].second] = 1;

  return min_val + offs;
}

/**********************************/

CardinalityDualFactorNode::CardinalityDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars,
                                                     const Math1D::Vector<float>& card_cost) :
  CardinalityDualFactorNodeBase(participating_vars), card_cost_(card_cost) {

  if (card_cost_.size() < participating_vars.size()+1) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
    exit(1);
  }
}

/*virtual*/ CardinalityDualFactorNode::~CardinalityDualFactorNode() {}

/*virtual*/ double CardinalityDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  
  uint sum = labels.sum();
  return card_cost_[sum];
}

/*virtual*/ void CardinalityDualFactorNode::update_duals(DualBCAMode mode) {

  CardinalityDualFactorNodeBase::update_duals(mode,card_cost_);
}

/*virtual*/ double CardinalityDualFactorNode::dual_value() const {

  return CardinalityDualFactorNodeBase::dual_value(card_cost_);
}

/*virtual*/ double CardinalityDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {

  return CardinalityDualFactorNodeBase::compute_minimizer(min_labels,card_cost_);
}

/*virtual*/ double CardinalityDualFactorRefNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {

  return CardinalityDualFactorNodeBase::compute_minimizer(min_labels,card_cost_);
}

/**********************************/

CardinalityDualFactorRefNode::CardinalityDualFactorRefNode(const Storage1D<DualVariableNode*>& participating_vars,
                                                           const Math1D::Vector<float>& card_cost) :
  CardinalityDualFactorNodeBase(participating_vars), card_cost_(card_cost) {
}

/*virtual*/ CardinalityDualFactorRefNode::~CardinalityDualFactorRefNode() {}

/*virtual*/ double CardinalityDualFactorRefNode::cost(const Math1D::Vector<uint>& labels) const {
  
  uint sum = labels.sum();
  return card_cost_[sum];
}

/*virtual*/ void CardinalityDualFactorRefNode::update_duals(DualBCAMode mode) {

  assert(card_cost_.size() >= participating_var_.size()+1);

  CardinalityDualFactorNodeBase::update_duals(mode,card_cost_);
}

/*virtual*/ double CardinalityDualFactorRefNode::dual_value() const {

  return CardinalityDualFactorNodeBase::dual_value(card_cost_);
}

/**********************************/

AllPosBILPConstraintDualFactorNode::AllPosBILPConstraintDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars,
                                                                       int rhs_lower, int rhs_upper) 
  : DualFactorNode(participating_vars), rhs_lower_(std::max(0,rhs_lower)), rhs_upper_(std::min<int>(participating_vars.size(),rhs_upper)) {

  if (rhs_lower_ > rhs_upper_ || rhs_upper < 0) {
    INTERNAL_ERROR << "constraint is unsatisfiable, so inference is pointless. Exiting." << std::endl;
    exit(1);
  }

  for (uint v=0; v < participating_vars.size(); v++) {
    if (participating_vars[v]->nLabels() != 2) {
      INTERNAL_ERROR << "instantiation of an AllPosBILP factor with non-binary variables. Exiting." << std::endl;
      exit(1);
    }
  }
}

/*virtual*/ AllPosBILPConstraintDualFactorNode::~AllPosBILPConstraintDualFactorNode() {}

/*virtual*/ double AllPosBILPConstraintDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {
  int sum = labels.sum();

  return (sum >= rhs_lower_ && sum <= rhs_upper_) ? 0.0 : 1e20;
}

/*virtual*/ 
void AllPosBILPConstraintDualFactorNode::update_duals(DualBCAMode mode) {

  const uint nVars = participating_var_.size();

  NamedStorage1D<Math1D::Vector<double> > msg(nVars, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    if (mode == DUAL_BCA_MODE_MPLP)     
      participating_var_[v]->compute_message(this, msg[v]);
    else
      msg[v].resize(2);
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }


  if (mode == DUAL_BCA_MODE_MPLP) {

    const double inv_nVars = 1.0 / nVars;

    Math1D::Vector<std::pair<double,uint> > rel_msg(nVars);

    double offs = 0.0;
    
    for (uint k=0; k < nVars; k++) {
      
      rel_msg[k].first = msg[k][1] - msg[k][0];
      rel_msg[k].second = k;
      offs += msg[k][0];
    }
    
    std::sort(rel_msg.direct_access(), rel_msg.direct_access() + nVars);
    
    Math1D::Vector<uint> order(nVars);
    for (uint k=0; k < nVars; k++) {
      order[rel_msg[k].second] = k + 1;
    }
    
    Math1D::Vector<double> cum_sum(nVars+1);
    cum_sum[0] = rel_msg[0].first;
    for (uint k=1; k < nVars; k++) {
      assert(rel_msg[k].first >= rel_msg[k-1].first);
      cum_sum[k] = cum_sum[k-1] + rel_msg[k].first;
    }

    Math1D::Vector<double> cum_best(nVars + 1,1e20);
    cum_best[0] = (rhs_lower_ == 0) ? 0.0 : 1e20;

    for (int k=1; k <= int(nVars); k++) {

      if (rhs_lower_ <= k && rhs_upper_ >= k) {
        double hyp = cum_sum[k-1];
        cum_best[k] = std::min(hyp, cum_best[k-1]);
      }
      else
        cum_best[k] = cum_best[k-1];
    }
    
    Math1D::Vector<double> cum_best_m1(nVars+1);
    cum_best_m1[0] = 1e300;
    cum_best_m1[1] = (rhs_lower_ <= 1 && rhs_upper_ >= 1) ? 0.0 : 1e20;
    
    for (int k=2; k <= int(nVars); k++) {
      if (k >= rhs_lower_ && k <= rhs_upper_)
        cum_best_m1[k] = std::min(cum_best_m1[k-1], cum_sum[k-2]  ); 
      else
        cum_best_m1[k] = cum_best_m1[k-1];
    }
    
    Math1D::Vector<double> rev_cum_best(nVars + 1,1e20);
    rev_cum_best[nVars] = (rhs_upper_ >= int(nVars)) ? cum_sum[nVars-1] : 1e20;
    for (int k=nVars-1; k >= 1; k--) {
      
      double hyp = (k >= rhs_lower_ && k <= rhs_upper_) ? cum_sum[k-1] : 0.0;
      rev_cum_best[k] = std::min(hyp, rev_cum_best[k+1]);
    }
    rev_cum_best[0] = 1e300;
    
    Math1D::Vector<double> rev_cum_best_p1(nVars+1);
    rev_cum_best_p1[nVars] = 1e300;
    for (int k=nVars-1; k >= 0; k--) {
      if (k >= rhs_lower_ && k <= rhs_upper_)
        rev_cum_best_p1[k] = std::min(rev_cum_best_p1[k+1], cum_sum[k]);
      else
        rev_cum_best_p1[k] = rev_cum_best_p1[k+1];
    }
    
    for (uint k=0; k < nVars; k++) {

      double* cur_dp = dual_ptr[k];
      const Math1D::Vector<double>& cur_msg = msg[k];
      
      const uint cur_order = order[k];
      
      const double cur_rel_msg = cur_msg[1] - cur_msg[0];
      
      const double val0 = std::min(cum_best[cur_order-1], rev_cum_best_p1[cur_order] - cur_rel_msg);
      cur_dp[0] = val0 + offs;

      const double val1 = std::min(rev_cum_best[cur_order], cum_best_m1[cur_order-1] + cur_rel_msg );
      cur_dp[1] = val1 + offs;
      

      for (uint i=0; i < 2; i++) {
        cur_dp[i] *= inv_nVars;
        cur_dp[i] -= cur_msg[i];
      }
    }


  }
  else { //MSD

    assert(nVars >= 2);

    Storage1D<std::pair<double,uint> > values(nVars);

    double offs = 0.0;

    for (uint v=0; v < nVars; v++) {
      const double* cur_dp = dual_ptr[v];
      values[v] = std::make_pair(cur_dp[0]-cur_dp[1],v);
      offs -= cur_dp[0];
    }

    std::sort(values.direct_access(),values.direct_access()+nVars);

    Storage1D<uint> order(nVars);
    for (uint v=0; v < nVars; v++) {
      order[values[v].second] = v;
    }

    for (uint v=0; v < nVars; v++) {

      Math1D::Vector<double>& cur_msg = msg[v];
      double* cur_dp = dual_ptr[v];      

      participating_var_[v]->compute_message(this, cur_msg);

      offs += cur_dp[0]; //the current dual vars don't enter!

      uint cur_order = order[v];

      double best_zero = (rhs_lower_ == 0) ? 0.0 : 1e20;
      double cum_sum = 0.0;
      for (uint c=1; c < std::min<uint>(nVars,rhs_upper_+1); c++) {

        if (cur_order < c)
          cum_sum += values[c].first;
        else
          cum_sum += values[c-1].first;

        if (c >= uint(rhs_lower_)) {
          if (cum_sum < best_zero)
            best_zero  = cum_sum;
        }
      }


      cum_sum = 0.0; // the current var doesn't enter!! 
      double best_one = (rhs_lower_ <= 1 && rhs_upper_ >= 1) ? 0.0 : 1e20;
      for (uint c=2; c <= std::min<uint>(rhs_upper_,nVars); c++) {
        if (cur_order < c-1)
          cum_sum += values[c-1].first;
        else
          cum_sum += values[c-2].first;
        
        if (uint(rhs_lower_) <= c) {
          if (cum_sum < best_one)
            best_one  = cum_sum;
        }
      }

      cur_dp[1] = 0.5 * (best_one + offs - cur_msg[1]);
      cur_dp[0] = 0.5 * (best_zero + offs - cur_msg[0]);
        
      offs -= cur_dp[0];

      //update values and order
      values[cur_order].first = cur_dp[0]-cur_dp[1];
      double val = values[cur_order].first;

      uint k=cur_order;
      while (k>0 && values[k-1].first > val)
        k--;
      if (k != cur_order) {
        for (uint kk=cur_order; kk > k; kk--) {
          values[kk] = values[kk-1];
          order[values[kk].second] = kk;
        }
      }
      else {
        while (k+1 < nVars && values[k+1].first < val)
          k++;
        for (uint kk=cur_order; kk < k; kk++) {
          values[kk] = values[kk+1];
          order[values[kk].second] = kk;
        }
      }
      values[k].first = val;
      values[k].second = v;
      //order[v] = cur_order; //not that we will ever need this again....
    }

  }
}

/*virtual*/ 
double AllPosBILPConstraintDualFactorNode::dual_value() const {

  const uint nVars = participating_var_.size();

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  Math1D::Vector<double> rel_msg(nVars);

  double offs = 0.0;

  for (uint k=0; k < nVars; k++) {

    const double val0 = - dual_ptr[k][0];
    const double val1 = - dual_ptr[k][1];

    rel_msg[k] = val1 - val0;
    offs += val0;
  }

  std::sort(rel_msg.direct_access(), rel_msg.direct_access() + nVars);

  double min_val = (rhs_lower_ == 0) ? 0.0 : 1e20;

  double cum_sum = 0.0;
  for (int c=1; c <= rhs_upper_; c++) {
    cum_sum += rel_msg[c-1];

    if (c >= rhs_lower_) {
      const double hyp = cum_sum;

      if (hyp < min_val) {
        min_val = hyp;
      }
    }
  }

  return min_val + offs;
}

/*virtual*/ 
double AllPosBILPConstraintDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {

  const uint nVars = participating_var_.size();

  min_labels.resize_dirty(nVars);

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  Storage1D<std::pair<double,uint> > rel_msg(nVars);

  double offs = 0.0;

  for (uint k=0; k < nVars; k++) {

    const double val0 = - dual_ptr[k][0];
    const double val1 = - dual_ptr[k][1];

    rel_msg[k] = std::make_pair(val1 - val0,k);
    offs += val0;
  }

  std::sort(rel_msg.direct_access(), rel_msg.direct_access() + nVars);

  double min_val = (rhs_lower_ == 0) ? 0.0 : 1e20;
  uint arg_min = 0;
  
  double cum_sum = 0.0;
  for (int c=1; c <= rhs_upper_; c++) {
    cum_sum += rel_msg[c-1].first;
    
    if (c >= rhs_lower_) {
      const double hyp = cum_sum;
      
      if (hyp < min_val) {
        min_val = hyp;
        arg_min = c;
      }
    }
  }

  min_labels.set_constant(0);
  for (uint k=0; k < arg_min; k++)
    min_labels[rel_msg[k].second] = 1;

  return min_val + offs;
}

/**********************************/


BILPConstraintDualFactorNode::BILPConstraintDualFactorNode(const Storage1D<DualVariableNode*>& participating_vars,
                                                           const Storage1D<bool>& positive, int rhs_lower, int rhs_upper) :
  DualFactorNode(participating_vars), rhs_lower_(rhs_lower), rhs_upper_(rhs_upper) {

  for (uint v=0; v < participating_vars.size(); v++) {
    if (participating_vars[v]->nLabels() != 2) {
      INTERNAL_ERROR << "instantiation of a BILP factor with non-binary variables. Exiting." << std::endl;
      exit(1);
    }
  }

  if (positive.size() < participating_vars.size()) {
    INTERNAL_ERROR << "dimension mismatch" << std::endl;
    exit(1);
  }

  if (rhs_lower_ > rhs_upper_) {
    INTERNAL_ERROR << "constraint is unsatisfiable, so inference is pointless. Exiting." << std::endl;
    exit(1);
  }


  Storage1D<DualVariableNode*> sorted_involved_vars(participating_vars.size());
  uint next_pos = 0;

  //pass 1 - find all positive
  for (uint v=0; v < participating_vars.size(); v++) {
    if (positive[v]) {
      sorted_involved_vars[next_pos] = participating_vars[v];
      next_pos++;
    } 
  }

  nPos_ = next_pos;

  //pass 2 - find all negative
  for (uint v=0; v < participating_vars.size(); v++) {
    if (!positive[v]) {
      sorted_involved_vars[next_pos] = participating_vars[v];
      next_pos++;
    }
  }

  participating_var_ = sorted_involved_vars;

  int nPositive = nPos_;
  int nNegative = participating_var_.size()-nPositive;

  int lower_bound = -nNegative;
  int upper_bound = nPositive;

  lower_bound = std::max(lower_bound, rhs_lower_ - nPositive);
  upper_bound = std::min(upper_bound, rhs_upper_ + nNegative);

  const int range = upper_bound - lower_bound + 1;
  const int zero_offset = -lower_bound;

  if (rhs_upper_ + zero_offset < 0 || rhs_lower + zero_offset >= range) {
    INTERNAL_ERROR << "constraint is unsatisfiable. Exiting..." << std::endl;
    exit(1);
  }

  //adjust the bounds to the actually possible range of values
  if (rhs_lower_ + zero_offset < 0) {
    rhs_lower_ -= (rhs_lower_ + zero_offset);
  }
  if (rhs_upper_ + zero_offset >= range) {
    rhs_upper_ -= (rhs_upper_ + zero_offset - range +1);
  }
}

/*virtual*/ BILPConstraintDualFactorNode::~BILPConstraintDualFactorNode() {}

/*virtual*/ double BILPConstraintDualFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  int sum = 0;

  const uint nVars = participating_var_.size();

  for (uint k=0; k < nPos_; k++)
    sum += labels[k];
  for (uint k=nPos_; k < nVars; k++)
    sum -= labels[k];

  return (sum >= rhs_lower_ && sum <= rhs_upper_) ? 0.0 : 1e15;
}

/*virtual*/ void BILPConstraintDualFactorNode::update_duals(DualBCAMode mode) {

  const uint nVars = participating_var_.size();

  NamedStorage1D<Math1D::Vector<double> > msg(nVars, MAKENAME(msg));

  NamedStorage1D<double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {

    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);

    if (mode == DUAL_BCA_MODE_MPLP)
      participating_var_[v]->compute_message(this, msg[v]);
    else
      msg[v].resize(participating_var_[v]->nLabels());
  }

  int nPositive = nPos_;
  int nNegative = nVars - nPositive;

  int lower_bound = -nNegative;
  int upper_bound = nPositive;

  //some intermediate values can never result in a satisfied constraint 
  // => prune the range to the useful values
  lower_bound = std::max(lower_bound, rhs_lower_ - nPositive);
  upper_bound = std::min(upper_bound, rhs_upper_ + nNegative);

  const int range = upper_bound - lower_bound + 1;
  const int zero_offset = -lower_bound;

  Math2D::Matrix<double> backward_light(range,nVars);

  const uint last_var = nVars-1;

  if (mode == DUAL_BCA_MODE_MPLP) {

    //MPLP

    const double inv_nVars = 1.0 / nVars;

    /**** backward ****/

    //init
    for (int sum=0; sum < range; sum++) 
      backward_light(sum,last_var) = 1e100;
    
    backward_light(zero_offset,last_var) = msg[last_var][0];
    const int end_val = zero_offset + ((last_var < nPos_) ? 1 : -1);
    
    if (end_val >= 0 && end_val < range ) {
      backward_light(end_val,last_var) = msg[last_var][1];
    }
    
    //proceed
    for (int v=last_var-1; v >= int(nPos_); v--) {
	
      const Math1D::Vector<double>& cur_msg = msg[v];

      for (int sum=zero_offset+1; sum < range; sum++)
	backward_light(sum,v) = 1e100;
      
      for (int sum=0; sum <= zero_offset; sum++) {
          
	double best_prev = 1e75;
          
	for (int l=0; l < 2; l++) {
	  
	  const int dest = sum +l;
	  double hyp = 1e75;
	  if (dest < range) {
	    hyp = backward_light(dest,v+1) + cur_msg[l];
	  }
	  
	  if (hyp < best_prev)
	    best_prev = hyp;
	}
	
	backward_light(sum,v) = best_prev;
      }

    }
    for (int v=std::min<int>(last_var-1,int(nPos_)-1); v >= 1; v--) {
      
      const Math1D::Vector<double>& cur_msg = msg[v];
      
      for (int sum=0; sum < range; sum++) {
	
	double best_prev = 1e75;
        
	for (int l=0; l < 2; l++) {
            
	  const int dest = sum - l;
	  double hyp = 1e75;
	  if (dest >= 0) {
	    hyp = backward_light(dest,v+1) + cur_msg[l];
	  }
	  
	  if (hyp < best_prev)
	    best_prev = hyp;
	}
	
	backward_light(sum,v) = best_prev;
      }
    }
    
    /**** jointly calculate forward and update the messages ****/

    Math2D::Matrix<double> temp_forward(2,range,1e100);
    Math1D::Vector<double> forward_light[2];
    forward_light[0].resize(range,1e100);
    forward_light[1].resize(range);
    
    uint cur_idx = 0;

    for (uint v=0; v < nVars; v++) {
        
      double* cur_dp = dual_ptr[v];
      const Math1D::Vector<double>& cur_msg = msg[v];
      
      if (v == 0) {
	
	temp_forward.set_constant(1e100);
	
	forward_light[0][zero_offset] = cur_msg[0];
	temp_forward(0,zero_offset) = cur_msg[0];
	
	if (nPos_ > 0) {
	  forward_light[0][zero_offset+1] = cur_msg[1];
	  temp_forward(1,zero_offset+1) = cur_msg[1];
	}
	else {
	  forward_light[0][zero_offset-1] = cur_msg[1];
	  temp_forward(1,zero_offset-1) = cur_msg[1];
	}
      }
      else {
	
	const Math1D::Vector<double>& prev_forward_light = forward_light[cur_idx];
	cur_idx = 1-cur_idx;
	
	Math1D::Vector<double>& cur_forward_light = forward_light[cur_idx];

	for (int k=0; k < range; k++) {
	  
	  for (uint l=0; l < 2; l++) {
	    
	    double best = 1e75;
	    
	    int dest = k + ((v < nPos_) ? (-l) : l);
	    if (dest >= 0 && dest < range)
	      best = prev_forward_light[dest] + cur_msg[l];
	    
	    temp_forward(l,k) = best;
	  }
	  
	  cur_forward_light[k] = std::min(temp_forward(0,k),temp_forward(1,k));
	}
	
      }
      
      for (uint l=0; l < 2; l++) {

	double min_msg = 1e300;
        
	for (int s=0; s < (int) range; s++) {
            
	  double hyp = temp_forward(l,s);

	  if (v+1 < nVars) {
	    
	    double best_bwd = 1e300;
            
	    const int diff = (s - zero_offset);
            
	    for (int r=rhs_lower_; r <= rhs_upper_; r++) {
	      const int other = r + zero_offset - diff; 
              
	      if (other >= 0 && other < (int) range) {
		
		best_bwd = std::min(best_bwd,backward_light(other,v+1));
	      }
	    }
            
	    hyp += best_bwd;
	  }
	  else {
	    if (s < int(rhs_lower_ + zero_offset) || s > int(rhs_upper_ + zero_offset)) 
	      hyp = 1e300;
	  }
	  
	  if (hyp < min_msg)
	    min_msg = hyp;
	}
          
	assert(!isnan(min_msg));
	
	cur_dp[l] = min_msg * inv_nVars - cur_msg[l];
      }
    }
  }
  else {
    //min sum diffusion (MSD) mode

    /**** solution based on forward-backward ****/
    
    //a) compute backward completely (and just once)
    //init
    for (int sum=0; sum < range; sum++) 
      backward_light(sum,last_var) = 1e100;
    
    backward_light(zero_offset,last_var) = -dual_ptr[last_var][0];
    const int end_val = zero_offset + ((last_var < nPos_) ? 1 : -1);
      
    if (end_val >= 0 && end_val < range ) 
      backward_light(end_val,last_var) = -dual_ptr[last_var][1];
    
    //proceed
    for (int v=last_var-1; v >= int(nPos_); v--) {
      
      const double* cur_dp = dual_ptr[v];
      
      for (int sum=zero_offset+1; sum < range; sum++)
	backward_light(sum,v) = 1e100;
      
      for (int sum=0; sum <= zero_offset; sum++) {
	
	double best_prev = 1e75;
        
	for (int l=0; l < 2; l++) {
	  
	  const int dest = sum + l;
	  double hyp = 1e75;
	  if (dest < range) 
	    hyp = backward_light(dest,v+1) - cur_dp[l]; 
	  
	  if (hyp < best_prev)
	    best_prev = hyp;
	}
	
	backward_light(sum,v) = best_prev;
      }
    }
    
    for (int v=std::min<int>(nPos_-1,last_var-1); v >= 0; v--) {

      const double* cur_dp = dual_ptr[v];
        
      for (int sum=0; sum < range; sum++) {
          
	double best_prev = 1e75;
        
	for (int l=0; l < 2; l++) {
            
	  const int dest = sum - l;
	  double hyp = 1e75;
	  if (dest >= 0) 
	    hyp = backward_light(dest,v+1) - cur_dp[l]; 
	  
	  if (hyp < best_prev)
	    best_prev = hyp;
	}
	
	backward_light(sum,v) = best_prev;
      }
    }
    
    //b) compute forward incrementally and derive messages

    Math2D::Matrix<double> temp_forward(range,2,1e100);
    Math1D::Vector<double> forward_light[2];
    forward_light[0].resize(range,1e100);
    forward_light[1].resize(range);
    
    uint cur_idx = 0;
      
    for (uint v=0; v < nVars; v++) {
      
      const double* cur_dp = dual_ptr[v];
      
      Math1D::Vector<double>& cur_msg = msg[v];
      
      participating_var_[v]->compute_message(this, cur_msg);
      
      if (v == 0) {

	temp_forward(zero_offset,0) = 0.0;
	const int init_val = zero_offset + ( (nPos_ > 0) ? 1 : -1);
	if (init_val >= 0 && init_val < range) { 
	  temp_forward(init_val,1) = 0.0;
	}
      }
      else {
	
	const Math1D::Vector<double>& prev_forward_light = forward_light[cur_idx];
	cur_idx = 1-cur_idx;
	
	for (int sum=0; sum < range; sum++) {
	    
	  for (int l=0; l < 2; l++) {
	    
	    double best_prev = 1e75;
            
	    int move = l;
	    if (v < nPos_)  //since we are tracing backward here
	      move *= -1;
	    
	    const int dest = sum + move;
	    if (dest >= 0 && dest < range) {
	      
	      best_prev = prev_forward_light[dest]; 
	    }
            
	    temp_forward(sum,l) = best_prev;
	  }
	}
      }

      //now compute new duals
      for (uint l=0; l < 2; l++) {
	
	double min_msg = 1e300;
	  
	for (int s=0; s < (int) range; s++) {
	  
	  double hyp = temp_forward(s,l);
          
	  if (v+1 < nVars) {
	    
	    double best_bwd = 1e300;
            
	    const int diff = (s - zero_offset);
            
	    for (int r=rhs_lower_; r <= rhs_upper_; r++) {
	      const int other = r + zero_offset - diff; 
              
	      if (other >= 0 && other < (int) range) {
		
		best_bwd = std::min(best_bwd,backward_light(other,v+1));
	      }
	    }
	    
	    hyp += best_bwd;
	  }
	  else {
	    if (s < int(rhs_lower_ + zero_offset) || s > int(rhs_upper_ + zero_offset)) 
	      hyp = 1e300;
	  }
	  
	  if (hyp < min_msg)
	    min_msg = hyp;
	}
          
	assert(!isnan(min_msg));
        

	dual_ptr[v][l] = 0.5 * (min_msg - cur_msg[l]);
      }
      
      Math1D::Vector<double>& cur_forward_light = forward_light[cur_idx];

      //correct the freshly computed forward term
      for (int sum=0; sum < range; sum++) {
	cur_forward_light[sum] = std::min(temp_forward(sum,0)-cur_dp[0], temp_forward(sum,1)-cur_dp[1]);
      }
    }
  }
}

/*virtual*/ double BILPConstraintDualFactorNode::dual_value() const {

  const uint nVars = participating_var_.size();

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  int nPositive = nPos_;
  int nNegative = nVars - nPositive;

  int lower_bound = -nNegative;
  int upper_bound = nPositive;

  //some intermediate values can never result in a satisfied constraint 
  // => prune the range to the useful values
  lower_bound = std::max(lower_bound, rhs_lower_ - nPositive);
  upper_bound = std::min(upper_bound, rhs_upper_ + nNegative);

  const int range = upper_bound - lower_bound + 1;
  const int zero_offset = -lower_bound;

  /**** forward ****/

  Math1D::Vector<double> forward_light[2];
  forward_light[0].resize(range,1e100);
  forward_light[1].resize(range);

  //init
  forward_light[0][zero_offset] = - dual_ptr[0][0];
  int init_mul = ((nPos_ > 0) ? 1 : -1) + zero_offset;
  if (init_mul >= 0 && init_mul < range)
    forward_light[0][zero_offset+init_mul] = - dual_ptr[0][1];

  uint cur_idx = 0;

  //proceed
  for (uint v=1; v < nVars; v++) {

    const Math1D::Vector<double>& prev_forward_light = forward_light[cur_idx];

    cur_idx = 1-cur_idx;

    Math1D::Vector<double>& cur_forward_light = forward_light[cur_idx];

    const double* cur_dp = dual_ptr[v];

    for (int sum=0; sum < range; sum++) {

      double best_val = 1e75;

      for (int l=0; l < 2; l++) {
	
	double hyp = 1e75;

        int move = l;
	if (v < nPos_)
          move *= -1; //since we are tracing backward here

        const int dest = sum + move;
        if (dest >= 0 && dest < range) {
	  
	  hyp = prev_forward_light[dest] - cur_dp[l]; 
        }

        if (hyp < best_val)
          best_val = hyp;
      }
      cur_forward_light[sum] = best_val;
    }    
  }

  const Math1D::Vector<double>& last_forward_light = forward_light[cur_idx];
    
  double min_val = 1e300;

  for (int k=rhs_lower_; k <= rhs_upper_; k++) {

    const double min_msg = last_forward_light[zero_offset + k];
    
    if (min_msg < min_val)
      min_val = min_msg;
  }

  return min_val;
}

/*virtual*/ double BILPConstraintDualFactorNode::compute_minimizer(Math1D::Vector<uint>& min_labels) const {

  const uint nVars = participating_var_.size();

  NamedStorage1D<const double*> dual_ptr(nVars, MAKENAME(dual_ptr));

  for (uint v=0; v < nVars; v++) {
    dual_ptr[v] = participating_var_[v]->get_dual_vars(this);
  }

  int nPositive = nPos_;
  int nNegative = nVars - nPositive;

  int lower_bound = -nNegative;
  int upper_bound = nPositive;

  //some intermediate values can never result in a satisfied constraint 
  // => prune the range to the useful values
  lower_bound = std::max(lower_bound, rhs_lower_ - nPositive);
  upper_bound = std::min(upper_bound, rhs_upper_ + nNegative);

  const int range = upper_bound - lower_bound + 1;
  const int zero_offset = -lower_bound;

  /**** forward ****/

  Math3D::NamedTensor<double> forward(range,2,participating_var_.size(),1e100,MAKENAME(forward));
  Math3D::NamedTensor<uint> trace(range,2,participating_var_.size(),MAX_UINT,MAKENAME(trace));

  //init
  forward(zero_offset,0,0) = - dual_ptr[0][0];
  int init_mul = (nPos_ > 0) ? 1 : -1;
  forward(zero_offset+init_mul,1,0) = - dual_ptr[0][1];

  //proceed
  for (uint v=1; v < nVars; v++) {

    for (int sum=0; sum < range; sum++) {

      for (int l=0; l < 2; l++) {
	
        double best_prev = 1e75;
	
        int move = l;
	if (v < nPos_) //since we are tracing backward here
          move *= -1;

        const int dest = sum + move;
        if (dest >= 0 && dest < range) {
          for (int l_prev = 0; l_prev < 2; l_prev ++) {
            const double hyp = forward(dest,l_prev,v-1);
            if (hyp < best_prev) {
              best_prev = hyp;
              trace(sum,l,v) = l_prev;
            }
          }
        }

        forward(sum,l,v) = best_prev - dual_ptr[v][l];
      }
    }
  }

  double min_val = 1e300;

  uint best_l = MAX_UINT;
  uint best_k = MAX_UINT;

  for (uint l=0; l < 2; l++) {

    for (int k=rhs_lower_; k <= rhs_upper_; k++) {

      double min_msg = forward(zero_offset + k,l,nVars-1); 
    
      if (min_msg < min_val) {
        min_val = min_msg;
        best_k = k + zero_offset;
        best_l = l;
      }
    }
  }

  // traceback
  min_labels[nVars-1] = best_l;

  for (int v=nVars-2; v >= 0; v--) {

    uint prev_k = best_k;

    if (v+1 < int(nPos_))
      best_k -= best_l;
    else
      best_k += best_l;

    best_l = trace(prev_k,best_l,v+1);

    assert(best_l != MAX_UINT);

    min_labels[v] = best_l;
  }

  return min_val;
}

/**********************************/
 
FactorDualOpt::FactorDualOpt(uint nVars, uint nFactors) : nUsedVars_(0), nUsedFactors_(0) {

  var_.resize(nVars,0);
  factor_.resize(nFactors,0);
}

//delete all allocated owned var. and factor nodes
FactorDualOpt::~FactorDualOpt() {
  for (uint i=0; i < nUsedVars_; i++)
    delete var_[i];

  for (uint i=0; i < nUsedFactors_; i++)
    delete factor_[i];
}

/***** Setting up a factor graph ******/

uint FactorDualOpt::add_factor(DualFactorNode* node) {
      
  nUsedFactors_++;

  uint prev_size = factor_.size();
  if (nUsedFactors_ > prev_size)
    factor_.resize(size_t(1.2*prev_size)+1,0);

  factor_[nUsedFactors_-1] = node;

  return nUsedFactors_-1;
}

//return var. number
uint FactorDualOpt::add_var(const Math1D::Vector<float>& unary_cost) {

  DualVariableNode* ptr = new DualVariableNode(unary_cost);

  uint nPrevVars = nUsedVars_;

  nUsedVars_++;

  uint prev_size = var_.size();
  if (nUsedVars_ > prev_size) {
    var_.resize(size_t(1.2*prev_size)+1,0);
  }
  var_[nPrevVars] = ptr;

  return nPrevVars;
}

DualVariableNode* FactorDualOpt::get_var(uint var_num) {

  if (var_num >= nUsedVars_) {
    INTERNAL_ERROR << "variable index out of bounds. Exiting" << std::endl;
    exit(1);
  }

  return var_[var_num];
}

uint FactorDualOpt::add_generic_factor(const Math1D::Vector<uint> var, const VarDimStorage<float>& cost) {

  assert(var.size() == cost.nDims());

  Storage1D<DualVariableNode*> var_nodes(var.size());
  for (uint k=0; k < var.size(); k++) {
    if (var[k] >= nUsedVars_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    var_nodes[k] = var_[var[k]];
  }

  GenericDualFactorNode* ptr = new GenericDualFactorNode(var_nodes,cost);
  return add_factor(ptr);
}

uint FactorDualOpt::add_generic_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost, bool ref) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  Storage1D<DualVariableNode*> var_nodes(2);
  var_nodes[0] = var_[var1];
  var_nodes[1] = var_[var2];

  if (ref) {
    BinaryDualRefFactorNode* ptr = new BinaryDualRefFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
  else {
    BinaryDualFactorNode* ptr = new BinaryDualFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
}
  
uint FactorDualOpt::add_potts_factor(uint var1, uint var2, double lambda) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  Storage1D<DualVariableNode*> var_nodes(2);
  var_nodes[0] = var_[var1];
  var_nodes[1] = var_[var2];

  PottsDualFactorNode* ptr = new PottsDualFactorNode(var_nodes, lambda);
  return add_factor(ptr);
}

//return factor number
uint FactorDualOpt::add_generic_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_ || var3 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  Storage1D<DualVariableNode*> var_nodes(3);
  var_nodes[0] = var_[var1];
  var_nodes[1] = var_[var2];
  var_nodes[2] = var_[var3];

  if (ref) {
    TernaryDualRefFactorNode* ptr = new TernaryDualRefFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
  else {
    TernaryDualFactorNode* ptr = new TernaryDualFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
}

uint FactorDualOpt::add_second_diff_factor(uint var1, uint var2, uint var3, float lambda) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_ || var3 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  Storage1D<DualVariableNode*> var_nodes(3);
  var_nodes[0] = var_[var1];
  var_nodes[1] = var_[var2];
  var_nodes[2] = var_[var3];

  SecondDiffDualFactorNode* ptr = new SecondDiffDualFactorNode(var_nodes,lambda);
  return add_factor(ptr);
}


//return factor number
uint FactorDualOpt::add_generic_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                                                    const Storage1D<Math3D::Tensor<float> >& cost, bool ref) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_ || var3 >= nUsedVars_ || var4 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  Storage1D<DualVariableNode*> var_nodes(4);
  var_nodes[0] = var_[var1];
  var_nodes[1] = var_[var2];
  var_nodes[2] = var_[var3];
  var_nodes[3] = var_[var4];

  if (ref) {
    FourthOrderDualRefFactorNode* ptr = new FourthOrderDualRefFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
  else {
    FourthOrderDualFactorNode* ptr = new FourthOrderDualFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }

}


//all participating variables must be binary
uint FactorDualOpt::add_one_of_N_factor(const Math1D::Vector<uint>& var) {

  Storage1D<DualVariableNode*> var_nodes(var.size());

  for (uint k=0; k < var.size(); k++) {

    if (var[k] >= nUsedVars_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    var_nodes[k] = var_[var[k]];

    if (var_[var[k]]->nLabels() != 2) {
      INTERNAL_ERROR << " variables of 1-of-N nodes must be binary. Exiting..." << std::endl;
      exit(1);
    }
  }

  if (var.size() == 1) {
    
    Math1D::Vector<float> cost(2);
    cost[0] = 1e20;
    cost[1] = 0.0;

    var_nodes[0]->add_cost(cost);

    return MAX_UINT;
  }

  OneOfNDualFactorNode* ptr = new OneOfNDualFactorNode(var_nodes);

  return add_factor(ptr);
}

//all participating variables must be binary
uint FactorDualOpt::add_cardinality_factor(const Math1D::Vector<uint>& var, const Math1D::Vector<float>& card_cost, bool ref) {

  if (card_cost.size() <= var.size()) {
    INTERNAL_ERROR << " dimension mismatch. Exiting..." << std::endl;
    exit(1);
  }


  if (var.size() == 1) {
    if (var[0] >= nUsedVars_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    var_[var[0]]->add_cost(card_cost);

    return MAX_UINT;
  }
  else {

    Storage1D<DualVariableNode*> var_nodes(var.size());

    for (uint k=0; k < var.size(); k++) {
      if (var[k] >= nUsedVars_) {
        INTERNAL_ERROR << "out of range. Exiting." << std::endl;
        exit(1);
      }
      
      var_nodes[k] = var_[var[k]];
      
      if (var_[var[k]]->nLabels() != 2) {
        INTERNAL_ERROR << " variables of Cardinality nodes must be binary. Exiting..." << std::endl;
        exit(1);
      }
    }

    DualFactorNode* ptr;
    if (!ref)
      ptr = new CardinalityDualFactorNode(var_nodes, card_cost);
    else
      ptr = new CardinalityDualFactorRefNode(var_nodes, card_cost);

    return add_factor(ptr);
  }
}

//all participating variables must be binary
uint FactorDualOpt::add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive,
                                          int rhs_lower, int rhs_upper) {


  if (rhs_lower > rhs_upper) {
    INTERNAL_ERROR << " INFEASIBLE CONSTRAINT" << std::endl;
    exit(1);
  }


  uint nUseful = 0;

  int nPos = 0;
  int nNeg = 0;

  //check for variables whose value is essentially fixed due to the cost vector
  for (uint k=0; k < var.size(); k++) {

    if (var[k] >= nUsedVars_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    if (var_[var[k]]->nLabels() != 2) {
      INTERNAL_ERROR << " variables of BILP nodes must be binary. Exiting..." << std::endl;
      exit(1);
    }
    
    const Math1D::Vector<float>& cur_cost = var_[var[k]]->cost();

    if (fabs(cur_cost[0] - cur_cost[1]) < 1e10) {
      nUseful++;
      if (positive[k])
        nPos++;
      else
        nNeg++;
    }
    else {
      if (cur_cost[0] > cur_cost[1]) {
        if (positive[k]) {
          rhs_lower--;
          rhs_upper--;
        }
        else {
          rhs_lower++;
          rhs_upper++;
        }
      }
    }
  }


  if (nUseful != 0 && rhs_lower <= -nNeg && rhs_upper >= nPos)
    nUseful = 0; //constraint is always true


  if (nUseful == 0) {

    std::cerr << "WARNING: removed superfluous constraint factor" << std::endl;
    
    return MAX_UINT;
  }
  

  if (rhs_upper < -nNeg || rhs_lower > nPos) {

    INTERNAL_ERROR << " INFEASIBLE CONSTRAINT" << std::endl;
    exit(1);
  }



  // if (nUseful != var.size())
  //   std::cerr << "removed " << (var.size() - nUseful) << " / " << var.size() << " vars for BILP node" << std::endl;

  // if (nUseful < 2) {
  //   std::cerr << "only " << nUseful << " out of " << var.size() << " variables are actually not fixed" << std::endl;
  
  //   for (uint k=0; k < var.size(); k++)
  //     std::cerr << "cost: " << var_[var[k]]->cost() << std::endl;
  
  //   std::cerr << "var: " << var << std::endl;
  // }
  
  Storage1D<DualVariableNode*> vars(nUseful);
  Storage1D<bool> reduced_positive(nUseful);
  
  uint next = 0;
  
  for (uint k=0; k < var.size(); k++) {
    
    if (fabs(var_[var[k]]->cost()[0] - var_[var[k]]->cost()[1]) < 1e10) {
      vars[next] = var_[var[k]];
      reduced_positive[next] = positive[k];
      next++;
    }
  }
  
  assert(next == nUseful);
  assert(nUsedFactors_ < factor_.size());
  

  if (nUseful == 1) {
    
    std::cerr << "happens" << std::endl;
    
    if (nPos == 0) {
      
      //invert the constraint
      
      double temp_lower = rhs_lower;
      rhs_lower = -rhs_upper;
      rhs_upper = -temp_lower;
    }
    
    Math1D::Vector<float> add_cost(2,0.0);
    
    if (rhs_lower == 1) {
      add_cost[0] = 1e30;
    }
    else if (rhs_upper == 0) {
      add_cost[1] = 1e30;
    }
    else {
      INTERNAL_ERROR << "STRANGE CONSTRAINT. exiting" << std::endl;
      exit(1);
    }
    
    vars[0]->add_cost(add_cost);
    return MAX_UINT;
  }
  else {

  
    DualFactorNode* ptr;
    
    if (nNeg == 0) 
    //if (false)
      ptr = new AllPosBILPConstraintDualFactorNode(vars, rhs_lower, rhs_upper);
    else
      ptr = new BILPConstraintDualFactorNode(vars, reduced_positive, rhs_lower, rhs_upper);
    
    return add_factor(ptr);    
  }
}

uint FactorDualOpt::best_of_n(uint fac_num) {

  if (fac_num >= nUsedFactors_) {
    INTERNAL_ERROR << "factor index out of bounds. Exiting." << std::endl;
    exit(1);
  }

  OneOfNDualFactorNode* fac = dynamic_cast<OneOfNDualFactorNode*>(factor_[fac_num]);
  
  if (fac == 0) {
    INTERNAL_ERROR << "not a 1-of-N factor. Exiting. " << std::endl;
    exit(1);
  }
  
  Math1D::Vector<uint> min_labels(fac->participating_nodes().size());
  fac->compute_minimizer(min_labels);

  for (uint k=0; k < min_labels.size(); k++) {
    if (min_labels[k] == 1)
      return k;
  }

  return MAX_UINT;
}

//NOTE: after calling this routine, owned factors can no longer be created
uint FactorDualOpt::pass_in_factor_node(DualFactorNode* factor) {

  return add_factor(factor);
}

/**** run inference ***/

double FactorDualOpt::dual_bca(uint nIter, DualBCAMode mode, bool init, bool quiet) {


  labeling_.resize(nUsedVars_,0);

  uint arg_min;

  std::cerr.precision(10);

  if (init) {
    for (uint v=0; v < nUsedVars_; v++) {
      var_[v]->init_dual_vars();
    }
  }

  size_t effort_per_iter = 0;
  for (uint f=0; f < nUsedFactors_; f++) {
    effort_per_iter += factor_[f]->participating_nodes().size() * factor_[f]->participating_nodes().size();
  }

  //DEBUG (for ICML)
  //std::string filename = (mode == DUAL_BCA_MODE_MPLP) ? "mplp.dat" : "msd.dat";
  //std::ofstream of(filename.c_str());
  //END_DEBUG

  for (uint iter = 1; iter <= nIter; iter++) {

    std::cerr << "******** iteration " << iter << " ************" << std::endl;

    for (uint f=0; f < nUsedFactors_; f++) {

      //DEBUG
      // std::cerr << "f: " << f << "/" << nUsedFactors_ << std::endl;

      // if (f > 10)
      // 	break;
      //END_DEBUG

      //std::cerr << "f: " << f << std::endl;
      factor_[f]->update_duals(mode);
    }
    
    double lower_bound = 0.0;

    if (!quiet && (iter % 1) == 0) {
      for (uint v=0; v < nUsedVars_; v++) {
        lower_bound += var_[v]->dual_value(arg_min);
        labeling_[v] = arg_min;
      }
      
      //std::cerr << "lb part for vars: " << lower_bound << std::endl;

      for (uint f=0; f < nUsedFactors_; f++)
        lower_bound += factor_[f]->dual_value();
      
      std::cerr << "lower bound: " << lower_bound << std::endl;
      
      //std::cerr << "primal energy: " << labeling_energy() << std::endl;

      //DEBUG (for ICML)
      //of << (effort_per_iter * iter) << " " << lower_bound << std::endl;
      //END_DEBUG
    }
  }

  size_t message_effort = 0;
  for (uint f=0; f < nUsedFactors_; f++) {
    message_effort += factor_[f]->participating_nodes().size() * factor_[f]->participating_nodes().size();
  }
  message_effort *= nIter;
  std::cerr << "message efffort " << message_effort << std::endl;

  double lower_bound = 0.0;

  for (uint v=0; v < nUsedVars_; v++) {
    lower_bound += var_[v]->dual_value(arg_min);
    labeling_[v] = arg_min;
  }
  
  for (uint f=0; f < nUsedFactors_; f++)
    lower_bound += factor_[f]->dual_value();

  std::cerr << "lower bound: " << lower_bound << std::endl;

  return lower_bound;
}

double FactorDualOpt::subgradient_opt(uint nIter, double start_step_size) {

  //TODO: implement projective version

  double best_dual = -1e300;

  Math1D::Vector<uint> var_label(nUsedVars_);

  Storage1D<Math1D::Vector<uint> > factor_label(nUsedFactors_);
  
  for (uint f=0; f < nUsedFactors_; f++) {
    factor_label[f].resize(factor_[f]->participating_nodes().size());
  }

  for (uint iter=1; iter <= nIter; iter++) {

    std::cerr << "iteration #" << iter << std::endl;

    double step_size = start_step_size / iter;

    double cur_bound = 0.0;

    for (uint v=0; v < nUsedVars_; v++) {
      uint cur_label;
      cur_bound += var_[v]->dual_value(cur_label);
      var_label[v] = cur_label;
    }
    
    std::cerr << "A" << std::endl;

    for (uint f=0; f < nUsedFactors_; f++) {
      cur_bound += factor_[f]->compute_minimizer(factor_label[f]);
    }

    if (cur_bound > best_dual)
      best_dual = cur_bound;

    std::cerr << "cur bound: " << cur_bound << ", best ever: " << best_dual << std::endl;

    for (uint v=0; v < nUsedVars_; v++) {
      const uint cur_label = var_label[v];

      const uint nCurFactors = var_[v]->neighboring_factor().size();
      const uint nLabels = var_[v]->nLabels();

      double* ptr = var_[v]->get_dual_var_start();

      for (uint k=0; k < nCurFactors; k++) {
        ptr[cur_label] += step_size;
        ptr += nLabels;
      }
    }

    for (uint f=0; f < nUsedFactors_; f++) {

      const Storage1D<DualVariableNode*>& cur_nodes = factor_[f]->participating_nodes();

      for (uint k=0; k < cur_nodes.size(); k++) {

        uint cur_label = factor_label[f][k];

        double* ptr = cur_nodes[k]->get_dual_vars(factor_[f]);
        ptr[cur_label] -= step_size;
      }
    }
  }

  labeling_ = var_label;

  return best_dual;
}


/**** get the state of the solver ****/

double FactorDualOpt::labeling_energy() {

  double energy = 0.0;

  std::map<DualVariableNode*,uint> label;

  for (uint k=0; k < nUsedVars_; k++) {
    assert(label.find(var_[k]) == label.end());
    label[var_[k]] = labeling_[k];
    energy += var_[k]->cost(labeling_[k]);
  }

  //std::cerr << "unary cost: " << energy << std::endl;

  //  std::cerr << "sum of labeling: " << labeling_.sum() << std::endl;
  //std::cerr << "nFactors: " << nUsedFactors_ << std::endl;

  for (uint k=0; k < nUsedFactors_; k++) {

    //std::cerr << "k: " << k << std::endl;

    const Storage1D<DualVariableNode*>& nodes = factor_[k]->participating_nodes();

    Math1D::Vector<uint> sub_labeling(nodes.size());
    for (uint i=0; i < nodes.size(); i++) {
      assert(label.find(nodes[i]) != label.end());
      sub_labeling[i] = label[nodes[i]];
    }

    energy += factor_[k]->cost(sub_labeling);
  }
  
  return energy;
}

const Math1D::Vector<uint>& FactorDualOpt::labeling() {

  return labeling_;
}


void FactorDualOpt::icm(uint nIter) {

  std::map<DualVariableNode*,uint> label;
  std::map<DualVariableNode*,uint> num;

  for (uint k=0; k < nUsedVars_; k++) {
    assert(label.find(var_[k]) == label.end());
    label[var_[k]] = labeling_[k];
    num[var_[k]] = k;
  }

  for (uint iter = 1; iter <= nIter; iter++) {

    std::cerr << "ICM iter " << iter << std::endl;

    for (uint v=0; v < nUsedVars_; v++) {

      const Storage1D<DualFactorNode*>& factors = var_[v]->neighboring_factor();

      Storage1D<Math1D::Vector<uint> > sublabeling(factors.size());
      
      Math1D::Vector<uint> var_index(factors.size(),MAX_UINT);

      for (uint f=0; f < factors.size(); f++) {

        const Storage1D<DualVariableNode*>& factor_var = factors[f]->participating_nodes();

        sublabeling[f].resize(factor_var.size());

        for (uint k=0; k < factor_var.size(); k++) {

          sublabeling[f][k] = label[factor_var[k]];
          if (factor_var[k] == var_[v])
            var_index[f] = k;
        }
      }

      /****** now search for the conditional mode *****/
      double mode_val = 1e300;
      uint mode = MAX_UINT;

      for (uint l=0; l < var_[v]->nLabels(); l++) {

        double hyp = var_[v]->cost(l);
        for (uint f=0; f < factors.size(); f++) {
          sublabeling[f][var_index[f]] = l;
          hyp += factors[f]->cost(sublabeling[f]);
        }

        if (hyp < mode_val) {
          mode_val = hyp;
          mode = l;
        }
      }

      labeling_[v] = mode;
      label[var_[v]] = mode;
    }
  }

}
