/******* written by Thomas Schoenemann as a private Person without employment, July 2011 *****/


/***** implements max-product belief propagation in negative log-space ****/

#include "factorMPBP.hh"
#include <map>

VariableNode::VariableNode(const Math1D::Vector<float>& cost) : cost_(cost) {

  message_matrix_.resize(cost_.size(),0); 
}

void VariableNode::add_factor(FactorNode* node) {

  uint nPrevFactors = neighboring_factor_.size();
  neighboring_factor_.resize(nPrevFactors+1,node);
  message_matrix_.resize(nLabels(),nPrevFactors+1,0.0); 
}

uint VariableNode::nLabels() const {
  return message_matrix_.xDim();
}

void VariableNode::init_messages() {
  message_matrix_.set_constant(0.0);
}

double* VariableNode::get_message(FactorNode* node)  {

  double* ptr = message_matrix_.direct_access();

  bool found = false;
  for (uint i=0; i < neighboring_factor_.size(); i++) {
    if (neighboring_factor_[i] == node) {
      found = true;
      break;
    }
    else
      ptr += message_matrix_.xDim();
  }

  if (!found) {
    INTERNAL_ERROR << "node not found" << std::endl;
    exit(1);
  }

  return ptr;
}

const Storage1D<FactorNode*>& VariableNode::neighboring_factor() const {
  return neighboring_factor_;
}


/*virtual*/ void VariableNode::compute_beliefs(Math1D::Vector<double>& beliefs) {

  Storage1D<const double*> factor_message(neighboring_factor_.size());

  for (uint i=0; i < neighboring_factor_.size(); i++) {
    factor_message[i] = neighboring_factor_[i]->get_message(this);
  }

  uint labels = nLabels();

  beliefs.resize(labels);

  for (uint k = 0; k < labels; k++) {
    
    double cur_belief = cost_[k];
    
    for (uint f=0; f < neighboring_factor_.size(); f++)
      cur_belief += factor_message[f][k];

    beliefs[k] = cur_belief;
  }

  //NOTE: we do not normalize as these are really log-beliefs
}

void VariableNode::compute_messages() {

  uint nLabels = cost_.size();

  Storage1D<const double*> factor_message(neighboring_factor_.size());

  for (uint i=0; i < neighboring_factor_.size(); i++) {
    factor_message[i] = neighboring_factor_[i]->get_message(this);
  }

  for (uint i=0; i < neighboring_factor_.size(); i++) {

    double min_message = 1e300;

    for (uint l=0; l < nLabels; l++) {
      
      assert(!isnan(cost_[l]));

      double cur_cost = cost_[l];

      for (uint j=0; j < neighboring_factor_.size(); j++) {

        if (j != i) {
          cur_cost += factor_message[j][l];
	  
          if (isnan(factor_message[j][l])) {
            std::cerr << "j: " << j << ", l: " << l << std::endl;
          }

          assert(!isnan(factor_message[j][l]));
          assert(factor_message[j][l] >= -1e-3);
        }
      }

      assert(!isnan(cur_cost));

      message_matrix_(l,i) = cur_cost;

      if (cur_cost < min_message) {
        min_message = cur_cost;
      }
    }

    for (uint l=0; l < nLabels; l++) {
      message_matrix_(l,i) -= min_message;

      if (isinf(message_matrix_(l,i))) {
	
        std::cerr << "min_message: " << min_message << std::endl;
        std::cerr << "cost: " << cost_ << std::endl;
      }

      assert(!isinf(message_matrix_(l,i)));
      assert(!isnan(message_matrix_(l,i)));
    }

  }
}

double VariableNode::cost(uint label) const {
  return cost_[label];
}

/***********************************/
  
FactorNode::FactorNode(const Storage1D<VariableNode*>& participating_vars) :
  participating_var_(participating_vars) {

  for (uint i=0; i < participating_var_.size(); i++)
    participating_var_[i]->add_factor(this);
}

/*virtual*/ FactorNode::~FactorNode() {}

const Storage1D<VariableNode*>& FactorNode::participating_nodes() const {

  return participating_var_;
}

/*virtual*/ bool FactorNode::process_labeling(Math1D::Vector<uint>& /*labels*/) const {
  //by default no action is taken
  return false;
}

/***********************************/

GenericFactorNode::GenericFactorNode(const Storage1D<VariableNode*>& participating_vars, 
                                     const VarDimStorage<float>& cost) :
  FactorNode(participating_vars), cost_(cost) {

  const uint size = participating_var_.size();

  message_.resize(size);
  for (uint k=0; k < size; k++) {
    message_[k].resize(participating_var_[k]->nLabels(),0.0);
  }
}

/*virtual*/ GenericFactorNode::~GenericFactorNode() {}

/*virtual*/ void GenericFactorNode::compute_messages() {


  Storage1D<const double*> var_message(participating_var_.size());

  for (uint k=0; k < participating_var_.size(); k++)
    var_message[k] = participating_var_[k]->get_message(this);

  const uint size = participating_var_.size();

  for (uint k=0; k < size; k++) {
    message_[k].set_constant(1e300);
  }

  Math1D::Vector<size_t> labeling(size,0);

  while (true) {

    double cost = cost_(labeling);
    for (uint k=0; k < size; k++) {
      cost += var_message[k][labeling[k]];
    }

    for (uint k=0; k < size; k++) {

      if (cost < message_[k][labeling[k]])
        message_[k][labeling[k]] = cost;
    }

    //increase labeling
    uint l;
    for (l=0; l < size; l++) {

      labeling[l] = (labeling[l] + 1) % message_[l].size(); //could be more efficient
      if (labeling[l] != 0)
        break;
    }

    if (l == size) //all zero after increase => cycle completed
      break;
  }

  //subtract the incorrectly added messages
  for (uint k=0; k < size; k++) {
    double min_msg = 1e300;
    for (uint l=0; l < message_[k].size(); l++) {
      message_[k][l] -= var_message[k][l];
      if (message_[k][l] < min_msg)
        min_msg = message_[k][l];
    }
    for (uint l=0; l < message_[k].size(); l++) {
      message_[k][l] -= min_msg;
      assert(message_[k][l] >= -1e-3);
    }
  }

}
  
/*virtual*/ double* GenericFactorNode::get_message(VariableNode* node) {

  const uint size = participating_var_.size();

  for (uint k=0; k < size; k++) {
    if (participating_var_[k] == node)
      return message_[k].direct_access();
  }

  INTERNAL_ERROR << "node not found" << std::endl;    
  exit(1);
  
  return 0;
}
  
/*virtual*/ double GenericFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  Math1D::Vector<size_t> size_t_labels(labels.size());

  for (uint k=0; k < labels.size(); k++)
    size_t_labels[k] = labels[k];

  return cost_(size_t_labels);
}

/*virtual*/ void GenericFactorNode::init_messages() {

  const uint size = participating_var_.size();

  for (uint k=0; k < size; k++) {
    message_[k].set_constant(0.0);
  }
}

/***********************************/

BinaryFactorNodeBase::BinaryFactorNodeBase(const Storage1D<VariableNode*>& participating_vars) : 
  FactorNode(participating_vars) {
}

/*virtual*/ void BinaryFactorNodeBase::init_messages() {

  message1_.set_constant(0.0);
  message2_.set_constant(0.0);
}

/*virtual*/ double* BinaryFactorNodeBase::get_message(VariableNode* node) {

  if (node == participating_var_[0])
    return message1_.direct_access();
  else if (node == participating_var_[1])
    return message2_.direct_access();
  else {
    INTERNAL_ERROR << "node not found" << std::endl;    
    exit(1);
  }
}

void BinaryFactorNodeBase::compute_messages(const Math2D::Matrix<float>& cost) {

  assert( participating_var_.size() == 2);

  Storage1D<const double*> var_message(participating_var_.size());

  for (uint k=0; k < participating_var_.size(); k++)
    var_message[k] = participating_var_[k]->get_message(this);

  //Message 1
  uint nLabels1 = participating_var_[0]->nLabels();
  uint nLabels2 = participating_var_[1]->nLabels();  

  for (uint l=0; l < nLabels1; l++) {

    double min_cost = 1e300;

    for (uint k=0; k < nLabels2; k++) {
      
      double hyp_cost = cost(l,k) + var_message[1][k];
      if (hyp_cost < min_cost)
        min_cost = hyp_cost;
    }

    message1_[l] = min_cost;
  }
  double m1 = message1_.min();
  for (uint k=0; k < message1_.size(); k++)
    message1_[k] -= m1;

  //Message 2
  for (uint k=0; k < nLabels2; k++) {

    double min_cost = 1e300;

    for (uint l=0; l < nLabels1; l++) {

      double hyp_cost = cost(l,k) + var_message[0][l];
      if (hyp_cost < min_cost)
        min_cost = hyp_cost;
    }

    message2_[k] = min_cost;
  }
  double m2 = message2_.min();
  for (uint k=0; k < message2_.size(); k++)
    message2_[k] -= m2;
}

/***********************************/

BinaryFactorNode::BinaryFactorNode(const Storage1D<VariableNode*>& participating_vars, const Math2D::Matrix<float>& cost) :
  BinaryFactorNodeBase(participating_vars), cost_(cost) {

  assert(participating_var_[0]->nLabels() == cost.xDim());
  assert(participating_var_[1]->nLabels() == cost.yDim());

  message1_.resize(cost.xDim(),0);
  message2_.resize(cost.yDim(),0);
}

/*virtual*/ BinaryFactorNode::~BinaryFactorNode(){}

/*virtual*/ void BinaryFactorNode::compute_messages() {

  BinaryFactorNodeBase::compute_messages(cost_);
}

/*virtual*/ double BinaryFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  return cost_(labels[0],labels[1]);
}


/***********************************/

BinaryRefFactorNode::BinaryRefFactorNode(const Storage1D<VariableNode*>& participating_vars, const Math2D::Matrix<float>& cost) :
  BinaryFactorNodeBase(participating_vars), cost_(cost) {

  assert(participating_var_[0]->nLabels() == cost.xDim());
  assert(participating_var_[1]->nLabels() == cost.yDim());

  message1_.resize(cost.xDim(),0);
  message2_.resize(cost.yDim(),0);
}

/*virtual*/ BinaryRefFactorNode::~BinaryRefFactorNode() {}

/*virtual*/ void BinaryRefFactorNode::compute_messages() {

  BinaryFactorNodeBase::compute_messages(cost_);
}

/*virtual*/ double BinaryRefFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  return cost_(labels[0],labels[1]);
}


/***********************************/

PottsFactorNode::PottsFactorNode(const Storage1D<VariableNode*>& participating_vars, float lambda) :
  BinaryFactorNodeBase(participating_vars), lambda_(lambda) {

  assert(participating_vars.size() == 2);
  assert(participating_vars[0]->nLabels() == participating_vars[1]->nLabels());

  message1_.resize(participating_vars[0]->nLabels(),0);
  message2_.resize(participating_vars[1]->nLabels(),0);  
}

/*virtual*/ PottsFactorNode::~PottsFactorNode() {}

/*virtual*/ void PottsFactorNode::compute_messages() {

  Storage1D<const double*> var_message(participating_var_.size());

  Math1D::Vector<double> msg_min(2);

  for (uint k=0; k < participating_var_.size(); k++) {
    var_message[k] = participating_var_[k]->get_message(this);

    double cur_min = 1e300;
    for (uint l=0; l < participating_var_[k]->nLabels(); l++) {
      if (var_message[k][l] < cur_min)
        cur_min = var_message[k][l];
    }

    msg_min[k] = cur_min;
  }

  //Message 1
  uint nLabels1 = participating_var_[0]->nLabels();
  uint nLabels2 = participating_var_[1]->nLabels();  

  for (uint l=0; l < nLabels1; l++) {

    message1_[l] = std::min(msg_min[1] + lambda_, var_message[1][l]);
  }

  //Message 2
  for (uint k=0; k < nLabels2; k++) {

    message2_[k] = std::min(msg_min[0] + lambda_, var_message[0][k]);
  }
}

/*virtual*/ double PottsFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  return (labels[0] != labels[1]) ? lambda_ : 0.0;
}


/***********************************/

TernaryFactorNodeBase::TernaryFactorNodeBase(const Storage1D<VariableNode*>& participating_vars) :
  FactorNode(participating_vars) {
}

/*virtual*/ double* TernaryFactorNodeBase::get_message(VariableNode* node) {

  if (node == participating_var_[0])
    return message1_.direct_access();
  else if (node == participating_var_[1])
    return message2_.direct_access();
  else if (node == participating_var_[2])
    return message3_.direct_access();
  else {
    INTERNAL_ERROR << "node not found" << std::endl;    
    exit(1);
    return 0;
  }
}

/*virtual*/ void TernaryFactorNodeBase::init_messages() {

  message1_.set_constant(0.0);
  message2_.set_constant(0.0);
  message3_.set_constant(0.0);
}

void TernaryFactorNodeBase::compute_messages(const Math3D::Tensor<float>& cost) {

  assert( participating_var_.size() == 3);

  Storage1D<const double*> var_message(participating_var_.size());

  for (uint k=0; k < participating_var_.size(); k++)
    var_message[k] = participating_var_[k]->get_message(this);

  //Message 1
  uint nLabels1 = participating_var_[0]->nLabels();
  uint nLabels2 = participating_var_[1]->nLabels();  
  uint nLabels3 = participating_var_[2]->nLabels();  
  
  for (uint l=0; l < nLabels1; l++) {

    double min_cost = 1e300;

    for (uint k=0; k < nLabels2; k++) {

      for (uint m=0; m < nLabels3; m++) {
    
        double hyp_cost = cost(l,k,m) + var_message[1][k] + var_message[2][m];
        if (hyp_cost < min_cost)
          min_cost = hyp_cost;
      }
    }

    message1_[l] = min_cost;
  }

  double msg_min = message1_.min();
  for (uint k=0; k < message1_.size(); k++)
    message1_[k] -= msg_min;
  
  //Message 2
  for (uint k=0; k < nLabels2; k++) {

    double min_cost = 1e300;

    for (uint l=0; l < nLabels1; l++) {

      for (uint m=0; m < nLabels3; m++) {

        double hyp_cost = cost(l,k,m) + var_message[0][l] + var_message[2][m];
        if (hyp_cost < min_cost)
          min_cost = hyp_cost;
      }
    }

    message2_[k] = min_cost;
  }

  msg_min = message2_.min();
  for (uint k=0; k < message2_.size(); k++)
    message2_[k] -= msg_min;

  //Message 3
  for (uint m=0; m < nLabels3; m++) {

    double min_cost = 1e300;

    for (uint l=0; l < nLabels1; l++) {

      for (uint k=0; k < nLabels2; k++) {

        double hyp_cost = cost(l,k,m) + var_message[0][l] + var_message[1][k];
        if (hyp_cost < min_cost)
          min_cost = hyp_cost;
      }
    }

    message3_[m] = min_cost;
  }

  msg_min = message3_.min();
  for (uint k=0; k < message3_.size(); k++)
    message3_[k] -= msg_min;
}


/***********************************/

TernaryFactorNode::TernaryFactorNode(const Storage1D<VariableNode*>& participating_vars, 
                                     const Math3D::Tensor<float>& cost) : 
  TernaryFactorNodeBase(participating_vars), cost_(cost) {

  assert(participating_var_[0]->nLabels() == cost.xDim());
  assert(participating_var_[1]->nLabels() == cost.yDim());
  assert(participating_var_[2]->nLabels() == cost.zDim());

  message1_.resize(cost.xDim(),0);
  message2_.resize(cost.yDim(),0);
  message3_.resize(cost.zDim(),0);
}

/*virtual*/ TernaryFactorNode::~TernaryFactorNode() {}

/*virtual*/ void TernaryFactorNode::compute_messages() {

  TernaryFactorNodeBase::compute_messages(cost_);
}

/*virtual*/ double TernaryFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  return cost_(labels[0],labels[1], labels[2]);
}

/***********************************/

TernaryRefFactorNode::TernaryRefFactorNode(const Storage1D<VariableNode*>& participating_vars, 
                                           const Math3D::Tensor<float>& cost) : 
  TernaryFactorNodeBase(participating_vars), cost_(cost) {

  assert(participating_var_[0]->nLabels() == cost.xDim());
  assert(participating_var_[1]->nLabels() == cost.yDim());
  assert(participating_var_[2]->nLabels() == cost.zDim());

  message1_.resize(cost.xDim(),0);
  message2_.resize(cost.yDim(),0);
  message3_.resize(cost.zDim(),0);
}

/*virtual*/ TernaryRefFactorNode::~TernaryRefFactorNode() {}

/*virtual*/ void TernaryRefFactorNode::compute_messages() {

  TernaryFactorNodeBase::compute_messages(cost_);
}

/*virtual*/ double TernaryRefFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  return cost_(labels[0],labels[1], labels[2]);
}


/***********************************/

FourthOrderFactorNodeBase::FourthOrderFactorNodeBase(const Storage1D<VariableNode*>& participating_vars) :
  FactorNode(participating_vars) {}

void FourthOrderFactorNodeBase::compute_messages(const Storage1D<Math3D::Tensor<float> >& cost) {

  assert( participating_var_.size() == 4);

  Storage1D<const double*> var_message(participating_var_.size());

  for (uint k=0; k < participating_var_.size(); k++)
    var_message[k] = participating_var_[k]->get_message(this);

  uint nLabels1 = participating_var_[0]->nLabels();
  uint nLabels2 = participating_var_[1]->nLabels();  
  uint nLabels3 = participating_var_[2]->nLabels();  
  uint nLabels4 = participating_var_[3]->nLabels();  

  //Message 1
  for (uint l1 = 0; l1 < nLabels1; l1++) {

    double min_cost = 1e300;

    for (uint l2 = 0; l2 < nLabels2; l2++) {

      for (uint l3 = 0; l3 < nLabels3; l3++) {

        for (uint l4 = 0; l4 < nLabels4; l4++) {

          double hyp =  cost[l1](l2,l3,l4) + var_message[1][l2] + var_message[2][l3] + var_message[3][l4];

          if (hyp < min_cost)
            min_cost = hyp;
        }
      }
    }

    message1_[l1] = min_cost;
  }
  double msg_min = message1_.min();
  for (uint k=0; k < message1_.size(); k++)
    message1_[k] -= msg_min;

  //Message 2
  for (uint l2 = 0; l2 < nLabels2; l2++) {

    double min_cost = 1e300;

    for (uint l1 = 0; l1 < nLabels1; l1++) {

      for (uint l3 = 0; l3 < nLabels3; l3++) {

        for (uint l4 = 0; l4 < nLabels4; l4++) {

          double hyp =  cost[l1](l2,l3,l4) + var_message[0][l1] + var_message[2][l3] + var_message[3][l4];

          if (hyp < min_cost)
            min_cost = hyp;
        }
      }
    }

    message2_[l2] = min_cost;
  }

  msg_min = message2_.min();
  for (uint k=0; k < message2_.size(); k++)
    message2_[k] -= msg_min;

  //Message 3
  for (uint l3 = 0; l3 < nLabels3; l3++) {
	
    double min_cost = 1e300;

    for (uint l1 = 0; l1 < nLabels1; l1++) {

      for (uint l2 = 0; l2 < nLabels2; l2++) {

        for (uint l4 = 0; l4 < nLabels4; l4++) {

          double hyp =  cost[l1](l2,l3,l4) + var_message[0][l1] + var_message[1][l2] + var_message[3][l4];

          if (hyp < min_cost)
            min_cost = hyp;
        }
      }
    }

    message3_[l3] = min_cost;
  }

  msg_min = message3_.min();
  for (uint k=0; k < message3_.size(); k++)
    message3_[k] -= msg_min;

  //Message 4
  for (uint l4 = 0; l4 < nLabels4; l4++) {
	
    double min_cost = 1e300;

    for (uint l1 = 0; l1 < nLabels1; l1++) {

      for (uint l2 = 0; l2 < nLabels2; l2++) {

        for (uint l3 = 0; l3 < nLabels3; l3++) {

          double hyp =  cost[l1](l2,l3,l4) + var_message[0][l1] + var_message[1][l2] + var_message[2][l3];

          if (hyp < min_cost)
            min_cost = hyp;
        }
      }
    }

    message4_[l4] = min_cost;
  }

  msg_min = message4_.min();
  for (uint k=0; k < message4_.size(); k++)
    message4_[k] -= msg_min;
}

/*virtual*/ void FourthOrderFactorNodeBase::init_messages() {
  message1_.set_constant(0.0);
  message2_.set_constant(0.0);
  message3_.set_constant(0.0);
  message4_.set_constant(0.0);
}

/*virtual*/ double* FourthOrderFactorNodeBase::get_message(VariableNode* node) {

  if (node == participating_var_[0])
    return message1_.direct_access();
  else if (node == participating_var_[1])
    return message2_.direct_access();
  else if (node == participating_var_[2])
    return message3_.direct_access();
  else if (node == participating_var_[3])
    return message4_.direct_access();
  else {
    INTERNAL_ERROR << "node not found" << std::endl;    
    exit(1);
    return 0;
  }  
}

/***********************************/

FourthOrderFactorNode::FourthOrderFactorNode(const Storage1D<VariableNode*>& participating_vars, 
                                             const Storage1D<Math3D::Tensor<float> >& cost) :
  FourthOrderFactorNodeBase(participating_vars), cost_(cost) {

  message1_.resize(cost.size());
  message2_.resize(cost[0].xDim());
  message3_.resize(cost[0].yDim());
  message4_.resize(cost[0].zDim());
}

/*virtual*/ FourthOrderFactorNode::~FourthOrderFactorNode() {}

/*virtual*/ void FourthOrderFactorNode::compute_messages() {

  FourthOrderFactorNodeBase::compute_messages(cost_);
}

/*virtual*/ double FourthOrderFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  assert(labels.size() == 4);

  uint l1 = labels[0];
  uint l2 = labels[1];
  uint l3 = labels[2];
  uint l4 = labels[3];

  return cost_[l1](l2,l3,l4);
}

/***********************************/

FourthOrderRefFactorNode::FourthOrderRefFactorNode(const Storage1D<VariableNode*>& participating_vars, 
                                                   const Storage1D<Math3D::Tensor<float> >& cost) :
  FourthOrderFactorNodeBase(participating_vars), cost_(cost) {
}

/*virtual*/ FourthOrderRefFactorNode::~FourthOrderRefFactorNode() {}

/*virtual*/ void FourthOrderRefFactorNode::compute_messages() {

  FourthOrderFactorNodeBase::compute_messages(cost_);
}

/*virtual*/ double FourthOrderRefFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  assert(labels.size() == 4);

  uint l1 = labels[0];
  uint l2 = labels[1];
  uint l3 = labels[2];
  uint l4 = labels[3];

  return cost_[l1](l2,l3,l4);
}


/***********************************/

OneOfNFactorNode::OneOfNFactorNode(const Storage1D<VariableNode*>& participating_vars)
  : FactorNode(participating_vars), message_matrix_(2,participating_vars.size(),0.0) {

  for (uint v=0; v < participating_vars.size(); v++) {
    if (participating_vars[v]->nLabels() != 2) {
      INTERNAL_ERROR << " for a 1-of-N potential all variables need to be binary. Exiting." << std::endl;
      exit(1);
    }
  }

  assert(participating_vars.size() >= 1);
}

/*virtual*/ OneOfNFactorNode::~OneOfNFactorNode() {}

/*virtual*/ void OneOfNFactorNode::init_messages() {
  message_matrix_.set_constant(0.0);
}

/*virtual*/ bool OneOfNFactorNode::process_labeling(Math1D::Vector<uint>& labels) const {

  if (labels.sum() != 1) {

    Math1D::Vector<double> one_beliefs(participating_var_.size());

    Math1D::Vector<double> belief(2);
    
    for (uint k=0; k < participating_var_.size(); k++) {
      participating_var_[k]->compute_beliefs(belief);
      one_beliefs[k] = belief[1];
    }

    double best_belief = one_beliefs[0];
    uint arg_best = 0;
    for (uint k=1; k < participating_var_.size(); k++) {

      if (one_beliefs[k] < best_belief) {
        best_belief = one_beliefs[k];
        arg_best = k;
      }
    }
    
    labels.set_constant(0);
    labels[arg_best] = 1;
    return true;
  }

  return false;
}

/*virtual*/ double* OneOfNFactorNode::get_message(VariableNode* node) {

  double* ptr = message_matrix_.direct_access();

  for (uint i=0; i < participating_var_.size(); i++) {

    if (participating_var_[i] == node)
      break;
    else
      ptr += 2;
  }

  return ptr;
}

/*virtual*/ double OneOfNFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  assert(labels.size() == participating_var_.size());

  uint sum = labels.sum();

  return (sum == 1) ? 0.0 : 1e15;
}

/*virtual*/ void OneOfNFactorNode::compute_messages() {

  const uint size = participating_var_.size();

  Math1D::Vector<double> rel_msg(size);

  for (uint k=0; k < size; k++) {

    const double* msg = participating_var_[k]->get_message(this);

    rel_msg[k] = msg[1] - msg[0];
  }

  double best = 1e50;
  uint arg_best = 0;
  double second_best = 1e50;
  
  for (uint k=0; k < size; k++) {

    if (rel_msg[k] < best) {

      second_best = best;

      best = rel_msg[k];
      arg_best = k;
    }
    else if (rel_msg[k] < second_best) {

      second_best = rel_msg[k];
    }
  }
    
  for (uint k=0; k < size; k++) {

    //msg0
    if (k == arg_best)
      message_matrix_(0,k) = second_best;
    else
      message_matrix_(0,k) = best;
    
    //msg 1
    message_matrix_(1,k) = 0.0;
 
    if (message_matrix_(0,k) < 0.0) {
      message_matrix_(1,k) -= message_matrix_(0,k);
      message_matrix_(0,k) = 0.0;
    }
  }
}


/***********************************/

CardinalityFactorNode::CardinalityFactorNode(const Storage1D<VariableNode*>& participating_vars, 
                                             const Math1D::Vector<float>& card_cost) :
  FactorNode(participating_vars), card_cost_(card_cost), message_matrix_(2,participating_vars.size(),0.0) {

  for (uint v=0; v < participating_vars.size(); v++) {
    if (participating_vars[v]->nLabels() != 2) {
      INTERNAL_ERROR << " for a cardinality potential all variables need to be binary. Exiting." << std::endl;
      exit(1);
    }
  }

  if (card_cost_.size() != participating_vars.size() + 1) {

    INTERNAL_ERROR << " Cardinality vector size does not match with the number of variables. Exiting." << std::endl;
    exit(1);
  }
}

/*virtual*/ CardinalityFactorNode::~CardinalityFactorNode() {}

/*virtual*/ void CardinalityFactorNode::init_messages() {
  message_matrix_.set_constant(0.0);
}

/*virtual*/ double CardinalityFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  assert(labels.size() == participating_var_.size());

  uint sum = labels.sum();

  return card_cost_[sum];
}

/*virtual*/ double* CardinalityFactorNode::get_message(VariableNode* node) {

  double* ptr = message_matrix_.direct_access();

  for (uint i=0; i < participating_var_.size(); i++) {

    if (participating_var_[i] == node)
      break;
    else
      ptr += 2;
  }

  return ptr;
}

/*virtual*/ void CardinalityFactorNode::compute_messages() {

  // according to the HOP-MAP paper [Tarlow, Givoni, Zemel 2010]

  const uint size = participating_var_.size();

  Storage1D<std::pair<double,uint> > rel_msg(size);

  for (uint k=0; k < size; k++) {

    const double* msg = participating_var_[k]->get_message(this);

    rel_msg[k] = std::make_pair(msg[1] - msg[0], k);
  }

  std::sort(rel_msg.direct_access(), rel_msg.direct_access() + size);

  Math1D::Vector<uint> order(size);
  for (uint k=0; k < size; k++) {
    order[rel_msg[k].second] = k + 1;
  }
  
  Math1D::Vector<double> cum_sum(size+1);
  cum_sum[0] = rel_msg[0].first;
  for (uint k=1; k < size; k++) {
    cum_sum[k] = cum_sum[k-1] + rel_msg[k].first;
  }

  Math1D::Vector<double> cum_best(size + 1);
  cum_best[0] = card_cost_[0]; 

  for (uint k=1; k <= size; k++) {

    double hyp = card_cost_[k] + cum_sum[k-1];
    cum_best[k] = std::min(hyp, cum_best[k-1]);
  }


  Math1D::Vector<double> cum_best_m1(size+1);
  cum_best_m1[0] = 1e300;
  cum_best_m1[1] = card_cost_[1];

  for (uint k=2; k <= size; k++) {
    cum_best_m1[k] = std::min(cum_best_m1[k-1], card_cost_[k] + cum_sum[k-2]  ); 
  }

  Math1D::Vector<double> rev_cum_best(size + 1);
  rev_cum_best[size] = card_cost_[size] + cum_sum[size-1];
  for (int k=size-1; k >= 1; k--) {

    double hyp = card_cost_[k] + cum_sum[k-1];
    rev_cum_best[k] = std::min(hyp, rev_cum_best[k+1]);
  }
  rev_cum_best[0] = 1e300;

  Math1D::Vector<double> rev_cum_best_p1(size+1);
  rev_cum_best_p1[size] = 1e300;
  for (int k=size-1; k >= 1; k--) {
    rev_cum_best_p1[k] = std::min(rev_cum_best_p1[k+1], card_cost_[k] + cum_sum[k]);
  }

  for (uint k=0; k < size; k++) {

    const uint cur_order = order[k];

    const double* msg = participating_var_[k]->get_message(this);
    double cur_rel_msg = msg[1] - msg[0];

    const double val0 = std::min(cum_best[cur_order-1], rev_cum_best_p1[cur_order] - cur_rel_msg);
    message_matrix_(0,k) = val0;

    const double val1 = std::min(rev_cum_best[cur_order], cum_best_m1[cur_order-1] + cur_rel_msg );
    message_matrix_(1,k) = val1;

    //DEBUG
#if 0
    //a) check msg0
    double best0 = card_cost_[0];
    uint arg_best0 = 0;
    for (uint l=1; l < size; l++) {

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
                << ", order: " << cur_order << ", factor size: " << size << std::endl;

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
    for (uint l=1; l <= size; l++) {

      double hyp = card_cost_[l];
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
                << ", factor size: " << size << std::endl;
      std::cerr << "rev_cum_best: " << rev_cum_best[cur_order] << std::endl;
      std::cerr << "cum_best_m1 value: " << cum_best_m1[cur_order-1] << std::endl;


      std::cerr << "cur_order: " << cur_order << std::endl;
      std::cerr << "cur_rel_msg: " << cur_rel_msg << std::endl;
    }

    assert(fabs(best1-val1) < 1e-2);
#endif
    //END_DEBUG

    if (val0 < val1) {
      message_matrix_(0,k) = 0.0;
      message_matrix_(1,k) -= val0;
    }
    else {
      message_matrix_(1,k) = 0.0;
      message_matrix_(0,k) -= val1;
    }
  }
}

/***********************************/

BILPConstraintFactorNode::BILPConstraintFactorNode(const Storage1D<VariableNode*>& participating_vars,
                                                   const Storage1D<bool>& positive, int rhs_lower, int rhs_upper) :
  FactorNode(participating_vars), rhs_lower_(rhs_lower), rhs_upper_(rhs_upper), positive_(positive),
  message_matrix_(2,participating_vars.size(),0.0) {


  int nPositive = 0;
  int nNegative = 0;

  for (uint k=0; k < positive_.size(); k++) {
    if (positive_[k])
      nPositive++;
    else
      nNegative++;
  }

  int lower_bound = -nNegative;
  int upper_bound = nPositive;

  lower_bound = std::max(lower_bound, rhs_lower_ - nPositive);
  upper_bound = std::min(upper_bound, rhs_upper_ + nNegative);

  const int range = upper_bound - lower_bound + 1;
  const int zero_offset = -lower_bound;

  //adjust bounds
  if (rhs_lower_ + zero_offset < 0) {
    rhs_lower_ -= (rhs_lower_ + zero_offset);
  }
  if (rhs_upper_ + zero_offset >= range) {
    rhs_upper_ -= (rhs_upper_ + zero_offset - range +1);
  }

  if (rhs_upper_ < rhs_lower_ ) {
    INTERNAL_ERROR << "constraint is unsatisfiable" << std::endl;
    exit(1);
  }

  assert(participating_vars.size() >= 1);
}

/*virtual*/ BILPConstraintFactorNode::~BILPConstraintFactorNode() {}

/*virtual*/ void BILPConstraintFactorNode::compute_messages() {

  //following [Potetz & Lee CVIU 2007]

  uint nPositive = 0;
  uint nNegative = 0;

  for (uint k=0; k < positive_.size(); k++) {
    if (positive_[k])
      nPositive++;
    else
      nNegative++;
  }

  uint range = positive_.size() + 1;
  uint zero_offset = nNegative;

  Storage1D<const double*> var_message(participating_var_.size());

  for (uint k=0; k < participating_var_.size(); k++)
    var_message[k] = participating_var_[k]->get_message(this);

  /**** forward ****/

  Math3D::Tensor<double> forward(range,2,participating_var_.size(),1e100);
  Math2D::Matrix<double> forward_light(range,participating_var_.size(),1e100);

  //init
  forward(zero_offset,0,0) = var_message[0][0];
  forward_light(zero_offset,0) = var_message[0][0];
  int init_mul = (positive_[0]) ? 1 : -1;
  if (int(zero_offset)+init_mul >= 0
      && int(zero_offset)+init_mul < int(range)) {
    forward(zero_offset+init_mul,1,0) = var_message[0][1];
    forward_light(zero_offset+init_mul,0) = var_message[0][1];
  }

  //proceed
  for (uint v=1; v < participating_var_.size(); v++) {

    for (int sum=0; sum < (int) range; sum++) {

      for (int l=0; l < 2; l++) {
	
        double best_prev = 1e75;
	
        int move = l;
        if (positive_[v]) //since we are tracing backward here
          move *= -1;

        int dest = sum + move;
        if (dest >= 0 && dest < (int) range) {

          double hyp = forward_light(dest,v-1);
          if (hyp < best_prev)
            best_prev = hyp;
        }

        forward(sum,l,v) = best_prev + var_message[v][l];
      }
      forward_light(sum,v) = std::min(forward(sum,0,v),forward(sum,1,v));
    }
  }

  Math3D::Tensor<double> backward(range,2,participating_var_.size(),1e100);
  Math2D::Matrix<double> backward_light(range,participating_var_.size(),1e100);

  /**** backward ****/

  uint last_var = participating_var_.size()-1;

  //init
  backward_light(zero_offset,last_var) = var_message[last_var][0];
  int end_mul = (positive_[last_var]) ? 1 : -1;
  backward_light(zero_offset + end_mul,last_var) = var_message[last_var][1];

  //proceed
  for (int v=last_var-1; v >= 1; v--) {

    for (int sum=0; sum < (int) range; sum++) {

      double best_light = 1e300;

      for (int l=0; l < 2; l++) {
	
        double best_prev = 1e75;
	
        int move = l;
        if (positive_[v]) //since we are tracing backward here
          move *= -1;

        int dest = sum + move;
        if (dest >= 0 && dest < (int) range) {
          double hyp = backward_light(dest,v+1);
          if (hyp < best_prev)
            best_prev = hyp;
        }

        double hyp = best_prev + var_message[v][l];
        if (hyp < best_light)
          best_light = hyp;
      }
      backward_light(sum,v) = best_light; 
    }
  }

  /*** derive messages ***/

  for (uint k=0; k < positive_.size(); k++) {

    double sub = 1e75;

    for (int l=0; l < 2; l++) {

      double min_msg = 1e300;

      for (int s=0; s < (int) range; s++) {

        double hyp = forward(s,l,k);
		
        double best_bwd = 1e300;
        if (k+1 < positive_.size()) {
          for (int r=rhs_lower_; r <= rhs_upper_; r++) {
            int other = r + zero_offset - (s - zero_offset);

            if (other >= 0 && other < (int) range) {
              double cur = backward_light(other,k+1);
              if (cur < best_bwd)
                best_bwd = cur;
            }
          }
        }
        else {
          if (s >= rhs_lower_ + int(zero_offset) && s <= rhs_upper_ + int(zero_offset))
            best_bwd = 0.0;
        }

        hyp += best_bwd;

        if (hyp < min_msg)
          min_msg = hyp;
      }

      assert(!isnan(min_msg));
      assert(!isnan(var_message[k][l]));
	
      message_matrix_(l,k) = min_msg -  var_message[k][l];
      if (message_matrix_(l,k) < sub)
        sub = message_matrix_(l,k);

      if (isnan(message_matrix_(l,k))) {
        std::cerr << "nan = " << min_msg << " - " <<  var_message[k][l] << std::endl;
      }

      assert(!isnan(message_matrix_(l,k)) );
    }

    for (uint l=0; l < 2; l++) {
      message_matrix_(l,k) -= sub;
    }
  }
}

/*virtual*/ bool BILPConstraintFactorNode::process_labeling(Math1D::Vector<uint>& labels) const {

  int eval = 0;

  for (uint k=0; k < labels.size(); k++) {

    if (positive_[k])
      eval += labels[k];
    else
      eval -= int(labels[k]);
  }

  if (eval < rhs_lower_ || eval > rhs_upper_) {

    uint nPositive = 0;
    uint nNegative = 0;
    
    for (uint k=0; k < positive_.size(); k++) {
      if (positive_[k])
        nPositive++;
      else
        nNegative++;
    }

    uint range = positive_.size() + 1;
    int zero_offset = nNegative;


    Math2D::Matrix<double> beliefs(participating_var_.size(),2);

    Math1D::Vector<double> belief(2);
    
    for (uint k=0; k < participating_var_.size(); k++) {
      participating_var_[k]->compute_beliefs(belief);
      beliefs(k,0) = belief[0];
      beliefs(k,1) = belief[1];      
    }


    Math3D::NamedTensor<double> forward(range,2,participating_var_.size(),1e250,MAKENAME(forward));
    Math3D::NamedTensor<uint> trace(range,2,participating_var_.size(),MAX_UINT,MAKENAME(trace));

    //init
    forward(zero_offset,0,0) = beliefs(0,0);
    int init_mul = (positive_[0]) ? 1 : -1;
    forward(zero_offset+init_mul,1,0) = beliefs(0,1);
    
    //proceed
    for (uint v=1; v < participating_var_.size(); v++) {
      
      for (int sum=0; sum < (int) range; sum++) {
	
        for (int l=0; l < 2; l++) {
	  
          double best_prev = 1e300;
          uint arg_best = MAX_UINT;
	
          int move = l;
          if (positive_[v]) //since we are tracing backward here
            move *= -1;
	  
          int dest = sum + move;
          if (dest >= 0 && dest < (int) range) {
            for (int l_prev = 0; l_prev < 2; l_prev ++) {
              double hyp = forward(dest,l_prev,v-1);
              if (hyp < best_prev) {
                best_prev = hyp;
                arg_best = l_prev;
              }
            }
          }
	  
          forward(sum,l,v) = best_prev + beliefs(v,l);
          trace(sum,l,v) = arg_best;
        }
      }
    }


    /**** traceback ***/

    double total_best = 1e300;
    int label = 32000;
    int sum = 32000;

    for (int r = rhs_lower_; r <= rhs_upper_; r++) {
      for (uint l=0; l < 2; l++) {
	
        double hyp = forward(r+zero_offset, l, participating_var_.size()-1);
	
        assert(!isnan(hyp));
	
        if (hyp < total_best) {
          total_best = hyp;
          label = l;
          sum = r;
        }
      }
    }
  
    labels[participating_var_.size()-1] = label;

    for (int v =  int(participating_var_.size())-2; v >= 0; v--) {

      int old_label = label;
      label = trace(sum + zero_offset,label,v+1);
      
      if (positive_[v+1])
        sum -= old_label;
      else
        sum += old_label;

      labels[v] = label;
    }
    

    return true;
  }

  return false;

}

/*virtual*/ double* BILPConstraintFactorNode::get_message(VariableNode* node) {

  double* ptr = message_matrix_.direct_access();

  for (uint i=0; i < participating_var_.size(); i++) {

    if (participating_var_[i] == node)
      break;
    else
      ptr += 2;
  }

  return ptr;
}

/*virtual*/ double BILPConstraintFactorNode::cost(const Math1D::Vector<uint>& labels) const {

  int sum = 0;

  for (uint k=0; k < participating_var_.size(); k++) {

    int label = labels[k];
    if (positive_[k])
      sum += label;
    else
      sum -= label;
  }

  return (sum >= rhs_lower_ && sum <= rhs_upper_) ? 0.0 : 1e15;
}

/*virtual*/ void BILPConstraintFactorNode::init_messages() {
  message_matrix_.set_constant(0.0);
}




/***********************************/

FactorMPBP::FactorMPBP(uint nVars, uint nFactors) : nUsedVars_(0), nUsedFactors_(0) {

  var_.resize(nVars,0);
  factor_.resize(nFactors,0);
}

//delete all allocated owned var. and factor nodes
FactorMPBP::~FactorMPBP() {
  for (uint i=0; i < nUsedVars_; i++)
    delete var_[i];

  for (uint i=0; i < nUsedFactors_; i++)
    delete factor_[i];
}

uint FactorMPBP::add_var(const Math1D::Vector<float>& unary_cost) {

  VariableNode* ptr = new VariableNode(unary_cost);

  uint nPrevVars = nUsedVars_;

  nUsedVars_++;

  uint prev_size = var_.size();
  if (nUsedVars_ > prev_size) {
    var_.resize(size_t(1.2*prev_size)+1,0);
  }
  var_[nPrevVars] = ptr;

  return nPrevVars;
}

VariableNode* FactorMPBP::get_var_node(uint v) {

  if (v >= nUsedVars_) {
    INTERNAL_ERROR << "variable index out of bounds. Exiting..." << std::endl;
  }

  return var_[v];
}

uint FactorMPBP::add_factor(FactorNode* node) {

  nUsedFactors_++;

  uint prev_size = factor_.size();
  if (nUsedFactors_ > prev_size)
    factor_.resize(size_t(1.2*prev_size)+1,0);

  factor_[nUsedFactors_-1] = node;

  return nUsedFactors_-1;
}

uint FactorMPBP::add_generic_factor(const Math1D::Vector<uint> var, VarDimStorage<float>& cost) {

  assert(var.size() == cost.nDims());

  Storage1D<VariableNode*> var_nodes(var.size());
  for (uint k=0; k < var.size(); k++) {
    if (var[k] >= nUsedVars_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }
    var_nodes[k] = var_[var[k]];
  }

  GenericFactorNode* ptr = new GenericFactorNode(var_nodes,cost);
  return add_factor(ptr);
}


uint FactorMPBP::add_generic_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost, bool ref) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  assert(var1 < var_.size());
  assert(var2 < var_.size());

  Storage1D<VariableNode*> var_nodes(2);
  var_nodes[0] = var_[var1];
  var_nodes[1] = var_[var2];

  if (ref) {
    BinaryRefFactorNode* ptr = new BinaryRefFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
  else {
    BinaryFactorNode* ptr = new BinaryFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
}

uint FactorMPBP::add_potts_factor(uint var1, uint var2, double lambda) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  assert(var1 < var_.size());
  assert(var2 < var_.size());

  Storage1D<VariableNode*> var_nodes(2);
  var_nodes[0] = var_[var1];
  var_nodes[1] = var_[var2];

  PottsFactorNode* ptr = new PottsFactorNode(var_nodes, lambda);
  return add_factor(ptr);
}


uint FactorMPBP::add_generic_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_ || var3 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  assert(var1 < var_.size());
  assert(var2 < var_.size());
  assert(var3 < var_.size());

  Storage1D<VariableNode*> var_nodes(3);
  var_nodes[0] = var_[var1];
  var_nodes[1] = var_[var2];
  var_nodes[2] = var_[var3];

  if (!ref) {
    TernaryFactorNode* ptr = new TernaryFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
  else {
    TernaryRefFactorNode* ptr = new TernaryRefFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
}

uint FactorMPBP::add_generic_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                                                 const Storage1D<Math3D::Tensor<float> >& cost, bool ref) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_ || var3 >= nUsedVars_ || var4 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  assert(var1 < var_.size());
  assert(var2 < var_.size());
  assert(var3 < var_.size());
  assert(var4 < var_.size());

  Storage1D<VariableNode*> var_nodes(4);
  var_nodes[0] = var_[var1];
  var_nodes[1] = var_[var2];
  var_nodes[2] = var_[var3];
  var_nodes[3] = var_[var4];

  if (!ref) {
    FourthOrderFactorNode* ptr = new FourthOrderFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
  else {
    FourthOrderRefFactorNode* ptr = new FourthOrderRefFactorNode(var_nodes, cost);
    return add_factor(ptr);
  }
}


uint FactorMPBP::add_one_of_N_factor(const Math1D::Vector<uint>& var) {

  Storage1D<VariableNode*> var_nodes(var.size());

  for (uint k=0; k < var.size(); k++) {

    assert(var[k] < var_.size());
    var_nodes[k] = var_[var[k]];
  }

  OneOfNFactorNode* ptr = new OneOfNFactorNode(var_nodes);

  return add_factor(ptr);
}

uint FactorMPBP::add_cardinality_factor(const Math1D::Vector<uint>& var, const Math1D::Vector<float>& card_cost) {

  assert(card_cost.size() == var.size()+1);

  Storage1D<VariableNode*> var_nodes(var.size());

  for (uint k=0; k < var.size(); k++) {

    if (var[k] >= nUsedVars_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    assert(var[k] < var_.size());
    var_nodes[k] = var_[var[k]];

    if (var_nodes[k]->nLabels() != 2) {
      INTERNAL_ERROR << " variables of 1-of-N nodes must be binary. Exiting..." << std::endl;
      exit(1);
    }
  }

  CardinalityFactorNode* ptr = new CardinalityFactorNode(var_nodes, card_cost);

  return add_factor(ptr);
}

uint FactorMPBP::add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive, 
                                       int rhs_lower, int rhs_upper) {

  assert(var.size() == positive.size());

  Storage1D<VariableNode*> var_nodes(var.size());

  for (uint k=0; k < var.size(); k++) {

    if (var[k] >= nUsedVars_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    assert(var[k] < var_.size());
    var_nodes[k] = var_[var[k]];

    if (var_nodes[k]->nLabels() != 2) {
      INTERNAL_ERROR << " variables of 1-of-N nodes must be binary. Exiting..." << std::endl;
      exit(1);
    }
  }

  BILPConstraintFactorNode* ptr = new BILPConstraintFactorNode(var_nodes,positive,rhs_lower,rhs_upper);

  return add_factor(ptr);
}

//NOTE: after calling this routine, owned factors can no longer be created
uint FactorMPBP::pass_in_factor_node(FactorNode* factor) {

  return add_factor(factor);
}

FactorNode* FactorMPBP::get_factor(uint f) {

  if (f >= nUsedFactors_) {
    INTERNAL_ERROR << "factor index out of bounds. Exiting..." << std::endl;
    exit(1);
  }

  return factor_[f];
}

const Math1D::Vector<uint>& FactorMPBP::labeling() {

  return labeling_;
}

double FactorMPBP::labeling_energy() {

  double energy = 0.0;

  std::map<VariableNode*,uint> label;

  for (uint k=0; k < nUsedVars_; k++) {
    assert(label.find(var_[k]) == label.end());
    label[var_[k]] = labeling_[k];
    energy += var_[k]->cost(labeling_[k]);
  }

  std::cerr << "unary cost: " << energy << std::endl;

  std::cerr << "sum of labeling: " << labeling_.sum() << std::endl;

  for (uint k=0; k < nUsedFactors_; k++) {

    const Storage1D<VariableNode*>& nodes = factor_[k]->participating_nodes();

    Math1D::Vector<uint> sub_labeling(nodes.size());
    for (uint i=0; i < nodes.size(); i++) {
      assert(label.find(nodes[i]) != label.end());
      sub_labeling[i] = label[nodes[i]];
    }

    energy += factor_[k]->cost(sub_labeling);
  }
  
  return energy;
}

void FactorMPBP::process_labeling() {

  std::map<VariableNode*,uint> label;
  std::map<VariableNode*,uint> num;

  for (uint k=0; k < nUsedVars_; k++) {
    assert(label.find(var_[k]) == label.end());
    label[var_[k]] = labeling_[k];
    num[var_[k]] = k;
  }

  for (uint k=0; k < nUsedFactors_; k++) {

    const Storage1D<VariableNode*>& nodes = factor_[k]->participating_nodes();

    Math1D::Vector<uint> sub_labeling(nodes.size());
    for (uint i=0; i < nodes.size(); i++) {
      assert(label.find(nodes[i]) != label.end());
      sub_labeling[i] = label[nodes[i]];
    }

    bool changed = factor_[k]->process_labeling(sub_labeling);

    if (changed) {
      for (uint i=0; i < nodes.size(); i++) {
        labeling_[num[nodes[i]]] = sub_labeling[i];
        label[nodes[i]] = sub_labeling[i];
      }    
    }
  }
}

void FactorMPBP::icm(uint nIter) {

  std::map<VariableNode*,uint> label;
  std::map<VariableNode*,uint> num;

  for (uint k=0; k < nUsedVars_; k++) {
    assert(label.find(var_[k]) == label.end());
    label[var_[k]] = labeling_[k];
    num[var_[k]] = k;
  }

  for (uint iter = 1; iter <= nIter; iter++) {

    std::cerr << "ICM iter " << iter << std::endl;

    for (uint v=0; v < nUsedVars_; v++) {

      const Storage1D<FactorNode*>& factors = var_[v]->neighboring_factor();

      Storage1D<Math1D::Vector<uint> > sublabeling(factors.size());
      
      Math1D::Vector<uint> var_index(factors.size(),MAX_UINT);

      for (uint f=0; f < factors.size(); f++) {

        const Storage1D<VariableNode*>& factor_var = factors[f]->participating_nodes();

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

/**** run inference ***/
void FactorMPBP::mpbp(uint nIter, bool quiet) {

  Math1D::Vector<double> belief;

  labeling_.resize(nUsedVars_,0);

  //init messages
  for (uint v=0; v < nUsedVars_; v++)
    var_[v]->init_messages();

  for (uint f=0; f < nUsedFactors_; f++)
    factor_[f]->init_messages();

  for (uint iter = 1; iter <= nIter; iter++) {

    std::cerr << "***** iteration " << iter << " ******" << std::endl;

    for (uint v=0; v < nUsedVars_; v++) {
      var_[v]->compute_messages();
    }

    for (uint f=0; f < nUsedFactors_; f++)
      factor_[f]->compute_messages();

    //compute labeling (minimal neg. log-beliefs)
  
    double min_best_belief = 1e300;
    double max_best_belief = -1e300;

    for (uint v=0; v < nUsedVars_; v++) {
      var_[v]->compute_beliefs(belief);
      
      assert(belief.size() == var_[v]->nLabels());

      double min_cost = 1e300;
      uint arg_min = MAX_UINT;
      for (uint k=0; k < belief.size(); k++) {
        if (belief[k] < min_cost) {
	  
          min_cost = belief[k];
          arg_min = k;
        }
      }

      min_best_belief = std::min(min_best_belief,min_cost);
      max_best_belief = std::max(max_best_belief,min_cost);

      assert(arg_min < belief.size());

      labeling_[v] = arg_min;
    }

    if (!quiet)
      std::cerr << "belief range: [" << min_best_belief << "," << max_best_belief << "]" << std::endl;
    
    process_labeling();

    if (!quiet) {
      double label_energy = labeling_energy();
      std::cerr << "labeling cost: " << label_energy << std::endl;

      //icm(5);
      //std::cerr << "after ICM: " << labeling_energy() << std::endl;
    }
  }

}

