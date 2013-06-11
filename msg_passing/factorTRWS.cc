/**** written by Thomas Schoenemann. Started as a private person without employment, August 2011 *****/
/**** continued at the University of Pisa, Italy, October - December 2011 ****/
/**** and at the Univeristy of Düsseldorf, Germany, 2012 ****/


#include "factorTRWS.hh"

#include <set>
#include "stl_out.hh"


CumTRWSVar::CumTRWSVar(const Math1D::Vector<float>& cost, uint rank) : cost_(cost), rank_(rank), nChains_(1) {
  cum_cost_.resize(cost.size());
  for (uint k=0; k < cum_cost_.size(); k++) {
    //convert float->double
    cum_cost_[k] = cost[k];
  }
}

void CumTRWSVar::add_cost(const Math1D::Vector<float>& add_cost) {

  if (add_cost.size() != cost_.size()) {
    INTERNAL_ERROR << "cannot add cost due to incompatible vector sizes: " << cost_.size() << " and " << add_cost.size() << std::endl;
    exit(1);
  }

  for (uint i=0; i < cost_.size(); i++) {
    cost_[i] += add_cost[i];
    cum_cost_[i] += add_cost[i] / float(nChains_);
  }
}

void CumTRWSVar::set_up_chains() {

  nChains_ = 0;

  {
    uint nOutgoing = 0;
    uint nMiddle = 0;
    uint nIncoming = 0;

    for (uint k=0; k < adjacent_factor_.size(); k++) {

      if (adjacent_factor_[k]->min_rank() == rank_)
        nOutgoing++;
      else if (adjacent_factor_[k]->max_rank() == rank_)
        nIncoming++;
      else
        nMiddle++;
    }

    uint nChainlessChains = nMiddle + std::max(nOutgoing,nIncoming);

    nChains_ = nChainlessChains;
  }


  //initialize cost
  for (uint l=0; l < cost_.size(); l++)
    cum_cost_[l] = cost_[l] / nChains_;
}

void CumTRWSVar::add_factor(CumTRWSFactor* factor) {

  uint size = adjacent_factor_.size();
  
  adjacent_factor_.resize(size+1);
  adjacent_factor_[size] = factor;
}

double CumTRWSVar::average(uint& arg_min) {


  const uint nLabels = cost_.size();
  const uint nFactors = adjacent_factor_.size();

  for (uint l=0; l < nLabels; l++)
    cum_cost_[l] = cost_[l]; //convert float -> double

  for (uint k=0; k < nFactors; k++)
    cum_cost_ += adjacent_factor_[k]->reparameterization(this);

  double offs = 1e300;

  for (uint l=0; l < nLabels; l++) {
    if (cum_cost_[l] < offs) {
      offs = cum_cost_[l];
      arg_min = l;
    }
  }

  for (uint l=0; l < nLabels; l++) {
    cum_cost_[l] -= offs;
  }

  cum_cost_ *= 1.0 / nChains_;

  return offs;
}

uint CumTRWSVar::rank() const {
  return rank_;
}

void CumTRWSVar::set_rank(uint rank) {
  rank_ = rank;
}

size_t CumTRWSVar::nLabels() const {
  return cost_.size();
}

const Storage1D<CumTRWSFactor*>& CumTRWSVar::adjacent_factor() const {
  return adjacent_factor_;
}

const Math1D::Vector<double>& CumTRWSVar::cost() const {
  return cum_cost_;
}

/**************************/

CumTRWSFactor::CumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars) :
  involved_var_(involved_vars)  {

  uint nVars = involved_vars.size();
  assert(nVars >= 2);

  reparameterization_.resize(nVars);

  for (uint v=0; v < nVars; v++) {

    uint nLabels = involved_vars[v]->nLabels();
    reparameterization_[v].resize(nLabels,0.0);

    involved_vars[v]->add_factor(this);
  }

  compute_rank_range();
}

/*virtual*/ CumTRWSFactor::~CumTRWSFactor() {}

/*virtual*/ void CumTRWSFactor::init() {
  //by default no action
}
  
uint CumTRWSFactor::min_rank() const {
  return min_rank_;
}

uint CumTRWSFactor::max_rank() const {
  return max_rank_;
}

void CumTRWSFactor::compute_rank_range() {

  min_rank_ = MAX_UINT;
  max_rank_ = 0;

  const uint nVars = involved_var_.size();

  for (uint v=0; v < nVars; v++) {

    min_rank_ = std::min(min_rank_,involved_var_[v]->rank());
    max_rank_ = std::max(max_rank_,involved_var_[v]->rank());
  }
}

void CumTRWSFactor::sort_by_rank() {

  const uint nVars = involved_var_.size();

  for (uint k=0; k < nVars-1; k++) {
    for (uint l=0; l < nVars-1-k; l++) {
      if (involved_var_[l]->rank() > involved_var_[l+1]->rank())
        std::swap(involved_var_[l],involved_var_[l+1]);
    }
  }

}

/*virtual*/ uint CumTRWSFactor::best_of_n() {

  TODO("best_of_n");
}

const Storage1D<CumTRWSVar*>& CumTRWSFactor::involved_vars() const {
  return involved_var_;
}

const Math1D::Vector<double>& CumTRWSFactor::reparameterization(const CumTRWSVar* var) const {

  for (uint k=0; k < involved_var_.size(); k++) {

    if (involved_var_[k] == var)
      return reparameterization_[k];
  }

  assert(false);
  return reparameterization_[0];
}

/********************/

GenericCumTRWSFactor::GenericCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, const VarDimStorage<float>& cost)
  : CumTRWSFactor(involved_vars), cost_(cost) {

  if (cost.nDims() != involved_vars.size()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
    exit(1);
  }

  for (uint v = 0; v < involved_vars.size(); v++) {
    if (cost.dim(v) < involved_vars[v]->nLabels()) {
      INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
      exit(1);
    }
  }
}

/*virtual*/ GenericCumTRWSFactor::~GenericCumTRWSFactor() {}

/*virtual*/ 
double GenericCumTRWSFactor::compute_reparameterization(const CumTRWSVar* var) {

  const uint nVars = involved_var_.size();

  Math1D::NamedVector<uint> nLabels(nVars, MAKENAME(nLabels));

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] == var)
      idx = k;

    param[k] -= involved_var_[k]->cost();
    nLabels[k] = involved_var_[k]->nLabels();
  }

  Math1D::Vector<double>& cur_repar = reparameterization_[idx];
  cur_repar.set_constant(1e300);

  Math1D::Vector<size_t> labeling(nVars,0);

  while (true) {
    
    double cost = cost_(labeling);
    for (uint v=0; v < nVars; v++) {
      if (v != idx)
        cost -= param[v][labeling[v]];
    }

    if (cost < cur_repar[labeling[idx]]) {
      cur_repar[labeling[idx]] = cost;
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


  double msg_offs = cur_repar.min();

  for (uint l=0; l < cur_repar.size(); l++)
    cur_repar[l] -= msg_offs;
  
  return msg_offs;
}

/********************/


BinaryCumTRWSFactorBase::BinaryCumTRWSFactorBase(const Storage1D<CumTRWSVar*>& involved_vars)
  : CumTRWSFactor(involved_vars) {

  if (involved_vars.size() != 2 ){
    INTERNAL_ERROR << "attempt to instantiate a binary factor with " << involved_vars.size() << " variables. Exiting." << std::endl;
    exit(1);
  }
}


double BinaryCumTRWSFactorBase::compute_reparameterization(const CumTRWSVar* var, const Math2D::Matrix<float>& cost) {

  double offs = 0.0;

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();

  //this routine also updates reparameterization_
  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++) {
    param[k] -= involved_var_[k]->cost();
  }

  uint idx = 0;

  if (var == involved_var_[0]) {

    idx = 0;

    for (uint l1 = 0; l1 < nLabels1; l1++) {

      double best = 1e300;

      for (uint l2 = 0; l2 < nLabels2; l2++) {

        double hyp = cost(l1,l2) - param[1][l2];

        if (hyp < best)
          best = hyp;
      }

      reparameterization_[idx][l1] = best;
    }
  }
  else {
    assert(var == involved_var_[1]);

    idx = 1;

    for (uint l2 = 0; l2 < nLabels2; l2++) {

      double best = 1e300;

      for (uint l1 = 0; l1 < nLabels1; l1++) {

        double hyp = cost(l1,l2) - param[0][l1];

        if (hyp < best)
          best = hyp;
      }

      reparameterization_[idx][l2] = best;
    }
  }

  double msg_offs = reparameterization_[idx].min();
  offs += msg_offs;

  for (uint l=0; l < reparameterization_[idx].size(); l++)
    reparameterization_[idx][l] -= msg_offs;
  
  return offs;
}

/********************/

BinaryCumTRWSFactor::BinaryCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                                         const Math2D::Matrix<float>& cost) : 
  BinaryCumTRWSFactorBase(involved_vars), cost_(cost) {

  if (cost_.xDim() < involved_vars[0]->nLabels() || cost_.yDim() < involved_vars[1]->nLabels()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
  }
}

/*virtual*/ BinaryCumTRWSFactor::~BinaryCumTRWSFactor() {}

/*virtual*/ double BinaryCumTRWSFactor::compute_reparameterization(const CumTRWSVar* var) {
  
  return BinaryCumTRWSFactorBase::compute_reparameterization(var,cost_);
}

/********************/

BinaryCumTRWSRefFactor::BinaryCumTRWSRefFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                                               const Math2D::Matrix<float>& cost) : 
  BinaryCumTRWSFactorBase(involved_vars), cost_(cost) {
}

/*virtual*/ BinaryCumTRWSRefFactor::~BinaryCumTRWSRefFactor() {}

/*virtual*/ double BinaryCumTRWSRefFactor::compute_reparameterization(const CumTRWSVar* var) {
  
  assert(cost_.xDim() >= involved_var_[0]->nLabels());
  assert(cost_.yDim() >= involved_var_[1]->nLabels());

  return BinaryCumTRWSFactorBase::compute_reparameterization(var,cost_);
}

/********************/

TernaryCumTRWSFactorBase::TernaryCumTRWSFactorBase(const Storage1D<CumTRWSVar*>& involved_vars) :
  CumTRWSFactor(involved_vars) {

  if (involved_vars.size() != 3 ){
    INTERNAL_ERROR << "attempt to instantiate a ternary factor with " << involved_vars.size() << " variables. Exiting." << std::endl;
    exit(1);
  }
}
  
double TernaryCumTRWSFactorBase::compute_reparameterization(const CumTRWSVar* var, const Math3D::Tensor<float>& cost) {

  double offs = 0.0;

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  
  //this routine also updates reparameterization_
  uint idx = 0;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++)
    param[k] -= involved_var_[k]->cost();
  
#ifndef NDEBUG
  const Math1D::Vector<double>& param0 = param[0];
  const Math1D::Vector<double>& param1 = param[1];
  const Math1D::Vector<double>& param2 = param[2];
#else
  const double* param0 = param[0].direct_access();
  const double* param1 = param[1].direct_access();
  const double* param2 = param[2].direct_access();
#endif

  if (var == involved_var_[0]) {
    
    idx = 0;

    Math1D::Vector<double>& message = reparameterization_[idx];

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;

      for (uint l2 = 0; l2 < nLabels2; l2++) {

        const double w2 = param1[l2];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
          double hyp = cost(l1,l2,l3) - w2 - param2[l3];
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      message[l1] = best;
    }
  }
  else if (var == involved_var_[1]) {

    idx = 1;

    Math1D::Vector<double>& message = reparameterization_[idx];

    for (uint l2 = 0; l2 < nLabels2; l2++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const double w1 = param0[l1];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
          double hyp = cost(l1,l2,l3) - w1 - param2[l3];
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      message[l2] = best;
    }    
  }
  else {
    assert(var == involved_var_[2]);
    
    idx = 2;

    Math1D::Vector<double>& message = reparameterization_[idx];

    for (uint l3 = 0; l3 < nLabels3; l3++) {
      
      double best = 1e300;
      
      for (uint l2 = 0; l2 < nLabels2; l2++) {

        const double w2 = param1[l2];

        for (uint l1 = 0; l1 < nLabels1; l1++) {

          double hyp = cost(l1,l2,l3) - w2  - param0[l1];
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      message[l3] = best;
    }
  }

  double msg_offs = reparameterization_[idx].min();
  offs += msg_offs;

  for (uint l=0; l < reparameterization_[idx].size(); l++)
    reparameterization_[idx][l] -= msg_offs;

  return offs;  
}

/********************/

TernaryCumTRWSFactor::TernaryCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, const Math3D::Tensor<float>& cost) :
  TernaryCumTRWSFactorBase(involved_vars), cost_(cost) {

  if (cost_.xDim() < involved_vars[0]->nLabels() || cost_.yDim() < involved_vars[1]->nLabels()
      || cost_.zDim() < involved_vars[2]->nLabels()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
  }
}

/*virtual*/ TernaryCumTRWSFactor::~TernaryCumTRWSFactor() {}

/*virtual*/ double TernaryCumTRWSFactor::compute_reparameterization(const CumTRWSVar* var) { 

  return TernaryCumTRWSFactorBase::compute_reparameterization(var,cost_);
}

/********************/

TernaryCumTRWSRefFactor::TernaryCumTRWSRefFactor(const Storage1D<CumTRWSVar*>& involved_vars, const Math3D::Tensor<float>& cost) :
  TernaryCumTRWSFactorBase(involved_vars), cost_(cost) {}

/*virtual*/ TernaryCumTRWSRefFactor::~TernaryCumTRWSRefFactor() {}

/*virtual*/ double TernaryCumTRWSRefFactor::compute_reparameterization(const CumTRWSVar* var) { 

  assert(cost_.xDim() >= involved_var_[0]->nLabels());
  assert(cost_.yDim() >= involved_var_[1]->nLabels());
  assert(cost_.zDim() >= involved_var_[2]->nLabels());

  return TernaryCumTRWSFactorBase::compute_reparameterization(var,cost_);
}

/********************/

SecondDiffCumTRWSFactor::SecondDiffCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, float lambda) :
  CumTRWSFactor(involved_vars), lambda_(lambda) {

  if (involved_vars.size() != 3) {
    INTERNAL_ERROR << "attempt to instantiate a second difference factor with " << involved_vars.size() << " variables" << std::endl;
    exit(1);
  }
}

/*virtual*/ SecondDiffCumTRWSFactor::~SecondDiffCumTRWSFactor() {}

/*virtual*/ double SecondDiffCumTRWSFactor::compute_reparameterization(const CumTRWSVar* var) {

  double offs = 0.0;

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  
  //this routine also updates reparameterization_
  uint idx = 0;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++)
    param[k] -= involved_var_[k]->cost();
  
  if (var == involved_var_[0]) {

    idx = 0;

    Math1D::Vector<double>& message = reparameterization_[idx];

    double base_cost = -param[1].max() - param[2].max() + 3*lambda_;

    for (int l1=0; l1 < int(nLabels1); l1++) {
      
      double best = base_cost;

      for (int l2 = std::max(0,l1-1); l2 <= std::min<int>(nLabels2-1,l1+1); l2++) {

        const double w2 = param[1][l2];

        const int part = - 2*l2 + l1;

        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {
	
          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = l3 + part; 
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -w2 - param[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = -w2 - param[2][l3] + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      message[l1] = best;
    }
  }
  else if (var == involved_var_[1]) {

    idx = 1;

    Math1D::Vector<double>& message = reparameterization_[idx];

    double base_cost = -param[0].max() - param[2].max() + 3*lambda_;

    for (int l2=0; l2 < int(nLabels2); l2++) {
      
      double best = base_cost;

      for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {

        const double w1 = param[0][l1];

        const int part = - 2*l2 + l1;

        for (int l3 = std::max(0,l2-1); l3 <= std::min<int>(nLabels3-1,l2+1); l3++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = l3 + part; 
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -w1 - param[2][l3];
          }
          else if (abs(so_diff) <= 1) {
            hyp = -w1 - param[2][l3] + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      message[l2] = best;
    }    

  }
  else {
    assert(var == involved_var_[2]);

    idx = 2;

    Math1D::Vector<double>& message = reparameterization_[idx];

    double base_cost = -param[0].max() - param[1].max() + 3*lambda_;

    for (int l3=0; l3 < int(nLabels3); l3++) {
      
      double best = base_cost;

      for (int l2 = std::max(0,l3-1); l2 <= std::min<int>(nLabels2-1,l3+1); l2++) {

        const double w2 = param[1][l2];

        const int part = l3 - 2*l2;

        for (int l1 = std::max(0,l2-1); l1 <= std::min<int>(nLabels1-1,l2+1); l1++) {

          assert(abs(l2-l1) <= 1);
          assert(abs(l3-l2) <= 1);

          const int so_diff = part + l1; 
	  
          double hyp = 1e300;

          if (so_diff == 0) {
            hyp = -param[0][l1] - w2;
          }
          else if (abs(so_diff) <= 1) {
            hyp = -param[0][l1] - w2 + lambda_;
          }

          if (hyp < best)
            best = hyp;
        }
      }

      message[l3] = best;
    }
  }

  double msg_offs = reparameterization_[idx].min();
  offs += msg_offs;

  for (uint l=0; l < reparameterization_[idx].size(); l++)
    reparameterization_[idx][l] -= msg_offs;

  return offs;  
}


/********************/

FourthOrderCumTRWSFactorBase::FourthOrderCumTRWSFactorBase(const Storage1D<CumTRWSVar*>& involved_vars)
  : CumTRWSFactor(involved_vars) {

  if (involved_vars.size() != 4) {
    INTERNAL_ERROR << "attempt to instantiate a 4th order factor with " << involved_vars.size() << " variables" << std::endl;
    exit(1);    
  }
}

double FourthOrderCumTRWSFactorBase::compute_reparameterization(const CumTRWSVar* var, const Storage1D<Math3D::Tensor<float> >& cost) {

  double offs = 0.0;
  
  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  const uint nLabels4 = involved_var_[3]->nLabels();
  
  //this routine also updates reparameterization_

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < involved_var_.size(); k++)
    param[k] -= involved_var_[k]->cost();

#ifndef NDEBUG
  const Math1D::Vector<double>& param0 = param[0];
  const Math1D::Vector<double>& param1 = param[1];
  const Math1D::Vector<double>& param2 = param[2];
  const Math1D::Vector<double>& param3 = param[3];
#else
  const double* param0 = param[0].direct_access();
  const double* param1 = param[1].direct_access();
  const double* param2 = param[2].direct_access();
  const double* param3 = param[3].direct_access();
#endif


  uint idx = 0;

  if (var == involved_var_[0]) {

    idx = 0;

    Math1D::Vector<double>& message = reparameterization_[idx];

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;

      const Math3D::Tensor<float>& cur_cost = cost[l1];

      for (uint l3 = 0; l3 < nLabels3; l3++) {

        const double w3 = param2[l3];

        for (uint l2 = 0; l2 < nLabels2; l2++) {

          const double sum2 = w3 + param1[l2];

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cur_cost(l2,l3,l4) 
              - sum2 - param3[l4];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }

      
      message[l1] = best;
    }
  }
  else if (var == involved_var_[1]) {

    idx = 1;

    Math1D::Vector<double>& message = reparameterization_[idx];
    
    for (uint l2 = 0; l2 < nLabels2; l2++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost[l1];
	
        const double w1 = param[0][l1];

        for (uint l3 = 0; l3 < nLabels3; l3++) {

          const double sum3 = w1 + param2[l3];

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cur_cost(l2,l3,l4) 
              - sum3 - param[3][l4];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      message[l2] = best;
    }    
  }
  else if (var == involved_var_[2]) {
    
    idx = 2;

    Math1D::Vector<double>& message = reparameterization_[idx];

    for (uint l3 = 0; l3 < nLabels3; l3++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost[l1];
	
        const double w1 = param0[l1];

        for (uint l2 = 0; l2 < nLabels2; l2++) {

          const double sum2 = w1 + param1[l2];

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cur_cost(l2,l3,l4) 
              - sum2 - param3[l4];
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      message[l3] = best;
    }
  }
  else {
    
    assert(var == involved_var_[3]);

    idx = 3;

    Math1D::Vector<double>& message = reparameterization_[idx];

    for (uint l4 = 0; l4 < nLabels4; l4++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost[l1];

        const double w1 = param0[l1];

        for (uint l3 = 0; l3 < nLabels3; l3++) {

          const double sum3 = w1 + param2[l3];
  
          for (uint l2 = 0; l2 < nLabels2; l2++) {

            double hyp = cur_cost(l2,l3,l4) 
              - sum3 - param1[l2];
	  
            if (hyp < best)
              best = hyp;
          }
        }

      }
      
      message[l4] = best;
    }    
  }

  double msg_offs = reparameterization_[idx].min();
  offs += msg_offs;

  for (uint l=0; l < reparameterization_[idx].size(); l++)
    reparameterization_[idx][l] -= msg_offs;
  
  return offs;
}

/********************/


FourthOrderCumTRWSFactor::FourthOrderCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                                                   const Storage1D<Math3D::Tensor<float> >& cost) :
  FourthOrderCumTRWSFactorBase(involved_vars), cost_(cost) {

  if (cost_.size() < involved_vars[0]->nLabels() || cost_[0].xDim() < involved_vars[1]->nLabels()
      || cost_[0].yDim() < involved_vars[2]->nLabels() || cost_[0].zDim() < involved_vars[3]->nLabels()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
  }
}

/*virtual*/ FourthOrderCumTRWSFactor::~FourthOrderCumTRWSFactor() {}

/*virtual*/ double FourthOrderCumTRWSFactor::compute_reparameterization(const CumTRWSVar* var) {

  return  FourthOrderCumTRWSFactorBase::compute_reparameterization(var,cost_);
}

/********************/

FourthOrderCumTRWSRefFactor::FourthOrderCumTRWSRefFactor(const Storage1D<CumTRWSVar*>& involved_vars,
                                                         const Storage1D<Math3D::Tensor<float> >& cost) :
  FourthOrderCumTRWSFactorBase(involved_vars), cost_(cost) {
}

/*virtual*/ FourthOrderCumTRWSRefFactor::~FourthOrderCumTRWSRefFactor() {}

/*virtual*/ double FourthOrderCumTRWSRefFactor::compute_reparameterization(const CumTRWSVar* var) {

  assert(cost_.size() >= involved_var_[0]->nLabels());
  assert(cost_[0].xDim() >= involved_var_[1]->nLabels());
  assert(cost_[0].yDim() >= involved_var_[2]->nLabels());
  assert(cost_[0].zDim() >= involved_var_[3]->nLabels());

  return  FourthOrderCumTRWSFactorBase::compute_reparameterization(var,cost_);
}

/********************/

OneOfNCumTRWSFactor::OneOfNCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars) :
  CumTRWSFactor(involved_vars) {

  for (uint v=0; v < involved_vars.size(); v++) {
    if (involved_vars[v]->nLabels() != 2) {
      INTERNAL_ERROR << "instantiation of a 1-of-N factor with non-binary variables. Exiting." << std::endl;
      exit(1);
    }
  }
}

/*virtual*/ OneOfNCumTRWSFactor::~OneOfNCumTRWSFactor() {}

/*virtual*/ uint OneOfNCumTRWSFactor::best_of_n() const {

  const uint nVars = involved_var_.size();

  double best = 1e300;
  uint arg_best = MAX_UINT;

  for (uint k=0; k < nVars; k++) {

    const double hyp = involved_var_[k]->cost()[1] - reparameterization_[k][1]
      - involved_var_[k]->cost()[0] + reparameterization_[k][0];

    if (hyp < best) {
      best = hyp;
      arg_best = k;
    }
  }

  return arg_best;
}

/*virtual*/ double OneOfNCumTRWSFactor::compute_reparameterization(const CumTRWSVar* var) {

  const uint nVars = involved_var_.size();

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->cost();
    else {
      idx = k;
    }
  }

  Math1D::Vector<double>& message = reparameterization_[idx];

  assert(idx < nVars);

  double best_gain = 1e300;

  double sum = 0.0;

  for (uint i=0; i < nVars; i++) {

    if (i != idx) {
      const double hyp = -param[i][1] + param[i][0];

      if (hyp < best_gain) 
        best_gain = hyp;
      
      sum -= param[i][0];
    }
  }

  message[0] = sum + best_gain;
  message[1] = sum;

  const double offs = message.min();

  for (uint k=0; k < 2; k++)
    message[k] -= offs;

  return offs;
}

/********************/

OneOfNCumTRWSFactorWithReuse::OneOfNCumTRWSFactorWithReuse(const Storage1D<CumTRWSVar*>& involved_vars) :
  CumTRWSFactor(involved_vars), to_update_(MAX_UINT) {

  for (uint v=0; v < involved_vars.size(); v++) {
    if (involved_vars[v]->nLabels() != 2) {
      INTERNAL_ERROR << "instantiation of a 1-of-N factor with non-binary variables. Exiting." << std::endl;
      exit(1);
    }
  }
}

/*virtual*/ OneOfNCumTRWSFactorWithReuse::~OneOfNCumTRWSFactorWithReuse() {}

/*virtual*/
void OneOfNCumTRWSFactorWithReuse::init() {
  sort_by_rank();
  to_update_ = MAX_UINT;
}

/*virtual*/ uint OneOfNCumTRWSFactorWithReuse::best_of_n() const {

  const uint nVars = involved_var_.size();

  double best = 1e300;
  uint arg_best = MAX_UINT;

  for (uint k=0; k < nVars; k++) {

    const double hyp = involved_var_[k]->cost()[1] - reparameterization_[k][1]
      - involved_var_[k]->cost()[0] + reparameterization_[k][0];

    if (hyp < best) {
      best = hyp;
      arg_best = k;
    }
  }

  return arg_best;
}


/*virtual*/ 
double OneOfNCumTRWSFactorWithReuse::compute_reparameterization(const CumTRWSVar* var) {

  uint idx = MAX_UINT;

  const uint nVars = involved_var_.size();

  if (to_update_ != MAX_UINT) {

    if (var == involved_var_[to_update_]) {
      idx = to_update_;
      TODO("should not happen if we check correctly for minimum (forward) and maximum (backward)");
    }
    else {
      if (to_update_+1 < nVars && involved_var_[to_update_+1] == var)
        idx = to_update_+1;
      else if (to_update_ > 0 && involved_var_[to_update_-1] == var)
        idx = to_update_-1;
      else {
        TODO("should not happen");
      }

      sum_ -= involved_var_[idx]->cost()[0] - reparameterization_[idx][0];
      const double update_val = involved_var_[to_update_]->cost()[0] - reparameterization_[to_update_][0];
      sum_ += update_val;

      enum {DidUpdate, Update, Recompute} status = Update;

      if (to_update_ == arg_best_) {

        const double hyp = involved_var_[to_update_]->cost()[1] - reparameterization_[to_update_][1]
          - update_val;

        if (hyp <= second_best_) {

          best_ = hyp;
          status = DidUpdate;
        }
        else
          status = Recompute;
      }
      else if (to_update_ == arg_second_best_) {

        const double hyp = involved_var_[to_update_]->cost()[1] - reparameterization_[to_update_][1]
          - update_val;
        
        if (hyp <= second_best_ && hyp >= best_) {
          second_best_ = hyp;
          status = DidUpdate;
        }
        else if (hyp < best_) {
          second_best_ = best_;
          arg_second_best_ = arg_best_;
          best_ = hyp;
          arg_best_ = to_update_;
          status = DidUpdate;
        }
        else
          status = Recompute;
      }

      if (status == Recompute) {

        //NOTE: for numerical stability we might want to recompute sum_ here as well

        Storage1D<Math1D::Vector<double> > param = reparameterization_;
        for (uint k=0; k < nVars; k++) {
          if (k != idx)
            param[k] -= involved_var_[k]->cost();
        }

        best_ = 1e300;
        second_best_ = 1e300;
        arg_best_ = MAX_UINT;
        arg_second_best_ = MAX_UINT;
        
        for (uint i=0; i < nVars; i++) {
          
          if (i != idx) {
            const double hyp = -param[i][1] + param[i][0];
            
            if (hyp < best_) {
              second_best_ = best_;
              arg_second_best_ = arg_best_;
              
              best_ = hyp;
              arg_best_ = i;
            }
            else if (hyp < second_best_) {
              second_best_ = hyp;
              arg_second_best_ = i;
            }
          }
        }
      }
      else if (status == Update) {

        assert(to_update_ != arg_second_best_);
        assert(to_update_ != arg_best_);

        const double hyp = involved_var_[to_update_]->cost()[1] - reparameterization_[to_update_][1]
          - update_val;
        
        if (hyp < best_) {
          second_best_ = best_;
          arg_second_best_ = arg_best_;
          
          best_ = hyp;
          arg_best_ = to_update_;
        }
        else if (hyp < second_best_) {
          second_best_ = hyp;
          arg_second_best_ = to_update_;
        }
      }
    }
  }
  else {

    Storage1D<Math1D::Vector<double> > param = reparameterization_;
    for (uint k=0; k < nVars; k++) {
      if (involved_var_[k] != var)
        param[k] -= involved_var_[k]->cost();
      else {
        idx = k;
      }
    }

    sum_ = 0.0;
    best_ = 1e300;
    second_best_ = 1e300;
    arg_best_ = MAX_UINT;
    arg_second_best_ = MAX_UINT;

    for (uint i=0; i < nVars; i++) {

      if (i != idx) {
        const double hyp = -param[i][1] + param[i][0];
        
        if (hyp < best_) {
          second_best_ = best_;
          arg_second_best_ = arg_best_;

          best_ = hyp;
          arg_best_ = i;
        }
        else if (hyp < second_best_) {
          second_best_ = hyp;
          arg_second_best_ = i;
        }
      
        sum_ -= param[i][0];
      }
    }

    assert(arg_best_ != idx && arg_best_ != MAX_UINT);
  }

  Math1D::Vector<double>& message = reparameterization_[idx];
    
  message[0] = sum_ + ((idx != arg_best_) ? best_ : second_best_);
  message[1] = sum_;

  const double offs = message.min();

  for (uint k=0; k < 2; k++)
    message[k] -= offs;

  to_update_ = idx;

  return offs;
}

/********************/

CardinalityCumTRWSFactorBase::CardinalityCumTRWSFactorBase(const Storage1D<CumTRWSVar*>& involved_vars)
  : CumTRWSFactor(involved_vars) {

  for (uint v=0; v < involved_vars.size(); v++) {
    if (involved_vars[v]->nLabels() != 2) {
      INTERNAL_ERROR << "instantiation of a cardinality factor with non-binary variables. Exiting." << std::endl;
      exit(1);
    }
  }
}

double CardinalityCumTRWSFactorBase::compute_reparameterization(const CumTRWSVar* var, const Math1D::Vector<float>& cost) {

  assert(var->nLabels() == 2);

  const uint nVars = involved_var_.size();

  assert(cost.size() >= nVars+1);

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->cost();
    else {
      idx = k;
    }
  }

  assert(idx < nVars);

  Math1D::Vector<double>& message = reparameterization_[idx];
  message.set_constant(1e300);

  Math1D::NamedVector<double> rel_param(nVars-1,MAKENAME(rel_param));

  double offs = 0.0;

  uint next = 0;
  for (uint k=0; k < nVars; k++) {

    if (k != idx) {
      rel_param[next] =  param[k][0] - param[k][1]; 
      offs += -param[k][0];

      next++;
    }
  }

  std::sort(rel_param.direct_access(), rel_param.direct_access() + nVars-1);

  double cum_sum = 0.0;

  for (uint c=0; c < nVars; c++) {

    const double hyp0 = cum_sum + cost[c];
    if (hyp0 < message[0])
      message[0] = hyp0;

    const double hyp1 = cum_sum + cost[c+1];
    if (hyp1 < message[1])
      message[1] = hyp1;

    if (c+1 < nVars) 
      cum_sum += rel_param[c];
  }

  for (uint l=0; l < 2; l++)
    message[l] += offs;

  const double total_offs = message.min();

  for (uint k=0; k < 2; k++)
    message[k] -= total_offs;

  return total_offs;
}

/********************/

CardinalityCumTRWSFactor::CardinalityCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, const Math1D::Vector<float>& cost) :
  CardinalityCumTRWSFactorBase(involved_vars), cost_(cost) {

  if (cost_.size() < involved_vars.size()+1) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
    exit(1);
  }
}

/*virtual*/ CardinalityCumTRWSFactor::~CardinalityCumTRWSFactor() {}
  
/*virtual*/ double CardinalityCumTRWSFactor::compute_reparameterization(const CumTRWSVar* var) {

  return CardinalityCumTRWSFactorBase::compute_reparameterization(var,cost_);
}

/********************/


CardinalityCumTRWSRefFactor::CardinalityCumTRWSRefFactor(const Storage1D<CumTRWSVar*>& involved_vars, const Math1D::Vector<float>& cost) :
  CardinalityCumTRWSFactorBase(involved_vars), cost_(cost) {
}

/*virtual*/ CardinalityCumTRWSRefFactor::~CardinalityCumTRWSRefFactor() {}
  
/*virtual*/ double CardinalityCumTRWSRefFactor::compute_reparameterization(const CumTRWSVar* var) {

  assert(cost_.size() >= involved_var_.size()+1);

  return CardinalityCumTRWSFactorBase::compute_reparameterization(var,cost_);
}

/********************/

CardinalityCumTRWSFactorBaseWithReuse::CardinalityCumTRWSFactorBaseWithReuse(const Storage1D<CumTRWSVar*>& involved_vars)
  : CumTRWSFactor(involved_vars), order_(involved_var_.size(),MAX_UINT), value_(involved_var_.size()) {


  for (uint v=0; v < involved_vars.size(); v++) {
    if (involved_vars[v]->nLabels() != 2) {
      INTERNAL_ERROR << "instantiation of a cardinality factor with non-binary variables. Exiting." << std::endl;
      exit(1);
    }
  }

  to_update_ = MAX_UINT;
}

/*virtual*/ CardinalityCumTRWSFactorBaseWithReuse::~CardinalityCumTRWSFactorBaseWithReuse() {}

/*virtual*/ 
void CardinalityCumTRWSFactorBaseWithReuse::init() {

  sort_by_rank();
  to_update_ = MAX_UINT;
}

double CardinalityCumTRWSFactorBaseWithReuse::compute_reparameterization(const CumTRWSVar* var, const Math1D::Vector<float>& cost) {

  assert(var->nLabels() == 2);

  const uint nVars = involved_var_.size();

  assert(cost.size() >= nVars+1);
  
  uint idx = MAX_UINT;

  if (to_update_ != MAX_UINT) {

    if (var == involved_var_[to_update_]) {
      idx = to_update_;
      TODO("should not happen if we check correctly for minimum (forward) and maximum (backward)");
    }
    else {
      if (to_update_+1 < nVars && involved_var_[to_update_+1] == var)
        idx = to_update_+1;
      else if (to_update_ > 0 && involved_var_[to_update_-1] == var)
        idx = to_update_-1;
      else {
        TODO("should not happen");
      }

      offs_ -= involved_var_[idx]->cost()[0] - reparameterization_[idx][0];
      const double update_val = involved_var_[to_update_]->cost()[0] - reparameterization_[to_update_][0];
      offs_ += update_val;


      const double new_val = (involved_var_[to_update_]->cost()[1] - reparameterization_[to_update_][1]) - update_val ;
      
      value_[order_[to_update_]].first = new_val;

      while (order_[to_update_] > 0 && value_[order_[to_update_]-1].first > new_val) {
        uint var = value_[order_[to_update_]-1].second;
        std::swap(value_[order_[to_update_]],value_[order_[to_update_]-1]);
        order_[var]++;
        order_[to_update_]--;
      }
     
      while (order_[to_update_]+1 < nVars && value_[order_[to_update_]+1].first < new_val) {
        uint var = value_[order_[to_update_]+1].second;
        std::swap(value_[order_[to_update_]],value_[order_[to_update_]+1]);
        order_[var]--;
        order_[to_update_]++;
      }
    }
  }
  else {

    Storage1D<Math1D::Vector<double> > param = reparameterization_;
    for (uint k=0; k < nVars; k++) {
      if (involved_var_[k] != var)
        param[k] -= involved_var_[k]->cost();
      else {
        idx = k;
      }
    }

    assert(idx < nVars);
    
    offs_ = 0.0;

    for (uint k=0; k < nVars; k++) {
      
      value_[k] =  std::make_pair(param[k][0] - param[k][1],k);
      if (k != idx) 
        offs_ += -param[k][0];
    }

    std::sort(value_.direct_access(),value_.direct_access()+nVars);
    for (uint k=0; k < nVars; k++) {
      order_[value_[k].second] = k;
    }
  }

  
  Math1D::Vector<double>& message = reparameterization_[idx];
  message[0] = cost[0];
  message[1] = cost[1];

  int true_k = -1;
  double cum_sum = 0.0;
  for (uint k=0; k < nVars; k++) {

    if (order_[idx] == k) 
      continue;

    true_k++;
    cum_sum += value_[k].first;

    const double hyp0 = cum_sum + cost[true_k+1];
    if (hyp0 < message[0])
      message[0] = hyp0;
    
    const double hyp1 = cum_sum + cost[true_k+2];
    if (hyp1 < message[1])
      message[1] = hyp1;
  }

  for (uint l=0; l < 2; l++)
    message[l] += offs_;

  to_update_ = idx;

  const double total_offs = message.min();


  for (uint k=0; k < 2; k++)
    message[k] -= total_offs;

  return total_offs;
}

/********************/

CardinalityCumTRWSFactorWithReuse::CardinalityCumTRWSFactorWithReuse(const Storage1D<CumTRWSVar*>& involved_vars,
                                                                     const Math1D::Vector<float>& cost) 
  : CardinalityCumTRWSFactorBaseWithReuse(involved_vars), cost_(cost) {}

/*virtual*/ CardinalityCumTRWSFactorWithReuse::~CardinalityCumTRWSFactorWithReuse() {}

/*virtual*/ 
double CardinalityCumTRWSFactorWithReuse::compute_reparameterization(const CumTRWSVar* var) {

  return CardinalityCumTRWSFactorBaseWithReuse::compute_reparameterization(var,cost_);
}

/********************/

AllPosBILPCumTRWSFactor::AllPosBILPCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, 
                                                 int rhs_lower, int rhs_upper)
  : CumTRWSFactor(involved_vars), rhs_lower_(std::max(0,rhs_lower)), rhs_upper_(std::min<int>(involved_vars.size(),rhs_upper)) {

  if (rhs_lower_ > rhs_upper_ || rhs_upper_ < 0) {
    INTERNAL_ERROR << "constraint is unsatisfiable, so inference is pointless. Exiting." << std::endl;
    exit(1);
  }

  for (uint v=0; v < involved_vars.size(); v++) {
    if (involved_vars[v]->nLabels() != 2) {
      INTERNAL_ERROR << "instantiation of an AllPosBILP factor with non-binary variables. Exiting." << std::endl;
      exit(1);
    }
  }
}

/*virtual*/ AllPosBILPCumTRWSFactor::~AllPosBILPCumTRWSFactor() {}

/*virtual*/ 
double AllPosBILPCumTRWSFactor::compute_reparameterization(const CumTRWSVar* var) {


  assert(var->nLabels() == 2);

  const uint nVars = involved_var_.size();

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->cost();
    else {
      idx = k;
    }
  }

  assert(idx < nVars);

  Math1D::Vector<double>& message = reparameterization_[idx];
  message.set_constant(1e300);

  Math1D::NamedVector<double> rel_param(nVars-1,MAKENAME(rel_param));

  double offs = 0.0;

  uint next = 0;
  for (uint k=0; k < nVars; k++) {

    if (k != idx) {
      rel_param[next] =  param[k][0] - param[k][1]; 
      offs += -param[k][0];

      next++;
    }
  }

  //TODO: think about partial_sort
  std::sort(rel_param.direct_access(), rel_param.direct_access() + nVars-1);

  double cum_sum = 0.0;

  for (int c=0; c <= rhs_upper_; c++) {

    if (c >= rhs_lower_) {
      if (cum_sum < message[0])
        message[0] = cum_sum;
    }

    if (c+1 >= rhs_lower_ && c+1 <= rhs_upper_) {
      if (cum_sum < message[1])
        message[1] = cum_sum;
    }

    if (c+1 < int(nVars))
      cum_sum += rel_param[c];
  }

  for (uint l=0; l < 2; l++)
    message[l] += offs;

  const double total_offs = message.min();

  for (uint k=0; k < 2; k++)
    message[k] -= total_offs;

  return total_offs;
}

/********************/

BILPCumTRWSFactor::BILPCumTRWSFactor(const Storage1D<CumTRWSVar*>& involved_vars, const Storage1D<bool>& positive,
                                     int rhs_lower, int rhs_upper) :
  CumTRWSFactor(involved_vars), rhs_lower_(rhs_lower), rhs_upper_(rhs_upper) {


  for (uint v=0; v < involved_vars.size(); v++) {
    if (involved_vars[v]->nLabels() != 2) {
      INTERNAL_ERROR << "instantiation of an BILP factor with non-binary variables. Exiting." << std::endl;
      exit(1);
    }
  }

  if (positive.size() < involved_vars.size()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
    exit(1);
  }

  if (rhs_lower_ > rhs_upper_) {
    INTERNAL_ERROR << "constraint is unsatisfiable, so inference is pointless. Exiting." << std::endl;
    exit(1);
  }

  assert(involved_vars.size() >= 2);


  Storage1D<CumTRWSVar*> sorted_involved_vars(involved_vars.size());
  uint next_pos = 0;

  //pass 1 - find all positive
  for (uint v=0; v < involved_vars.size(); v++) {
    if (positive[v]) {
      sorted_involved_vars[next_pos] = involved_vars[v];
      next_pos++;
    } 
  }

  nPos_ = next_pos;

  //pass 2 - find all negative
  for (uint v=0; v < involved_vars.size(); v++) {
    if (!positive[v]) {
      sorted_involved_vars[next_pos] = involved_vars[v];
      next_pos++;
    }
  }

  involved_var_ = sorted_involved_vars;

  const int nPositive = nPos_;
  const int nNegative = involved_var_.size()-nPositive;

  int lower_bound = -nNegative;
  
  /*** since we process the positive vars first, we need not compute entries below rhs_lower_-1 ***/
  /*** the offset of -1 is because we always leave one variable (the target of the message) out 
       in the forward computation, and this variable can have a positive sign ***/
  lower_bound = std::max(lower_bound, rhs_lower_ - 1);

  int upper_bound = nPositive;
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

  range_ = range;
  zero_offset_ = zero_offset;
}

/*virtual*/ BILPCumTRWSFactor::~BILPCumTRWSFactor() {}

/*virtual*/ double BILPCumTRWSFactor::compute_reparameterization(const CumTRWSVar* var) {

  const uint nVars = involved_var_.size();

  uint idx = MAX_UINT;

  Storage1D<Math1D::Vector<double> > param = reparameterization_;
  for (uint k=0; k < nVars; k++) {
    if (involved_var_[k] != var)
      param[k] -= involved_var_[k]->cost();
    else {
      idx = k;
      param[k].set_constant(0.0);
    }
  }

  Math1D::Vector<double>& message = reparameterization_[idx];

  assert(idx < nVars);

  /*** solution based on forward ****/


  /**** forward ****/
  assert(nVars >= 2);
  
  Math1D::Vector<double> forward_vec[2];
  forward_vec[0].resize(range_,1e100);
  forward_vec[1].resize(range_,1e100);
  
  Math1D::Vector<double>& start_forward_vec = forward_vec[0];
  
  const uint start_idx = (idx != 0) ? 0 : 1;
  
  const Math1D::Vector<double>& start_param = param[start_idx];
  
  start_forward_vec[zero_offset_] = -start_param[0];
  const int init_val = zero_offset_ + ((start_idx < nPos_) ? 1 : -1);
  if (init_val >= 0 && init_val < range_) {
    start_forward_vec[init_val] = -start_param[1];
  }
  
  uint cur_idx = 0;

  for (uint v = start_idx+1; v < nPos_; v++) {

    if (v != idx) {
      
      const Math1D::Vector<double>& cur_param = param[v];
      
      const Math1D::Vector<double>& prev_forward_vec = forward_vec[cur_idx];
      
      cur_idx = 1 - cur_idx;
      
      Math1D::Vector<double>& cur_forward_vec = forward_vec[cur_idx];
      
      for (int sum=zero_offset_; sum < std::min<int>(range_,zero_offset_+v+2); sum++) {
        
	double best = 1e300;

	for (int l=0; l < 2; l++) {
	  
	  double best_prev = 1e75;
          
            const int dest = sum - l;
            if (dest >= 0) {
	      
              best_prev = prev_forward_vec[dest]; 
            }
            
            const double hyp = best_prev - cur_param[l];
            if (hyp < best)
              best = hyp;
	}
	cur_forward_vec[sum] = best;
      }
    }
  }
  
  for (uint v = std::max(start_idx+1,nPos_); v < nVars; v++) {

    if (v != idx) {
      
      const Math1D::Vector<double>& cur_param = param[v];
	
      const Math1D::Vector<double>& prev_forward_vec = forward_vec[cur_idx];
      
      cur_idx = 1 - cur_idx;
      
      Math1D::Vector<double>& cur_forward_vec = forward_vec[cur_idx];
      
      for (int sum=0; sum < range_; sum++) {
	
	double best = 1e300;

	for (int l=0; l < 2; l++) {
	  
	  double best_prev = 1e75;
          
	  const int dest = sum + l;
	  if (dest < range_) {
	    
	    best_prev = prev_forward_vec[dest]; 
	  }
          
	  const double hyp = best_prev - cur_param[l];
	  if (hyp < best)
	    best = hyp;
	}
	cur_forward_vec[sum] = best;
      }
    }
  }

  const Math1D::Vector<double>& last_forward_vec = forward_vec[cur_idx];

  //now derive the message
  for (uint l=0; l < 2; l++) {
      
    double min_msg = 1e300;

    for (int s = rhs_lower_ + zero_offset_; s <= rhs_upper_ + zero_offset_; s++) { 
      
      double best_prev = 1e75;
      
      int move = l;
      if (idx < nPos_) 
	move *= -1; //since we are tracing backward here
      
      const int dest = s + move;
      if (dest >= 0 && dest < range_) {
	
	best_prev = last_forward_vec[dest]; 
      }
      
      double hyp = best_prev - param[idx][l];
      
      if (hyp < min_msg)
	min_msg = hyp;
    }      
    message[l] = min_msg;
  }

  const double offs = message.min();

  if (fabs(offs) > 1e10) {

    std::cerr << "idx: " << idx << "/" << nVars << std::endl;
    std::cerr << "message: " << message << std::endl;
    for (uint v=0; v < nVars; v++) {
      std::cerr << "combined parameters[" << v <<  "]: " << param[v] << std::endl;
      std::cerr << "reparameterization[" << v << "]: " << reparameterization_[v] << std::endl;
    }

    std::cerr << "factor: ";
    for (uint v=0; v < involved_var_.size(); v++)
      std::cerr << involved_var_[v]->rank() << " ";
    std::cerr << std::endl;
    std::cerr << "var params:" << std::endl;
    for (uint v=0; v < involved_var_.size(); v++)
      std::cerr << involved_var_[v]->cost() << std::endl;

    //DEBUG
    exit(1);
    //END_DEBUG
  }

  for (uint k=0; k < 2; k++)
    message[k] -= offs;

  return offs;
}

/********************/


BILPCumTRWSFactorWithReuse::BILPCumTRWSFactorWithReuse(const Storage1D<CumTRWSVar*>& involved_vars, const Storage1D<bool>& positive,
                                                       int rhs_lower, int rhs_upper) :
  CumTRWSFactor(involved_vars), positive_(positive), rhs_lower_(rhs_lower), rhs_upper_(rhs_upper) {

  const uint nVars = involved_var_.size();

  for (uint v=0; v < nVars; v++)
    assert(involved_vars[v]->nLabels() == 2);

  assert(rhs_lower_ <= rhs_upper_);

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

  zero_offset_ = zero_offset;

  forward_light_.resize(range,nVars-1);
  backward_light_.resize(range,nVars-1);

  to_update_ = MAX_UINT;
}

/*virtual*/ BILPCumTRWSFactorWithReuse::~BILPCumTRWSFactorWithReuse() {}

/*virtual*/ 
void BILPCumTRWSFactorWithReuse::init() {
  sort_by_rank();
}


/*virtual*/ 
double BILPCumTRWSFactorWithReuse::compute_reparameterization(const CumTRWSVar* var) {

  const uint nVars = involved_var_.size();

  uint idx = MAX_UINT;

  const int range = forward_light_.xDim();

  if (to_update_ != MAX_UINT) {
    Math1D::Vector<double> prev_param = reparameterization_[to_update_];
    prev_param -= involved_var_[to_update_]->cost();


    if (to_update_+1 < nVars && involved_var_[to_update_+1] == var) {
      idx = to_update_+1;

      //need to update forward

      if (to_update_ == 0) {

        forward_light_(zero_offset_,0) = -prev_param[0];
        const int init_mul = (positive_[0]) ? 1 : -1;
        if (int(zero_offset_)+init_mul >= 0
            && int(zero_offset_)+init_mul < range) {
          forward_light_(zero_offset_+init_mul,0) = -prev_param[1];
        }        
      }
      else {
                
        for (int sum=0; sum < range; sum++) {

          double best = 1e300;
        
          for (int l=0; l < 2; l++) {
            
            double best_prev = 1e75;
            
            int move = l;
            if (positive_[to_update_]) //since we are tracing backward here
              move *= -1;
            
            const int dest = sum + move;
            if (dest >= 0 && dest < range) {
              
              best_prev = forward_light_(dest,to_update_-1);
            }
            
            double forward = best_prev - prev_param[l];
            if (forward < best)
              best = forward;
          }

          forward_light_(sum,to_update_) = best;
        }
      }
    }
    else if (to_update_ > 0 && involved_var_[to_update_-1] == var) {
      idx = to_update_-1;

      //need to update backward_light
      if (to_update_ == nVars-1) {

        backward_light_(zero_offset_,to_update_-1) = -prev_param[0];
        const int end_mul = (positive_[to_update_]) ? 1 : -1;
        if (int(zero_offset_) + end_mul >= 0
            && int(zero_offset_) + end_mul < range)
          backward_light_(zero_offset_ + end_mul,to_update_-1) = -prev_param[1];
      }
      else {

        for (int sum=0; sum < range; sum++) {
          
          double best_prev = 1e75;
          
          for (int l=0; l < 2; l++) {
            
            int move = l;
            if (positive_[to_update_]) //since we are tracing backward here
              move *= -1;
            
            const int dest = sum + move;
            double hyp = 1e75;
            if (dest >= 0 && dest < range) {
              hyp = backward_light_(dest,to_update_) - prev_param[l];
            }
            
            if (hyp < best_prev)
              best_prev = hyp;
          }
          
          backward_light_(sum,to_update_-1) = best_prev;
        }
      }
    }
    else {
      TODO("should not happen when vars are properly ordered");
    }
  }
  else {

    Storage1D<Math1D::Vector<double> > param = reparameterization_;
    for (uint k=0; k < nVars; k++) {
      if (involved_var_[k] != var)
        param[k] -= involved_var_[k]->cost();
      else {
        idx = k;
        param[k].set_constant(0.0);
      }
    }
    
    assert(idx < nVars);

    //init
    for (int sum=0; sum < range; sum++) {
      
      forward_light_(sum,0) = 1e100;
    }

    forward_light_(zero_offset_,0) = -param[0][0];
    const int init_mul = (positive_[0]) ? 1 : -1;
    if (int(zero_offset_)+init_mul >= 0
        && int(zero_offset_)+init_mul < range) {
      forward_light_(zero_offset_+init_mul,0) = -param[0][1];
    }
    
    //proceed
    for (uint v=1; v <= idx; v++) {
    
      for (int sum=0; sum < range; sum++) {
        
        double best = 1e300;
        
        for (int l=0; l < 2; l++) {
          
          double best_prev = 1e75;
          
          int move = l;
          if (positive_[v]) //since we are tracing backward here
            move *= -1;

          const int dest = sum + move;
          if (dest >= 0 && dest < range) {
            
            best_prev = forward_light_(dest,v-1);
          }
          
          double forward = best_prev - param[v][l];

          if (forward < best)
            best = forward;
        }
        forward_light_(sum,v) = best;
      }
    }
    
    const uint last_var = nVars-1;
    
    //init
    for (int sum=0; sum < range; sum++) 
      backward_light_(sum,last_var-1) = 1e100;
    
    backward_light_(zero_offset_,last_var-1) = -param[last_var][0];
    const int end_mul = (positive_[last_var]) ? 1 : -1;
    if (int(zero_offset_) + end_mul >= 0
        && int(zero_offset_) + end_mul < range)
      backward_light_(zero_offset_ + end_mul,last_var-1) = -param[last_var][1];

    //proceed
    for (int v=last_var-1; v > int(idx); v--) {
      
      for (int sum=0; sum < range; sum++) {
      
        double best_prev = 1e75;
        
        for (int l=0; l < 2; l++) {
          
          int move = l;
          if (positive_[v]) //since we are tracing backward here
            move *= -1;

          const int dest = sum + move;
          double hyp = 1e75;
          if (dest >= 0 && dest < range) {
            hyp = backward_light_(dest,v) - param[v][l];
          }
          
          if (hyp < best_prev)
            best_prev = hyp;
        }
        
        backward_light_(sum,v-1) = best_prev;
      }
    }
  }

  Math1D::Vector<double>& message = reparameterization_[idx];

  for (uint l=0; l < 2; l++) {
    
    double min_msg = 1e300;
    
    for (int s=0; s < (int) range; s++) {
      
      double forward = 1e300;
      if (idx == 0) {

        if (l == 0 && s == zero_offset_)
          forward = 0.0;
        if (l == 1) {
          const int init_mul = (positive_[0]) ? 1 : -1;
          if (s == int(zero_offset_)+init_mul)
            forward = 0.0;
        }
      }
      else {
        
        int move = l;
        if (positive_[idx]) //since we are tracing backward here
          move *= -1;
        
        const int dest = s + move;
        if (dest >= 0 && dest < range) {
          forward = forward_light_(dest,idx-1);
        }
      }

      double hyp = forward; 
      
      if (idx+1 < nVars) {
        
        double best_bwd = 1e300;
        
        const int diff = (s - zero_offset_);
        
        for (int r=rhs_lower_; r <= rhs_upper_; r++) {
          const int other = r + zero_offset_ - diff; 
          
          if (other >= 0 && other < (int) range) {
            
            best_bwd = std::min(best_bwd,backward_light_(other,idx));
          }
        }
        
        hyp += best_bwd;
      }
      else {
        if (s < int(rhs_lower_ + zero_offset_) || s > int(rhs_upper_ + zero_offset_)) 
          hyp = 1e300;
      }
      
      if (hyp < min_msg)
        min_msg = hyp;
    }
    
    message[l] = min_msg;
  }

  double offs = message.min();

  for (uint k=0; k < 2; k++)
    message[k] -= offs;

  to_update_ = idx;

  return offs;
}


/********************/

CumFactorTRWS::CumFactorTRWS(uint nVars, uint nFactors) : var_(nVars), factor_(nFactors), labeling_(nVars,0), rank2var_(nVars,MAX_UINT),
                                                          nUsedVars_(0), nUsedFactors_(0), optimize_called_(false) {}

CumFactorTRWS::~CumFactorTRWS() {

  for (uint v=0; v < nUsedVars_; v++)
    delete var_[v];

  for (uint f=0; f < nUsedFactors_; f++)
    delete factor_[f];
}

void CumFactorTRWS::set_ranks(const Math1D::Vector<uint>& ranks) {

  if (ranks.size() != nUsedVars_) {
    std::cerr << "setting ranks is currently only possible for all variables at once." << std::endl;
    return;
  }

  Storage1D<bool> hit(nUsedVars_,false);

  for (uint i=0; i < nUsedVars_; i++) {
    const uint rank = ranks[i];
    if (rank >= nUsedVars_) {
      INTERNAL_ERROR << "rank #" << rank << " is out of range" << std::endl;
      exit(1);      
    }
    if (hit[rank]) {
      INTERNAL_ERROR << "rank #" << rank << " was doubly specified" << std::endl;
      exit(1);
    }
    hit[rank] = true;
    var_[i]->set_rank(ranks[i]);
    rank2var_[ranks[i]] = i;
  }
}

uint CumFactorTRWS::add_var(const Math1D::Vector<float>& cost) {

  assert(!optimize_called_);

  if (nUsedVars_ == var_.size())
    var_.resize(uint(nUsedVars_*1.2)+4);

  assert(nUsedVars_ < var_.size());
  var_[nUsedVars_] = new CumTRWSVar(cost,nUsedVars_);
  rank2var_[nUsedVars_] = nUsedVars_;

  nUsedVars_++;

  return nUsedVars_-1;
}

CumTRWSVar* CumFactorTRWS::get_variable(uint v) {

  if (v >= nUsedVars_) {
    INTERNAL_ERROR << "variable index out of bounds. Exiting. " << std::endl;
    exit(1);
  }

  return var_[v];
}

CumTRWSFactor* CumFactorTRWS::get_factor(uint f) {

  if (f >= nUsedFactors_) {
    INTERNAL_ERROR << "factor index out of bounds. Exiting. " << std::endl;
    exit(1);
  }

  return factor_[f];
}

//CAUTION: the passed factor will be deleted together with all truly owned factors
uint CumFactorTRWS::pass_in_factor(CumTRWSFactor* new_fac) {

  return add_factor(new_fac);
}

uint CumFactorTRWS::add_factor(CumTRWSFactor* fac) {

  assert(!optimize_called_);

  if (factor_.size() == nUsedFactors_)
    factor_.resize(uint(nUsedFactors_*1.2)+4);

  factor_[nUsedFactors_] = fac;
  nUsedFactors_++;

  return nUsedFactors_-1;
}
  
uint CumFactorTRWS::add_generic_factor(const Math1D::Vector<uint>& var, const VarDimStorage<float>& cost) {

  Storage1D<CumTRWSVar*> vars(var.size());

  for (uint k=0; k < var.size(); k++) {
    if (var[k] >= nUsedVars_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    vars[k] = var_[var[k]];
  }

  return add_factor(new GenericCumTRWSFactor(vars,cost));
}

uint CumFactorTRWS::add_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost, bool ref) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }


  Storage1D<CumTRWSVar*> vars(2);
  vars[0] = var_[var1];
  vars[1] = var_[var2];

  CumTRWSFactor* newFac;
  if (ref)
    newFac = new BinaryCumTRWSRefFactor(vars,cost);
  else
    newFac = new BinaryCumTRWSFactor(vars,cost);
      
  return add_factor(newFac);
}

uint CumFactorTRWS::add_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_ || var3 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  Storage1D<CumTRWSVar*> vars(3);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];

  CumTRWSFactor* new_fac = 0;

  if (!ref)
    new_fac = new TernaryCumTRWSFactor(vars,cost);
  else
    new_fac = new TernaryCumTRWSRefFactor(vars,cost);

  return add_factor(new_fac);
}

uint CumFactorTRWS::add_second_diff_factor(uint var1, uint var2, uint var3, float lambda) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_ || var3 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  Storage1D<CumTRWSVar*> vars(3);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];

  return add_factor(new SecondDiffCumTRWSFactor(vars,lambda));
}

uint CumFactorTRWS::add_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                                            const Storage1D<Math3D::Tensor<float> >& cost, bool ref) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_ || var3 >= nUsedVars_ || var4 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  Storage1D<CumTRWSVar*> vars(4);
  vars[0] = var_[var1];
  vars[1] = var_[var2];
  vars[2] = var_[var3];
  vars[3] = var_[var4];

  CumTRWSFactor* newFac;

  if (ref)
    newFac = new FourthOrderCumTRWSRefFactor(vars,cost);
  else
    newFac = new FourthOrderCumTRWSFactor(vars,cost);

  return add_factor(newFac);
}

uint CumFactorTRWS::add_one_of_n_factor(const Math1D::Vector<uint>& var, bool reuse) {

  Storage1D<CumTRWSVar*> vars(var.size());
  
  for (uint k=0; k < var.size(); k++) {
   
    if (var[k] >= nUsedVars_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    vars[k] = var_[var[k]];

    if (vars[k]->nLabels() != 2) {
      INTERNAL_ERROR << " variables of 1-of-N nodes must be binary. Exiting..." << std::endl;
      exit(1);
    }
  }


  if (var.size() == 0) {
    std::cerr << "WARNING: ignoring empty 1-of-N factor" << std::endl;

    return MAX_UINT;
  }
  else if (var.size() == 1) {

    Math1D::Vector<float> cost(2);
    cost[0] = 10000.0;
    cost[1] = 0.0;

    vars[0]->add_cost(cost);

    return MAX_UINT;
  }
  else {
    if (!reuse)
      return add_factor(new OneOfNCumTRWSFactor(vars));
    else
      return add_factor(new OneOfNCumTRWSFactorWithReuse(vars));
  }
}

uint CumFactorTRWS::add_cardinality_factor(const Math1D::Vector<uint>& var, const Math1D::Vector<float>& cost, bool ref, bool reuse) {

  if (var.size() == 0) {
    std::cerr << "WARNING: ignoring empty cardinality factor" << std::endl;

    return MAX_UINT;
  }
  else if (var.size() == 1) {
    if (var[0] >= nUsedVars_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }
    if (var_[var[0]]->nLabels() != 2) {
      INTERNAL_ERROR << " variables of cardinality nodes must be binary. Exiting..." << std::endl;
      exit(1);
    }

    var_[var[0]]->add_cost(cost);

    return MAX_UINT;
  }
  else {

    Storage1D<CumTRWSVar*> vars(var.size());
  
    for (uint k=0; k < var.size(); k++) {
      if (var[k] >= nUsedVars_) {
        INTERNAL_ERROR << "out of range. Exiting." << std::endl;
        exit(1);
      }

      vars[k] = var_[var[k]];

      if (vars[k]->nLabels() != 2) {
        INTERNAL_ERROR << " variables of cardinality nodes must be binary. Exiting..." << std::endl;
        exit(1);
      }
    }

    if (!reuse) {
      if (!ref)
        return add_factor(new CardinalityCumTRWSFactor(vars,cost));
      else
        return add_factor(new CardinalityCumTRWSRefFactor(vars,cost));
    }
    else
      return add_factor(new CardinalityCumTRWSFactorWithReuse(vars,cost));
  }
}

uint CumFactorTRWS::add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive,
                                          int rhs_lower, int rhs_upper, bool reuse) {

  uint nUseful = 0;
  int nPos = 0;
  int nNeg = 0;
  for (uint k=0; k < var.size(); k++) {
    
    if (var[k] >= nUsedVars_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    if (var_[var[k]]->nLabels() != 2) {
      INTERNAL_ERROR << " variables of BILP nodes must be binary. Exiting..." << std::endl;
      exit(1);
    }

    const Math1D::Vector<double>& cur_cost = var_[var[k]]->cost();

    if (fabs(cur_cost[0] - cur_cost[1]) < 1e10) {
      nUseful++;
      if (positive[k])
        nPos++;
      else
        nNeg++;
    }
    else {
      if (cur_cost[0] > cur_cost[1]) {
        //var is fixed to 1 => adjust rhs_lower and rhs_upper
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

  if (nUseful != 0) {

    assert(nUseful >= 2);
    
    Storage1D<CumTRWSVar*> vars(nUseful);
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

    if (nNeg == 0 && !reuse) 
      return add_factor(new AllPosBILPCumTRWSFactor(vars,rhs_lower,rhs_upper));
    else {
      if (reuse)
        return add_factor(new BILPCumTRWSFactorWithReuse(vars,reduced_positive,rhs_lower,rhs_upper));
      else
        return add_factor(new BILPCumTRWSFactor(vars,reduced_positive,rhs_lower,rhs_upper));
    }
  }
  else {
    std::cerr << "WARNING: removed superfluous constraint factor" << std::endl;

    return MAX_UINT;
  }
}

uint CumFactorTRWS::best_of_n(uint fac_num) const {

  if (fac_num >= nUsedFactors_) {

    INTERNAL_ERROR << "factor index out of bounds. Exiting..." << std::endl;
    exit(1);
  }

  return factor_[fac_num]->best_of_n();
}

double CumFactorTRWS::optimize(uint nIter, bool quiet) {

  if (!quiet)
    std::cerr << "efficient Cum-scheme" << std::endl;

  if (!optimize_called_) {
    size_t max_order = 0;

    for (uint f=0; f < nUsedFactors_; f++) {
      factor_[f]->compute_rank_range();
      max_order = std::max(max_order,factor_[f]->involved_vars().size());
    }
    
    if (!quiet)
      std::cerr << "maximal order of all factors: " << max_order << std::endl;
    
    for (uint v=0; v < nUsedVars_; v++) {
      var_[v]->set_up_chains();
    }
  }
  optimize_called_ = true;

  uint arg_min;

  double bound = 1e300;

  if (!quiet) {
    size_t message_effort = 0;
    for (uint f=0; f < nUsedFactors_; f++) {
      message_effort += 2* (factor_[f]->involved_vars().size()-1) * (factor_[f]->involved_vars().size());
    } 
    message_effort *= nIter;
    std::cerr << "predicted message effort: " << message_effort << std::endl;
  }

  for (uint f=0; f < nUsedFactors_; f++)
    factor_[f]->init();

  for (uint iter=1; iter <= nIter; iter++) {

    if (!quiet)
      std::cerr << "******* iteration " << iter << " *********" << std::endl;

    //forward
    double forward_lower = 0.0;

    for (uint i=0; i < nUsedVars_; i++) {

      if (rank2var_[i] == MAX_UINT)
        continue; //can happen when variables do not have adjacent factors
      
      CumTRWSVar* cur_var = var_[rank2var_[i]];
      
      Storage1D<CumTRWSFactor*> adj_factors = cur_var->adjacent_factor();
      
      for (uint k=0; k < adj_factors.size(); k++) {
	
        if (i != adj_factors[k]->min_rank()) {

          double cur_offs = adj_factors[k]->compute_reparameterization(cur_var);
	
	  if (i == adj_factors[k]->max_rank())
	    forward_lower += cur_offs;
        }
      }

      forward_lower += cur_var->average(arg_min);
      labeling_[i] = arg_min;
    }

    //backward
    double backward_lower = 0.0;
    for (int i = nUsedVars_-1; i >= 0; i--) {

      if (rank2var_[i] == MAX_UINT)
        continue; //can happen when variables do not have adjacent factors
      
      CumTRWSVar* cur_var = var_[rank2var_[i]];

      Storage1D<CumTRWSFactor*> adj_factors = cur_var->adjacent_factor();
      
      for (uint k=0; k < adj_factors.size(); k++) {

        if (i != int(adj_factors[k]->max_rank())) {	
          double cur_offs = adj_factors[k]->compute_reparameterization(cur_var);
	
          if (i == int(adj_factors[k]->min_rank()))
            backward_lower += cur_offs;
        }
      }

      backward_lower += cur_var->average(arg_min);
      labeling_[i] = arg_min;
    }

    if (!quiet)
      std::cerr << "bwd bound: " << backward_lower << std::endl;

    bound = backward_lower;
  }

  if (!quiet) {

    size_t message_effort = 0;
    for (uint f=0; f < nUsedFactors_; f++) {
      message_effort += 2* (factor_[f]->involved_vars().size()-1) * (factor_[f]->involved_vars().size());
    } 
    message_effort *= nIter;
    std::cerr << "message effort: " << message_effort << std::endl;
  }

  return bound;
}

const Math1D::Vector<uint>& CumFactorTRWS::labeling() {
  return labeling_;
}
