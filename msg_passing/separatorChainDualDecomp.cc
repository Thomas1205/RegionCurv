/******* written by Thomas Schoenemann as an employee of the University of Düsseldorf, Germany, 2012 *****/

#include "separatorChainDualDecomp.hh"

#include <set>
#include <map>

SepChainDDVar::SepChainDDVar(const Math1D::Vector<float>& cost) :
  cost_(cost) {}

void SepChainDDVar::add_cost(const Math1D::Vector<float>& cost) {
  cost_ += cost;
}
  
void SepChainDDVar::add_factor(SepChainDDFactor* factor) {

  uint size = neighboring_factor_.size();

  neighboring_factor_.resize(size+1);
  neighboring_factor_[size] = factor;
}
  
void SepChainDDVar::add_separator(SepChainDDPairSeparator* sep) {

  uint size = neighboring_separator_.size();
  
  neighboring_separator_.resize(size+1);
  neighboring_separator_[size] = sep;
}

const Math1D::Vector<float>& SepChainDDVar::cost() const {
  return cost_;
}
  
uint SepChainDDVar::nLabels() const {
  return cost_.size();
}

const Storage1D<SepChainDDFactor*>& SepChainDDVar::neighboring_factor() const {
  return neighboring_factor_;
}

uint SepChainDDVar::nChains() const {
  TODO("");
}

void SepChainDDVar::set_up_chains() {

  uint nChains = 0;

  for (uint f=0; f < neighboring_factor_.size(); f++) {
    
    SepChainDDPairSeparator* prev_sep = neighboring_factor_[f]->prev_sep();

    if (neighboring_factor_[f]->prev_var() != this
        && (prev_sep == 0 || (prev_sep->var1() != this && prev_sep->var2() != this) ) ) {

      nChains++;
    }
  }

  cost_ *= 1.0 / nChains;
}

/***********************************/

SepChainDDPairSeparator::SepChainDDPairSeparator(SepChainDDVar* var1, SepChainDDVar* var2) :
  var1_(var1), var2_(var2) {
  var1->add_separator(this);
  var2->add_separator(this);
}

SepChainDDVar* SepChainDDPairSeparator::var1() {
  return var1_;
}

SepChainDDVar* SepChainDDPairSeparator::var2() {
  return var2_;
}

const SepChainDDVar* SepChainDDPairSeparator::var1() const {
  return var1_;
}

const SepChainDDVar* SepChainDDPairSeparator::var2() const {
  return var2_;
}

void SepChainDDPairSeparator::add_factor(SepChainDDFactor* factor) {

  uint size = neighboring_factor_.size();
  neighboring_factor_.resize(size+1);
  neighboring_factor_[size] = factor;
}

Storage1D<SepChainDDFactor*> SepChainDDPairSeparator::neighboring_factor() const {
  return neighboring_factor_;
}

/***********************************/

SepChainDDFactor::SepChainDDFactor(const Storage1D<SepChainDDVar*>& involved_vars, 
                                   const Storage1D<SepChainDDPairSeparator*>& separators) :
  prev_var_(0), next_var_(0), prev_sep_(0), next_sep_(0), prev_factor_(0), next_factor_(0), 
  involved_var_(involved_vars), involved_separator_(separators) {

  dual_var_.resize(involved_vars.size());

  for (uint v=0; v < involved_vars.size(); v++) {
    involved_vars[v]->add_factor(this);
    dual_var_[v].resize(involved_vars[v]->nLabels(),0.0);
  }

  dual_pair_var_.resize(involved_vars.size());

  for (uint s=0; s < separators.size(); s++) {
    separators[s]->add_factor(this);

    dual_pair_var_[s].resize(separators[s]->var1()->nLabels(),separators[s]->var2()->nLabels(),0.0);
  }
}

/*virtual*/ SepChainDDFactor::~SepChainDDFactor() {}

SepChainDDVar* SepChainDDFactor::prev_var() const {
  return prev_var_;
}

SepChainDDVar* SepChainDDFactor::next_var() const {
  return next_var_;
}

SepChainDDPairSeparator* SepChainDDFactor::prev_sep() const {
  return prev_sep_;
}

SepChainDDPairSeparator* SepChainDDFactor::next_sep() const {
  return next_sep_;
}

void SepChainDDFactor::set_prev_var(SepChainDDVar* var) {
  prev_var_ = var;
}

void SepChainDDFactor::set_next_var(SepChainDDVar* var) {
  next_var_ = var;
}

void SepChainDDFactor::set_prev_sep(SepChainDDPairSeparator* sep) {
  prev_sep_ = sep;
}

void SepChainDDFactor::set_next_sep(SepChainDDPairSeparator* sep) {
  next_sep_ = sep;
}

Math1D::Vector<double>& SepChainDDFactor::get_duals(const SepChainDDVar* var) {

  for (uint k=0; k < involved_var_.size(); k++) {

    if (involved_var_[k] == var)
      return dual_var_[k];
  }

  assert(false);
  return dual_var_[0];
}

Math2D::Matrix<double>& SepChainDDFactor::get_pair_duals(const SepChainDDPairSeparator* sep) {

  for (uint s=0; s < involved_separator_.size(); s++) {
    
    if (involved_separator_[s] == sep)
      return dual_pair_var_[s];
  }

  assert(false);
  return dual_pair_var_[0];
}

SepChainDDFactor* SepChainDDFactor::prev_factor() const {
  return prev_factor_;
}

SepChainDDFactor* SepChainDDFactor::next_factor() const {
  return next_factor_;
}

void SepChainDDFactor::set_prev_factor(SepChainDDFactor* factor) {
  prev_factor_ = factor;
}

void SepChainDDFactor::set_next_factor(SepChainDDFactor* factor) {
  next_factor_ = factor;
}

const Storage1D<SepChainDDVar*>& SepChainDDFactor::involved_vars() const {
  return involved_var_;
}

const Storage1D<SepChainDDPairSeparator*>& SepChainDDFactor::involved_separators() const {
  return involved_separator_;
}

/***********************************/

BinarySepChainDDFactor::BinarySepChainDDFactor(const Storage1D<SepChainDDVar*>& involved_vars,
                                               const Math2D::Matrix<float>& cost) :
  SepChainDDFactor(involved_vars,Storage1D<SepChainDDPairSeparator*>()), cost_(cost) {

  if (involved_vars.size() != 2) {
    INTERNAL_ERROR << " attempt to instantiate binary factor with " << involved_vars.size() 
                   << " variables. Exiting." << std::endl;
    exit(1);
  }

  if (cost.xDim() < involved_vars[0]->nLabels() || cost.yDim() < involved_vars[1]->nLabels()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
    exit(1);
  }
}

/*virtual*/ BinarySepChainDDFactor::~BinarySepChainDDFactor() {}

/*virtual */
double BinarySepChainDDFactor::compute_forward(const SepChainDDPairSeparator* incoming_sep, const SepChainDDVar* incoming_var, 
                                               const SepChainDDVar* outgoing_var,
                                               const Math2D::Matrix<double>& /*prev_pair_forward*/, 
					       const Math1D::Vector<double>& prev_var_forward, 
                                               Math1D::Vector<double>& forward, 
                                               Math2D::Matrix<uint>& trace) const {

  assert(incoming_var != outgoing_var);

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();

  uint var_idx = 0;

  Storage1D<Math1D::Vector<double> > param = dual_var_;
  for (uint k=0; k < involved_var_.size(); k++) {
    if (involved_var_[k] == incoming_var) {
      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= prev_var_forward[l];
      }
    }
    else {

      if (incoming_sep != 0 && (incoming_sep->var1() == involved_var_[k] || incoming_sep->var2() == involved_var_[k] ) ) {
        //in this case we only want the dual vars of the factor

        assert(involved_var_[k] != outgoing_var);
      }
      else {

        if (involved_var_[k] == outgoing_var)
          var_idx = k;

        for (uint l=0; l < param[k].size(); l++) {
          param[k][l] -= involved_var_[k]->cost()[l];
        }
      }
    }
  }

  assert(var_idx != MAX_UINT);
  
  if (var_idx == 0) {
    
    forward.resize_dirty(nLabels1);
    trace.resize_dirty(nLabels1,2);

    for (uint l1=0; l1 < nLabels1; l1++) {

      double best = 1e300;
      uint best2 = 0;

      for (uint l2=0; l2 < nLabels2; l2++) {
      
        double hyp = cost_(l1,l2) - param[1][l2];

        if (hyp < best) {
          best = hyp;
          best2 = l2;
        }
      }

      forward[l1] = best - param[0][l1];
      trace(l1,0) = l1;
      trace(l1,1) = best2;
    }
  }
  else  {

    assert(var_idx == 1);

    forward.resize_dirty(nLabels2);
    trace.resize_dirty(nLabels2,2);

    for (uint l2=0; l2 < nLabels2; l2++) {

      double best = 1e300;
      uint best1 = 0;

      for (uint l1=0; l1 < nLabels1; l1++) {

        double hyp = cost_(l1,l2) - param[0][l1];
        
        if (hyp < best) {
          best = hyp;
          best1 = l1;
        }
      }

      forward[l2] = best - param[1][l2];
      trace(l2,0) = best1;
      trace(l2,1) = l2;
    }
  }

  return 0.0;
}

/*virtual*/
double BinarySepChainDDFactor::compute_forward(const SepChainDDPairSeparator* /*incoming_sep*/, const SepChainDDVar* /*incoming_var*/, 
                                               const SepChainDDPairSeparator* /*outgoing_sep*/,
                                               const Math2D::Matrix<double>& /*prev_pair_forward*/, 
					       const Math1D::Vector<double>& /*prev_var_forward*/, 
                                               Math2D::Matrix<double>& /*forward*/, 
                                               Math3D::Tensor<uint>& /*trace*/) const {

  //this routine should never be called
  assert(false);
  return 0.0;
}


/***********************************/

TernarySepChainDDFactor::TernarySepChainDDFactor(const Storage1D<SepChainDDVar*>& involved_vars, 
                                                 Storage1D<SepChainDDPairSeparator*>& separators,
                                                 const Math3D::Tensor<float>& cost) :
  SepChainDDFactor(involved_vars,separators), cost_(cost) {

  if (involved_vars.size() != 3) {
    INTERNAL_ERROR << " attempt to instantiate ternary factor with " 
                   << involved_vars.size() << " variables. Exiting." << std::endl;
    exit(1);
  }

  if (cost.xDim() < involved_vars[0]->nLabels() || cost.yDim() < involved_vars[1]->nLabels() 
      || cost.zDim() < involved_vars[2]->nLabels()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
    exit(1);
  }
}


/*virtual*/ TernarySepChainDDFactor::~TernarySepChainDDFactor() {}

/*virtual*/ 
double TernarySepChainDDFactor::compute_forward(const SepChainDDPairSeparator* incoming_sep, const SepChainDDVar* incoming_var, 
                                                const SepChainDDVar* outgoing_var,
                                                const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                                Math1D::Vector<double>& forward, 
                                                Math2D::Matrix<uint>& trace) const {

  assert(incoming_var != outgoing_var);

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();

  const uint nSeps = involved_separator_.size();

  uint var_idx = 0;

  Storage1D<Math1D::Vector<double> > param = dual_var_;
  for (uint k=0; k < involved_var_.size(); k++) {
    if (involved_var_[k] == incoming_var) {
      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= prev_var_forward[l];
      }
    }
    else {

      if (incoming_sep != 0 && (incoming_sep->var1() == involved_var_[k] || incoming_sep->var2() == involved_var_[k] ) ) {
        //in this case we only want the dual vars of the factor

        assert(involved_var_[k] != outgoing_var);
      }
      else {

        if (involved_var_[k] == outgoing_var)
          var_idx = k;

        for (uint l=0; l < param[k].size(); l++) {
          param[k][l] -= involved_var_[k]->cost()[l];
        }
      }
    }
  }

  Storage1D<Math2D::Matrix<double> > pair_param = dual_pair_var_;
  for (uint s=0; s < nSeps; s++) {

    if (involved_separator_[s] == incoming_sep) {
      pair_param[s] -= prev_pair_forward; 
    }
  }

  assert(var_idx != MAX_UINT);
  
  if (var_idx == 0) {
    
    forward.resize_dirty(nLabels1);
    trace.resize_dirty(nLabels1,3);

    for (uint l1=0; l1 < nLabels1; l1++) {

      double best = 1e300;
      uint best2 = 0;
      uint best3 = 0;

      for (uint l2=0; l2 < nLabels2; l2++) {

        const double w2 = param[1][l2];

        for (uint l3=0; l3 < nLabels3; l3++) {
      
          double hyp = cost_(l1,l2,l3) - w2 - param[2][l3];

          for (uint s=0; s < nSeps; s++)
            hyp -= eval_pair(s,l1,l2,l3,pair_param);
          
          if (hyp < best) {
            best = hyp;
            best2 = l2;
            best3 = l3;
          }
        }
      }

      forward[l1] = best - param[0][l1];
      trace(l1,0) = l1;
      trace(l1,1) = best2;
      trace(l1,2) = best3;
    }
  }
  else if (var_idx == 1) {

    forward.resize_dirty(nLabels2);
    trace.resize_dirty(nLabels2,3);

    for (uint l2=0; l2 < nLabels2; l2++) {

      double best = 1e300;
      uint best1 = 0;
      uint best3 = 0;

      for (uint l1=0; l1 < nLabels1; l1++) {

        const double w1 = param[0][l1];

        for (uint l3=0; l3 < nLabels3; l3++) {

          double hyp = cost_(l1,l2,l3) - w1 - param[2][l3];

          for (uint s=0; s < nSeps; s++)
            hyp -= eval_pair(s,l1,l2,l3,pair_param);
          
          if (hyp < best) {
            best = hyp;
            best1 = l1;
            best3 = l3;
          }
        }
      }

      forward[l2] = best - param[1][l2];
      trace(l2,0) = best1;
      trace(l2,1) = l2;
      trace(l2,2) = best3;
    }
  }
  else {
    assert(var_idx == 2);

    forward.resize_dirty(nLabels3);
    trace.resize_dirty(nLabels3,3);

    for (uint l3=0; l3 < nLabels3; l3++) {

      double best = 1e300;
      uint best1 = 0;
      uint best2 = 0;

      for (uint l1=0; l1 < nLabels1; l1++) {

        const double w1 = param[0][l1];

        for (uint l2=0; l2 < nLabels2; l2++) {

          double hyp = cost_(l1,l2,l3) - w1 - param[1][l2];

          for (uint s=0; s < nSeps; s++)
            hyp -= eval_pair(s,l1,l2,l3,pair_param);
          
          if (hyp < best) {
            best = hyp;
            best1 = l1;
            best2 = l2;
          }
        }
      }

      forward[l3] = best - param[2][l3];
      trace(l3,0) = best1;
      trace(l3,1) = best2;
      trace(l3,2) = l3;
    }
  }

  return 0.0;
}


/*virtual*/ 
double TernarySepChainDDFactor::compute_forward(const SepChainDDPairSeparator* incoming_sep, const SepChainDDVar* incoming_var, 
                                                const SepChainDDPairSeparator* outgoing_sep,
                                                const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                                Math2D::Matrix<double>& forward, 
                                                Math3D::Tensor<uint>& trace) const {

  assert(incoming_sep != outgoing_sep);

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();

  const uint nSeps = involved_separator_.size();

  uint pair_idx = 0;

  Storage1D<Math1D::Vector<double> > param = dual_var_;
  for (uint k=0; k < involved_var_.size(); k++) {
    if (involved_var_[k] == incoming_var) {
      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= prev_var_forward[l];
      }
    }
    else {

      if (incoming_sep != 0 && (incoming_sep->var1() == involved_var_[k] || incoming_sep->var2() == involved_var_[k] ) ) {
        //in this case we only want the dual vars of the factor
      }
      else {

        for (uint l=0; l < param[k].size(); l++) {
          param[k][l] -= involved_var_[k]->cost()[l];
        }
      }
    }
  }

  Storage1D<Math2D::Matrix<double> > pair_param = dual_pair_var_;
  for (uint s=0; s < nSeps; s++) {

    if (involved_separator_[s] == incoming_sep) {
      pair_param[s] -= prev_pair_forward; 
    }
    else if (involved_separator_[s] == outgoing_sep) {
      pair_idx = s;
    }
  }

  assert(pair_idx != MAX_UINT);

  SepChainDDVar* v1 = involved_separator_[pair_idx]->var1();
  SepChainDDVar* v2 = involved_separator_[pair_idx]->var2();

  if (v1 == involved_var_[0]) {

    if (v2 == involved_var_[1]) {

      forward.resize_dirty(nLabels1,nLabels2);
      trace.resize_dirty(nLabels1,nLabels2,3);

      for (uint l1=0; l1 < nLabels1; l1++) {

        for (uint l2=0; l2 < nLabels2; l2++) {

          double best = 1e300;
          uint arg_best = 0;
          
          for (uint l3 = 0; l3 < nLabels3; l3++) {

            double hyp = cost_(l1,l2,l3) - param[2][l3];
            
            for (uint s=0; s < nSeps; s++)
              hyp -= eval_pair(s,l1,l2,l3,pair_param);

            if (hyp < best) {
              best = hyp;
              arg_best = l3;
            }
          }

          forward(l1,l2) = best - param[0][l1] - param[1][l2];
          trace(l1,l2,0) = l1;
          trace(l1,l2,1) = l2;
          trace(l1,l2,2) = arg_best;
        }
      }

    }
    else {
      assert(v2 == involved_var_[2]);

      forward.resize_dirty(nLabels1,nLabels3);
      trace.resize_dirty(nLabels1,nLabels3,3);

      for (uint l1=0; l1 < nLabels1; l1++) {

        for (uint l3 = 0; l3 < nLabels3; l3++) {

          double best = 1e300;
          uint arg_best = 0;
          
          for (uint l2 = 0; l2 < nLabels2; l2++) {

            double hyp = cost_(l1,l2,l3) - param[1][l2];
            
            for (uint s=0; s < nSeps; s++)
              hyp -= eval_pair(s,l1,l2,l3,pair_param);

            if (hyp < best) {
              best = hyp;
              arg_best = l2;
            }
          }

          forward(l1,l3) = best - param[0][l1] - param[2][l3];
          trace(l1,l3,0) = l1;
          trace(l1,l3,1) = arg_best;
          trace(l1,l3,2) = l3;
        }
      }
    }
  }
  else {
    assert(v1 == involved_var_[1]);
    
    forward.resize_dirty(nLabels2,nLabels3);
    trace.resize_dirty(nLabels2,nLabels3,3);
    
    for (uint l2=0; l2 < nLabels2; l2++) {
      
      for (uint l3 = 0; l3 < nLabels3; l3++) {
        
        double best = 1e300;
        uint arg_best = 0;
        
        for (uint l1 = 0; l1 < nLabels1; l1++) {
          
          double hyp = cost_(l1,l2,l3) - param[0][l1];
          
          for (uint s=0; s < nSeps; s++)
            hyp -= eval_pair(s,l1,l2,l3,pair_param);

          if (hyp < best) {
            best = hyp;
            arg_best = l1;
          }
        }
        
        forward(l2,l3) = best - param[1][l2] - param[2][l3];
        trace(l2,l3,0) = arg_best;
        trace(l2,l3,1) = l2;
        trace(l2,l3,2) = l3;
      }
    }
  }
  
  return 0.0; //presently not subtracting an offset
}

double TernarySepChainDDFactor::eval_pair(uint s, uint x, uint y, uint z, const Storage1D< Math2D::Matrix<double> >& pair_param) const {

  SepChainDDVar* v1 = involved_separator_[s]->var1();
  SepChainDDVar* v2 = involved_separator_[s]->var2();
  
  uint a=MAX_UINT;
  uint b=MAX_UINT;
  
  if (involved_var_[0] == v1)
    a = x;
  else { 
    assert(involved_var_[1] == v1);
    a = y;
  }

  if (involved_var_[1] == v2)
    b = y;
  else {
    assert(involved_var_[2] == v2);
    b = z;
  }
  
  return pair_param[s](a,b);
}

/***********************************/

FourthOrderSepChainDDFactor::FourthOrderSepChainDDFactor(const Storage1D<SepChainDDVar*>& involved_vars, 
                                                         Storage1D<SepChainDDPairSeparator*>& separators,
                                                         const Storage1D<Math3D::Tensor<float> >& cost) 
  : SepChainDDFactor(involved_vars,separators), cost_(cost) {

  if (involved_vars.size() != 4) {
    INTERNAL_ERROR << " attempt to instantiate ternary factor with " 
                   << involved_vars.size() << " variables. Exiting." << std::endl;
    exit(1);
  }

  if (cost.size() < involved_vars[0]->nLabels() || cost[0].xDim() < involved_vars[1]->nLabels() 
      || cost[0].yDim() < involved_vars[2]->nLabels() || cost[0].zDim() < involved_vars[3]->nLabels()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
    exit(1);
  }

}

/*virtual*/ FourthOrderSepChainDDFactor::~FourthOrderSepChainDDFactor() {}

/*virtual */
double FourthOrderSepChainDDFactor::compute_forward(const SepChainDDPairSeparator* incoming_sep, const SepChainDDVar* incoming_var, 
                                                    const SepChainDDVar* outgoing_var,
                                                    const Math2D::Matrix<double>& prev_pair_forward, 
                                                    const Math1D::Vector<double>& prev_var_forward, 
                                                    Math1D::Vector<double>& forward, 
                                                    Math2D::Matrix<uint>& trace) const {

  assert(incoming_var != outgoing_var);

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  const uint nLabels4 = involved_var_[3]->nLabels();

  const uint nSeps = involved_separator_.size();

  uint var_idx = 0;

  Storage1D<Math1D::Vector<double> > param = dual_var_;
  for (uint k=0; k < involved_var_.size(); k++) {
    if (involved_var_[k] == incoming_var) {
      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= prev_var_forward[l];
      }
    }
    else {

      if (incoming_sep != 0 && (incoming_sep->var1() == involved_var_[k] || incoming_sep->var2() == involved_var_[k] ) ) {
        //in this case we only want the dual vars of the factor

        assert(involved_var_[k] != outgoing_var);
      }
      else {

        if (involved_var_[k] == outgoing_var)
          var_idx = k;

        for (uint l=0; l < param[k].size(); l++) {
          param[k][l] -= involved_var_[k]->cost()[l];
        }
      }
    }
  }

  Storage1D<Math2D::Matrix<double> > pair_param = dual_pair_var_;
  for (uint s=0; s < nSeps; s++) {

    if (involved_separator_[s] == incoming_sep) {
      pair_param[s] -= prev_pair_forward; 
    }
  }

  assert(var_idx != MAX_UINT);
  
  if (var_idx == 0) {
    
    forward.resize_dirty(nLabels1);
    trace.resize_dirty(nLabels1,4);

    for (uint l1=0; l1 < nLabels1; l1++) {

      double best = 1e300;
      uint best2 = 0;
      uint best3 = 0;
      uint best4 = 0;

      const Math3D::Tensor<float>& cur_cost = cost_[l1];

      for (uint l2=0; l2 < nLabels2; l2++) {

        const double w2 = param[1][l2];
        
        for (uint l3=0; l3 < nLabels3; l3++) {

          const double w3 = w2 + param[2][l3];

          for (uint l4=0; l4 < nLabels4; l4++) {
          
            double hyp = cur_cost(l2,l3,l4) - w3 - param[3][l4];

            for (uint s=0; s < nSeps; s++)
              hyp -= eval_pair(s,l1,l2,l3,l4,pair_param);
          
            if (hyp < best) {
              best = hyp;
              best2 = l2;
              best3 = l3;
              best4 = l4;
            }
          }
        }

        forward[l1] = best - param[0][l1];
        trace(l1,0) = l1;
        trace(l1,1) = best2;
        trace(l1,2) = best3;
        trace(l1,3) = best4;
      }
    }
  }
  else if (var_idx == 1) {

    forward.resize_dirty(nLabels2);
    trace.resize_dirty(nLabels2,4);

    for (uint l2=0; l2 < nLabels2; l2++) {

      double best = 1e300;
      uint best1 = 0;
      uint best3 = 0;
      uint best4 = 0;

      for (uint l1=0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost_[l1];

        const double w1 = param[0][l1];

        for (uint l3=0; l3 < nLabels3; l3++) {

          const double w3 = w1 + param[2][l3];

          for (uint l4=0; l4 < nLabels4; l4++) {

            double hyp = cur_cost(l2,l3,l4) - w3 - param[3][l4];

            for (uint s=0; s < nSeps; s++)
              hyp -= eval_pair(s,l1,l2,l3,l4,pair_param);
            
            if (hyp < best) {
              best = hyp;
              best1 = l1;
              best3 = l3;
              best4 = l4;
            }
          }
        }
      }

      forward[l2] = best - param[1][l2];
      trace(l2,0) = best1;
      trace(l2,1) = l2;
      trace(l2,2) = best3;
      trace(l2,3) = best4;
    }
  }
  else if (var_idx == 2) {

    forward.resize_dirty(nLabels3);
    trace.resize_dirty(nLabels3,4);

    for (uint l3=0; l3 < nLabels3; l3++) {

      double best = 1e300;
      uint best1 = 0;
      uint best2 = 0;
      uint best4 = 0;

      for (uint l1=0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost_[l1];

        const double w1 = param[0][l1];

        for (uint l2=0; l2 < nLabels2; l2++) {

          const double w2 = w1 + param[1][l2];

          for (uint l4=0; l4 < nLabels4; l4++) {

            double hyp = cur_cost(l2,l3,l4) - w2 - param[3][l4];

            for (uint s=0; s < nSeps; s++)
              hyp -= eval_pair(s,l1,l2,l3,l4,pair_param);
            
            if (hyp < best) {
              best = hyp;
              best1 = l1;
              best2 = l2;
              best4 = l4;
            }
          }
        }
      }

      forward[l3] = best - param[2][l3];
      trace(l3,0) = best1;
      trace(l3,1) = best2;
      trace(l3,2) = l3;
      trace(l3,3) = best4;
    }
  }
  else {

    forward.resize_dirty(nLabels4);
    trace.resize_dirty(nLabels4,4);

    for (uint l4=0; l4 < nLabels4; l4++) {

      double best = 1e300;
      uint best1 = 0;
      uint best2 = 0;
      uint best3 = 0;

      for (uint l1=0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost_[l1];

        const double w1 = param[0][l1];

        for (uint l2=0; l2 < nLabels2; l2++) {

          const double w2 = w1 + param[1][l2];

          for (uint l3=0; l3 < nLabels3; l3++) {

            double hyp = cur_cost(l2,l3,l4) - w2 - param[2][l3];

            for (uint s=0; s < nSeps; s++)
              hyp -= eval_pair(s,l1,l2,l3,l4,pair_param);
            
            if (hyp < best) {
              best = hyp;
              best1 = l1;
              best2 = l2;
              best3 = l3;
            }            
          }
        }
      }

      forward[l4] = best - param[3][l4];
      trace(l4,0) = best1;
      trace(l4,1) = best2;
      trace(l4,2) = best3;
      trace(l4,3) = l4;      
    }
  }

  return 0.0;
}

/*virtual*/ 
double FourthOrderSepChainDDFactor::compute_forward(const SepChainDDPairSeparator* incoming_sep, const SepChainDDVar* incoming_var, 
                                                    const SepChainDDPairSeparator* outgoing_sep,
                                                    const Math2D::Matrix<double>& prev_pair_forward, 
                                                    const Math1D::Vector<double>& prev_var_forward, 
                                                    Math2D::Matrix<double>& forward, 
                                                    Math3D::Tensor<uint>& trace) const {

  assert(incoming_sep != outgoing_sep);

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  const uint nLabels4 = involved_var_[3]->nLabels();

  const uint nSeps = involved_separator_.size();

  uint pair_idx = 0;

  Storage1D<Math1D::Vector<double> > param = dual_var_;
  for (uint k=0; k < involved_var_.size(); k++) {
    if (involved_var_[k] == incoming_var) {
      for (uint l=0; l < param[k].size(); l++) {
        param[k][l] -= prev_var_forward[l];
      }
    }
    else {

      if (incoming_sep != 0 && (incoming_sep->var1() == involved_var_[k] || incoming_sep->var2() == involved_var_[k] ) ) {
        //in this case we only want the dual vars of the factor
      }
      else {

        for (uint l=0; l < param[k].size(); l++) {
          param[k][l] -= involved_var_[k]->cost()[l];
        }
      }
    }
  }

  Storage1D<Math2D::Matrix<double> > pair_param = dual_pair_var_;
  for (uint s=0; s < nSeps; s++) {

    if (involved_separator_[s] == incoming_sep) {
      pair_param[s] -= prev_pair_forward; 
    }
    else if (involved_separator_[s] == outgoing_sep) {
      pair_idx = s;
    }
  }

  assert(pair_idx != MAX_UINT);

  SepChainDDVar* v1 = involved_separator_[pair_idx]->var1();
  SepChainDDVar* v2 = involved_separator_[pair_idx]->var2();


  if (v1 == involved_var_[0]) {

    if (v2 == involved_var_[1]) {

      forward.resize_dirty(nLabels1,nLabels2);
      trace.resize_dirty(nLabels1,nLabels2,4);

      for (uint l1=0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost_[l1];

        for (uint l2=0; l2 < nLabels2; l2++) {

          double best = 1e300;
          uint arg_best3 = 0;
          uint arg_best4 = 0;
          
          for (uint l3 = 0; l3 < nLabels3; l3++) {

            const double w3 = param[2][l3];

            for (uint l4 = 0; l4 < nLabels4; l4++) {

              double hyp = cur_cost(l2,l3,l4) - w3 - param[3][l4];
            
              for (uint s=0; s < nSeps; s++)
                hyp -= eval_pair(s,l1,l2,l3,l4,pair_param);
              
              if (hyp < best) {
                best = hyp;
                arg_best3 = l3;
                arg_best4 = l4;
              }
            }
          }
           
          forward(l1,l2) = best - param[0][l1] - param[1][l2];
          trace(l1,l2,0) = l1;
          trace(l1,l2,1) = l2;
          trace(l1,l2,2) = arg_best3;
          trace(l1,l2,3) = arg_best4;
        }
      }

    }
    else if (v2 == involved_var_[2]) {

      forward.resize_dirty(nLabels1,nLabels3);
      trace.resize_dirty(nLabels1,nLabels3,4);

      for (uint l1=0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost_[l1];

        for (uint l3 = 0; l3 < nLabels3; l3++) {

          double best = 1e300;
          uint arg_best2 = 0;
          uint arg_best4 = 0;
          
          for (uint l2 = 0; l2 < nLabels2; l2++) {

            const double w2 = param[1][l2];

            for (uint l4 = 0; l4 < nLabels4; l4++) {

              double hyp = cur_cost(l2,l3,l4) - w2 - param[3][l4];
              
              for (uint s=0; s < nSeps; s++)
                hyp -= eval_pair(s,l1,l2,l3,l4,pair_param);
              
              if (hyp < best) {
                best = hyp;
                arg_best2 = l2;
                arg_best4 = l4;
              }
            }
          }

          forward(l1,l3) = best - param[0][l1] - param[2][l3];
          trace(l1,l3,0) = l1;
          trace(l1,l3,1) = arg_best2;
          trace(l1,l3,2) = l3;
          trace(l1,l3,3) = arg_best4;
        }
      }
    }
    else {

      assert(v2 == involved_var_[3]);

      forward.resize_dirty(nLabels1,nLabels4);
      trace.resize_dirty(nLabels1,nLabels4,4);


      for (uint l1=0; l1 < nLabels1; l1++) {

        const Math3D::Tensor<float>& cur_cost = cost_[l1];

        for (uint l4 = 0; l4 < nLabels4; l4++) {

          double best = 1e300;
          uint arg_best2 = 0;
          uint arg_best3 = 0;

          for (uint l2 = 0; l2 < nLabels2; l2++) {

            const double w2 = param[1][l2];

            for (uint l3 = 0; l3 < nLabels3; l3++) {

              double hyp = cur_cost(l2,l3,l4) - w2 - param[2][l3];
              
              for (uint s=0; s < nSeps; s++)
                hyp -= eval_pair(s,l1,l2,l3,l4,pair_param);
              
              if (hyp < best) {
                best = hyp;
                arg_best2 = l2;
                arg_best3 = l3;
              }
            }
          }

          forward(l1,l4) = best - param[0][l1] - param[3][l4];
          trace(l1,l4,0) = l1;
          trace(l1,l4,1) = arg_best2;
          trace(l1,l4,2) = arg_best3;
          trace(l1,l4,3) = l4;          
        }
      }
    }
  }
  else if (v1 == involved_var_[1]) {
    
    if (v2 == involved_var_[2]) {

      forward.resize_dirty(nLabels2,nLabels3);
      trace.resize_dirty(nLabels2,nLabels3,4);
      
      for (uint l2=0; l2 < nLabels2; l2++) {
        
        for (uint l3 = 0; l3 < nLabels3; l3++) {
          
          double best = 1e300;
          uint arg_best1 = 0;
          uint arg_best4 = 0;
          
          for (uint l1 = 0; l1 < nLabels1; l1++) {
            
            const Math3D::Tensor<float>& cur_cost = cost_[l1];

            const double w1 = param[0][l1];

            for (uint l4 = 0; l4 < nLabels4; l4++) {
              
              double hyp = cur_cost(l2,l3,l4) - w1 - param[3][l4];
              
              for (uint s=0; s < nSeps; s++)
                hyp -= eval_pair(s,l1,l2,l3,l4,pair_param);
              
              if (hyp < best) {
                best = hyp;
                arg_best1 = l1;
                arg_best4 = l4;
              }
            }
          }
          
          forward(l2,l3) = best - param[1][l2] - param[2][l3];
          trace(l2,l3,0) = arg_best1;
          trace(l2,l3,1) = l2;
          trace(l2,l3,2) = l3;
          trace(l2,l3,3) = arg_best4;
        }
      }
    }
    else {

      assert(v2 == involved_var_[3]);

      forward.resize_dirty(nLabels2,nLabels4);
      trace.resize_dirty(nLabels2,nLabels4,4);

      for (uint l2=0; l2 < nLabels2; l2++) {
        
        for (uint l4 = 0; l4 < nLabels4; l4++) {
          
          double best = 1e300;
          uint arg_best1 = 0;
          uint arg_best3 = 0;

          for (uint l1 = 0; l1 < nLabels1; l1++) {
            
            const Math3D::Tensor<float>& cur_cost = cost_[l1];

            const double w1 = param[0][l1];

            for (uint l3 = 0; l3 < nLabels3; l3++) {
              
              double hyp = cur_cost(l2,l3,l4) - w1 - param[2][l3];
              
              for (uint s=0; s < nSeps; s++)
                hyp -= eval_pair(s,l1,l2,l3,l4,pair_param);
              
              if (hyp < best) {
                best = hyp;
                arg_best1 = l1;
                arg_best3 = l3;
              }
            }
          }

          forward(l2,l4) = best - param[1][l2] - param[3][l4];
          trace(l2,l4,0) = arg_best1;
          trace(l2,l4,1) = l2;
          trace(l2,l4,2) = arg_best3;
          trace(l2,l4,3) = l4;
        }
      }
    }
  }
  else {

    assert(v1 == involved_var_[2]);
    assert(v2 == involved_var_[3]);

    forward.resize_dirty(nLabels3,nLabels4);
    trace.resize_dirty(nLabels3,nLabels4,4);

    for (uint l3=0; l3 < nLabels3; l3++) {
        
      for (uint l4 = 0; l4 < nLabels4; l4++) {
        
        double best = 1e300;
        uint arg_best1 = 0;
        uint arg_best2 = 0;
        
        for (uint l1 = 0; l1 < nLabels1; l1++) {
            
          const Math3D::Tensor<float>& cur_cost = cost_[l1];

          const double w1 = param[0][l1];

          for (uint l2 = 0; l2 < nLabels2; l2++) {
              
            double hyp = cur_cost(l2,l3,l4) - w1 - param[1][l2];
              
            for (uint s=0; s < nSeps; s++)
              hyp -= eval_pair(s,l1,l2,l3,l4,pair_param);
            
            if (hyp < best) {
              best = hyp;
              arg_best1 = l1;
              arg_best2 = l2;
            }
          }
        }
        
        forward(l3,l4) = best - param[2][l3] - param[3][l4];
        trace(l3,l4,0) = arg_best1;
        trace(l3,l4,1) = arg_best2;
        trace(l3,l4,2) = l3;
        trace(l3,l4,3) = l4;
      }
    }
  }
    
  return 0.0; //presently not subtracting an offset
}

double FourthOrderSepChainDDFactor::eval_pair(uint s, uint x, uint y, uint z, uint w, 
                                              const Storage1D< Math2D::Matrix<double> >& pair_param) const {

  SepChainDDVar* v1 = involved_separator_[s]->var1();
  SepChainDDVar* v2 = involved_separator_[s]->var2();
  
  uint a=MAX_UINT;
  uint b=MAX_UINT;
  
  if (involved_var_[0] == v1)
    a = x;
  else if (involved_var_[1] == v1) {
    a = y;
  }
  else {
    assert(involved_var_[2] == v1);
    a = z;
  }

  if (involved_var_[1] == v2)
    b = y;
  else if (involved_var_[2] == v2) {
    b = z;
  }
  else {
    assert(involved_var_[3] == v2);
    b = w;
  }
  
  return pair_param[s](a,b);

}


/***********************************/

SeparatorChainDualDecomposition::SeparatorChainDualDecomposition(uint nVars, uint nSeparators, uint nFactors) :
  nUsedVars_(0), nUsedSeparators_(0), nUsedFactors_(0), optimize_called_(false) {
  
  var_.resize(nVars);
  separator_.resize(nSeparators);
  factor_.resize(nFactors);
}


SeparatorChainDualDecomposition::~SeparatorChainDualDecomposition() {
  for (uint v=0; v < nUsedVars_; v++)
    delete var_[v];
  for (uint s=0; s < nUsedSeparators_; s++)
    delete separator_[s];
  for (uint f=0; f < nUsedFactors_; f++)
    delete factor_[f];
}

uint SeparatorChainDualDecomposition::add_var(const Math1D::Vector<float>& cost) {

  assert(!optimize_called_);

  if (nUsedVars_ == var_.size())
    var_.resize(uint(1.2*nUsedVars_)+4);

  assert(nUsedVars_ < var_.size());
  var_[nUsedVars_] = new SepChainDDVar(cost);

  nUsedVars_++;
  return (nUsedVars_-1);
}

uint SeparatorChainDualDecomposition::add_pair_separator(uint var1, uint var2) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  assert(!optimize_called_);

  if (nUsedSeparators_ == separator_.size())
    separator_.resize(uint(1.2*nUsedSeparators_)+4);

  assert(nUsedSeparators_ < separator_.size());

  assert(var1 < var2);

  separator_[nUsedSeparators_] = new SepChainDDPairSeparator(var_[var1],var_[var2]);

  nUsedSeparators_++;
  return nUsedSeparators_-1;
}

void SeparatorChainDualDecomposition::add_factor(SepChainDDFactor* fac) {

  assert(!optimize_called_);
  
  if (nUsedFactors_ == factor_.size())
    factor_.resize(uint(1.2*nUsedFactors_)+4);

  factor_[nUsedFactors_] = fac;
  nUsedFactors_++;
}

void SeparatorChainDualDecomposition::add_binary_factor(uint v1, uint v2, const Math2D::Matrix<float>& cost) {

  if (v1 >= nUsedVars_ || v2 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  Storage1D<SepChainDDVar*> vars(2);
  vars[0] = var_[v1];
  vars[1] = var_[v2];

  add_factor(new BinarySepChainDDFactor(vars,cost));
}

void SeparatorChainDualDecomposition::add_ternary_factor(uint v1, uint v2, uint v3, 
                                                         const Storage1D<uint>& separators, 
                                                         const Math3D::Tensor<float>& cost) {

  if (v1 >= nUsedVars_ || v2 >= nUsedVars_ || v3 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  Math3D::Tensor<float> cost_copy = cost;

#if 1
  if (v1 > v2) {

    if (var_[v1]->nLabels() != var_[v2]->nLabels())
      TODO("non-standard variable order with heterogeneous number of labels");

    for (uint x=0; x < var_[v1]->nLabels(); x++)
      for (uint y=x+1; y < var_[v2]->nLabels(); y++)
        for (uint z=0; z < var_[v3]->nLabels(); z++) {
          std::swap(cost_copy(x,y,z),cost_copy(y,x,z));
        }

    std::swap(v1,v2);
  }
  if (v2 > v3) {

    if (var_[v2]->nLabels() != var_[v3]->nLabels())
      TODO("non-standard variable order with heterogeneous number of labels");

    for (uint x=0; x < var_[v1]->nLabels(); x++)
      for (uint y=0; y < var_[v2]->nLabels(); y++)
        for (uint z=y+1; z < var_[v3]->nLabels(); z++)
          std::swap(cost_copy(x,y,z),cost_copy(x,z,y));
    
    std::swap(v2,v3);
  }
  if (v1 > v2) {

    if (var_[v1]->nLabels() != var_[v2]->nLabels())
      TODO("non-standard variable order with heterogeneous number of labels");

    for (uint x=0; x < var_[v1]->nLabels(); x++)
      for (uint y=x+1; y < var_[v2]->nLabels(); y++)
        for (uint z=0; z < var_[v3]->nLabels(); z++)
          std::swap(cost_copy(x,y,z),cost_copy(y,x,z));

    std::swap(v1,v2);
  }

  assert(v1 < v2);
  assert(v2 < v3);
#endif


  Storage1D<SepChainDDVar*> vars(3);
  vars[0] = var_[v1];
  vars[1] = var_[v2];
  vars[2] = var_[v3];

  Storage1D<SepChainDDPairSeparator*> seps(separators.size());
  for (uint s=0; s < separators.size(); s++) {
    if (separators[s] >= nUsedSeparators_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    seps[s] = separator_[separators[s]];
  }

  add_factor(new TernarySepChainDDFactor(vars,seps,cost_copy));
}

void SeparatorChainDualDecomposition::add_fourth_order_factor(uint v1, uint v2, uint v3, uint v4, 
                                                              const Storage1D<uint>& separators, 
                                                              const Storage1D<Math3D::Tensor<float> >& cost) {

  if (v1 >= nUsedVars_ || v2 >= nUsedVars_ || v3 >= nUsedVars_ || v4 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  Storage1D<Math3D::Tensor<float> > cost_copy = cost;

#if 1
  bool changed = true;
  while (changed) {

    changed = false;

    if (v1 > v2) {
      changed = true;

      if (var_[v1]->nLabels() != var_[v2]->nLabels())
        TODO("non-standard variable order with heterogeneous number of labels");

      for (uint x=0; x < var_[v1]->nLabels(); x++)
        for (uint y=x+1; y < var_[v2]->nLabels(); y++)
          for (uint z=0; z < var_[v3]->nLabels(); z++)
            for (uint w=0; w < var_[v4]->nLabels(); w++)
              std::swap(cost_copy[x](y,z,w), cost_copy[y](x,z,w));

      std::swap(v1,v2);
    }
    if (v2 > v3) {
      changed = true;

      if (var_[v2]->nLabels() != var_[v3]->nLabels())
        TODO("non-standard variable order with heterogeneous number of labels");

      for (uint x=0; x < var_[v1]->nLabels(); x++)
        for (uint y=0; y < var_[v2]->nLabels(); y++)
          for (uint z=y+1; z < var_[v3]->nLabels(); z++)
            for (uint w=0; w < var_[v4]->nLabels(); w++)
              std::swap(cost_copy[x](y,z,w), cost_copy[x](z,y,w));

      std::swap(v2,v3);
    }
    if (v3 > v4) {
      changed = true;

      if (var_[v3]->nLabels() != var_[v4]->nLabels())
        TODO("non-standard variable order with heterogeneous number of labels");

      for (uint x=0; x < var_[v1]->nLabels(); x++)
        for (uint y=0; y < var_[v2]->nLabels(); y++)
          for (uint z=0; z < var_[v3]->nLabels(); z++)
            for (uint w=z+1; w < var_[v4]->nLabels(); w++)
              std::swap(cost_copy[x](y,z,w), cost_copy[x](y,w,z));

      std::swap(v3,v4);
    }
  }


  assert(v1 < v2);
  assert(v2 < v3);
  assert(v3 < v4);
#endif

  Storage1D<SepChainDDVar*> vars(4);
  vars[0] = var_[v1];
  vars[1] = var_[v2];
  vars[2] = var_[v3];
  vars[3] = var_[v4];

  Storage1D<SepChainDDPairSeparator*> seps(separators.size());
  for (uint s=0; s < separators.size(); s++) {

    if (separators[s] >= nUsedSeparators_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    seps[s] = separator_[separators[s]];
  }

  add_factor(new FourthOrderSepChainDDFactor(vars,seps,cost_copy));
}

const Math1D::Vector<uint>& SeparatorChainDualDecomposition::labeling() {
  return labeling_;
}

void SeparatorChainDualDecomposition::set_up_chains() {

  uint nChains = 0;
  uint nAtLeast5 = 0;
  uint nAtLeast10 = 0;
  uint nAtLeast25 = 0;

  for (uint pass=0; pass < 2; pass++) {

    for (uint f=0; f < nUsedFactors_; f++) {
      
      //std::cerr << "factor #" << f << "/" << nUsedFactors_ << std::endl;
      
      if (factor_[f]->prev_var() == 0 && factor_[f]->prev_sep() == 0
          && factor_[f]->next_var() == 0 && factor_[f]->next_sep() == 0) {
        
        uint length = 1;
        
        std::set<SepChainDDVar*> current_vars;
        
        bool extension_found = true;
        
        SepChainDDFactor* cur_factor = factor_[f];
        
        if (pass > 0) {

          SepChainDDFactor* temp_fac = cur_factor;

          while(temp_fac != 0) {

            const Storage1D<SepChainDDVar*>& involved_vars = temp_fac->involved_vars();

            for (uint k=0; k < involved_vars.size(); k++) {
              current_vars.insert(involved_vars[k]);
            }
            
            temp_fac = temp_fac->prev_factor();
          }

          temp_fac = cur_factor->next_factor();

          while(temp_fac != 0) {

            const Storage1D<SepChainDDVar*>& involved_vars = temp_fac->involved_vars();

            for (uint k=0; k < involved_vars.size(); k++) {
              current_vars.insert(involved_vars[k]);
            }
            
            temp_fac = temp_fac->next_factor();
          }
        }

#if 1
        //extend lower end
        while (extension_found) {
          
          extension_found = false;
          
          const Storage1D<SepChainDDVar*>& involved_vars = cur_factor->involved_vars();
          const Storage1D<SepChainDDPairSeparator*>& involved_separators = cur_factor->involved_separators();
          
          for (uint k=0; k < involved_vars.size(); k++) {
            current_vars.insert(involved_vars[k]);
          }

          if (pass == 0) {
            for (uint s=0; s < involved_separators.size(); s++) {
            
              SepChainDDPairSeparator* sep = involved_separators[s];
            
              if (sep == cur_factor->next_sep())
                continue;
              if (cur_factor->next_var() != 0 && (sep->var1() == cur_factor->next_var() || sep->var2() == cur_factor->next_var()))
                continue;
              
              const Storage1D<SepChainDDFactor*>& adjacent_factor = sep->neighboring_factor();
              
              for (uint l=0; l < adjacent_factor.size(); l++) {
                
                SepChainDDFactor* hyp_factor = adjacent_factor[l];
                const Storage1D<SepChainDDVar*>& hyp_involved_vars = hyp_factor->involved_vars();
                
                if (hyp_factor != cur_factor && hyp_factor->prev_var() == 0 && hyp_factor->prev_sep() == 0
                    && hyp_factor->next_var() == 0 && hyp_factor->next_sep() == 0) {
                  
                  bool is_valid_extension = true;
                  
                  for (uint v=0; v < hyp_involved_vars.size(); v++) {
                    
                    if (hyp_involved_vars[v] != sep->var1() && hyp_involved_vars[v] != sep->var2() 
                        && current_vars.find(hyp_involved_vars[v]) != current_vars.end())
                      is_valid_extension = false;
                  }
                  
                  if (is_valid_extension) {
                  
                    extension_found = true;
                    cur_factor->set_prev_sep(sep);
                    cur_factor->set_prev_factor(hyp_factor);
                    
                    hyp_factor->set_next_sep(sep);
                    hyp_factor->set_next_factor(cur_factor);
                    
                    cur_factor = hyp_factor;
                    
                    length++;
                    
                    break;
                  }
                }
                
                if (extension_found)
                  break;
              }
            }
          }
          else {
            //second pass

            for (double k=0; k < involved_vars.size(); k++) {
              
              SepChainDDVar* var = involved_vars[k];

              if (var == cur_factor->next_var())
                continue;
              
              const Storage1D<SepChainDDFactor*>& adjacent_factor = var->neighboring_factor();
              
              for (uint l=0; l < adjacent_factor.size(); l++) {
                
                SepChainDDFactor* hyp_factor = adjacent_factor[l];
                const Storage1D<SepChainDDVar*>& hyp_involved_vars = hyp_factor->involved_vars();
                
                bool is_valid_extension = false;
                
                if (hyp_factor != cur_factor
                    && hyp_factor->next_var() == 0 && hyp_factor->next_sep() == 0) {
                  
                  is_valid_extension = true;
                  
                  for (uint v=0; v < hyp_involved_vars.size(); v++) {
                    
                    if (hyp_involved_vars[v] != var && current_vars.find(hyp_involved_vars[v]) != current_vars.end())
                      is_valid_extension = false;
                  }
                  
                  SepChainDDFactor* temp_fac = hyp_factor->prev_factor();
                  while (is_valid_extension && temp_fac != 0) {

                    const Storage1D<SepChainDDVar*>& temp_involved_vars = temp_fac->involved_vars();
                    
                    for (uint v=0; v < temp_involved_vars.size(); v++) {
                    
                      if (current_vars.find(temp_involved_vars[v]) != current_vars.end()) {
                        is_valid_extension = false;
                        break;
                      }
                    }

                    temp_fac = temp_fac->prev_factor();
                  }
                  
                  if (is_valid_extension) {
                    
                    extension_found = true;
                    cur_factor->set_prev_var(var);
                    cur_factor->set_prev_factor(hyp_factor);
                    
                    hyp_factor->set_next_var(var);
                    hyp_factor->set_next_factor(cur_factor);
                    
                    cur_factor = hyp_factor;
                    
                    length++;
                    
                    break;
                  }
                }
              }
              
              if (extension_found)
                break;
            }
          }
        }
#endif

#if 1
        //extend upper end
        cur_factor = factor_[f];
        extension_found = true;
        
        while (extension_found) {
          
          extension_found = false;
        
          const Storage1D<SepChainDDVar*>& involved_vars = cur_factor->involved_vars();
          const Storage1D<SepChainDDPairSeparator*>& involved_separators = cur_factor->involved_separators();
          
          for (uint k=0; k < involved_vars.size(); k++) {
            current_vars.insert(involved_vars[k]);
          }

          if (pass == 0) {
            for (uint s=0; s < involved_separators.size(); s++) {
            
              SepChainDDPairSeparator* sep = involved_separators[s];

              if (sep == cur_factor->prev_sep())
                continue;

              const Storage1D<SepChainDDFactor*>& adjacent_factor = sep->neighboring_factor();
              
              for (uint l=0; l < adjacent_factor.size(); l++) {
                
                SepChainDDFactor* hyp_factor = adjacent_factor[l];

                bool is_valid_extension = false;

                if (hyp_factor != cur_factor && hyp_factor->prev_var() == 0 && hyp_factor->prev_sep() == 0
                    && hyp_factor->next_var() == 0 && hyp_factor->next_sep() == 0) {

                  const Storage1D<SepChainDDVar*>& hyp_involved_vars = hyp_factor->involved_vars();

                  for (uint v=0; v < hyp_involved_vars.size(); v++) {
                    
                    if (hyp_involved_vars[v] != sep->var1() && hyp_involved_vars[v] != sep->var2() 
                        && current_vars.find(hyp_involved_vars[v]) != current_vars.end())
                      is_valid_extension = false;
                  }

                  if (is_valid_extension) {

                    extension_found = true;
                    cur_factor->set_next_sep(sep);
                    cur_factor->set_next_factor(hyp_factor);
                    
                    hyp_factor->set_prev_sep(sep);
                    hyp_factor->set_prev_factor(cur_factor);
                    
                    cur_factor = hyp_factor;
                    
                    length++;
                    
                    break;
                  }
                }                

                if (extension_found)
                  break;
              }
            }
          }
          else {
            //second pass (use variables as links)
            
            for (double k=0; k < involved_vars.size(); k++) {
	
              SepChainDDVar* var = involved_vars[k];
              
              if (var == cur_factor->prev_var())
                continue;
              
              const Storage1D<SepChainDDFactor*>& adjacent_factor = var->neighboring_factor();
              
              for (uint l=0; l < adjacent_factor.size(); l++) {
                
                SepChainDDFactor* hyp_factor = adjacent_factor[l];
                
                bool is_valid_extension = false;

                if (hyp_factor != cur_factor 
                    && hyp_factor->prev_var() == 0 && hyp_factor->prev_sep() == 0) {

                  const Storage1D<SepChainDDVar*>& hyp_involved_vars = hyp_factor->involved_vars();
                  
                  is_valid_extension = true;
                  
                  for (uint v=0; v < hyp_involved_vars.size(); v++) {
                    
                    if (hyp_involved_vars[v] != var && current_vars.find(hyp_involved_vars[v]) != current_vars.end())
                      is_valid_extension = false;
                  }
	      
                  SepChainDDFactor* temp_fac = hyp_factor->next_factor();
                  while (is_valid_extension && temp_fac != 0) {

                    const Storage1D<SepChainDDVar*>& temp_involved_vars = temp_fac->involved_vars();
                    
                    for (uint v=0; v < temp_involved_vars.size(); v++) {
                    
                      if (current_vars.find(temp_involved_vars[v]) != current_vars.end()) {
                        is_valid_extension = false;
                        break;
                      }
                    }

                    temp_fac = temp_fac->next_factor();
                  }


                  if (is_valid_extension) {

                    extension_found = true;
                    cur_factor->set_next_var(var);
                    cur_factor->set_next_factor(hyp_factor);
                    
                    hyp_factor->set_prev_var(var);
                    hyp_factor->set_prev_factor(cur_factor);
                    
                    cur_factor = hyp_factor;
                    
                    length++;
                    
                    break;
                  }
                }	    

                if (extension_found)
                  break;
              }
            }
          }
        }
#endif
      }
    }
  }

}

double SeparatorChainDualDecomposition::optimize(uint nIter, double start_step_size) {

  std::cerr.precision(10);

  std::cerr << "subgradient optimization" << std::endl;
  std::cerr << nUsedFactors_ << " factors" << std::endl;

  if (!optimize_called_) {
    set_up_chains();
    std::cerr << "chains created" << std::endl;
    
    for (uint v=0; v < nUsedVars_; v++)
      var_[v]->set_up_chains();
  }
  optimize_called_ = true;

  double best_dual = -1e300;

  labeling_.resize(nUsedVars_,0);

  Storage1D<Math1D::Vector<uint> > factor_label(nUsedFactors_);
  for (uint f=0; f < nUsedFactors_; f++)
    factor_label[f].resize(factor_[f]->involved_vars().size());

  std::map<SepChainDDVar*,uint> var_num;
  for (uint v=0; v < nUsedVars_; v++)
    var_num[var_[v]] = v;

  std::map<SepChainDDFactor*,uint> factor_num;
  for (uint f=0; f < nUsedFactors_; f++)
    factor_num[factor_[f]] = f;

  uint nIncreases = 1;


  size_t effort_per_iteration = 0;
    
  for (uint f=0; f < nUsedFactors_; f++) {
    
    uint var_size = factor_[f]->involved_vars().size();
    uint sep_size = factor_[f]->involved_separators().size();
    
    effort_per_iteration += (var_size + sep_size);
  }

  for (uint iter=1; iter < nIter; iter++) {

    double step_size = start_step_size / nIncreases;

    double cur_bound = 0.0;

    for (uint f=0; f < nUsedFactors_; f++) {

      //std::cerr << "f: " << f << std::endl;

      if (factor_[f]->prev_var() == 0 && factor_[f]->prev_sep() == 0) {

        SepChainDDFactor* cur_factor = factor_[f];
        
        std::vector<SepChainDDFactor*> chain;
        std::vector<SepChainDDVar*> out_var;
        std::vector<SepChainDDPairSeparator*> out_sep;

        //find chain start
        SepChainDDVar* in_var = 0;
        for (uint k=0; k < cur_factor->involved_vars().size(); k++) {

          if (cur_factor->involved_vars()[k] != cur_factor->next_var()
              && (cur_factor->next_sep() == 0 
                  || (cur_factor->next_sep()->var1() != cur_factor->involved_vars()[k] 
                      && cur_factor->next_sep()->var2() != cur_factor->involved_vars()[k]) ) ) {
            in_var = cur_factor->involved_vars()[k];
            break;
          }
        }

        assert(in_var != 0);

        while (cur_factor != 0) {
          chain.push_back(cur_factor);
          assert(cur_factor->next_var() != in_var);
          out_var.push_back(cur_factor->next_var());
          out_sep.push_back(cur_factor->next_sep());
          cur_factor = cur_factor->next_factor();
        }

        uint chain_length = chain.size();

        //find chain end
        for (uint k=0; k < chain.back()->involved_vars().size(); k++) {

          if (chain.back()->involved_vars()[k] != chain.back()->prev_var() 
              && (chain.back()->prev_sep() == 0
                  || (chain.back()->prev_sep()->var1() != chain.back()->involved_vars()[k]
                      && chain.back()->prev_sep()->var2() != chain.back()->involved_vars()[k]) ) 
              && chain.back()->involved_vars()[k] != in_var) { //can happen for chains of length 1
            out_var.back() = chain.back()->involved_vars()[k];
            break;
          }
        }

        Math1D::NamedVector<double> forward1(MAKENAME(forward1));
        Math1D::NamedVector<double> forward2(in_var->nLabels(),MAKENAME(forward2));
        for (uint l=0; l < forward2.size(); l++)
          forward2[l] = in_var->cost()[l];

        Math2D::Matrix<double> pair_forward1;
        Math2D::Matrix<double> pair_forward2;
        
        NamedStorage1D<Math2D::Matrix<uint> > var_trace(chain_length,MAKENAME(var_trace));
        NamedStorage1D<Math3D::Tensor<uint> > pair_trace(chain_length,MAKENAME(pair_trace));

        //compute forward
        if (out_var[0] != 0)
          cur_bound += chain[0]->compute_forward(0,in_var,out_var[0],pair_forward2,forward2,forward1,var_trace[0]);
        else
          cur_bound += chain[0]->compute_forward(0,in_var,out_sep[0],pair_forward2,forward2,pair_forward1,pair_trace[0]);

        for (uint k=1; k < chain_length; k++) {

          Math1D::Vector<double>& last_forward = ((k % 2) == 1) ? forward1 : forward2;
          Math1D::Vector<double>& new_forward = ((k % 2) == 0) ? forward1 : forward2;

          Math2D::Matrix<double>& last_pair_forward = ((k % 2) == 1) ? pair_forward1 : pair_forward2;
          Math2D::Matrix<double>& new_pair_forward = ((k % 2) == 0) ? pair_forward1 : pair_forward2;

          if (out_var[k] != 0) 
            cur_bound += chain[k]->compute_forward(out_sep[k-1],out_var[k-1],out_var[k],
                                                   last_pair_forward,last_forward,new_forward,var_trace[k]);
          else
            cur_bound += chain[k]->compute_forward(out_sep[k-1],out_var[k-1],out_sep[k],
                                                   last_pair_forward,last_forward,new_pair_forward,pair_trace[k]);          
        }

        //traceback
        Math1D::Vector<double>& total_forward = ((chain_length-1) % 2 == 0) ? forward1 : forward2;

        assert(total_forward.size() == out_var[chain_length-1]->nLabels());

        double best = 1e300;
        uint arg_best1 = MAX_UINT;
        uint arg_best2 = MAX_UINT;
        for (uint l=0; l < total_forward.size(); l++) {
          if (total_forward[l] < best) {

            best = total_forward[l];
            arg_best1 = l;
          }
        }

        assert(arg_best1 < MAX_UINT);

        cur_bound += best;

        for (int k=chain_length-1; k >= 0; k--) {

          Math1D::Vector<uint>& labeling = factor_label[factor_num[chain[k]]];	  

          if (out_var[k] != 0) {
            for (uint v=0; v < labeling.size(); v++)
              labeling[v] = var_trace[k](arg_best1,v);
          }
          else {
            for (uint v=0; v < labeling.size(); v++)
              labeling[v] = pair_trace[k](arg_best1,arg_best2,v);
          }

          //update arg_best
          if (k > 0) {
            for (uint v=0; v < labeling.size(); v++) {
	      
              if (out_var[k-1] != 0) {
                if (chain[k]->involved_vars()[v] == out_var[k-1]) {
                  arg_best1 = labeling[v];
                }
              }
              else {
                if (chain[k]->involved_vars()[v] == out_sep[k-1]->var1())
                  arg_best1 = labeling[v];
                else if (chain[k]->involved_vars()[v] == out_sep[k-1]->var2())
                  arg_best2 = labeling[v];
              }
            }
          }

        }
      }
    }

    if (cur_bound > best_dual) {
      best_dual = cur_bound;
    }
    else {
      nIncreases++;
    }

    std::cerr << "iter " << iter << ", cur bound: " << cur_bound << ", best ever: " << best_dual << std::endl;

    //std::cerr << "go in grad. direction" << std::endl;

    //take the next step in subgradient direction

    for (uint f=0; f < nUsedFactors_; f++) {

      for (uint k=0; k < factor_label[f].size(); k++) {
        
        const SepChainDDVar* cur_var = factor_[f]->involved_vars()[k];
        
        uint cur_fac_label = factor_label[f][k];
        factor_[f]->get_duals(cur_var)[cur_fac_label] -= step_size;
        assert(!isinf(factor_[f]->get_duals(cur_var)[cur_fac_label]));
      }

      for (uint s=0; s < factor_[f]->involved_separators().size(); s++) {

        const SepChainDDPairSeparator* cur_sep = factor_[f]->involved_separators()[s];

        const SepChainDDVar* v1 = cur_sep->var1();
        const SepChainDDVar* v2 = cur_sep->var2();

        uint l1 = MAX_UINT;
        uint l2 = MAX_UINT;
        
        for (uint k=0; k < factor_label[f].size(); k++) {
        
          const SepChainDDVar* cur_var = factor_[f]->involved_vars()[k];

          if (cur_var == v1)
            l1 = factor_label[f][k];
          else if (cur_var == v2)
            l2 = factor_label[f][k];
        }

        factor_[f]->get_pair_duals(cur_sep)(l1,l2) -= step_size;
        assert(!isinf(factor_[f]->get_pair_duals(cur_sep)(l1,l2)));
      }
    }

    //std::cerr << "re-project" << std::endl;

    //re-project
    for (uint v=0; v < nUsedVars_; v++) {
      
      if (var_[v]->neighboring_factor().size() > 0) {

        Math1D::Vector<double> sum(var_[v]->nLabels(),0.0);
	  
        for (uint k=0; k < var_[v]->neighboring_factor().size(); k++) {
          sum += var_[v]->neighboring_factor()[k]->get_duals(var_[v]);
        }
        
        sum *= 1.0 / var_[v]->neighboring_factor().size();
	  
        for (uint k=0; k < var_[v]->neighboring_factor().size(); k++) {
	    
          var_[v]->neighboring_factor()[k]->get_duals(var_[v]) -= sum;
          for (uint l=0; l < var_[v]->nLabels(); l++)
            assert(!isinf(var_[v]->neighboring_factor()[k]->get_duals(var_[v])[l]));
        }
      }
    }
    for (uint s=0; s < nUsedSeparators_; s++) {

      SepChainDDPairSeparator* cur_sep = separator_[s];

      if (cur_sep->neighboring_factor().size() > 0) {

        Math2D::Matrix<double> sum(cur_sep->var1()->nLabels(),cur_sep->var2()->nLabels(),0.0);

        for (uint k=0; k < cur_sep->neighboring_factor().size(); k++) {
          sum += cur_sep->neighboring_factor()[k]->get_pair_duals(cur_sep);
        }
        
        sum *= 1.0 / cur_sep->neighboring_factor().size();

        for (uint k=0; k < cur_sep->neighboring_factor().size(); k++) {
          cur_sep->neighboring_factor()[k]->get_pair_duals(cur_sep) -= sum;
        }
      }
    }
  }

  for (uint f=0; f < nUsedFactors_; f++) {
      
    for (uint k=0; k < factor_label[f].size(); k++) {
      
      SepChainDDVar* cur_var = factor_[f]->involved_vars()[k];
      
      labeling_[var_num[cur_var]] = factor_label[f][k];
    }
  }


  size_t message_effort = 0;
    
  for (uint f=0; f < nUsedFactors_; f++) {
    
    uint var_size = factor_[f]->involved_vars().size();
    uint sep_size = factor_[f]->involved_separators().size();
    
    message_effort += (var_size + sep_size);
  }
  message_effort *= nIter;

  std::cerr << "message effort: " << message_effort << std::endl;

  return best_dual;
}

