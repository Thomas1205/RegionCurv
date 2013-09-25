/**** written by Thomas Schoenemann as an employee of the University of Pisa, Italy, Oct. 2011 ****/
/**** and continued at the University of Düsseldorf, Germany, 2012 ***/

#include "separatorDualOpt.hh"

#include <map>


sepDualOptVar::sepDualOptVar(const Math1D::Vector<float>& cost) : cost_(cost) {}

void sepDualOptVar::add_factor(sepDualOptFactor* adjacent_fac) {

  uint size = adjacent_factor_.size();
  adjacent_factor_.resize(size+1);
  adjacent_factor_[size] = adjacent_fac;
}

void sepDualOptVar::add_cost(const Math1D::Vector<float>& cost) {

  cost_ += cost;
}

void sepDualOptVar::add_pair_separator(sepDualOptPairSeparator* adjacent_sep) {

  std::set<sepDualOptVar*> vars = adjacent_sep->involved_vars();

  sepDualOptVar* other = 0;
  assert(vars.size() == 2);
  
  for (std::set<sepDualOptVar*>::iterator it = vars.begin(); it != vars.end(); it++) {

    if (*it != this)
      other = *it;
  }
  assert(other != 0);

  uint size = adjacent_separator_.size();
  adjacent_separator_.resize(size+1);
  adjacent_separator_[size] = std::make_pair(other,adjacent_sep);
}

double sepDualOptVar::dual_value(uint& arg_min) const {

  Math1D::Vector<double> sum(cost_.size());
  for (uint k=0; k < sum.size(); k++)
    sum[k] = cost_[k];

  for (uint k=0; k < adjacent_separator_.size(); k++)
    sum += adjacent_separator_[k].second->dual_var(this);
  
  for (uint k=0; k < adjacent_factor_.size(); k++) 
    sum += adjacent_factor_[k]->dual_var(this);

  double best = 1e300;
  for (uint k=0; k < sum.size(); k++) {

    if (sum[k] < best) {
      best = sum[k];
      arg_min = k;
    }
  }
  
  return best;
}

uint sepDualOptVar::nLabels() {
  return cost_.size();
}

void sepDualOptVar::compute_message(const sepDualOptFactor* factor, Math1D::Vector<double>& msg) {

  msg.resize_dirty(cost_.size());
  for (uint k=0; k < msg.size(); k++)
    msg[k] = cost_[k];

  for (uint s=0; s < adjacent_separator_.size(); s++)
    msg += adjacent_separator_[s].second->dual_var(this);

  for (uint f=0; f < adjacent_factor_.size(); f++) {
    if (adjacent_factor_[f] != factor)
      msg += adjacent_factor_[f]->dual_var(this);
  }
}

void sepDualOptVar::compute_message(const sepDualOptPairSeparator* sep, Math1D::Vector<double>& msg) {

  msg.resize_dirty(cost_.size());
  for (uint k=0; k < msg.size(); k++)
    msg[k] = cost_[k];

  for (uint s=0; s < adjacent_separator_.size(); s++) {
    if (adjacent_separator_[s].second != sep)
      msg += adjacent_separator_[s].second->dual_var(this);
  }

  for (uint f=0; f < adjacent_factor_.size(); f++) {
    msg += adjacent_factor_[f]->dual_var(this);
  }
}

const Math1D::Vector<float>& sepDualOptVar::cost() const {
  return cost_;
}

/**********************/

sepDualOptPairSeparator::sepDualOptPairSeparator(sepDualOptVar* var1, sepDualOptVar* var2) :
  var1_(var1), var2_(var2), dual_var_(2) {

  var1_->add_pair_separator(this);
  var2_->add_pair_separator(this);
 
  dual_var_[0].resize(var1_->nLabels(),0.0);
  dual_var_[1].resize(var2_->nLabels(),0.0);
}

/*virtual*/ sepDualOptPairSeparator::~sepDualOptPairSeparator() {}

const Math1D::Vector<double>& sepDualOptPairSeparator::dual_var(const sepDualOptVar* var) const {

  if (var == var1_)
    return dual_var_[0];
  else if (var == var2_)
    return dual_var_[1];
  else {
    assert(false);
    return dual_var_[0];
  }
}

/*virtual*/ void sepDualOptPairSeparator::update_duals(DualBCAMode /*mode*/) {

  const uint xDim = var1_->nLabels();
  const uint yDim = var2_->nLabels();

  Math2D::Matrix<double> fac_duals(xDim,yDim,0.0);

  for (uint f=0; f < adjacent_factor_.size(); f++)
    fac_duals += adjacent_factor_[f]->pair_dual(this);

  Math1D::Vector<double> msg;

  //handle var1
  var1_->compute_message(this,msg);

  for (uint x=0; x < xDim; x++) {
    
    double best = 1e300;

    for (uint y=0; y < yDim; y++) {

      double hyp = fac_duals(x,y) - dual_var_[1][y];

      if (hyp < best)
        best = hyp;
    }

    dual_var_[0][x] = 0.5 * (best - msg[x]);
  }

  //handle var2
  var2_->compute_message(this,msg);

  for (uint y=0; y < yDim; y++) {
  
    double best = 1e300;

    for (uint x=0; x < xDim; x++) {

      double hyp = fac_duals(x,y) - dual_var_[0][x];

      if (hyp < best)
        best = hyp;
    }

    dual_var_[1][y] = 0.5 * (best - msg[y]);
  }

}

void sepDualOptPairSeparator::add_factor(sepDualOptFactor* adjacent_fac) {
  uint size = adjacent_factor_.size();
  adjacent_factor_.resize(size+1);
  adjacent_factor_[size] = adjacent_fac;  
}

std::set<sepDualOptVar*> sepDualOptPairSeparator::involved_vars() {

  std::set<sepDualOptVar*> s;
  s.insert(var1_);
  s.insert(var2_);

  return s;
}

double sepDualOptPairSeparator::dual_value() const {

  uint nLabels1 = var1_->nLabels();
  uint nLabels2 = var2_->nLabels();

  assert(nLabels1 == dual_var_[0].size());
  assert(nLabels2 == dual_var_[1].size());

  Math2D::Matrix<double> msum(var1_->nLabels(),var2_->nLabels(),0.0);

  for (uint k=0; k < adjacent_factor_.size(); k++)
    msum += adjacent_factor_[k]->pair_dual(this);

  for (uint x=0; x < nLabels1; x++)
    for (uint y=0; y < nLabels2; y++)
      msum(x,y) -= dual_var_[0][x] + dual_var_[1][y];

  return msum.min();
}

void sepDualOptPairSeparator::compute_message(const sepDualOptFactor* factor, Math2D::Matrix<double>& msg) {

  msg.resize_dirty(var1_->nLabels(),var2_->nLabels());

  const uint xDim = var1_->nLabels();
  const uint yDim = var2_->nLabels();

  //this is separable - can we exploit that?
  for (uint x=0; x < xDim; x++)
    for (uint y=0; y < yDim; y++)
      msg(x,y) = - dual_var_[0][x] - dual_var_[1][y];

  const uint nFactors = adjacent_factor_.size();

  for (uint f=0; f < nFactors; f++) {

    if (adjacent_factor_[f] != factor)
      msg += adjacent_factor_[f]->pair_dual(this);
  }
}

sepDualOptVar* sepDualOptPairSeparator::var1() {
  return var1_;
}

sepDualOptVar* sepDualOptPairSeparator::var2() {
  return var2_;
}

/**********************/

sepDualOptFactor::sepDualOptFactor(const Storage1D<sepDualOptVar*>& vars, 
                                   const Storage1D<sepDualOptPairSeparator*>& separators, bool minimal_links) :
  var_(vars), separator_(separators), dual_var_(vars.size()), pair_dual_(separators.size()) {

  std::map<sepDualOptVar*,uint> var_num;

  for (uint v=0; v < vars.size(); v++) {
    var_num[vars[v]] = v;
  }

  std::set<sepDualOptVar*> sep_vars;
  for (uint k=0; k < separator_.size(); k++) {

    separator_[k]->add_factor(this);

    pair_dual_[k].resize(separator_[k]->var1()->nLabels(),separator_[k]->var2()->nLabels(),0.0);

    assert(var_num.find(separator_[k]->var2()) != var_num.end());
    assert(var_num.find(separator_[k]->var1()) != var_num.end());

    if (var_num[separator_[k]->var2()] < var_num[separator_[k]->var1()]) {
      std::cerr << "ERROR: wrong variable order" << std::endl;
      exit(1);
    }

    std::set<sepDualOptVar*> cur_set = separator_[k]->involved_vars();
    for (std::set<sepDualOptVar*>::iterator it = cur_set.begin(); it != cur_set.end(); it++) {
      sep_vars.insert(*it);
    }
  }

  for (uint v=0; v < vars.size(); v++) {
    if (!minimal_links || sep_vars.find(vars[v]) == sep_vars.end()) {
      vars[v]->add_factor(this);
      dual_var_[v].resize(vars[v]->nLabels(),0.0);
    }
  }
}

/*virtual*/ sepDualOptFactor::~sepDualOptFactor() {}

const Math1D::Vector<double>& sepDualOptFactor::dual_var(const sepDualOptVar* var) const {

  for (uint k=0; k < var_.size(); k++) {

    if (var_[k] == var)
      return dual_var_[k];
  }
  assert(false);
  return dual_var_[0];
}

const Math2D::Matrix<double>& sepDualOptFactor::pair_dual(const sepDualOptPairSeparator* pair_sep) const {

  for (uint k=0; k < separator_.size(); k++) {
    if (separator_[k] == pair_sep)
      return pair_dual_[k];
  }
  assert(false);
  return pair_dual_[0];
}

const Storage1D<sepDualOptVar*>& sepDualOptFactor::vars() const {
  return var_;
}

const Storage1D<sepDualOptPairSeparator*>& sepDualOptFactor::separators() const {
  return separator_;
}

/**********************/

BinarySepDualOptFactor::BinarySepDualOptFactor(const Storage1D<sepDualOptVar*>& vars, 
                                               const Math2D::Matrix<float>& cost, bool minimal_links) :
  sepDualOptFactor(vars,Storage1D<sepDualOptPairSeparator*>(),minimal_links), cost_(cost) {

  if (vars.size() != 2) {
    INTERNAL_ERROR << " attempt to instantiate binary factor with " << vars.size() << " variables. Exiting." << std::endl;
    exit(1);
  }

  if (cost.xDim() < vars[0]->nLabels() || cost.yDim() < vars[1]->nLabels()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
    exit(1);
  }
}

/*virtual*/ BinarySepDualOptFactor::~BinarySepDualOptFactor() {}

/*virtual*/ void BinarySepDualOptFactor::update_duals(DualBCAMode /*mode*/) {

  Math1D::Vector<double> msg;

  const uint nLabels1 = cost_.xDim();
  const uint nLabels2 = cost_.yDim();


  // handle var1
  var_[0]->compute_message(this,msg);

  for (uint x=0; x < nLabels1; x++) {

    double best = 1e300;

    for (uint y=0; y < nLabels2; y++) {

      double hyp = cost_(x,y) - dual_var_[1][y];

      if (hyp < best)
        best = hyp;
    }

    dual_var_[0][x] = 0.5 * (best - msg[x]);
  }
  
  // handle var2
  var_[1]->compute_message(this,msg);

  for (uint y=0; y < nLabels2; y++) {

    double best = 1e300;

    for (uint x=0; x < nLabels1; x++) {

      double hyp = cost_(x,y) - dual_var_[0][x];

      if (hyp < best)
        best = hyp;
    }

    dual_var_[1][y] = 0.5 * (best - msg[y]);
  }

}

/*virtual*/ void BinarySepDualOptFactor::write_cost(std::ostream& out, double factor) {

  for (uint x=0; x < cost_.xDim(); x++) 
    for (uint y=0; y < cost_.yDim(); y++) 
      out << (factor*cost_(x,y)) << " ";
}


/*virtual */
double BinarySepDualOptFactor::dual_value() const {

  uint nLabels1 = cost_.xDim();
  uint nLabels2 = cost_.yDim();

  double best = 1e300;

  for (uint y=0; y < nLabels2; y++) {
    for (uint x=0; x < nLabels1; x++) {

      double hyp = cost_(x,y) - dual_var_[0][x] - dual_var_[1][y];
      if (hyp < best)
        best = hyp;
    }
  }

  return best;
}

/**********************/

TernarySepDualOptFactor::TernarySepDualOptFactor(const Storage1D<sepDualOptVar*>& vars, 
                                                 const Storage1D<sepDualOptPairSeparator*>& separators,
                                                 const Math3D::Tensor<float>& cost, bool minimal_links) :
  sepDualOptFactor(vars,separators,minimal_links), cost_(cost) {

  if (vars.size() != 3) {
    INTERNAL_ERROR << " attempt to instantiate ternary factor with " << vars.size() << " variables. Exiting." << std::endl;
    exit(1);
  }

  if (cost.xDim() < vars[0]->nLabels() || cost.yDim() < vars[1]->nLabels() || cost.zDim() < vars[2]->nLabels()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
    exit(1);
  }
}

/*virtual*/ TernarySepDualOptFactor::~TernarySepDualOptFactor() {}

double TernarySepDualOptFactor::eval_pair(uint pair_num, uint x, uint y, uint z) const {

  sepDualOptVar* v1 = separator_[pair_num]->var1();
  sepDualOptVar* v2 = separator_[pair_num]->var2();

  uint a=MAX_UINT;
  uint b=MAX_UINT;

  if (var_[0] == v1)
    a = x;
  else {
    assert(var_[1] == v1);
    a = y;
  }

  if (var_[1] == v2)
    b = y;
  else 
    b = z;

  return pair_dual_[pair_num](a,b);  
}

/*virtual*/ void TernarySepDualOptFactor::update_duals(DualBCAMode /*mode*/) {

  const uint xDim = cost_.xDim();
  const uint yDim = cost_.yDim();
  const uint zDim = cost_.zDim();

  const uint nSeps = separator_.size();

  //1.) update pairwise separators
#if 1
  Math2D::Matrix<double> pair_msg;

  for (uint s=0; s < nSeps; s++) {
    separator_[s]->compute_message(this,pair_msg);

    if (separator_[s]->var1() == var_[0]) {

      if (separator_[s]->var2() == var_[1]) {

        for (uint y=0; y < yDim; y++) {
          for (uint x=0; x < xDim; x++) {

            double best = 1e300;

            for (uint z=0; z < zDim; z++) {

              double hyp = cost_(x,y,z);

              if (dual_var_[0].size() > 0)
                hyp -= dual_var_[0][x];
              if (dual_var_[1].size() > 0)
                hyp -= dual_var_[1][y];
              if (dual_var_[2].size() > 0)
                hyp -= dual_var_[2][z];
	      
              for (uint ss=0; ss < nSeps; ss++) {
	      
                if (ss != s) {

                  assert(separator_[ss]->var2() == var_[2]);

                  if (separator_[ss]->var1() == var_[0]) {
                    hyp -= pair_dual_[ss](x,z);
                  }
                  else {
                    assert(separator_[ss]->var1() == var_[1]);
                    hyp -= pair_dual_[ss](y,z);
                  }
                }
              }

              if (hyp < best)
                best = hyp;
            }

            pair_dual_[s](x,y) = 0.5 * (best - pair_msg(x,y));
          }
        }
      }
      else {

        for (uint z=0; z < zDim; z++) {
          for (uint x=0; x < xDim; x++) {

            double best = 1e300;

            for (uint y=0; y < yDim; y++) {

              double hyp = cost_(x,y,z);

              if (dual_var_[0].size() > 0)
                hyp -= dual_var_[0][x];
              if (dual_var_[1].size() > 0)
                hyp -= dual_var_[1][y];
              if (dual_var_[2].size() > 0)
                hyp -= dual_var_[2][z];

              for (uint ss=0; ss < nSeps; ss++) {
	      
                if (ss != s) {
		  
                  if (separator_[ss]->var1() == var_[0]) {
                    assert(separator_[ss]->var2() == var_[1]);
                    hyp -= pair_dual_[ss](x,y);
                  }
                  else {
                    assert(separator_[ss]->var1() == var_[1]);
                    assert(separator_[ss]->var2() == var_[2]);
                    hyp -= pair_dual_[ss](y,z);
                  }
                }
              }

              if (hyp < best)
                best = hyp;
            }
	    
            pair_dual_[s](x,z) = 0.5 * (best - pair_msg(x,z));
          }
        }

      }
    }
    else {

      for (uint z=0; z < zDim; z++) {
        for (uint y=0; y < yDim; y++) {

          double best = 1e300;

          for (uint x=0; x < xDim; x++) {

            double hyp = cost_(x,y,z);

            if (dual_var_[0].size() > 0)
              hyp -= dual_var_[0][x];
            if (dual_var_[1].size() > 0)
              hyp -= dual_var_[1][y];
            if (dual_var_[2].size() > 0)
              hyp -= dual_var_[2][z];
	    
            for (uint ss=0; ss < nSeps; ss++) {
	    
              if (ss != s) {

                assert(separator_[ss]->var1() != var_[1] || separator_[ss]->var2() != var_[2]);
                assert(separator_[ss]->var1() == var_[0]);

                if (separator_[ss]->var2() == var_[1]) {
                  hyp -= pair_dual_[ss](x,y);
                }
                else {
                  hyp -= pair_dual_[ss](x,z);
                }
              }
            }

            if (hyp < best)
              best = hyp;
          }

          pair_dual_[s](y,z) = 0.5 * (best - pair_msg(y,z));
        }
      }

    }
  }
#endif

  //2.) update single separators
  Math1D::Vector<double> msg;

  //update var1
  if (dual_var_[0].size() > 0) {
    var_[0]->compute_message(this,msg);

    for (uint x=0; x < xDim; x++) {

      double best = 1e300;

      for (uint z=0; z < zDim; z++) {
        for (uint y=0; y < yDim; y++) {

          double hyp = cost_(x,y,z);
	  
          if (dual_var_[1].size() > 0)
            hyp -= dual_var_[1][y];
          if (dual_var_[2].size() > 0)
            hyp -= dual_var_[2][z];
	  
          for (uint s=0; s < nSeps; s++) {

            if (separator_[s]->var1() == var_[0]) {

              if (separator_[s]->var2() == var_[1]) {
                hyp -= pair_dual_[s](x,y);
              }
              else {
                hyp -= pair_dual_[s](x,z);
              }
            }
            else
              hyp -= pair_dual_[s](y,z);
          }
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      dual_var_[0][x] = 0.5 * (best - msg[x]);
    }
  }

  //update var2
  if (dual_var_[1].size() > 0) {
    var_[1]->compute_message(this,msg);

    for (uint y=0; y < yDim; y++) {

      double best = 1e300;

      for (uint z=0; z < zDim; z++) {
        for (uint x=0; x < xDim; x++) {

          double hyp = cost_(x,y,z);

          if (dual_var_[0].size() > 0)
            hyp -= dual_var_[0][x];
          if (dual_var_[2].size() > 0)
            hyp -= dual_var_[2][z];
	  
          for (uint s=0; s < nSeps; s++) {

            if (separator_[s]->var1() == var_[0]) {

              if (separator_[s]->var2() == var_[1]) {
                hyp -= pair_dual_[s](x,y);
              }
              else {
                hyp -= pair_dual_[s](x,z);
              }
            }
            else
              hyp -= pair_dual_[s](y,z);
          }
	  
          if (hyp < best)
            best = hyp;
        }
      }

      dual_var_[1][y] = 0.5 * (best - msg[y]);
    }
  }

  //update var3
  if (dual_var_[2].size() > 0) {

    var_[2]->compute_message(this,msg);

    for (uint z=0; z < zDim; z++) {

      double best = 1e300;

      for (uint y=0; y < yDim; y++) {
        for (uint x=0; x < xDim; x++) {
	
          double hyp = cost_(x,y,z);

          if (dual_var_[0].size() > 0)
            hyp -= dual_var_[0][x];
          if (dual_var_[1].size() > 0)
            hyp -= dual_var_[1][y];
	  
          for (uint s=0; s < nSeps; s++) {

            if (separator_[s]->var1() == var_[0]) {

              if (separator_[s]->var2() == var_[1]) {
                hyp -= pair_dual_[s](x,y);
              }
              else {
                hyp -= pair_dual_[s](x,z);
              }
            }
            else
              hyp -= pair_dual_[s](y,z);
          }
	  
          if (hyp < best)
            best = hyp;	  
        }
      }

      dual_var_[2][z] = 0.5 * (best - msg[z]);
    }
  }

}

/*virtual*/ void TernarySepDualOptFactor::write_cost(std::ostream& out, double factor) {

  for (uint x=0; x < cost_.xDim(); x++) 
    for (uint y=0; y < cost_.yDim(); y++) 
      for (uint z=0; z < cost_.zDim(); z++) 
        out << (factor*cost_(x,y,z)) << " ";
}

/*virtual*/ double TernarySepDualOptFactor::dual_value() const {

  const uint xDim = cost_.xDim();
  const uint yDim = cost_.yDim();
  const uint zDim = cost_.zDim();

  Math3D::Tensor<double> sum(cost_.xDim(),cost_.yDim(),cost_.zDim());
  for (uint k=0; k < sum.size(); k++)
    sum.direct_access(k) = cost_.direct_access(k);

  if (dual_var_[0].size() + dual_var_[1].size() + dual_var_[2].size() >= 1) {

    for (uint z=0; z < zDim; z++) {
      for (uint y=0; y < yDim; y++) {
        for (uint x=0; x < xDim; x++) {
	  
          if (dual_var_[0].size() > 0)
            sum(x,y,z) -= dual_var_[0][x];
          if (dual_var_[1].size() > 0)
            sum(x,y,z) -= dual_var_[1][y];
          if (dual_var_[2].size() > 0)
            sum(x,y,z) -= dual_var_[2][z];
        }
      }
    }
  }

  //now add pairwise separators
  for (uint s=0; s < separator_.size(); s++) {

    if (separator_[s]->var1() == var_[0]) {

      if (separator_[s]->var2() == var_[1]) {

        for (uint z=0; z < zDim; z++) 
          for (uint y=0; y < yDim; y++) 
            for (uint x=0; x < xDim; x++) 
              sum(x,y,z) -= pair_dual_[s](x,y);
      }
      else {
	
        for (uint z=0; z < zDim; z++) 
          for (uint y=0; y < yDim; y++) 
            for (uint x=0; x < xDim; x++) 
              sum(x,y,z) -= pair_dual_[s](x,z);
      }
    }
    else {

      for (uint z=0; z < zDim; z++) 
        for (uint y=0; y < yDim; y++) 
          for (uint x=0; x < xDim; x++) 
            sum(x,y,z) -= pair_dual_[s](y,z);
    }
  }

  return sum.min();
}

/**********************/

FourthOrderSepDualOptFactor::FourthOrderSepDualOptFactor(const Storage1D<sepDualOptVar*>& vars, 
                                                         const Storage1D<sepDualOptPairSeparator*>& separators,
                                                         const Storage1D<Math3D::Tensor<float> >& cost, bool minimal_links)
  :   sepDualOptFactor(vars,separators,minimal_links), cost_(cost) {

  if (vars.size() != 4) {
    INTERNAL_ERROR << " attempt to instantiate ternary factor with " << vars.size() << " variables. Exiting." << std::endl;
    exit(1);
  }

  if (cost.size() < vars[0]->nLabels() || cost[0].xDim() < vars[1]->nLabels() || cost[0].yDim() < vars[2]->nLabels()
      || cost[0].zDim() < vars[3]->nLabels()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
    exit(1);
  }
}

/*virtual*/ FourthOrderSepDualOptFactor::~FourthOrderSepDualOptFactor() {}

double FourthOrderSepDualOptFactor::eval_pair(uint pair_num, 
                                              uint x, uint y, uint z, uint w) const {

  
  sepDualOptVar* v1 = separator_[pair_num]->var1();
  sepDualOptVar* v2 = separator_[pair_num]->var2();

  uint a=MAX_UINT;
  uint b=MAX_UINT;

  if (var_[0] == v1)
    a = x;
  else if (var_[1] == v1)
    a = y;
  else {
    assert (var_[2] == v1);
    a = z;
  }

  if (var_[1] == v2)
    b = y;
  else if (var_[2] == v2)
    b = z;
  else
    b = w;

  return pair_dual_[pair_num](a,b);
}


/*virtual*/ void FourthOrderSepDualOptFactor::update_duals(DualBCAMode /*mode*/) {

  const uint xDim = cost_.size();
  const uint yDim = cost_[0].xDim();
  const uint zDim = cost_[0].yDim();
  const uint wDim = cost_[0].zDim();

  const uint nSeps = separator_.size();

  //1.) update pair separators
  Math2D::Matrix<double> pair_msg;

  for (uint s=0; s < nSeps; s++) {
    separator_[s]->compute_message(this,pair_msg);

    if (separator_[s]->var1() == var_[0]) {

      if (separator_[s]->var2() == var_[1]) {

        for (uint x=0; x < xDim; x++) {

          const Math3D::Tensor<float>& cur_cost = cost_[x];

          for (uint y=0; y < yDim; y++) {

            double best = 1e300;

            for (uint z=0; z < zDim; z++) {
	      
              for (uint w=0; w < wDim; w++) {

                double hyp = cur_cost(y,z,w);
		
                if (dual_var_[0].size() > 0)
                  hyp -= dual_var_[0][x];
                if (dual_var_[1].size() > 0)
                  hyp -= dual_var_[1][y];
                if (dual_var_[2].size() > 0)
                  hyp -= dual_var_[2][z];
                if (dual_var_[3].size() > 0)
                  hyp -= dual_var_[3][w];

	      
                for (uint ss=0; ss < nSeps; ss++) {
		  
                  if (ss != s) {
                    hyp -= eval_pair(ss,x,y,z,w);
                  }
                }
		
                if (hyp < best)
                  best = hyp;		
              }
            }
            pair_dual_[s](x,y) = 0.5 * (best - pair_msg(x,y));
          }
        }
      }
      else if (separator_[s]->var2() == var_[2]) {

        for (uint x=0; x < xDim; x++) {

          const Math3D::Tensor<float>& cur_cost = cost_[x];

          for (uint z=0; z < zDim; z++) {

            double best = 1e300;

            for (uint y=0; y < yDim; y++) {	    

              for (uint w=0; w < wDim; w++) {

                double hyp = cur_cost(y,z,w);
		
                if (dual_var_[0].size() > 0)
                  hyp -= dual_var_[0][x];
                if (dual_var_[1].size() > 0)
                  hyp -= dual_var_[1][y];
                if (dual_var_[2].size() > 0)
                  hyp -= dual_var_[2][z];
                if (dual_var_[3].size() > 0)
                  hyp -= dual_var_[3][w];

	      
                for (uint ss=0; ss < nSeps; ss++) {
		  
                  if (ss != s) {
                    hyp -= eval_pair(ss,x,y,z,w);
                  }
                }
		
                if (hyp < best)
                  best = hyp;		
              }
            }

            pair_dual_[s](x,z) = 0.5 * (best - pair_msg(x,z));
          }
        }

      }
      else {

        for (uint x=0; x < xDim; x++) {

          const Math3D::Tensor<float>& cur_cost = cost_[x];

          for (uint w=0; w < wDim; w++) {

            double best = 1e300;
	    
            for (uint y=0; y < yDim; y++) {
              for (uint z=0; z < zDim; z++) {

                double hyp = cur_cost(y,z,w);
		
                if (dual_var_[0].size() > 0)
                  hyp -= dual_var_[0][x];
                if (dual_var_[1].size() > 0)
                  hyp -= dual_var_[1][y];
                if (dual_var_[2].size() > 0)
                  hyp -= dual_var_[2][z];
                if (dual_var_[3].size() > 0)
                  hyp -= dual_var_[3][w];

	      
                for (uint ss=0; ss < nSeps; ss++) {
		  
                  if (ss != s) {
                    hyp -= eval_pair(ss,x,y,z,w);
                  }
                }
		
                if (hyp < best)
                  best = hyp;
              }
            }

            pair_dual_[s](x,w) = 0.5 * (best - pair_msg(x,w));
          }
        }
      }
    }
    else if (separator_[s]->var1() == var_[1]) { 

      if (separator_[s]->var2() == var_[2]) {      

        for (uint y=0; y < yDim; y++) {
          for (uint z=0; z < zDim; z++) {

            double best = 1e300;

            for (uint x=0; x < xDim; x++) {

              const Math3D::Tensor<float>& cur_cost = cost_[x];
          
              for (uint w=0; w < wDim; w++) {
	    
                double hyp = cur_cost(y,z,w);
		
                if (dual_var_[0].size() > 0)
                  hyp -= dual_var_[0][x];
                if (dual_var_[1].size() > 0)
                  hyp -= dual_var_[1][y];
                if (dual_var_[2].size() > 0)
                  hyp -= dual_var_[2][z];
                if (dual_var_[3].size() > 0)
                  hyp -= dual_var_[3][w];

	      
                for (uint ss=0; ss < nSeps; ss++) {
		  
                  if (ss != s) {
                    hyp -= eval_pair(ss,x,y,z,w);
                  }
                }
		
                if (hyp < best)
                  best = hyp;

              }
            }

            pair_dual_[s](y,z) = 0.5 * (best - pair_msg(y,z));
          }
        }
      }
      else  {
     
        for (uint y=0; y < yDim; y++) {
          for (uint w=0; w < wDim; w++) {

            double best = 1e300;

            for (uint x=0; x < xDim; x++) {

              const Math3D::Tensor<float>& cur_cost = cost_[x];

              for (uint z=0; z < zDim; z++) {

                double hyp = cur_cost(y,z,w);
		
                if (dual_var_[0].size() > 0)
                  hyp -= dual_var_[0][x];
                if (dual_var_[1].size() > 0)
                  hyp -= dual_var_[1][y];
                if (dual_var_[2].size() > 0)
                  hyp -= dual_var_[2][z];
                if (dual_var_[3].size() > 0)
                  hyp -= dual_var_[3][w];

	      
                for (uint ss=0; ss < nSeps; ss++) {
		  
                  if (ss != s) {
                    hyp -= eval_pair(ss,x,y,z,w);
                  }
                }
		
                if (hyp < best)
                  best = hyp;
              }
            }

            pair_dual_[s](y,w) = 0.5 * (best - pair_msg(y,w));
          }
        }
 
      }
    }
    else {

      for (uint z=0; z < zDim; z++) {
        for (uint w=0; w < wDim; w++) {

          double best = 1e300;

          for (uint x=0; x < xDim; x++) {

            const Math3D::Tensor<float>& cur_cost = cost_[x];

            for (uint y=0; y < yDim; y++) {

              double hyp = cur_cost(y,z,w);
		
              if (dual_var_[0].size() > 0)
                hyp -= dual_var_[0][x];
              if (dual_var_[1].size() > 0)
                hyp -= dual_var_[1][y];
              if (dual_var_[2].size() > 0)
                hyp -= dual_var_[2][z];
              if (dual_var_[3].size() > 0)
                hyp -= dual_var_[3][w];
	      
	      
              for (uint ss=0; ss < nSeps; ss++) {
		
                if (ss != s) {
                  hyp -= eval_pair(ss,x,y,z,w);
                }
              }
	      
              if (hyp < best)
                best = hyp;
            }
          }	  

          pair_dual_[s](z,w) = 0.5 * (best - pair_msg(z,w));
        }
      }
    }

  }

  //2.) update single separators
  Math1D::Vector<double> msg;

  //update var1
  if (dual_var_[0].size() > 0) {
    var_[0]->compute_message(this,msg);
    
    for (uint x=0; x < xDim; x++) {

      double best = 1e300;

      const Math3D::Tensor<float>& cur_cost = cost_[x];

      for (uint w=0; w < wDim; w++) {
        for (uint z=0; z < zDim; z++) {
          for (uint y=0; y < yDim; y++) {
	    
            double hyp = cur_cost(y,z,w);
	    
            if (dual_var_[1].size() > 0)
              hyp -= dual_var_[1][y];
            if (dual_var_[2].size() > 0)
              hyp -= dual_var_[2][z];
            if (dual_var_[3].size() > 0)
              hyp -= dual_var_[3][w];
	    
	  
            for (uint s=0; s < nSeps; s++) {
	      
              if (separator_[s]->var1() == var_[0]) {
		
                if (separator_[s]->var2() == var_[1]) {
                  hyp -= pair_dual_[s](x,y);
                }
                else if (separator_[s]->var2() == var_[2]) {
                  hyp -= pair_dual_[s](x,z);
                }
                else {
                  hyp -= pair_dual_[s](x,w);
                }
              }
              else if (separator_[s]->var1() == var_[1]) {
                if (separator_[s]->var2() == var_[2]) 
                  hyp -= pair_dual_[s](y,z);
                else
                  hyp -= pair_dual_[s](y,w);
              }
              else {
                hyp -= pair_dual_[s](z,w);
              }
            }
	    
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      dual_var_[0][x] = 0.5 * (best - msg[x]);
    }
  }

  //update var2
  if (dual_var_[1].size() > 0) {
    var_[1]->compute_message(this,msg);

    for (uint y=0; y < yDim; y++) {

      double best = 1e300;

      for (uint x=0; x < xDim; x++) {

        const Math3D::Tensor<float>& cur_cost = cost_[x];

        for (uint w=0; w < wDim; w++) {
          for (uint z=0; z < zDim; z++) {

            double hyp = cur_cost(y,z,w);
	    
            if (dual_var_[0].size() > 0)
              hyp -= dual_var_[0][x];
            if (dual_var_[2].size() > 0)
              hyp -= dual_var_[2][z];
            if (dual_var_[3].size() > 0)
              hyp -= dual_var_[3][w];
	    
	  
            for (uint s=0; s < nSeps; s++) {
	      
              if (separator_[s]->var1() == var_[0]) {
		
                if (separator_[s]->var2() == var_[1]) {
                  hyp -= pair_dual_[s](x,y);
                }
                else if (separator_[s]->var2() == var_[2]) {
                  hyp -= pair_dual_[s](x,z);
                }
                else {
                  hyp -= pair_dual_[s](x,w);
                }
              }
              else if (separator_[s]->var1() == var_[1]) {
                if (separator_[s]->var2() == var_[2]) 
                  hyp -= pair_dual_[s](y,z);
                else
                  hyp -= pair_dual_[s](y,w);
              }
              else {
                hyp -= pair_dual_[s](z,w);
              }
            }
	    
            if (hyp < best)
              best = hyp;
          }
        }
      }

      dual_var_[1][y] = 0.5 * (best - msg[y]);
    }
  }

  //update var3
  if (dual_var_[2].size() > 0) {
    var_[2]->compute_message(this,msg);

    for (uint z=0; z < zDim; z++) {

      double best = 1e300;

      for (uint x=0; x < xDim; x++) {

        const Math3D::Tensor<float>& cur_cost = cost_[x];

        for (uint w=0; w < wDim; w++) {
          for (uint y=0; y < yDim; y++) {

            double hyp = cur_cost(y,z,w);

            if (dual_var_[0].size() > 0)
              hyp -= dual_var_[0][x];
            if (dual_var_[1].size() > 0)
              hyp -= dual_var_[1][y];
            if (dual_var_[3].size() > 0)
              hyp -= dual_var_[3][w];
	    	  
            for (uint s=0; s < nSeps; s++) {
	      
              if (separator_[s]->var1() == var_[0]) {
		
                if (separator_[s]->var2() == var_[1]) {
                  hyp -= pair_dual_[s](x,y);
                }
                else if (separator_[s]->var2() == var_[2]) {
                  hyp -= pair_dual_[s](x,z);
                }
                else {
                  hyp -= pair_dual_[s](x,w);
                }
              }
              else if (separator_[s]->var1() == var_[1]) {
                if (separator_[s]->var2() == var_[2]) 
                  hyp -= pair_dual_[s](y,z);
                else
                  hyp -= pair_dual_[s](y,w);
              }
              else {
                hyp -= pair_dual_[s](z,w);
              }
            }
	    
            if (hyp < best)
              best = hyp;

          }
        }
      }

      dual_var_[2][z] = 0.5 * (best - msg[z]);
    }
  }


  //update var4
  if (dual_var_[3].size() > 0) {
    var_[3]->compute_message(this,msg);

    for (uint w=0; w < wDim; w++) {

      double best = 1e300;

      for (uint x=0; x < xDim; x++) {

        const Math3D::Tensor<float>& cur_cost = cost_[x];

        for (uint z=0; z < zDim; z++) {
          for (uint y=0; y < yDim; y++) {

            double hyp = cur_cost(y,z,w);

            if (dual_var_[0].size() > 0)
              hyp -= dual_var_[0][x];
            if (dual_var_[1].size() > 0)
              hyp -= dual_var_[1][y];
            if (dual_var_[2].size() > 0)
              hyp -= dual_var_[2][z];
	  
            for (uint s=0; s < nSeps; s++) {
	      
              if (separator_[s]->var1() == var_[0]) {
		
                if (separator_[s]->var2() == var_[1]) {
                  hyp -= pair_dual_[s](x,y);
                }
                else if (separator_[s]->var2() == var_[2]) {
                  hyp -= pair_dual_[s](x,z);
                }
                else {
                  hyp -= pair_dual_[s](x,w);
                }
              }
              else if (separator_[s]->var1() == var_[1]) {
                if (separator_[s]->var2() == var_[2]) 
                  hyp -= pair_dual_[s](y,z);
                else
                  hyp -= pair_dual_[s](y,w);
              }
              else {
                hyp -= pair_dual_[s](z,w);
              }
            }
	    
            if (hyp < best)
              best = hyp;
	    
          }
        }
      }

      dual_var_[3][w] = 0.5 * (best - msg[w]);
    }
  }

}

/*virtual*/ void FourthOrderSepDualOptFactor::write_cost(std::ostream& out, double factor) {

  for (uint x=0; x < cost_.size(); x++) 
    for (uint y=0; y < cost_[x].xDim(); y++) 
      for (uint z=0; z < cost_[x].yDim(); z++) 
        for (uint w=0; w < cost_[x].zDim(); w++) 
          out << (factor*cost_[x](y,z,w)) << " ";
}


/*virtual*/ double FourthOrderSepDualOptFactor::dual_value() const {

  const uint xDim = cost_.size();
  const uint yDim = cost_[0].xDim();
  const uint zDim = cost_[0].yDim();
  const uint wDim = cost_[0].zDim();

  Storage1D<Math3D::Tensor<double> > sum(cost_.size());

  for (uint x=0; x < xDim; x++)  {

    sum[x].resize(cost_[x].xDim(), cost_[x].yDim(),cost_[x].zDim());
    for (uint k=0; k < sum[x].size(); k++)
      sum[x].direct_access(k) = cost_[x].direct_access(k);
  }

  if (dual_var_[0].size() + dual_var_[1].size() + dual_var_[2].size() + dual_var_[3].size() >= 1) {

    for (uint x=0; x < xDim; x++) {

      for (uint w=0; w < wDim; w++) {
        for (uint z=0; z < zDim; z++) {
          for (uint y=0; y < yDim; y++) {
	  
            if (dual_var_[0].size() > 0)
              sum[x](y,z,w) -= dual_var_[0][x];
            if (dual_var_[1].size() > 0)
              sum[x](y,z,w) -= dual_var_[1][y];
            if (dual_var_[2].size() > 0)
              sum[x](y,z,w) -= dual_var_[2][z];
            if (dual_var_[3].size() > 0)
              sum[x](y,z,w) -= dual_var_[3][w];
          }
        }
      }
    }
  }

  //now add pairwise separators
  for (uint s=0; s < separator_.size(); s++) {

    //std::cerr << "s: " << s << std::endl;

    sepDualOptVar* v1 = separator_[s]->var1();
    sepDualOptVar* v2 = separator_[s]->var2();

    if (v1 == var_[0]) {

      if (v2 == var_[1]) {

        for (uint x=0; x < xDim; x++) 
          for (uint w=0; w < wDim; w++) 
            for (uint z=0; z < zDim; z++) 
              for (uint y=0; y < yDim; y++) 
                sum[x](y,z,w) -= pair_dual_[s](x,y);
      }
      else if (v2 == var_[2]) {

        for (uint x=0; x < xDim; x++) 
          for (uint w=0; w < wDim; w++) 
            for (uint z=0; z < zDim; z++) 
              for (uint y=0; y < yDim; y++) 
                sum[x](y,z,w) -= pair_dual_[s](x,z);
      }
      else {

        for (uint x=0; x < xDim; x++) 
          for (uint w=0; w < wDim; w++) 
            for (uint z=0; z < zDim; z++) 
              for (uint y=0; y < yDim; y++) 
                sum[x](y,z,w) -= pair_dual_[s](x,w);
      }
    }
    else if (v1 == var_[1]) {

      if (v2 == var_[2]) {

        for (uint x=0; x < xDim; x++) 
          for (uint w=0; w < wDim; w++) 
            for (uint z=0; z < zDim; z++) 
              for (uint y=0; y < yDim; y++) 
                sum[x](y,z,w) -= pair_dual_[s](y,z);
      }
      else {

        for (uint x=0; x < xDim; x++) 
          for (uint w=0; w < wDim; w++) 
            for (uint z=0; z < zDim; z++) 
              for (uint y=0; y < yDim; y++) 
                sum[x](y,z,w) -= pair_dual_[s](y,w);
      }
    }
    else {

      for (uint x=0; x < xDim; x++) 
        for (uint w=0; w < wDim; w++) 
          for (uint z=0; z < zDim; z++) 
            for (uint y=0; y < yDim; y++) 
              sum[x](y,z,w) -= pair_dual_[s](z,w);
    }
  }
  
  double best = 1e300;

  for (uint x=0; x < xDim; x++) {
    for (uint w=0; w < wDim; w++) {
      for (uint z=0; z < zDim; z++) {
        for (uint y=0; y < yDim; y++) {

          best = std::min(best,sum[x](y,z,w));
        }
      }
    }
  }

  return best;
}



/**********************/

SeparatorDualOptimization::SeparatorDualOptimization(uint nVars, uint nSeparators, uint nFactors, bool minimal_links) :
  nUsedVars_(0), nUsedSeparators_(0), nUsedFactors_(0), minimal_links_(minimal_links) {

  var_.resize(nVars,0);
  separator_.resize(nSeparators,0);
  factor_.resize(nFactors,0);
}

SeparatorDualOptimization::~SeparatorDualOptimization() {

  for (uint k=0; k < nUsedVars_; k++)
    delete var_[k];
  for (uint k=0; k < nUsedSeparators_; k++)
    delete separator_[k];
  for (uint k=0; k < nUsedFactors_; k++)
    delete factor_[k];
}

uint SeparatorDualOptimization::add_var(const Math1D::Vector<float>& cost) {

  if (nUsedVars_ == var_.size())
    var_.resize(uint(1.2*nUsedVars_)+4);

  assert(nUsedVars_ < var_.size());

  var_[nUsedVars_] = new sepDualOptVar(cost);

  nUsedVars_++;
  return nUsedVars_-1;
}

void SeparatorDualOptimization::add_factor(sepDualOptFactor* fac) {

  if (nUsedFactors_ == factor_.size())
    factor_.resize(uint(1.2*nUsedFactors_)+4);

  factor_[nUsedFactors_] = fac;
  nUsedFactors_++;
}

uint SeparatorDualOptimization::add_separator(uint v1, uint v2) {

  if (v1 >= nUsedVars_ || v2 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  if (nUsedSeparators_ == separator_.size())
    separator_.resize(uint(1.2*nUsedSeparators_)+4);

  assert(v1 < v2);

  assert(nUsedSeparators_ < separator_.size());

  separator_[nUsedSeparators_] = new sepDualOptPairSeparator(var_[v1],var_[v2]);
  
  nUsedSeparators_++;
  return nUsedSeparators_-1;
}

void SeparatorDualOptimization::add_generic_binary_factor(uint v1, uint v2, const Math2D::Matrix<float>& cost) {

  if (v1 >= nUsedVars_ || v2 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }


  Storage1D<sepDualOptVar*> var(2);
  var[0] = var_[v1];
  var[1] = var_[v2];

  add_factor(new BinarySepDualOptFactor(var,cost,minimal_links_));
}

void SeparatorDualOptimization::add_generic_ternary_factor(uint v1, uint v2, uint v3, const Storage1D<uint>& separators,
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

  Storage1D<sepDualOptVar*> var(3);
  var[0] = var_[v1];
  var[1] = var_[v2];
  var[2] = var_[v3];

  Storage1D<sepDualOptPairSeparator*> sep(separators.size());
  for (uint k=0; k < separators.size(); k++) {
    if (separators[k] >= nUsedSeparators_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    sep[k] = separator_[separators[k]];
  }

  add_factor(new TernarySepDualOptFactor(var,sep,cost_copy,minimal_links_));
}

void SeparatorDualOptimization::add_fourth_order_factor(uint v1, uint v2, uint v3, uint v4,
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

  Storage1D<sepDualOptVar*> var(4);
  var[0] = var_[v1];
  var[1] = var_[v2];
  var[2] = var_[v3];
  var[3] = var_[v4];

  Storage1D<sepDualOptPairSeparator*> sep(separators.size());
  for (uint k=0; k < separators.size(); k++) {
    if (separators[k] >= nUsedSeparators_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    sep[k] = separator_[separators[k]];
  }

  add_factor(new FourthOrderSepDualOptFactor(var,sep,cost_copy,minimal_links_));
}

sepDualOptVar* SeparatorDualOptimization::get_variable(uint v) {

  if (v < nUsedVars_)
    return var_[v];
  return 0;
}

const Math1D::Vector<uint>& SeparatorDualOptimization::labeling() {
  return labeling_;
}


void SeparatorDualOptimization::save_problem() {

  std::map<sepDualOptVar*,uint> var_idx;
  for (uint v=0; v < nUsedVars_; v++) {
    var_idx[var_[v]] = v;
  }

  std::map<sepDualOptPairSeparator*,uint> sep_idx;
  for (uint s=0; s < nUsedSeparators_; s++) {
    sep_idx[separator_[s]] = s;
  }

  std::ofstream of("regions.txt");

  for (uint f=0; f < nUsedFactors_; f++) {

    for (uint v=0; v < factor_[f]->vars().size(); v++)
      of << (var_idx[factor_[f]->vars()[v]] + 1) << " ";
    of << std::endl;
  }

  //unaries
  for (uint v=0; v < nUsedVars_; v++) 
    of << (v+1) << std::endl;
  of.close();

  of.open("intersects.txt");
  for (uint v=0; v < nUsedVars_; v++) {
    of << (v+1) << std::endl;
  }
  // for (uint s=0; s < nUsedSeparators_; s++) {
  //   of << (var_idx[separator_[s]->var1()] + 1) << " "
  //      << (var_idx[separator_[s]->var2()] + 1) << std::endl;
  // }
  of.close();
  
  of.open("region_intersects.txt");
  for (uint f=0; f < nUsedFactors_; f++) {
    
    for (uint v=0; v < factor_[f]->vars().size(); v++)
      of << (var_idx[factor_[f]->vars()[v]] + 1) << " ";
    
    // for (uint s=0; s < factor_[f]->separators().size(); s++) {
    //   of << (sep_idx[factor_[f]->separators()[s]] + nUsedVars_ + 1) << " "; 
    // }
    of << std::endl;
  }
  for (uint v=0; v < nUsedVars_; v++) 
    of << (v+1) << std::endl;

  of.close();

  of.open("lambdas.txt");
  for (uint f=0; f < nUsedFactors_; f++) {
    factor_[f]->write_cost(of,-1.0);
    of << std::endl;
  }
  for (uint v=0; v < nUsedVars_; v++) {
    Math1D::Vector<float> cur_cost = var_[v]->cost();
    cur_cost *= -1.0;
    of << cur_cost << std::endl;
  }
  of.close();

  of.open("var_sizes.txt");
  for (uint v=0; v < nUsedVars_; v++)
    of << var_[v]->nLabels() << std::endl;
  of.close();
}

double SeparatorDualOptimization::optimize(uint nIter, DualBCAMode mode, bool quiet) {

  double bound = -1e300;

  labeling_.resize(nUsedVars_,0);

  std::cerr.precision(8);

  uint label;

  //DEBUG
  bound = 0.0;
  for (uint v=0; v < nUsedVars_; v++) {
    bound += var_[v]->dual_value(label);
    labeling_[v] = label;
  }
  
  for (uint s=0; s < nUsedSeparators_; s++) {
    bound += separator_[s]->dual_value();
  }
  
  for (uint f=0; f < nUsedFactors_; f++) {
    bound += factor_[f]->dual_value();
  }
  std::cerr << "first bound: " << bound << std::endl;
  //END_DEBUG

  size_t effort_per_iteration = 0;
    
  for (uint f=0; f < nUsedFactors_; f++) {
    
    uint var_size = factor_[f]->vars().size();
    uint sep_size = factor_[f]->separators().size();

    effort_per_iteration += (var_size + sep_size) * (var_size + sep_size);
  }

  for (uint iter = 1; iter <= nIter; iter++) {

    std::cerr << "******** iteration " << iter << " ***********" << std::endl;

    if (!minimal_links_) {
      for (uint s=0; s < nUsedSeparators_; s++) {
	separator_[s]->update_duals(mode);
      }
    }

#if 1
    for (uint f=0; f < nUsedFactors_; f++) {
      factor_[f]->update_duals(mode);
    }
#endif

    /*** compute bound ***/
    if (!quiet) {
      bound = 0.0;
      
      for (uint v=0; v < nUsedVars_; v++) {
	bound += var_[v]->dual_value(label);
	labeling_[v] = label;
      }
      
      for (uint s=0; s < nUsedSeparators_; s++)
	bound += separator_[s]->dual_value();
      
      for (uint f=0; f < nUsedFactors_; f++)
	bound += factor_[f]->dual_value();
      
      std::cerr << "lower bound: " << bound << std::endl;
    }
  }

  bound = 0.0;

  for (uint v=0; v < nUsedVars_; v++) {
    bound += var_[v]->dual_value(label);
    labeling_[v] = label;
  }

  for (uint s=0; s < nUsedSeparators_; s++)
    bound += separator_[s]->dual_value();
      
  for (uint f=0; f < nUsedFactors_; f++)
    bound += factor_[f]->dual_value();

  if (!minimal_links_) {
    size_t message_effort = 0;
    
    for (uint f=0; f < nUsedFactors_; f++) {

      uint var_size = factor_[f]->vars().size();
      uint sep_size = factor_[f]->separators().size();

      message_effort += (var_size + sep_size) * (var_size + sep_size);
    }
    message_effort *= nIter;

    std::cerr << "message effort: " << message_effort << std::endl;
  }

  return bound;
}

