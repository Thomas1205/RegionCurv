/***** written by Thomas Schoenemann as an employee of the University of Pisa, Italy, December 2011 *****/
/****        and as an employee of the University of Düsseldorf, Germany, January 2012 *******/


#include "sepTRWS.hh"

#include <map>

//DEBUG
//#include "stl_out.hh"
//END_DEBUG

AllInclusiveSepCumTRWSSepChainLink::AllInclusiveSepCumTRWSSepChainLink() : sep_(0), left_fac_(0), right_fac_(0) {}

AllInclusiveSepCumTRWSSepChainLink::AllInclusiveSepCumTRWSSepChainLink(AllInclusiveSepCumTRWSFactor* lfac, AllInclusiveSepCumTRWSFactor* rfac,
								       AllInclusiveSepCumTRWSPairSeparator* sep) :
  sep_(sep), left_fac_(lfac), right_fac_(rfac) {}

/***************/

AllInclusiveSepCumTRWSVarChainLink::AllInclusiveSepCumTRWSVarChainLink() {}

/***************/

AllInclusiveSepCumTRWSVariable::AllInclusiveSepCumTRWSVariable(const Math1D::Vector<float>& cost, uint rank) :
  cost_(cost), rank_(rank) {

  cum_cost_.resize(cost_.size());
  for (uint l=0; l < cost_.size(); l++)
    cum_cost_[l] = cost_[l];
}

//will delete chain links
AllInclusiveSepCumTRWSVariable::~AllInclusiveSepCumTRWSVariable() {

  for (uint k=0; k < chain_link_.size(); k++) {
    
    delete chain_link_[k];
  }
}
  
void AllInclusiveSepCumTRWSVariable::add_factor(AllInclusiveSepCumTRWSFactor* adjacent_fac) {

  uint nFac = adjacent_factor_.size();
  adjacent_factor_.resize(nFac+1);
  adjacent_factor_[nFac] = adjacent_fac;
}
  
void AllInclusiveSepCumTRWSVariable::add_pair_separator(AllInclusiveSepCumTRWSPairSeparator* adjacent_sep) {

  uint nSep = adjacent_separator_.size();
  adjacent_separator_.resize(nSep+1);
  adjacent_separator_[nSep] = adjacent_sep;

  //bubble sort to get the correct order
  for (uint k=0; k < nSep; k++) {
    for (uint l=0; l < nSep-1; l++) {
      AllInclusiveSepCumTRWSPairSeparator* sep1 = adjacent_separator_[l];
      AllInclusiveSepCumTRWSPairSeparator* sep2 = adjacent_separator_[l+1];

      if (sep1->var2() == sep2->var2()) {
        if (sep1->var1()->rank() > sep2->var1()->rank())
          std::swap(adjacent_separator_[l],adjacent_separator_[l+1]);
      }
      else {
        if (sep1->var2()->rank() > sep2->var2()->rank())
          std::swap(adjacent_separator_[l],adjacent_separator_[l+1]);
      }
    }
  }
}

uint AllInclusiveSepCumTRWSVariable::nLabels() const {
  return cost_.size();
}

uint AllInclusiveSepCumTRWSVariable::rank() const {
  return rank_;
}

const Storage1D<AllInclusiveSepCumTRWSFactor*>& AllInclusiveSepCumTRWSVariable::adjacent_factors() const {
  return adjacent_factor_;
}

const Storage1D<AllInclusiveSepCumTRWSPairSeparator*>& AllInclusiveSepCumTRWSVariable::adjacent_separators() const {
  return adjacent_separator_;
}

double AllInclusiveSepCumTRWSVariable::reparameterize_forward() {

  double offs = 0.0;

  for (uint c = 0; c < chain_link_.size(); c++) {

    AllInclusiveSepCumTRWSVarChainLink* cur_chain_link = chain_link_[c];

    AllInclusiveSepCumTRWSFactor* fac = cur_chain_link->involved_factor_[0];
    
    if (cur_chain_link->involved_factor_.size() > 1) {
      
      // we use the pair separator as the end of the first factor in a chain link.
      //  if the variable is the higher ranked in the separator, it comes after the separator,
      //  so we move to the second factor
      
      //NOTE: this requires that start- and end-separators are set to 0 
      //   if the chain link is a single variable
      
      uint f=0;
      while (fac->end_separator() != 0 && 
	     fac->end_separator()->var2() == this) {
	f++;
	if (f >= cur_chain_link->involved_factor_.size())
	  break;
	
	fac = cur_chain_link->involved_factor_[f];
      }
      
      //DEBUG
      // if ( !(f < chain_link_[c]->involved_factor_.size())) {
      //   std::cerr << "var " << rank_ << std::endl;
      //   std::cerr << "involved factors: " << std::endl;
      //   for (uint k=0; k < chain_link_[c]->involved_factor_.size(); k++)
      //     chain_link_[c]->involved_factor_[k]->print_factor();
      // }
      //END_DEBUG
      
      assert(f < chain_link_[c]->involved_factor_.size());
    }

    double temp = fac->compute_var_reparameterization(this);


    if (fac->max_rank() == rank_ && fac->end_separator() == 0) {

      //CAREFUL here: this depends heavily on whether pair-seps are precomputed
      // (we need to make sure that a chain/factor end is not counted twice)
      offs += temp;
    }
  }

  return offs;
}

double AllInclusiveSepCumTRWSVariable::reparameterize_backward() {

  double offs = 0.0;

  for (uint c = 0; c < chain_link_.size(); c++) {

    AllInclusiveSepCumTRWSVarChainLink* cur_chain_link = chain_link_[c];

    if (true) {

      uint nInvolved = cur_chain_link->involved_factor_.size();
      //assert(nInvolved <= 2);
      AllInclusiveSepCumTRWSFactor* fac = cur_chain_link->involved_factor_[nInvolved-1];

      if (cur_chain_link->involved_factor_.size() > 1) {

        //NOTE: this requires that start- and end-separators are set to 0 
        //   if the chain link is a single variable

        int f=nInvolved-1;
        while (fac->start_separator() != 0 && 
               fac->start_separator()->var1() == this) {
          f--;
          if (f < 0)
            break;

          fac = cur_chain_link->involved_factor_[f];
        }
        
        assert(f >= 0);
      }

      double temp = fac->compute_var_reparameterization(this);

      if (fac->min_rank() == rank_ && fac->start_separator() == 0) {

        //CAREFUL here: this depends heavily on whether pair-seps are precomputed
        // (we need to make sure that a chain/factor end is not counted twice)
        offs += temp;
      }
    }
  }

  return offs;
}

double AllInclusiveSepCumTRWSVariable::average(uint& arg_min) {

  const uint nLabels = cost_.size();

  for (uint l=0; l < nLabels; l++)
    cum_cost_[l] = cost_[l];

  //NOTE: I think we could also loop over the adjacent factors here
  for (uint c = 0; c < chain_link_.size(); c++) {

    for (uint f=0; f < chain_link_[c]->involved_factor_.size(); f++) {
      cum_cost_ += chain_link_[c]->involved_factor_[f]->var_reparameterization(this);
    }
  }

  double offs = 1e300;
  for (uint l=0; l < nLabels; l++) {
    if (cum_cost_[l] < offs) {

      offs = cum_cost_[l];
      arg_min = l;
    }
  }


  for (uint l=0; l < nLabels; l++)
    cum_cost_[l] = (cum_cost_[l] - offs) / chain_link_.size();

  return offs;
}


const Math1D::Vector<double>& AllInclusiveSepCumTRWSVariable::cost() const {
  return cum_cost_;
}

void AllInclusiveSepCumTRWSVariable::set_up_chains() {

  //check if the assumptions of the current implementation are satisfied
  for (uint f=0; f < adjacent_factor_.size(); f++) {
    
    AllInclusiveSepCumTRWSFactor* cur_fac = adjacent_factor_[f];
    AllInclusiveSepCumTRWSPairSeparator* start_sep = cur_fac->start_separator();
    AllInclusiveSepCumTRWSPairSeparator* end_sep = cur_fac->end_separator();
    
    if (start_sep == 0 || end_sep == 0)
      continue;
    
    bool in_start_sep = (this == start_sep->var1() || this == start_sep->var2());
    bool in_end_sep = (this == end_sep->var1() || this == end_sep->var2());
    
    if (in_start_sep && in_end_sep && (this != start_sep->var2() || this != end_sep->var1())) {
      
      INTERNAL_ERROR << " if a variable is part of both the start and end separator of a factor, " << std::endl
                     << "    it currently must be the end of the start sepator and the start of the end separator. Exiting.." 
                     << std::endl;

      exit(1);
    }
    
    if (start_sep->var2()->rank() > end_sep->var1()->rank() ) {
      
      INTERNAL_ERROR << " start and end separator cannot be processed sequentially. Exiting" << std::endl;
      exit(1);
    } 
  }
  
  std::set<AllInclusiveSepCumTRWSFactor*> processed;
  
  for (uint f=0; f < adjacent_factor_.size(); f++) {

    if (processed.find(adjacent_factor_[f]) != processed.end())
      continue;

    processed.insert(adjacent_factor_[f]);

    AllInclusiveSepCumTRWSVarChainLink* link = new AllInclusiveSepCumTRWSVarChainLink;
    
    link->involved_factor_.push_back(adjacent_factor_[f]);

    bool can_form_var_end_link = (rank_ == adjacent_factor_[f]->max_rank() && adjacent_factor_[f]->next_factor() == 0);
    bool can_form_var_start_link = (rank_ == adjacent_factor_[f]->min_rank() && adjacent_factor_[f]->prev_factor() == 0);

    bool is_var_end_link = false;
    bool is_var_start_link = false;

    if (can_form_var_end_link) {

      //try to find a following factor
      uint ff;
      for (ff=f+1; ff < adjacent_factor_.size(); ff++) {

        if (processed.find(adjacent_factor_[ff]) != processed.end())
          continue;

        if (adjacent_factor_[ff]->min_rank() == rank_ && adjacent_factor_[ff]->prev_factor() == 0) {
          assert(adjacent_factor_[ff]->start_separator() == 0 || adjacent_factor_[ff]->start_separator()->var1() == this);
          break;
        }
      }

      if (ff < adjacent_factor_.size()) {

        adjacent_factor_[f]->set_next_factor(adjacent_factor_[ff]);
        adjacent_factor_[f]->set_end_separator(0);
        adjacent_factor_[ff]->set_prev_factor(adjacent_factor_[f]);
        adjacent_factor_[ff]->set_start_separator(0);
        processed.insert(adjacent_factor_[ff]);
        is_var_end_link = true;

        //this assertion may fail if not all possibilities for pair-separator chains were exploited
        if (adjacent_factor_[f]->end_separator() == adjacent_factor_[ff]->start_separator())
          assert(adjacent_factor_[f]->end_separator() == 0);
      }
    }

    if (can_form_var_start_link) {

      //try to find a previous factor
      uint ff;
      for (ff=f+1; ff < adjacent_factor_.size(); ff++) {

        if (processed.find(adjacent_factor_[ff]) != processed.end())
          continue;

        if (adjacent_factor_[ff]->max_rank() == rank_ && adjacent_factor_[ff]->next_factor() == 0) {
          assert(adjacent_factor_[ff]->start_separator() == 0 || adjacent_factor_[ff]->start_separator()->var1() == this);
          break;
        }
      }

      if (ff < adjacent_factor_.size()) {

        adjacent_factor_[f]->set_prev_factor(adjacent_factor_[ff]);
        adjacent_factor_[f]->set_start_separator(0);
        adjacent_factor_[ff]->set_next_factor(adjacent_factor_[f]);
        adjacent_factor_[ff]->set_end_separator(0);
        processed.insert(adjacent_factor_[ff]);        

        is_var_start_link = true;

        //this assertion may fail if not all possibilities for pair-separator chains were exploited
        if (adjacent_factor_[f]->start_separator() == adjacent_factor_[ff]->end_separator())
          assert(adjacent_factor_[f]->start_separator() == 0);
      }
    }

    if (is_var_start_link) {

      AllInclusiveSepCumTRWSFactor* prev_fac = adjacent_factor_[f]->prev_factor();
      if (prev_fac->end_separator() != 0) {
        link->sep_.insert(link->sep_.begin(),prev_fac->end_separator()); 
      }
      if (prev_fac->start_separator() == 0 || prev_fac->start_separator()->var2() != this)
        link->involved_factor_.insert(link->involved_factor_.begin(),prev_fac);
    }
    
    AllInclusiveSepCumTRWSFactor* cur_fac = adjacent_factor_[f];

    while (true) {
      
      AllInclusiveSepCumTRWSPairSeparator* start_sep = cur_fac->start_separator();
      if (start_sep == 0 || (start_sep->var1() != this && start_sep->var2() != this))
        break;
      
      link->sep_.insert(link->sep_.begin(),start_sep);
      
      if (cur_fac->prev_factor() == 0)
        break;
      
      cur_fac = cur_fac->prev_factor();
      if (std::find(adjacent_factor_.direct_access(),adjacent_factor_.direct_access()+adjacent_factor_.size(),cur_fac)
          == adjacent_factor_.direct_access() + adjacent_factor_.size())
        break;

      processed.insert(cur_fac);
      link->involved_factor_.insert(link->involved_factor_.begin(),cur_fac);
    }

    cur_fac = adjacent_factor_[f];

    if (is_var_end_link) {

      AllInclusiveSepCumTRWSFactor* next_fac = adjacent_factor_[f]->next_factor();
      if (next_fac->start_separator() != 0) {
        link->sep_.push_back(next_fac->start_separator());
      }
      if (next_fac->end_separator() == 0 || next_fac->end_separator()->var1() != this)
        link->involved_factor_.push_back(next_fac);
    }

    while (true) {
      
      AllInclusiveSepCumTRWSPairSeparator* end_sep = cur_fac->end_separator();
      if (end_sep == 0 || (end_sep->var1() != this && end_sep->var2() != this)) {
        break;
      }
      link->sep_.push_back(end_sep);
      
      if (cur_fac->next_factor() == 0)
        break;
      
      cur_fac = cur_fac->next_factor();

      if (std::find(adjacent_factor_.direct_access(),adjacent_factor_.direct_access()+adjacent_factor_.size(),cur_fac)
          == adjacent_factor_.direct_access() + adjacent_factor_.size())
        break;
      
      processed.insert(cur_fac);
      link->involved_factor_.push_back(cur_fac);
    }
      
    assert(link->involved_factor_.size() <= 3);
    assert(link->sep_.size() <= 2);

    uint size = chain_link_.size();
    chain_link_.resize(size+1);
    chain_link_[size] = link;
  }

  for (uint l=0; l < cost_.size(); l++)
    cum_cost_[l] = cost_[l] / chain_link_.size();

}

/************************/

AllInclusiveSepCumTRWSPairSeparator::AllInclusiveSepCumTRWSPairSeparator(AllInclusiveSepCumTRWSVariable* var1, AllInclusiveSepCumTRWSVariable* var2)
  :  sep_rank_(MAX_UINT), var1_(var1), var2_(var2) {

  var1_->add_pair_separator(this);
  var2_->add_pair_separator(this);

  pair_parameters_.resize(var1_->nLabels(),var2_->nLabels(),0.0);
}

//will delete chain links
AllInclusiveSepCumTRWSPairSeparator::~AllInclusiveSepCumTRWSPairSeparator() {

  for (uint k=0; k < chain_link_.size(); k++)
    delete chain_link_[k];
}

void AllInclusiveSepCumTRWSPairSeparator::add_factor(AllInclusiveSepCumTRWSFactor* adjacent_fac) {

  uint nFac = adjacent_factor_.size();
  adjacent_factor_.resize(nFac+1);
  adjacent_factor_[nFac] = adjacent_fac;
}

AllInclusiveSepCumTRWSVariable* AllInclusiveSepCumTRWSPairSeparator::var1() const {
  return var1_;
}

AllInclusiveSepCumTRWSVariable* AllInclusiveSepCumTRWSPairSeparator::var2() const {
  return var2_;
}

const Storage1D<AllInclusiveSepCumTRWSFactor*>& AllInclusiveSepCumTRWSPairSeparator::adjacent_factors() const {
  return adjacent_factor_;
}

double AllInclusiveSepCumTRWSPairSeparator::reparameterize_forward() {

  double offs = 0.0;

  const uint nChains = chain_link_.size();

  for (uint c=0; c < nChains; c++) {

    AllInclusiveSepCumTRWSFactor* fac = chain_link_[c]->left_fac_;

    double temp;

    if (this != fac->valid_separator()) {
      temp = fac->compute_pair_reparameterization(this);
    }
    else
      temp = fac->valid_offs();

    if (fac->end_separator() == this) {
      offs += temp;
    }
  }

  return offs;
}

double AllInclusiveSepCumTRWSPairSeparator::reparameterize_backward() {

  double offs = 0.0;

  const uint nChains = chain_link_.size();

  for (uint c=0; c < nChains; c++) {

    AllInclusiveSepCumTRWSFactor* fac = (chain_link_[c]->right_fac_ != 0) ? chain_link_[c]->right_fac_ : chain_link_[c]->left_fac_;

    double temp;

    if (this != fac->valid_separator())
      temp = fac->compute_pair_reparameterization(this);
    else
      temp = fac->valid_offs();
    
    if (fac->start_separator() == this) {
      offs += temp;
    }
  }

  return offs;
}

double AllInclusiveSepCumTRWSPairSeparator::average() {

  pair_parameters_.set_constant(0.0);

  const uint nChains = chain_link_.size();

  for (uint c=0; c < nChains; c++) {
  
    pair_parameters_ += chain_link_[c]->left_fac_->pair_reparameterization(this);
    if (chain_link_[c]->right_fac_ != 0)
      pair_parameters_ += chain_link_[c]->right_fac_->pair_reparameterization(this);
  }

  double offs = pair_parameters_.min();

  double inv_nChains = 1.0 / chain_link_.size();

  const uint xDim = pair_parameters_.xDim();
  const uint yDim = pair_parameters_.yDim();

  for (uint l1=0; l1 < xDim; l1++)
    for (uint l2=0; l2 < yDim; l2++)
      pair_parameters_(l1,l2) = (pair_parameters_(l1,l2) - offs) * inv_nChains;

  return offs;
}

void AllInclusiveSepCumTRWSPairSeparator::set_up_chains() {

  std::vector<AllInclusiveSepCumTRWSFactor*> left_factor;
  std::vector<AllInclusiveSepCumTRWSFactor*> right_factor;

  for (uint f=0; f < adjacent_factor_.size(); f++) {

    AllInclusiveSepCumTRWSFactor* cur_fac = adjacent_factor_[f];

    if (cur_fac->start_separator() == this)
      right_factor.push_back(cur_fac);
    else if (cur_fac->end_separator() == this)
      left_factor.push_back(cur_fac);
    else {
      
      AllInclusiveSepCumTRWSSepChainLink* link = new AllInclusiveSepCumTRWSSepChainLink(cur_fac,0,this);
      uint size = chain_link_.size();
      chain_link_.resize(size+1);
      chain_link_[size] = link;
    }
  }

  uint nLinks = std::min(left_factor.size(),right_factor.size());

  for (uint l=0; l < nLinks; l++) {
    AllInclusiveSepCumTRWSSepChainLink* link = new AllInclusiveSepCumTRWSSepChainLink(left_factor[l],right_factor[l],this);

    left_factor[l]->set_next_factor(right_factor[l]);
    right_factor[l]->set_prev_factor(left_factor[l]);

    uint size = chain_link_.size();
    chain_link_.resize(size+1);
    chain_link_[size] = link;
  }

  for (uint l=nLinks; l < left_factor.size(); l++) {

    AllInclusiveSepCumTRWSSepChainLink* link = new AllInclusiveSepCumTRWSSepChainLink(left_factor[l],0,this);

    uint size = chain_link_.size();
    chain_link_.resize(size+1);
    chain_link_[size] = link;
  }

  for (uint l=nLinks; l < right_factor.size(); l++) {

    AllInclusiveSepCumTRWSSepChainLink* link = new AllInclusiveSepCumTRWSSepChainLink(right_factor[l],0,this);
    uint size = chain_link_.size();

    chain_link_.resize(size+1);
    chain_link_[size] = link;
  }

}

Math2D::Matrix<double>& AllInclusiveSepCumTRWSPairSeparator::pair_parameters() {
  return pair_parameters_;
}

const Math2D::Matrix<double>& AllInclusiveSepCumTRWSPairSeparator::pair_parameters() const {
  return pair_parameters_;
}

uint AllInclusiveSepCumTRWSPairSeparator::sep_rank() const {
  return sep_rank_;
}

void AllInclusiveSepCumTRWSPairSeparator::set_sep_rank(uint rank) {
  sep_rank_ = rank;
}

/********************/


AllInclusiveSepCumTRWSFactor::AllInclusiveSepCumTRWSFactor(const Storage1D<AllInclusiveSepCumTRWSVariable*>& vars, 
							   const Storage1D<AllInclusiveSepCumTRWSPairSeparator*>& separators)
  :  involved_var_(vars), adjacent_separator_(separators), 
     prev_factor_(0), next_factor_(0), valid_sep_(255), start_separator_(0), end_separator_(0) {

  assert(separators.size() <= 255);

  pair_reparameterization_.resize(adjacent_separator_.size());
  for (uint s=0; s < adjacent_separator_.size(); s++) {
    adjacent_separator_[s]->add_factor(this);
    pair_reparameterization_[s].resize(adjacent_separator_[s]->var1()->nLabels(),
                                       adjacent_separator_[s]->var2()->nLabels(),0.0);
  }

  var_reparameterization_.resize(involved_var_.size());
  for (uint v=0; v < involved_var_.size(); v++) {
    involved_var_[v]->add_factor(this);
    var_reparameterization_[v].resize(involved_var_[v]->nLabels(),0.0);
  }
}

/*virtual*/ AllInclusiveSepCumTRWSFactor::~AllInclusiveSepCumTRWSFactor() {}

/*virtual*/
 double AllInclusiveSepCumTRWSFactor::compute_forward(const AllInclusiveSepCumTRWSPairSeparator* /*incoming_sep*/, 
						      const AllInclusiveSepCumTRWSVariable* /*incoming_var*/, 
						      const AllInclusiveSepCumTRWSVariable* /*outgoing_var*/,
						      const Math2D::Matrix<double>& /*prev_pair_forward*/, 
						      const Math1D::Vector<double>& /*prev_var_forward*/, 
						      Math1D::Vector<double>& /*forward*/)
{

  //this routine is only needed for debugging and independent computations of the current bound
  // => we add an idle standard implementation
  assert(false);
  return 0.0;
}

/*virtual*/
double AllInclusiveSepCumTRWSFactor::compute_forward(const AllInclusiveSepCumTRWSPairSeparator* /*incoming_sep*/, 
						     const AllInclusiveSepCumTRWSVariable* /*incoming_var*/, 
						     const AllInclusiveSepCumTRWSPairSeparator* /*outgoing_sep*/,
						     const Math2D::Matrix<double>& /*prev_pair_forward*/, 
						     const Math1D::Vector<double>& /*prev_var_forward*/, 
						     Math2D::Matrix<double>& /*forward*/)
{

  //this routine is only needed for debugging and independent computations of the current bound
  // => we add an idle standard implementation
  assert(false);
  return 0.0;
}

void AllInclusiveSepCumTRWSFactor::compute_rank_range() {
  min_rank_ = MAX_UINT;
  max_rank_ = 0;

  for (uint v=0; v < involved_var_.size(); v++) {

    uint cur_rank = involved_var_[v]->rank();
    min_rank_ = std::min(min_rank_,cur_rank);
    max_rank_ = std::max(max_rank_,cur_rank);
  }

  uint min_sep_rank = MAX_UINT;
  uint max_sep_rank = 0;

  for (uint s=0; s < adjacent_separator_.size(); s++) {

    uint sep_rank = adjacent_separator_[s]->sep_rank();
    min_sep_rank = std::min(min_sep_rank,sep_rank);
    max_sep_rank = std::max(max_sep_rank,sep_rank);
  }

  //now set start and end separator
  for (uint s=0; s < adjacent_separator_.size(); s++) {

    uint sep_rank = adjacent_separator_[s]->sep_rank();
    if (sep_rank == min_sep_rank) {
      //note that this will fail in certain cases
      set_start_separator(adjacent_separator_[s]);
    }
    if (sep_rank == max_sep_rank) {
      //note that this will fail in certain cases
      set_end_separator(adjacent_separator_[s]);
    }
  }  
}

const Math1D::Vector<double>& AllInclusiveSepCumTRWSFactor::var_reparameterization(AllInclusiveSepCumTRWSVariable* var) const {

  for (uint v=0; v < involved_var_.size(); v++) {
    if (involved_var_[v] == var)
      return var_reparameterization_[v];
  }

  assert(false);
  return var_reparameterization_[0];
}
  
const Math2D::Matrix<double>& AllInclusiveSepCumTRWSFactor::pair_reparameterization(AllInclusiveSepCumTRWSPairSeparator* pair) const {

  for (uint s=0; s < adjacent_separator_.size(); s++) {
    
    if (adjacent_separator_[s] == pair)
      return pair_reparameterization_[s];
  }

  assert(false);
  return pair_reparameterization_[0];
}

uint AllInclusiveSepCumTRWSFactor::min_rank() const {
  return min_rank_;
}

uint AllInclusiveSepCumTRWSFactor::max_rank() const {
  return max_rank_;
}

AllInclusiveSepCumTRWSFactor* AllInclusiveSepCumTRWSFactor::prev_factor() const {
  return prev_factor_;
}

AllInclusiveSepCumTRWSFactor* AllInclusiveSepCumTRWSFactor::next_factor() const {
  return next_factor_;
}

double AllInclusiveSepCumTRWSFactor::valid_offs() const {
  return valid_offs_;
}

void AllInclusiveSepCumTRWSFactor::set_prev_factor(AllInclusiveSepCumTRWSFactor* prev_factor) {
  prev_factor_ = prev_factor;
}

void AllInclusiveSepCumTRWSFactor::set_next_factor(AllInclusiveSepCumTRWSFactor* next_factor) {
  next_factor_ = next_factor;
}

AllInclusiveSepCumTRWSPairSeparator* AllInclusiveSepCumTRWSFactor::start_separator() const {
  return start_separator_;
}

AllInclusiveSepCumTRWSPairSeparator* AllInclusiveSepCumTRWSFactor::end_separator() const {
  return end_separator_;
}

const Storage1D<AllInclusiveSepCumTRWSVariable*>& AllInclusiveSepCumTRWSFactor::involved_var() const {
  return involved_var_;
}

const Storage1D<AllInclusiveSepCumTRWSPairSeparator*>& AllInclusiveSepCumTRWSFactor::adjacent_separator() const {
  return adjacent_separator_;
}

//returns if successful
bool AllInclusiveSepCumTRWSFactor::set_start_separator(AllInclusiveSepCumTRWSPairSeparator* sep) {
  if (sep == 0) {
    start_separator_ = sep;
    return true;
  }

  for (uint v=0; v < involved_var_.size(); v++) {
    if (involved_var_[v]->rank() < sep->var2()->rank()
        && involved_var_[v] != sep->var1()) {
      start_separator_ = 0;
      return false;
    }
  }

  for (uint s=0; s < adjacent_separator_.size(); s++) {

    if (adjacent_separator_[s] == sep)
      break;
    else if (adjacent_separator_[s]->var1()->rank() < sep->var2()->rank()) {    

      start_separator_ = 0;
      return false;
    }
  }

  start_separator_ = sep;
  return true;
}

//returns if successful
bool AllInclusiveSepCumTRWSFactor::set_end_separator(AllInclusiveSepCumTRWSPairSeparator* sep) {

  if (sep == 0) {
    end_separator_ = sep;
    return true;
  }

  for (uint v=0; v < involved_var_.size(); v++) {

    if (involved_var_[v]->rank() > sep->var1()->rank()
        && involved_var_[v] != sep->var2()) {
      end_separator_ = 0;
      return false;
    }
  }
  for (int s= adjacent_separator_.size()-1; s >= 0; s--) {

    if (adjacent_separator_[s] == sep)
      break;
    else if (adjacent_separator_[s]->var2()->rank() > sep->var1()->rank()) {

      end_separator_ = 0;
      return false;
    }
  }

  end_separator_ = sep;
  return true;
}

bool AllInclusiveSepCumTRWSFactor::sep_is_interrupted(AllInclusiveSepCumTRWSPairSeparator* sep) const {

  uint rank1 = sep->var1()->rank();
  uint rank2 = sep->var2()->rank();

  //a) interrupted by a variable
  for (uint v=0; v < involved_var_.size(); v++) {

    uint cur_rank = involved_var_[v]->rank();
    if (cur_rank > rank1 && cur_rank < rank2)
      return true;
  }

  //b) interrupted by a pair-separator
  for (uint s=0; s < adjacent_separator_.size(); s++) {

    if (adjacent_separator_[s] != sep) {

      if (adjacent_separator_[s]->var2()->rank() == rank2)
        return true;
    }
  }
  
  return false;
}

void AllInclusiveSepCumTRWSFactor::print_factor(bool short_form) const {

  std::cerr << "factor " << this << " has variables ";
  for (uint v=0; v < involved_var_.size(); v++) {
    std::cerr << involved_var_[v]->rank() << " ";
  }
  std::cerr << std::endl << "     and separators ";
  for (uint s=0; s < adjacent_separator_.size(); s++) {
    std::cerr << adjacent_separator_[s] << " (" << adjacent_separator_[s]->var1()->rank() 
              << "," << adjacent_separator_[s]->var2()->rank() << ") ;";
  }
  std::cerr << std::endl;
  std::cerr << "     prev factor: " << prev_factor_ << ", next factor: " << next_factor_ << std::endl;  
  std::cerr << " start sep: " << start_separator_ << ", end sep: " << end_separator_ << std::endl;

  if (!short_form) {
    std::cerr << "var repars: " << std::endl;
    for (uint v=0; v < involved_var_.size(); v++) 
      std::cerr << "     " << var_reparameterization_[v] << std::endl;
    std::cerr << "var cost: " << std::endl;
    for (uint v=0; v < involved_var_.size(); v++) 
      std::cerr << "     " << involved_var_[v]->cost() << std::endl;
  

    std::cerr << "pair repars: " << std::endl;
    for (uint s=0; s < adjacent_separator_.size(); s++)
      std::cerr << pair_reparameterization_[s] << std::endl;
    std::cerr << "pair cost: " << std::endl;
    for (uint s=0; s < adjacent_separator_.size(); s++)
      std::cerr << adjacent_separator_[s]->pair_parameters() << std::endl;
  }
}

const AllInclusiveSepCumTRWSPairSeparator* AllInclusiveSepCumTRWSFactor::valid_separator() const {

  if (valid_sep_ < adjacent_separator_.size())
    return adjacent_separator_[valid_sep_];
  else
    return 0;
}

/***********************/


BinaryAllInclusiveSepCumTRWSFactor::BinaryAllInclusiveSepCumTRWSFactor(const Storage1D<AllInclusiveSepCumTRWSVariable*>& vars, 
                                                                       const Math2D::Matrix<float>& cost) :
  AllInclusiveSepCumTRWSFactor(vars,Storage1D<AllInclusiveSepCumTRWSPairSeparator*>()), cost_(cost) {

  if (vars.size() != 2) {
    INTERNAL_ERROR << " attempt to instantiate binary factor with " << vars.size() << " variables. Exiting." << std::endl;
    exit(1);
  }

  if (cost.xDim() < vars[0]->nLabels() || cost.yDim() < vars[1]->nLabels()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
    exit(1);
  }
}

/*virtual*/ BinaryAllInclusiveSepCumTRWSFactor::~BinaryAllInclusiveSepCumTRWSFactor() {}
  
/*virtual*/ 
double BinaryAllInclusiveSepCumTRWSFactor::compute_var_reparameterization(AllInclusiveSepCumTRWSVariable* var) {

  double offs = 0.0;

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();

  Storage1D<Math1D::Vector<double> > param = var_reparameterization_;
  for (uint v=0; v < 2; v++)
    param[v] -= involved_var_[v]->cost();

  //this routine also updates reparameterization_
  uint idx = 0;

  if (var == involved_var_[0]) {
    
    idx = 0;

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;

      for (uint l2 = 0; l2 < nLabels2; l2++) {

        const double w2 = param[1][l2];
	
        double hyp = cost_(l1,l2) - w2;
	  
        if (hyp < best)
          best = hyp;
      }
      
      var_reparameterization_[0][l1] = best;
    }
  }
  else {
    assert(var == involved_var_[1]);

    idx = 1;

    for (uint l2 = 0; l2 < nLabels2; l2++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const double w1 = param[0][l1];
	
        double hyp = cost_(l1,l2) - w1;
	  
        if (hyp < best)
          best = hyp;
      }
      
      var_reparameterization_[1][l2] = best;
    }    
  }

  const double msg_offs = var_reparameterization_[idx].min();
  offs += msg_offs;

  for (uint l=0; l < var_reparameterization_[idx].size(); l++)
    var_reparameterization_[idx][l] -= msg_offs;

  valid_sep_ = 255;
  
  return offs;
}

/*virtual*/ 
double BinaryAllInclusiveSepCumTRWSFactor::compute_pair_reparameterization(AllInclusiveSepCumTRWSPairSeparator* /*pair*/) {

  //should never be called
  assert(false);
  return 0.0;
}

/*virtual*/ 
double BinaryAllInclusiveSepCumTRWSFactor::best_value() {

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();

  double min_val = 1e300;

  for (uint l1 = 0; l1 < nLabels1; l1++) {
    
    for (uint l2 = 0; l2 < nLabels2; l2++) {
	
      double hyp = cost_(l1,l2);
      hyp += involved_var_[0]->cost()[l1] - var_reparameterization_[0][l1];
      hyp += involved_var_[1]->cost()[l2] - var_reparameterization_[1][l2];
      
      if (hyp < min_val)
        min_val = hyp;
    }
  }

  return min_val;
}


/*virtual*/ 
double BinaryAllInclusiveSepCumTRWSFactor::compute_forward(const AllInclusiveSepCumTRWSPairSeparator* /*incoming_sep*/, 
                                                           const AllInclusiveSepCumTRWSVariable* incoming_var, 
                                                           const AllInclusiveSepCumTRWSVariable* outgoing_var,
                                                           const Math2D::Matrix<double>& /*prev_pair_forward*/, 
                                                           const Math1D::Vector<double>& prev_var_forward, 
                                                           Math1D::Vector<double>& forward) {

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  
  Storage1D<Math1D::Vector<double> > param = var_reparameterization_;
  for (uint v=0; v < involved_var_.size(); v++) {
    //do NOT include cost of variables included in the pair-separator
    if (involved_var_[v] != outgoing_var)
      param[v] -= involved_var_[v]->cost();
    if (involved_var_[v] == incoming_var)
      param[v] -= prev_var_forward;
  }

  if (outgoing_var == involved_var_[0]) {
    
    forward.resize(nLabels1);

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;

      for (uint l2 = 0; l2 < nLabels2; l2++) {

        const double w2 = param[1][l2];
	
        double hyp = cost_(l1,l2) - w2;

        if (hyp < best)
          best = hyp;
      }
      
      forward[l1] = best - param[0][l1];
    }
  }
  else {
    assert(outgoing_var == involved_var_[1]);

    forward.resize(nLabels2);

    for (uint l2 = 0; l2 < nLabels2; l2++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const double w1 = param[0][l1];
	
        double hyp = cost_(l1,l2) - w1;
        
        if (hyp < best)
          best = hyp;
      }
      
      forward[l2] = best  - param[1][l2];
    }    
  }

  const double offs = forward.min();
  for (uint k=0; k < forward.size(); k++)
    forward[k] -= offs;

  return offs;

}

/***********************/

TernaryAllInclusiveSepCumTRWSFactor::TernaryAllInclusiveSepCumTRWSFactor(const Storage1D<AllInclusiveSepCumTRWSVariable*>& vars, 
									 const Storage1D<AllInclusiveSepCumTRWSPairSeparator*>& separators,
									 const Math3D::Tensor<float>& cost) :
  AllInclusiveSepCumTRWSFactor(vars,separators), cost_(cost) {

  if (vars.size() != 3) {
    INTERNAL_ERROR << " attempt to instantiate ternary factor with " << vars.size() << " variables. Exiting." << std::endl;
    exit(1);
  }

  if (cost.xDim() < vars[0]->nLabels() || cost.yDim() < vars[1]->nLabels() || cost.zDim() < vars[2]->nLabels()) {
    INTERNAL_ERROR << "dimension mismatch. Exiting." << std::endl;
    exit(1);
  }
}

/*virtual*/ TernaryAllInclusiveSepCumTRWSFactor::~TernaryAllInclusiveSepCumTRWSFactor() {}
  
/*virtual*/ 
double TernaryAllInclusiveSepCumTRWSFactor::compute_var_reparameterization(AllInclusiveSepCumTRWSVariable* var) {

  double offs = 0.0;

  uint idx = 0;

  const AllInclusiveSepCumTRWSPairSeparator* valid_sep = this->valid_separator();

  if (valid_sep != 0 && (valid_sep->var1() == var || valid_sep->var2() == var)) {
    //if (false) {

    offs += valid_offs_;

    for (uint v=0; v < involved_var_.size(); v++) {
      if (involved_var_[v] == var) {
	idx = v;
	break;
      }
    }
    assert(involved_var_[idx] == var);

    Math1D::Vector<double>& repar = var_reparameterization_[idx];

    const Math2D::Matrix<double>& pair_cost = valid_sep->pair_parameters();
    const Math1D::Vector<double>& otherpar = (valid_sep->var1() == var) ? valid_sep->var2()->cost() : valid_sep->var1()->cost();

    if (valid_sep->var1() == var) {

      for (uint l1=0; l1 < var->nLabels(); l1++) {
	
	double best = 1e300;
	
	for (uint l2=0; l2 < otherpar.size(); l2++) {

	  const double hyp = pair_cost(l1,l2) + otherpar[l2];
	  
	  if (hyp < best)
	    best = hyp;
	}

	//repar was used for the computation of pair_cost -> cancel this
	repar[l1] = best + repar[l1];
      }
    }
    else {

      for (uint l2=0; l2 < var->nLabels(); l2++) {

	double best = 1e300;

	for (uint l1=0; l1 < otherpar.size(); l1++) {

	  const double hyp = pair_cost(l1,l2) + otherpar[l1];

	  if (hyp < best)
	    best = hyp;
	}

	//repar was used for the computation of pair_cost -> cancel this
	repar[l2] = best + repar[l2];
      }
    }
  }
  else {

    const uint nLabels1 = involved_var_[0]->nLabels();
    const uint nLabels2 = involved_var_[1]->nLabels();
    const uint nLabels3 = involved_var_[2]->nLabels();
    
    const uint nSeps = adjacent_separator_.size();

    //std::cerr << "before" << std::endl;

    Storage1D<Math1D::Vector<double> > param = var_reparameterization_;
    for (uint v=0; v < involved_var_.size(); v++) {
      param[v] -= involved_var_[v]->cost();
    }

    if (var == involved_var_[0]) {
    
      idx = 0;

      Math1D::Vector<double>& repar = var_reparameterization_[0];
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
	double best = 1e300;
	
	for (uint l2 = 0; l2 < nLabels2; l2++) {
	  
	  const double w2 = param[1][l2];
	
	  for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
	    double hyp = cost_(l1,l2,l3) - w2;
	    hyp -= param[2][l3];
	    
	    for (uint s=0; s < nSeps; s++)
	      hyp += eval_pair(s,l1,l2,l3);
	  
	    if (hyp < best)
	      best = hyp;
	  }
	}
      
	repar[l1] = best;
      }
    }
    else if (var == involved_var_[1]) {

      idx = 1;

      Math1D::Vector<double>& repar = var_reparameterization_[1];

      for (uint l2 = 0; l2 < nLabels2; l2++) {
      
	double best = 1e300;
	
	for (uint l1 = 0; l1 < nLabels1; l1++) {
	
	  const double w1 = param[0][l1];
	
	  for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
	    double hyp = cost_(l1,l2,l3) - w1;
	    hyp -= param[2][l3];

	    for (uint s=0; s < nSeps; s++)
	      hyp += eval_pair(s,l1,l2,l3);
	    
	    if (hyp < best)
	      best = hyp;
	  }
	}

	repar[l2] = best;
      }    
    }
    else {
      assert(var == involved_var_[2]);
    
      idx = 2;

      Math1D::Vector<double>& repar = var_reparameterization_[2];

      for (uint l3 = 0; l3 < nLabels3; l3++) {
	
	double best = 1e300;
	
	for (uint l1 = 0; l1 < nLabels1; l1++) {
	  
	  const double w1 = param[0][l1];
	
	  for (uint l2 = 0; l2 < nLabels2; l2++) {
	  
	    double hyp = cost_(l1,l2,l3) - w1;
	    hyp -= param[1][l2];
	    
	    for (uint s=0; s < nSeps; s++)
	      hyp += eval_pair(s,l1,l2,l3);
	    
	    if (hyp < best)
	      best = hyp;
	  }
	}
	
	repar[l3] = best;
      }
    }    
  }

  Math1D::Vector<double>& repar = var_reparameterization_[idx];

  const double msg_offs = repar.min();

  for (uint l=0; l < repar.size(); l++)
    repar[l] -= msg_offs;

  offs += msg_offs;

  valid_sep_ = 255;

  return offs;

}

/*virtual*/ 
double TernaryAllInclusiveSepCumTRWSFactor::compute_pair_reparameterization(AllInclusiveSepCumTRWSPairSeparator* pair) {

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  
  Storage1D<Math1D::Vector<double> > param = var_reparameterization_;
  for (uint v=0; v < involved_var_.size(); v++) {
    //do NOT include cost of variables included in the pair-separator
    if (pair->var1() != involved_var_[v] && pair->var2() != involved_var_[v])
      param[v] -= involved_var_[v]->cost();
  }

  const uint nSeps = adjacent_separator_.size();

  uint s=0;
  while(adjacent_separator_[s] != pair) {
    s++;
    assert (s < nSeps);
  }

  Math2D::Matrix<double>& repar = pair_reparameterization_[s];

  if (pair->var1() == involved_var_[0]) {

    if (pair->var2() == involved_var_[1]) {

      for (uint l1=0; l1 < nLabels1; l1++) {
        
        const double w1 = param[0][l1];

        for (uint l2=0; l2 < nLabels2; l2++) {

          double best = 1e300;
	  
          for (uint l3=0; l3 < nLabels3; l3++) {

            double hyp = cost_(l1,l2,l3);

            hyp -= param[2][l3]; 

            for (uint ss=0; ss < nSeps; ss++) {

              if (ss != s)
                hyp += eval_pair(ss,l1,l2,l3);
            }

            if (hyp < best)
              best = hyp;
          }

          repar(l1,l2) = best - w1 - param[1][l2];
        }
      }
    }
    else {

      for (uint l1=0; l1 < nLabels1; l1++) {

        const double w1 = param[0][l1];

        for (uint l3=0; l3 < nLabels3; l3++) {

          double best = 1e300;

          for (uint l2=0; l2 < nLabels2; l2++) {

            double hyp = cost_(l1,l2,l3);

            hyp -= param[1][l2]; 

            for (uint ss=0; ss < nSeps; ss++) {

              if (ss != s)
                hyp += eval_pair(ss,l1,l2,l3);
            }

            if (hyp < best)
              best = hyp;
          }

          repar(l1,l3) = best - w1 - param[2][l3];
        }
      }
    }
  }
  else {

    for (uint l2=0; l2 < nLabels2; l2++) {

      const double w2 = param[1][l2];

      for (uint l3=0; l3 < nLabels3; l3++) {

        double best = 1e300;

        for (uint l1=0; l1 < nLabels1; l1++) {

          double hyp = cost_(l1,l2,l3);
	  
          hyp -= param[0][l1]; 
	  
          for (uint ss=0; ss < nSeps; ss++) {
	    
            if (ss != s)
              hyp += eval_pair(ss,l1,l2,l3);
          }
	  
          if (hyp < best)
            best = hyp;
        }
        repar(l2,l3) = best - w2 - param[2][l3];
      }
    }
  }

  const double offs = repar.min();
  for (uint k=0; k < repar.size(); k++)
    repar.direct_access(k) -= offs;

  valid_sep_ = s;
  valid_offs_ = offs;

  return offs;
}

/*virtual*/ 
double TernaryAllInclusiveSepCumTRWSFactor::best_value() {

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();

  double min_val = 1e300;

  for (uint l1 = 0; l1 < nLabels1; l1++) {
    
    for (uint l2 = 0; l2 < nLabels2; l2++) {
	
      for (uint l3 = 0; l3 < nLabels3; l3++) {

        double hyp = cost_(l1,l2,l3);
        hyp += involved_var_[0]->cost()[l1] - var_reparameterization_[0][l1];
        hyp += involved_var_[1]->cost()[l2] - var_reparameterization_[1][l2];
        hyp += involved_var_[2]->cost()[l3] - var_reparameterization_[2][l3];

        for (uint s=0; s < adjacent_separator_.size(); s++)
          hyp += eval_pair(s,l1,l2,l3);

        if (hyp < min_val)
          min_val = hyp;
      }
    }
  }

  return min_val;
}


/*virtual*/ 
double TernaryAllInclusiveSepCumTRWSFactor::compute_forward(const AllInclusiveSepCumTRWSPairSeparator* incoming_sep, 
							    const AllInclusiveSepCumTRWSVariable* incoming_var, 
							    const AllInclusiveSepCumTRWSVariable* outgoing_var,
							    const Math2D::Matrix<double>& prev_pair_forward, 
							    const Math1D::Vector<double>& prev_var_forward, 
							    Math1D::Vector<double>& forward) {

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  
  Storage1D<Math1D::Vector<double> > param = var_reparameterization_;
  for (uint v=0; v < involved_var_.size(); v++) {
    //do NOT include cost of variables included in the pair-separator
    if (involved_var_[v] != outgoing_var)
      param[v] -= involved_var_[v]->cost();
    if (involved_var_[v] == incoming_var)
      param[v] -= prev_var_forward;
  }

  Storage1D<Math2D::Matrix<double> > pair_param = pair_reparameterization_;
  for (uint s=0; s < pair_param.size(); s++) {
    
    pair_param[s] -= adjacent_separator_[s]->pair_parameters();
    if (adjacent_separator_[s] == incoming_sep)
      pair_param[s] -= prev_pair_forward;
  }

  const uint nSeps = adjacent_separator_.size();

  if (outgoing_var == involved_var_[0]) {
    
    forward.resize(nLabels1);

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;

      for (uint l2 = 0; l2 < nLabels2; l2++) {

        const double w2 = param[1][l2];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
          double hyp = cost_(l1,l2,l3) - w2;
          hyp -= param[2][l3];

          for (uint s=0; s < nSeps; s++)
            hyp -= eval_pair(s,l1,l2,l3,pair_param);
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      forward[l1] = best - param[0][l1];
    }
  }
  else if (outgoing_var == involved_var_[1]) {

    forward.resize(nLabels2);

    for (uint l2 = 0; l2 < nLabels2; l2++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const double w1 = param[0][l1];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {
	  
          double hyp = cost_(l1,l2,l3) - w1;
          hyp -= param[2][l3];

          for (uint s=0; s < nSeps; s++)
            hyp -= eval_pair(s,l1,l2,l3,pair_param);
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      forward[l2] = best  - param[1][l2];
    }    
  }
  else {
    assert(outgoing_var == involved_var_[2]);

    forward.resize(nLabels3);

    for (uint l3 = 0; l3 < nLabels3; l3++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const double w1 = param[0][l1];
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {
	  
          double hyp = cost_(l1,l2,l3) - w1;
          hyp -= param[1][l2];

          for (uint s=0; s < nSeps; s++)
            hyp -= eval_pair(s,l1,l2,l3,pair_param);
	  
          if (hyp < best)
            best = hyp;
        }
      }
      
      forward[l3] = best  - param[2][l3];
    }
  }

  double offs = forward.min();
  for (uint k=0; k < forward.size(); k++)
    forward[k] -= offs;

  return offs;
}

/*virtual*/ 
double TernaryAllInclusiveSepCumTRWSFactor::compute_forward(const AllInclusiveSepCumTRWSPairSeparator* incoming_sep, 
							    const AllInclusiveSepCumTRWSVariable* incoming_var, 
							    const AllInclusiveSepCumTRWSPairSeparator* outgoing_sep,
							    const Math2D::Matrix<double>& prev_pair_forward, 
							    const Math1D::Vector<double>& prev_var_forward, 
							    Math2D::Matrix<double>& forward) {


  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  
  Storage1D<Math1D::Vector<double> > param = var_reparameterization_;
  for (uint v=0; v < involved_var_.size(); v++) {
    //do NOT include cost of variables included in the pair-separator
    if (outgoing_sep->var1() != involved_var_[v] && outgoing_sep->var2() != involved_var_[v])
      param[v] -= involved_var_[v]->cost();
    if (involved_var_[v] == incoming_var)
      param[v] -= prev_var_forward;
  }

  Storage1D<Math2D::Matrix<double> > pair_param = pair_reparameterization_;
  for (uint s=0; s < pair_param.size(); s++) {

    if (adjacent_separator_[s] != outgoing_sep)
      pair_param[s] -= adjacent_separator_[s]->pair_parameters();
    if (adjacent_separator_[s] == incoming_sep)
      pair_param[s] -= prev_pair_forward;
  }

  const uint nSeps = adjacent_separator_.size();

  uint s=0;
  while(adjacent_separator_[s] != outgoing_sep) {
    s++;
    assert (s < nSeps);
  }

  if (outgoing_sep->var1() == involved_var_[0]) {

    if (outgoing_sep->var2() == involved_var_[1]) {

      forward.resize(nLabels1,nLabels2);

      for (uint l1=0; l1 < nLabels1; l1++) {
        for (uint l2=0; l2 < nLabels2; l2++) {

          double best = 1e300;
	  
          for (uint l3=0; l3 < nLabels3; l3++) {

            double hyp = cost_(l1,l2,l3);

            hyp -= param[2][l3]; 

            for (uint ss=0; ss < nSeps; ss++) {

              hyp -= eval_pair(ss,l1,l2,l3,pair_param);
            }

            if (hyp < best)
              best = hyp;
          }

          forward(l1,l2) = best - param[0][l1] - param[1][l2];
        }
      }
    }
    else {

      forward.resize(nLabels1,nLabels3);

      for (uint l1=0; l1 < nLabels1; l1++) {
        for (uint l3=0; l3 < nLabels3; l3++) {

          double best = 1e300;

          for (uint l2=0; l2 < nLabels2; l2++) {

            double hyp = cost_(l1,l2,l3);

            hyp -= param[1][l2]; 

            for (uint ss=0; ss < nSeps; ss++) {

              hyp -= eval_pair(ss,l1,l2,l3,pair_param);
            }

            if (hyp < best)
              best = hyp;
          }

          forward(l1,l3) = best - param[0][l1] - param[2][l3];
        }
      }
    }
  }
  else {

    forward.resize(nLabels2,nLabels3);

    for (uint l2=0; l2 < nLabels2; l2++) {
      for (uint l3=0; l3 < nLabels3; l3++) {

        double best = 1e300;

        for (uint l1=0; l1 < nLabels1; l1++) {

          double hyp = cost_(l1,l2,l3);
	  
          hyp -= param[0][l1]; 
	  
          for (uint ss=0; ss < nSeps; ss++) {
	    
            hyp -= eval_pair(ss,l1,l2,l3,pair_param);
          }
	  
          if (hyp < best)
            best = hyp;
        }
        forward(l2,l3) = best - param[1][l2] - param[2][l3];
      }
    }
  }

  return 0.0;
}

double TernaryAllInclusiveSepCumTRWSFactor::eval_pair(uint pair_num, uint x, uint y, uint z) const {

  AllInclusiveSepCumTRWSVariable* v1 = adjacent_separator_[pair_num]->var1();
  AllInclusiveSepCumTRWSVariable* v2 = adjacent_separator_[pair_num]->var2();
  
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
  
  return adjacent_separator_[pair_num]->pair_parameters()(a,b) - pair_reparameterization_[pair_num](a,b);
}

double TernaryAllInclusiveSepCumTRWSFactor::eval_pair(uint pair_num, uint x, uint y, uint z, 
						      const Storage1D< Math2D::Matrix<double> >& pair_param) const {
  
  AllInclusiveSepCumTRWSVariable* v1 = adjacent_separator_[pair_num]->var1();
  AllInclusiveSepCumTRWSVariable* v2 = adjacent_separator_[pair_num]->var2();
  
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
  
  return pair_param[pair_num](a,b);
}


/*******************/

FourthOrderAllInclusiveSepCumTRWSFactor::FourthOrderAllInclusiveSepCumTRWSFactor(const Storage1D<AllInclusiveSepCumTRWSVariable*>& vars, 
                                                                                 const Storage1D<AllInclusiveSepCumTRWSPairSeparator*>& separators,
                                                                                 const Storage1D<Math3D::Tensor<float> >& cost) :
  AllInclusiveSepCumTRWSFactor(vars,separators), cost_(cost) {

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

/*virtual*/ FourthOrderAllInclusiveSepCumTRWSFactor::~FourthOrderAllInclusiveSepCumTRWSFactor() {}

/*virtual*/ 
double FourthOrderAllInclusiveSepCumTRWSFactor::compute_var_reparameterization(AllInclusiveSepCumTRWSVariable* var) {

  double offs = 0.0;

  uint idx = 0;

  const AllInclusiveSepCumTRWSPairSeparator* valid_sep = this->valid_separator();

  if (valid_sep != 0 && (valid_sep->var1() == var || valid_sep->var2() == var)) {
    //if (false) {

    offs += valid_offs_;

    for (uint v=0; v < involved_var_.size(); v++) {
      if (involved_var_[v] == var) {
	idx = v;
	break;
      }
    }
    assert(involved_var_[idx] == var);

    Math1D::Vector<double>& repar = var_reparameterization_[idx];

    const Math2D::Matrix<double>& pair_cost = valid_sep->pair_parameters();
    const Math1D::Vector<double>& otherpar = (valid_sep->var1() == var) ? valid_sep->var2()->cost() : valid_sep->var1()->cost();

    if (valid_sep->var1() == var) {

      for (uint l1=0; l1 < var->nLabels(); l1++) {
	
	double best = 1e300;
	
	for (uint l2=0; l2 < otherpar.size(); l2++) {

	  const double hyp = pair_cost(l1,l2) + otherpar[l2];
	  
	  if (hyp < best)
	    best = hyp;
	}

	//repar was used for the computation of pair_cost -> cancel this
	repar[l1] = best + repar[l1];
      }
    }
    else {

      for (uint l2=0; l2 < var->nLabels(); l2++) {

	double best = 1e300;

	for (uint l1=0; l1 < otherpar.size(); l1++) {

	  const double hyp = pair_cost(l1,l2) + otherpar[l1];

	  if (hyp < best)
	    best = hyp;
	}

	//repar was used for the computation of pair_cost -> cancel this
	repar[l2] = best + repar[l2];
      }
    }
  }
  else {

    const uint nLabels1 = involved_var_[0]->nLabels();
    const uint nLabels2 = involved_var_[1]->nLabels();
    const uint nLabels3 = involved_var_[2]->nLabels();
    const uint nLabels4 = involved_var_[3]->nLabels();
    
    const uint nSeps = adjacent_separator_.size();

    Storage1D<Math1D::Vector<double> > param = var_reparameterization_;
    for (uint v=0; v < involved_var_.size(); v++) {
      param[v] -= involved_var_[v]->cost();
    }

    if (var == involved_var_[0]) {
    
      idx = 0;

      Math1D::Vector<double>& repar = var_reparameterization_[0];
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {
	
	double best = 1e300;

        const Math3D::Tensor<float>& cur_cost = cost_[l1];
	
	for (uint l2 = 0; l2 < nLabels2; l2++) {
	  
	  const double w2 = param[1][l2];
	
	  for (uint l3 = 0; l3 < nLabels3; l3++) {

            const double w3 =  w2 + param[2][l3];

            for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
              double hyp = cur_cost(l2,l3,l4) - w3;
              hyp -= param[3][l4];
              
              for (uint s=0; s < nSeps; s++)
                hyp += eval_pair(s,l1,l2,l3,l4);
              
              if (hyp < best)
                best = hyp;
            }
	  }
	}
      
	repar[l1] = best;
      }
    }
    else if (var == involved_var_[1]) {

      idx = 1;

      Math1D::Vector<double>& repar = var_reparameterization_[1];

      for (uint l2 = 0; l2 < nLabels2; l2++) {
      
	double best = 1e300;
	
	for (uint l1 = 0; l1 < nLabels1; l1++) {
	
	  const double w1 = param[0][l1];

          const Math3D::Tensor<float>& cur_cost = cost_[l1];
	
	  for (uint l3 = 0; l3 < nLabels3; l3++) {

            const double w3 = w1 + param[2][l3];

            for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
              double hyp = cur_cost(l2,l3,l4) - w3;
              hyp -= param[3][l4];
              
              for (uint s=0; s < nSeps; s++)
                hyp += eval_pair(s,l1,l2,l3,l4);
              
              if (hyp < best)
                best = hyp;
            }
	  }
	}

	repar[l2] = best;
      }    
    }
    else if (var == involved_var_[2]) {
    
      idx = 2;

      Math1D::Vector<double>& repar = var_reparameterization_[2];

      for (uint l3 = 0; l3 < nLabels3; l3++) {
	
	double best = 1e300;
	
	for (uint l1 = 0; l1 < nLabels1; l1++) {
	  
	  const double w1 = param[0][l1];
	
          const Math3D::Tensor<float>& cur_cost = cost_[l1];

	  for (uint l2 = 0; l2 < nLabels2; l2++) {

            const double w2 = w1 + param[1][l2];

            for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
              double hyp = cur_cost(l2,l3,l4) - w2;
              hyp -= param[3][l4];
              
              for (uint s=0; s < nSeps; s++)
                hyp += eval_pair(s,l1,l2,l3,l4);
              
              if (hyp < best)
                best = hyp;
            }
	  }
	}
	
	repar[l3] = best;
      }
    }    
    else {
      assert(var == involved_var_[3]);

      idx = 3;

      Math1D::Vector<double>& repar = var_reparameterization_[3];

      for (uint l4 = 0; l4 < nLabels4; l4++) {

	double best = 1e300;

	for (uint l1 = 0; l1 < nLabels1; l1++) {
	  
	  const double w1 = param[0][l1];
	
          const Math3D::Tensor<float>& cur_cost = cost_[l1];

	  for (uint l2 = 0; l2 < nLabels2; l2++) {

            const double w2 = w1 + param[1][l2];

            for (uint l3 = 0; l3 < nLabels3; l3++) {

              double hyp = cur_cost(l2,l3,l4) - w2;
              hyp -= param[2][l3];
              
              for (uint s=0; s < adjacent_separator_.size(); s++)
                hyp += eval_pair(s,l1,l2,l3,l4);
              
              if (hyp < best)
                best = hyp;
            }
          }
        }
        
        repar[l4] = best;
      }
    }
  }


  Math1D::Vector<double>& repar = var_reparameterization_[idx];

  const double msg_offs = repar.min();

  for (uint l=0; l < repar.size(); l++)
    repar[l] -= msg_offs;

  offs += msg_offs;

  valid_sep_ = 255;

  return offs;
}

/*virtual*/ 
double FourthOrderAllInclusiveSepCumTRWSFactor::compute_pair_reparameterization(AllInclusiveSepCumTRWSPairSeparator* pair) {

  double offs = 0.0;

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  const uint nLabels4 = involved_var_[3]->nLabels();

  Storage1D<Math1D::Vector<double> > param = var_reparameterization_;
  for (uint v=0; v < involved_var_.size(); v++) {
    //do NOT include cost of variables included in the pair-separator
    if (pair->var1() != involved_var_[v] && pair->var2() != involved_var_[v])
      param[v] -= involved_var_[v]->cost();
  }

  const uint nSeps = adjacent_separator_.size();

  uint s=0;
  while(adjacent_separator_[s] != pair) {
    s++;
    assert (s < nSeps);
  }

  Math2D::Matrix<double>& repar = pair_reparameterization_[s];

  if (pair->var1() == involved_var_[0]) {

    if (pair->var2() == involved_var_[1]) {

      repar.resize(nLabels1,nLabels2);

      for (uint l1=0; l1 < nLabels1; l1++) {
        
        const double w1 = param[0][l1];

        for (uint l2=0; l2 < nLabels2; l2++) {

          double best = 1e300;
	  
          for (uint l3=0; l3 < nLabels3; l3++) {
            
            const double w3 = param[2][l3];
            
            for (uint l4=0; l4 < nLabels4; l4++) {

              double hyp = cost_[l1](l2,l3,l4) - w3;
              hyp -= param[3][l4]; 
              
              for (uint ss=0; ss < nSeps; ss++) {
                
                if (ss != s)
                  hyp += eval_pair(ss,l1,l2,l3,l4);
              }
              
              if (hyp < best)
                best = hyp;
            }
          }

          repar(l1,l2) = best - w1 - param[1][l2];
        }
      }
    }
    else if (pair->var2() == involved_var_[2]) {

      repar.resize(nLabels1,nLabels3);

      for (uint l1=0; l1 < nLabels1; l1++) {

        const double w1 = param[0][l1];

        for (uint l3=0; l3 < nLabels3; l3++) {

          double best = 1e300;

          for (uint l2=0; l2 < nLabels2; l2++) {

            const double w2 = param[1][l2];

            for (uint l4=0; l4 < nLabels4; l4++) {

              double hyp = cost_[l1](l2,l3,l4) - w2;
              hyp -= param[3][l4]; 
              
              for (uint ss=0; ss < nSeps; ss++) {
                
                if (ss != s)
                  hyp += eval_pair(ss,l1,l2,l3,l4);
              }
              
              if (hyp < best)
                best = hyp;
            }
          }

          repar(l1,l3) = best - w1 - param[2][l3];
        }
      }
    }
    else {

      assert(pair->var2() == involved_var_[3]);

      repar.resize(nLabels1,nLabels4);

      for (uint l1=0; l1 < nLabels1; l1++) {

        const double w1 = param[0][l1];

        for (uint l4=0; l4 < nLabels4; l4++) {

          double best = 1e300;

          for (uint l2=0; l2 < nLabels2; l2++) {
            
            const double w2 = param[1][l2];

            for (uint l3=0; l3 < nLabels3; l3++) {

              double hyp = cost_[l1](l2,l3,l4) - w2;
              hyp -= param[2][l3]; 
              
              for (uint ss=0; ss < nSeps; ss++) {
                
                if (ss != s)
                  hyp += eval_pair(ss,l1,l2,l3,l4);
              }
              
              if (hyp < best)
                best = hyp;
            }
          }

          repar(l1,l4) = best - w1 - param[3][l4];
        }
      }
    }
  }
  else if (pair->var1() == involved_var_[1]) {
    
    if (pair->var2() == involved_var_[2]) {

      repar.resize(nLabels2,nLabels3);

      for (uint l2=0; l2 < nLabels2; l2++) {

        const double w2 = param[1][l2];

        for (uint l3=0; l3 < nLabels3; l3++) {

          double best = 1e300;

          for (uint l1=0; l1 < nLabels1; l1++) {

            const double w1 = param[0][l1];

            for (uint l4=0; l4 < nLabels4; l4++) {

              double hyp = cost_[l1](l2,l3,l4) - w1;
              hyp -= param[3][l4]; 
              
              for (uint ss=0; ss < nSeps; ss++) {
                
                if (ss != s)
                  hyp += eval_pair(ss,l1,l2,l3,l4);
              }
              
              if (hyp < best)
                best = hyp;
            }
          }

          repar(l2,l3) = best - w2 - param[2][l3];
        }
      }
    }
    else {

      assert(pair->var2() == involved_var_[3]);

      repar.resize(nLabels2,nLabels4);

      for (uint l2=0; l2 < nLabels2; l2++) {

        const double w2 = param[1][l2];

        for (uint l4=0; l4 < nLabels4; l4++) {

          double best = 1e300;

          for (uint l1=0; l1 < nLabels1; l1++) {

            const double w1 = param[0][l1];

            for (uint l3=0; l3 < nLabels3; l3++) {

              double hyp = cost_[l1](l2,l3,l4) - w1;
              hyp -= param[2][l3]; 
              
              for (uint ss=0; ss < nSeps; ss++) {
                
                if (ss != s)
                  hyp += eval_pair(ss,l1,l2,l3,l4);
              }
              
              if (hyp < best)
                best = hyp;
            }
          }

          repar(l2,l4) = best - w2 - param[3][l4];
        }
      }
    }
  }
  else {

    assert (pair->var1() == involved_var_[2]);

    repar.resize(nLabels3,nLabels4);

    for (uint l3=0; l3 < nLabels3; l3++) {

      const double w3 = param[2][l3];

      for (uint l4=0; l4 < nLabels4; l4++) {

        double best = 1e300;
        
        for (uint l1=0; l1 < nLabels1; l1++) {

          const double w1 = param[0][l1];

          for (uint l2=0; l2 < nLabels2; l2++) {

            double hyp = cost_[l1](l2,l3,l4) - w1;
            hyp -= param[1][l2]; 
              
            for (uint ss=0; ss < nSeps; ss++) {
              
              if (ss != s)
                hyp += eval_pair(ss,l1,l2,l3,l4);
            }
              
            if (hyp < best)
              best = hyp;
          }
        }
        
        repar(l3,l4) = best - w3 - param[3][l4];
      }
    }
  }

  offs = repar.min();
  for (uint k=0; k < repar.size(); k++)
    repar.direct_access(k) -= offs;

  valid_sep_ = s;
  valid_offs_ = offs;

  return offs;
}

/*virtual*/ 
double FourthOrderAllInclusiveSepCumTRWSFactor::best_value() {

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  const uint nLabels4 = involved_var_[2]->nLabels();

  double min_val = 1e300;

  const uint nSeps = adjacent_separator_.size();

  for (uint l1 = 0; l1 < nLabels1; l1++) {
    
    for (uint l2 = 0; l2 < nLabels2; l2++) {
	
      for (uint l3 = 0; l3 < nLabels3; l3++) {

        for (uint l4 = 0; l4 < nLabels4; l4++) {

          double hyp = cost_[l1](l2,l3,l4);
          hyp += involved_var_[0]->cost()[l1] - var_reparameterization_[0][l1];
          hyp += involved_var_[1]->cost()[l2] - var_reparameterization_[1][l2];
          hyp += involved_var_[2]->cost()[l3] - var_reparameterization_[2][l3];
          hyp += involved_var_[3]->cost()[l4] - var_reparameterization_[3][l4];

          for (uint s=0; s < nSeps; s++)
            hyp += eval_pair(s,l1,l2,l3,l4);
  
          if (hyp < min_val)
            min_val = hyp;
        }
      }
    }
  }

  return min_val;
}

/*virtual*/ 
double FourthOrderAllInclusiveSepCumTRWSFactor::compute_forward(const AllInclusiveSepCumTRWSPairSeparator* incoming_sep, 
                                                                const AllInclusiveSepCumTRWSVariable* incoming_var, 
                                                                const AllInclusiveSepCumTRWSVariable* outgoing_var,
                                                                const Math2D::Matrix<double>& prev_pair_forward, 
                                                                const Math1D::Vector<double>& prev_var_forward, 
                                                                Math1D::Vector<double>& forward) {

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  const uint nLabels4 = involved_var_[3]->nLabels();
  
  Storage1D<Math1D::Vector<double> > param = var_reparameterization_;
  for (uint v=0; v < involved_var_.size(); v++) {
    //do NOT include cost of variables included in the pair-separator
    if (involved_var_[v] != outgoing_var)
      param[v] -= involved_var_[v]->cost();
    if (involved_var_[v] == incoming_var)
      param[v] -= prev_var_forward;
  }

  const uint nSeps = adjacent_separator_.size();

  Storage1D<Math2D::Matrix<double> > pair_param = pair_reparameterization_;
  for (uint s=0; s < pair_param.size(); s++) {
    
    pair_param[s] -= adjacent_separator_[s]->pair_parameters();
    if (adjacent_separator_[s] == incoming_sep)
      pair_param[s] -= prev_pair_forward;
  }

  if (outgoing_var == involved_var_[0]) {

    forward.resize(nLabels1);

    for (uint l1 = 0; l1 < nLabels1; l1++) {
      
      double best = 1e300;

      for (uint l2 = 0; l2 < nLabels2; l2++) {

        const double w2 = param[1][l2];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cost_[l1](l2,l3,l4) - w2;
            hyp -= param[2][l3];
            hyp -= param[3][l4];
            
            for (uint s=0; s < nSeps; s++)
              hyp -= eval_pair(s,l1,l2,l3,l4,pair_param);
            
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      forward[l1] = best - param[0][l1];
    }
  }
  else if (outgoing_var == involved_var_[1]) {

    forward.resize(nLabels2);

    for (uint l2 = 0; l2 < nLabels2; l2++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const double w1 = param[0][l1];
	
        for (uint l3 = 0; l3 < nLabels3; l3++) {

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cost_[l1](l2,l3,l4) - w1;
            hyp -= param[2][l3];
            hyp -= param[3][l4];
            
            for (uint s=0; s < nSeps; s++)
              hyp -= eval_pair(s,l1,l2,l3,l4,pair_param);
            
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      forward[l2] = best  - param[1][l2];
    }    
  }
  else if (outgoing_var == involved_var_[2]) {

    forward.resize(nLabels3);

    for (uint l3 = 0; l3 < nLabels3; l3++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const double w1 = param[0][l1];
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {

          for (uint l4 = 0; l4 < nLabels4; l4++) {
	  
            double hyp = cost_[l1](l2,l3,l4) - w1;
            hyp -= param[1][l2];
            hyp -= param[3][l4];

            for (uint s=0; s < nSeps; s++)
              hyp -= eval_pair(s,l1,l2,l3,l4,pair_param);
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }
      
      forward[l3] = best  - param[2][l3];
    }
  }
  else {

    forward.resize(nLabels4);

    for (uint l4 = 0; l4 < nLabels4; l4++) {
      
      double best = 1e300;
      
      for (uint l1 = 0; l1 < nLabels1; l1++) {

        const double w1 = param[0][l1];
	
        for (uint l2 = 0; l2 < nLabels2; l2++) {

          for (uint l3 = 0; l3 < nLabels3; l3++) {

            double hyp = cost_[l1](l2,l3,l4) - w1;
            hyp -= param[1][l2];
            hyp -= param[2][l3];

            for (uint s=0; s < nSeps; s++)
              hyp -= eval_pair(s,l1,l2,l3,l4,pair_param);
	  
            if (hyp < best)
              best = hyp;
          }
        }
      }

      forward[l4] = best  - param[3][l4];
    }
  }

  double offs = forward.min();
  for (uint k=0; k < forward.size(); k++)
    forward[k] -= offs;

  return offs;
}

/*virtual*/ 
double FourthOrderAllInclusiveSepCumTRWSFactor::compute_forward(const AllInclusiveSepCumTRWSPairSeparator* incoming_sep, 
                                                                const AllInclusiveSepCumTRWSVariable* incoming_var, 
                                                                const AllInclusiveSepCumTRWSPairSeparator* outgoing_sep,
                                                                const Math2D::Matrix<double>& prev_pair_forward, 
                                                                const Math1D::Vector<double>& prev_var_forward, 
                                                                Math2D::Matrix<double>& forward) {

  const uint nLabels1 = involved_var_[0]->nLabels();
  const uint nLabels2 = involved_var_[1]->nLabels();
  const uint nLabels3 = involved_var_[2]->nLabels();
  const uint nLabels4 = involved_var_[3]->nLabels();
  
  Storage1D<Math1D::Vector<double> > param = var_reparameterization_;
  for (uint v=0; v < involved_var_.size(); v++) {
    //do NOT include cost of variables included in the pair-separator
    if (outgoing_sep->var1() != involved_var_[v] && outgoing_sep->var2() != involved_var_[v])
      param[v] -= involved_var_[v]->cost();
    if (involved_var_[v] == incoming_var)
      param[v] -= prev_var_forward;
  }

  const uint nSeps = adjacent_separator_.size();

  Storage1D<Math2D::Matrix<double> > pair_param = pair_reparameterization_;
  for (uint s=0; s < pair_param.size(); s++) {

    if (adjacent_separator_[s] != outgoing_sep)
      pair_param[s] -= adjacent_separator_[s]->pair_parameters();
    if (adjacent_separator_[s] == incoming_sep)
      pair_param[s] -= prev_pair_forward;
  }

  uint s=0;
  while(adjacent_separator_[s] != outgoing_sep) {
    s++;
    assert (s < nSeps);
  }


  if (outgoing_sep->var1() == involved_var_[0]) {

    if (outgoing_sep->var2() == involved_var_[1]) {

      forward.resize(nLabels1,nLabels2);

      for (uint l1=0; l1 < nLabels1; l1++) {
        for (uint l2=0; l2 < nLabels2; l2++) {

          double best = 1e300;
	  
          for (uint l3=0; l3 < nLabels3; l3++) {
            for (uint l4=0; l4 < nLabels4; l4++) {

              double hyp = cost_[l1](l2,l3,l4);

              hyp -= param[2][l3]; 
              hyp -= param[3][l4]; 
              
              for (uint ss=0; ss < nSeps; ss++) {
                
                hyp -= eval_pair(ss,l1,l2,l3,l4,pair_param);
              }
              
              if (hyp < best)
                best = hyp;
            }
          }

          forward(l1,l2) = best - param[0][l1] - param[1][l2];
        }
      }
    }
    else if (outgoing_sep->var2() == involved_var_[2]) {

      forward.resize(nLabels1,nLabels3);

      for (uint l1=0; l1 < nLabels1; l1++) {
        for (uint l3=0; l3 < nLabels3; l3++) {

          double best = 1e300;

          for (uint l2=0; l2 < nLabels2; l2++) {
            for (uint l4=0; l4 < nLabels4; l4++) {

              double hyp = cost_[l1](l2,l3,l4);

              hyp -= param[1][l2]; 
              hyp -= param[3][l4]; 
              
              for (uint ss=0; ss < nSeps; ss++) {
                
                hyp -= eval_pair(ss,l1,l2,l3,l4,pair_param);
              }
              
              if (hyp < best)
                best = hyp;
            }
          }

          forward(l1,l3) = best - param[0][l1] - param[2][l3];
        }
      }
    }
    else {

      assert(outgoing_sep->var2() == involved_var_[3]);

      forward.resize(nLabels1,nLabels4);

      for (uint l1=0; l1 < nLabels1; l1++) {
        for (uint l4=0; l4 < nLabels4; l4++) {

          double best = 1e300;

          for (uint l2=0; l2 < nLabels2; l2++) {
            for (uint l3=0; l3 < nLabels3; l3++) {

              double hyp = cost_[l1](l2,l3,l4);

              hyp -= param[1][l2]; 
              hyp -= param[2][l3]; 
              
              for (uint ss=0; ss < nSeps; ss++) {
                
                hyp -= eval_pair(ss,l1,l2,l3,l4,pair_param);
              }
              
              if (hyp < best)
                best = hyp;
            }
          }

          forward(l1,l4) = best - param[0][l1] - param[3][l4];
        }
      }
    }
  }
  else if (outgoing_sep->var1() == involved_var_[1]) {
    
    if (outgoing_sep->var2() == involved_var_[2]) {

      forward.resize(nLabels2,nLabels3);

      for (uint l2=0; l2 < nLabels2; l2++) {
        for (uint l3=0; l3 < nLabels3; l3++) {

          double best = 1e300;

          for (uint l1=0; l1 < nLabels1; l1++) {
            for (uint l4=0; l4 < nLabels4; l4++) {

              double hyp = cost_[l1](l2,l3,l4);

              hyp -= param[0][l1]; 
              hyp -= param[3][l4]; 
              
              for (uint ss=0; ss < nSeps; ss++) {
                
                hyp -= eval_pair(ss,l1,l2,l3,l4,pair_param);
              }
              
              if (hyp < best)
                best = hyp;
            }
          }

          forward(l2,l3) = best - param[1][l2] - param[2][l3];
        }
      }
    }
    else {

      assert(outgoing_sep->var2() == involved_var_[3]);

      forward.resize(nLabels2,nLabels4);

      for (uint l2=0; l2 < nLabels2; l2++) {
        for (uint l4=0; l4 < nLabels4; l4++) {

          double best = 1e300;

          for (uint l1=0; l1 < nLabels1; l1++) {
            for (uint l3=0; l3 < nLabels3; l3++) {

              double hyp = cost_[l1](l2,l3,l4);

              hyp -= param[0][l1]; 
              hyp -= param[2][l3]; 
              
              for (uint ss=0; ss < nSeps; ss++) {
                
                hyp -= eval_pair(ss,l1,l2,l3,l4,pair_param);
              }
              
              if (hyp < best)
                best = hyp;
            }
          }

          forward(l2,l4) = best - param[1][l2] - param[3][l4];
        }
      }
    }
  }
  else {

    assert (outgoing_sep->var1() == involved_var_[2]);

    forward.resize(nLabels3,nLabels4);

    for (uint l3=0; l3 < nLabels3; l3++) {
      for (uint l4=0; l4 < nLabels4; l4++) {

        double best = 1e300;
        
        for (uint l1=0; l1 < nLabels1; l1++) {
          for (uint l2=0; l2 < nLabels2; l2++) {

            double hyp = cost_[l1](l2,l3,l4);

            hyp -= param[0][l1]; 
            hyp -= param[1][l2]; 
              
            for (uint ss=0; ss < nSeps; ss++) {
              
              hyp -= eval_pair(ss,l1,l2,l3,l4,pair_param);
            }
              
            if (hyp < best)
              best = hyp;
          }
        }
        
        forward(l3,l4) = best - param[2][l3] - param[3][l4];
      }
    }
  }
  
  return 0.0;
}


double FourthOrderAllInclusiveSepCumTRWSFactor::eval_pair(uint pair_num, uint x, uint y, uint z, uint w) const {

  AllInclusiveSepCumTRWSVariable* v1 = adjacent_separator_[pair_num]->var1();
  AllInclusiveSepCumTRWSVariable* v2 = adjacent_separator_[pair_num]->var2();
  
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

  return adjacent_separator_[pair_num]->pair_parameters()(a,b) - pair_reparameterization_[pair_num](a,b);
}

double FourthOrderAllInclusiveSepCumTRWSFactor::eval_pair(uint pair_num, uint x, uint y, uint z, uint w,
                                                          const Storage1D< Math2D::Matrix<double> >& pair_param) const {

  AllInclusiveSepCumTRWSVariable* v1 = adjacent_separator_[pair_num]->var1();
  AllInclusiveSepCumTRWSVariable* v2 = adjacent_separator_[pair_num]->var2();
  
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

  return pair_param[pair_num](a,b);
}

/*******************/

AllInclusiveSepCumTRWS::AllInclusiveSepCumTRWS(uint nVars, uint nSeparators, uint nFactors) :
  nUsedVars_(0), nUsedSeparators_(0), nUsedFactors_(0), optimize_called_(false) {

  var_.resize(nVars,0);
  separator_.resize(nSeparators);
  factor_.resize(nFactors);
}

AllInclusiveSepCumTRWS::~AllInclusiveSepCumTRWS() {
  
  for (uint v=0; v < nUsedVars_; v++)
    delete var_[v];

  for (uint s=0; s < nUsedSeparators_; s++) 
    delete separator_[s];

  for (uint f=0; f < nUsedFactors_; f++)
    delete factor_[f];
}

uint AllInclusiveSepCumTRWS::add_var(const Math1D::Vector<float>& cost) {

  assert(!optimize_called_);

  if (nUsedVars_ == var_.size())
    var_.resize(uint(1.2*nUsedVars_)+4);

  assert(nUsedVars_ < var_.size());

  var_[nUsedVars_] = new AllInclusiveSepCumTRWSVariable(cost,nUsedVars_);

  nUsedVars_++;
  return nUsedVars_-1;
}

uint AllInclusiveSepCumTRWS::add_pair_separator(uint var1, uint var2) {

  if (var1 >= nUsedVars_ || var2 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  assert(!optimize_called_);

  if (nUsedSeparators_ == separator_.size())
    separator_.resize(uint(1.2*nUsedSeparators_)+4);

  assert(nUsedSeparators_ < separator_.size());

  assert(var1 < var2);

  separator_[nUsedSeparators_] = new AllInclusiveSepCumTRWSPairSeparator(var_[var1],var_[var2]);

  nUsedSeparators_++;
  return nUsedSeparators_-1;
}

void AllInclusiveSepCumTRWS::add_factor(AllInclusiveSepCumTRWSFactor* fac) {

  assert(!optimize_called_);

  if (nUsedFactors_ == factor_.size())
    factor_.resize(uint(1.2*nUsedFactors_)+4);

  factor_[nUsedFactors_] = fac;
  nUsedFactors_++;
}

void AllInclusiveSepCumTRWS::add_binary_factor(uint v1, uint v2, const Math2D::Matrix<float>& cost)
{
  if (v1 >= nUsedVars_ || v2 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  Storage1D<AllInclusiveSepCumTRWSVariable*> vars(2);
  vars[0] = var_[v1];
  vars[1] = var_[v2];

  add_factor(new BinaryAllInclusiveSepCumTRWSFactor(vars,cost));
}

void AllInclusiveSepCumTRWS::add_ternary_factor(uint v1, uint v2, uint v3, const Storage1D<uint>& separators,
						const Math3D::Tensor<float>& cost) {

  if (v1 >= nUsedVars_ || v2 >= nUsedVars_ || v3 >= nUsedVars_) {
    INTERNAL_ERROR << "out of range. Exiting." << std::endl;
    exit(1);
  }

  Math3D::Tensor<float> cost_copy = cost;

#if 1
  if (v1 > v2) {

    //std::cerr << "v1-v2 1." << std::endl;

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

    //std::cerr << "v2-v3" << std::endl;

    if (var_[v2]->nLabels() != var_[v3]->nLabels())
      TODO("non-standard variable order with heterogeneous number of labels");

    for (uint x=0; x < var_[v1]->nLabels(); x++)
      for (uint y=0; y < var_[v2]->nLabels(); y++)
        for (uint z=y+1; z < var_[v3]->nLabels(); z++)
          std::swap(cost_copy(x,y,z),cost_copy(x,z,y));
    
    std::swap(v2,v3);
  }
  if (v1 > v2) {

    //std::cerr << "v1-v2 2." << std::endl;

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

  assert(nUsedFactors_ < factor_.size());

  Storage1D<AllInclusiveSepCumTRWSVariable*> vars(3);
  vars[0] = var_[v1];
  vars[1] = var_[v2];
  vars[2] = var_[v3];

  Storage1D<AllInclusiveSepCumTRWSPairSeparator*> seps(separators.size());
  for (uint s=0; s < separators.size(); s++) {

    if (separators[s] >= nUsedSeparators_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    seps[s] = separator_[separators[s]];
  }


  add_factor(new TernaryAllInclusiveSepCumTRWSFactor(vars,seps,cost_copy));
}

void AllInclusiveSepCumTRWS::add_fourth_order_factor(uint v1, uint v2, uint v3, uint v4, const Storage1D<uint>& separators,
						     const Storage1D<Math3D::Tensor<float> >& cost)
{

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

  assert(nUsedFactors_ < factor_.size());
  
  Storage1D<AllInclusiveSepCumTRWSVariable*> vars(4);
  vars[0] = var_[v1];
  vars[1] = var_[v2];
  vars[2] = var_[v3];
  vars[3] = var_[v4];
  
  Storage1D<AllInclusiveSepCumTRWSPairSeparator*> seps(separators.size());
  for (uint s=0; s < separators.size(); s++) {

    if (separators[s] >= nUsedSeparators_) {
      INTERNAL_ERROR << "out of range. Exiting." << std::endl;
      exit(1);
    }

    seps[s] = separator_[separators[s]];
  }

  add_factor(new FourthOrderAllInclusiveSepCumTRWSFactor(vars,seps,cost_copy));
}

double AllInclusiveSepCumTRWS::optimize(uint nIter)
{

  labeling_.resize(nUsedVars_,0);

  //dry-run to set the separator ranks and start- and end pair-separators of each factor

  double bound = -1e300;


  if (!optimize_called_) {

    uint next_sep_rank = 0;
    
    for (uint v=0; v < nUsedVars_; v++) {
      
      const Storage1D<AllInclusiveSepCumTRWSPairSeparator*>& adj_sep = var_[v]->adjacent_separators();
      
      for (uint s=0; s < adj_sep.size(); s++) {
        
        if (var_[v] == adj_sep[s]->var2()) {
          
          adj_sep[s]->set_sep_rank(next_sep_rank);
          next_sep_rank++;
        }
      }
    }
    
    for (uint f=0; f < nUsedFactors_; f++)
      factor_[f]->compute_rank_range();
    
    for (uint s=0; s < nUsedSeparators_; s++) {
      separator_[s]->set_up_chains();
    }
    
    for (uint v=0; v < nUsedVars_; v++) {
      var_[v]->set_up_chains();
    }
    
    for (uint f=0; f < nUsedFactors_; f++) {
      if (factor_[f]->prev_factor() == 0)
        factor_[f]->set_start_separator(0);
      if (factor_[f]->next_factor() == 0)
        factor_[f]->set_end_separator(0);      
    }
  }

  uint arg_min;
  
  std::cerr << "starting main loop of memory-efficient scheme" << std::endl;

  optimize_called_ = true;

  std::cerr.precision(8);

  //std::cerr << "start energy: " << cur_bound() << std::endl;

  size_t effort_per_iteration = 0;

  {
    size_t message_effort = 0;
    
    for (uint f=0; f < nUsedFactors_; f++) {
      
      uint var_size = factor_[f]->involved_var().size();
      uint sep_size = factor_[f]->adjacent_separator().size();
      
      effort_per_iteration += (var_size + sep_size-1) * (var_size + sep_size);
    }
    message_effort = effort_per_iteration * nIter;

    std::cerr << "predicted message effort: " << message_effort << std::endl;
  }

  for (uint iter=1; iter <= nIter; iter++) {

    /*** forward ***/

    
    double fwd_bound = 0.0;
    for (uint v=0; v < nUsedVars_; v++) {

      //std::cerr << "v: " << v << std::endl;

      //average the separators
      const Storage1D<AllInclusiveSepCumTRWSPairSeparator*>& adj_sep = var_[v]->adjacent_separators();

      for (uint s=0; s < adj_sep.size(); s++) {
        
        if (var_[v] == adj_sep[s]->var2()) {

	  fwd_bound += adj_sep[s]->reparameterize_forward();

	  double temp = adj_sep[s]->average();
	  fwd_bound += temp;
	}
      }

      //std::cerr << "now averaging var" << std::endl;

      //now average the var
      double temp = var_[v]->reparameterize_forward();
      fwd_bound += temp;

      temp = var_[v]->average(arg_min);
      fwd_bound += temp;
      labeling_[v] = arg_min;
    }

    std::cerr << "iteration " << iter << ", forward bound: " << fwd_bound << std::endl;

    /*** backward ***/

    double bwd_bound = 0.0;
    for (int v=nUsedVars_-1; v >= 0; v--) {

      //average the var
      
      bwd_bound += var_[v]->reparameterize_backward();

      double temp = var_[v]->average(arg_min);
      bwd_bound += temp;
      labeling_[v] = arg_min;

      //now average the pair separators
      const Storage1D<AllInclusiveSepCumTRWSPairSeparator*>& adj_sep = var_[v]->adjacent_separators();

      for (int s= adj_sep.size() - 1; s >= 0; s--) {
        
        if (var_[v] == adj_sep[s]->var2()) {

	  bwd_bound += adj_sep[s]->reparameterize_backward();

	  double temp = adj_sep[s]->average();
	  bwd_bound += temp;
	}
      }

    }


    std::cerr << "iteration " << iter << ", backward bound: " << bwd_bound << std::endl;

    bound = bwd_bound;
  }

  return bound;
}

const Math1D::Vector<uint>& AllInclusiveSepCumTRWS::labeling() const {
  return labeling_;
}

double AllInclusiveSepCumTRWS::cur_bound(bool backward) {

  double bound = 0.0; 

  for (uint f=0; f < nUsedFactors_; f++) {

    if (factor_[f]->prev_factor() == 0) {

      AllInclusiveSepCumTRWSFactor* cur_factor = factor_[f];
        
      std::vector<AllInclusiveSepCumTRWSFactor*> chain;
      std::vector<AllInclusiveSepCumTRWSVariable*> out_var;
      std::vector<AllInclusiveSepCumTRWSPairSeparator*> out_sep;

      //find chain start
      AllInclusiveSepCumTRWSVariable* in_var = 0;
      for (uint k=0; k < cur_factor->involved_var().size(); k++) {
        
        if (cur_factor->involved_var()[k]->rank() == cur_factor->min_rank()) {
          in_var = cur_factor->involved_var()[k];
          break;
        }
      }
      
      assert(in_var != 0);
      
      while (cur_factor != 0) {

        chain.push_back(cur_factor);
        if (cur_factor->end_separator() != 0) {
          out_sep.push_back(cur_factor->end_separator());
          out_var.push_back(0);
        }
        else {
          out_sep.push_back(0);
          for (uint k=0; k < cur_factor->involved_var().size(); k++) {
            
            if (cur_factor->involved_var()[k]->rank() == cur_factor->max_rank()) {
              out_var.push_back(cur_factor->involved_var()[k]);
              break;
            }
          }
        }
        cur_factor = cur_factor->next_factor();
      }
      
      uint chain_length = chain.size();
      
      //find chain end
      for (uint k=0; k < chain.back()->involved_var().size(); k++) {
        
        if (chain.back()->involved_var()[k]->rank() == chain.back()->max_rank()) {
          out_var.back() = chain.back()->involved_var()[k];
          break;
        }
      }

      if (backward) {

	Math1D::NamedVector<double> backward1(MAKENAME(backward1));
	Math1D::NamedVector<double> backward2(in_var->nLabels(),0.0,MAKENAME(backward2));

	Math2D::Matrix<double> pair_backward1;
	Math2D::Matrix<double> pair_backward2;

	if (chain_length > 1) {
	  if (out_var[chain_length-2] != 0) {
	    double offs = chain.back()->compute_forward(0,out_var.back(),out_var[chain_length-2],pair_backward2,backward2,backward1);
	    bound += offs;
	  }
	  else {
	    double offs = chain.back()->compute_forward(0,out_var.back(),out_sep[chain_length-2],pair_backward2,backward2,pair_backward1);
	    bound += offs;

	    assert(offs == 0.0);
	  }
	}
	else {
	  backward1.resize(out_var[0]->nLabels(),0.0);
	}

	for (int k=chain_length-2; k >= 0; k--) {

	  Math1D::Vector<double>& last_backward = (( (int(chain_length)-1-k) % 2) == 1) ? backward1 : backward2;
	  Math1D::Vector<double>& new_backward = (((int(chain_length)-1-k) % 2) == 0) ? backward1 : backward2;
        
	  Math2D::Matrix<double>& last_pair_backward = (((int(chain_length)-1-k) % 2) == 1) ? pair_backward1 : pair_backward2;
	  Math2D::Matrix<double>& new_pair_backward = (((int(chain_length)-1-k) % 2) == 0) ? pair_backward1 : pair_backward2;
	
	  if (k == 0) {

	    double offs = chain[0]->compute_forward(out_sep[k],out_var[k],in_var,last_pair_backward,last_backward,new_backward);
	    bound += offs;
	  }
	  else {

	    if (out_var[k-1] != 0) {
	      double offs = chain[k]->compute_forward(out_sep[k],out_var[k],out_var[k-1],last_pair_backward,last_backward,new_backward);
	      bound += offs;
	    }
	    else {
	      double offs = chain[k]->compute_forward(out_sep[k],out_var[k],out_sep[k-1],last_pair_backward,last_backward,new_pair_backward);
	      bound += offs;
	    }
	  }
	}
	

	Math1D::Vector<double>& new_backward = (((chain_length-1) % 2) == 0) ? backward1 : backward2;

	new_backward += in_var->cost();

	bound += new_backward.min();
      }
      else {


	Math1D::NamedVector<double> forward1(MAKENAME(forward1));
	Math1D::NamedVector<double> forward2(in_var->nLabels(),0.0,MAKENAME(forward2));
	
	Math2D::Matrix<double> pair_forward1;
	Math2D::Matrix<double> pair_forward2;

	//compute forward
	if (out_var[0] != 0) {
	  double offs = chain[0]->compute_forward(0,in_var,out_var[0],pair_forward2,forward2,forward1);
	  bound += offs;
	}
	else {
	  bound += chain[0]->compute_forward(0,in_var,out_sep[0],pair_forward2,forward2,pair_forward1);
	}
	
	for (uint k=1; k < chain_length; k++) {
	  
	  Math1D::Vector<double>& last_forward = ((k % 2) == 1) ? forward1 : forward2;
	  Math1D::Vector<double>& new_forward = ((k % 2) == 0) ? forward1 : forward2;
        
	  Math2D::Matrix<double>& last_pair_forward = ((k % 2) == 1) ? pair_forward1 : pair_forward2;
	  Math2D::Matrix<double>& new_pair_forward = ((k % 2) == 0) ? pair_forward1 : pair_forward2;
	  
	  if (out_var[k] != 0) {
	    double offs = chain[k]->compute_forward(out_sep[k-1],out_var[k-1],out_var[k],
						    last_pair_forward,last_forward,new_forward);
	    bound += offs;
	  }
	  else {
	    bound += chain[k]->compute_forward(out_sep[k-1],out_var[k-1],out_sep[k],
					       last_pair_forward,last_forward,new_pair_forward);
	  }
	}

	Math1D::Vector<double>& new_forward = (((chain_length-1) % 2) == 0) ? forward1 : forward2;
      
	assert(out_var[chain_length-1] != 0);

	new_forward += out_var[chain_length-1]->cost();

	bound += new_forward.min();
      }
    }
  }

  return bound;
}
