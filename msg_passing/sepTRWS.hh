/***** written by Thomas Schoenemann as an employee of the University of Pisa, Italy, December 2011 *****/
/****        and as an employee of the University of Düsseldorf, Germany, January - August 2012 *******/

/***** implements the solving of dual relaxations with pairwise separators *****/
/**** via TRWS with factor (junction) chains *****/


#ifndef SEP_TRWS_HH
#define SEP_TRWS_HH


#include "vector.hh"
#include "matrix.hh"
#include "tensor.hh"

#include <set>


enum VarLinkStatus {PairReuse, PairRecompute, PairIgnore};


class AllInclusiveSepCumTRWSFactor;
class AllInclusiveSepCumTRWSPairSeparator;
class AllInclusiveSepCumTRWSVariable;

//here a pairwise separator links two outer factors of a chain
class AllInclusiveSepCumTRWSSepChainLink {
public:

  AllInclusiveSepCumTRWSSepChainLink();

  AllInclusiveSepCumTRWSSepChainLink(AllInclusiveSepCumTRWSFactor* lfac, AllInclusiveSepCumTRWSFactor* rfac,
				     AllInclusiveSepCumTRWSPairSeparator* sep);

  AllInclusiveSepCumTRWSPairSeparator* sep_;

  AllInclusiveSepCumTRWSFactor* left_fac_; //the factor that is processed first in forward
  AllInclusiveSepCumTRWSFactor* right_fac_;
};

//here a single variable links two outer factors of a chain
class AllInclusiveSepCumTRWSVarChainLink {
public:

  AllInclusiveSepCumTRWSVarChainLink();

  std::vector<AllInclusiveSepCumTRWSPairSeparator*> sep_;
  std::vector<AllInclusiveSepCumTRWSFactor*> involved_factor_;
};

//variable class
class AllInclusiveSepCumTRWSVariable {
public:

  AllInclusiveSepCumTRWSVariable(const Math1D::Vector<float>& cost, uint rank);

  //will delete chain links
  ~AllInclusiveSepCumTRWSVariable();
  
  void add_factor(AllInclusiveSepCumTRWSFactor* adjacent_fac);
  
  void add_pair_separator(AllInclusiveSepCumTRWSPairSeparator* adjacent_sep);

  uint nLabels() const;

  uint rank() const;

  const Storage1D<AllInclusiveSepCumTRWSFactor*>& adjacent_factors() const;

  const Storage1D<AllInclusiveSepCumTRWSPairSeparator*>& adjacent_separators() const;

  double average(uint& arg_min);

  const Math1D::Vector<double>& cost() const;

  void set_up_chains();

  double reparameterize_forward();

  double reparameterize_backward();

protected:  

  Storage1D<AllInclusiveSepCumTRWSFactor*> adjacent_factor_;

  Storage1D<AllInclusiveSepCumTRWSPairSeparator*> adjacent_separator_;

  //current assumption: lists only variables that are not included in a separator
  //NOTE: right now, all listed entries are dead ends
  Storage1D<AllInclusiveSepCumTRWSVarChainLink*> chain_link_;

  const Math1D::Vector<float> cost_;

  Math1D::Vector<double> cum_cost_;

  uint rank_;
};

//pairwise separator class
class AllInclusiveSepCumTRWSPairSeparator {
public:

  AllInclusiveSepCumTRWSPairSeparator(AllInclusiveSepCumTRWSVariable* var1, AllInclusiveSepCumTRWSVariable* var2);

  //will delete chain links
  ~AllInclusiveSepCumTRWSPairSeparator();

  void add_factor(AllInclusiveSepCumTRWSFactor* adjacent_fac);

  AllInclusiveSepCumTRWSVariable* var1() const;

  AllInclusiveSepCumTRWSVariable* var2() const;

  const Storage1D<AllInclusiveSepCumTRWSFactor*>& adjacent_factors() const;

  double average();

  void set_up_chains();

  Math2D::Matrix<double>& pair_parameters();

  const Math2D::Matrix<double>& pair_parameters() const;

  Storage1D<AllInclusiveSepCumTRWSSepChainLink*> chain_link_;

  uint sep_rank() const;

  void set_sep_rank(uint rank);

  double reparameterize_forward();

  double reparameterize_backward();

protected:

  uint sep_rank_;

  AllInclusiveSepCumTRWSVariable* var1_;
  AllInclusiveSepCumTRWSVariable* var2_;

  Storage1D<AllInclusiveSepCumTRWSFactor*> adjacent_factor_;

  Math2D::Matrix<double> pair_parameters_;
};

/******************/

//abstract base class for a factor
/*abstract*/ class AllInclusiveSepCumTRWSFactor {
public:

  AllInclusiveSepCumTRWSFactor(const Storage1D<AllInclusiveSepCumTRWSVariable*>& vars, 
			       const Storage1D<AllInclusiveSepCumTRWSPairSeparator*>& separators);

  virtual ~AllInclusiveSepCumTRWSFactor();
  
  const Math1D::Vector<double>& var_reparameterization(AllInclusiveSepCumTRWSVariable* var) const;
  
  const Math2D::Matrix<double>& pair_reparameterization(AllInclusiveSepCumTRWSPairSeparator* pair) const;

  virtual double compute_var_reparameterization(AllInclusiveSepCumTRWSVariable* var) = 0;

  virtual double compute_pair_reparameterization(AllInclusiveSepCumTRWSPairSeparator* pair) = 0;


  virtual double compute_forward(const AllInclusiveSepCumTRWSPairSeparator* incoming_sep, const AllInclusiveSepCumTRWSVariable* incoming_var, 
                                 const AllInclusiveSepCumTRWSVariable* outgoing_var,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math1D::Vector<double>& forward);

  virtual double compute_forward(const AllInclusiveSepCumTRWSPairSeparator* incoming_sep, const AllInclusiveSepCumTRWSVariable* incoming_var, 
                                 const AllInclusiveSepCumTRWSPairSeparator* outgoing_sep,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math2D::Matrix<double>& forward);

  virtual double best_value() = 0;

  void compute_rank_range();

  uint min_rank() const;

  uint max_rank() const;

  AllInclusiveSepCumTRWSFactor* prev_factor() const;

  AllInclusiveSepCumTRWSFactor* next_factor() const;

  void set_prev_factor(AllInclusiveSepCumTRWSFactor* prev_factor);

  void set_next_factor(AllInclusiveSepCumTRWSFactor* next_factor);

  double valid_offs() const;

  AllInclusiveSepCumTRWSPairSeparator* start_separator() const;

  AllInclusiveSepCumTRWSPairSeparator* end_separator() const;

  //returns if successful
  bool set_start_separator(AllInclusiveSepCumTRWSPairSeparator* sep);

  //returns if successful
  bool set_end_separator(AllInclusiveSepCumTRWSPairSeparator* sep);

  bool sep_is_interrupted(AllInclusiveSepCumTRWSPairSeparator* sep) const;

  void print_factor(bool short_form=false) const;

  const AllInclusiveSepCumTRWSPairSeparator* valid_separator() const;

  const Storage1D<AllInclusiveSepCumTRWSVariable*>& involved_var() const;

  const Storage1D<AllInclusiveSepCumTRWSPairSeparator*>& adjacent_separator() const;

protected:

  uint min_rank_;
  uint max_rank_;

  Storage1D<AllInclusiveSepCumTRWSVariable*> involved_var_;
  Storage1D<Math1D::Vector<double> > var_reparameterization_;

  Storage1D<AllInclusiveSepCumTRWSPairSeparator*> adjacent_separator_;
  Storage1D<Math2D::Matrix<double> > pair_reparameterization_;

  //DO WE REALLY NEED TO KEEP TRACK OF THESE??
  AllInclusiveSepCumTRWSFactor* prev_factor_;
  AllInclusiveSepCumTRWSFactor* next_factor_;

  uchar valid_sep_;
  double valid_offs_;

  Math1D::Vector<uchar> parent_sep_;
  
protected:
  AllInclusiveSepCumTRWSPairSeparator* start_separator_;
  AllInclusiveSepCumTRWSPairSeparator* end_separator_;
};

//binary factor
class BinaryAllInclusiveSepCumTRWSFactor : public AllInclusiveSepCumTRWSFactor {
public:

  BinaryAllInclusiveSepCumTRWSFactor(const Storage1D<AllInclusiveSepCumTRWSVariable*>& vars, 
                                     const Math2D::Matrix<float>& cost);

  virtual ~BinaryAllInclusiveSepCumTRWSFactor();
  
  virtual double compute_var_reparameterization(AllInclusiveSepCumTRWSVariable* var);

  virtual double compute_pair_reparameterization(AllInclusiveSepCumTRWSPairSeparator* pair);

  virtual double best_value();


  virtual double compute_forward(const AllInclusiveSepCumTRWSPairSeparator* incoming_sep, 
                                 const AllInclusiveSepCumTRWSVariable* incoming_var, 
                                 const AllInclusiveSepCumTRWSVariable* outgoing_var,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math1D::Vector<double>& forward);

protected:

  const Math2D::Matrix<float> cost_;
};

//ternary factor
class TernaryAllInclusiveSepCumTRWSFactor : public AllInclusiveSepCumTRWSFactor {
public:

  TernaryAllInclusiveSepCumTRWSFactor(const Storage1D<AllInclusiveSepCumTRWSVariable*>& vars, 
				      const Storage1D<AllInclusiveSepCumTRWSPairSeparator*>& separators,
				      const Math3D::Tensor<float>& cost);

  virtual ~TernaryAllInclusiveSepCumTRWSFactor();
  
  virtual double compute_var_reparameterization(AllInclusiveSepCumTRWSVariable* var);

  virtual double compute_pair_reparameterization(AllInclusiveSepCumTRWSPairSeparator* pair);

  virtual double best_value();


  virtual double compute_forward(const AllInclusiveSepCumTRWSPairSeparator* incoming_sep, const AllInclusiveSepCumTRWSVariable* incoming_var, 
                                 const AllInclusiveSepCumTRWSVariable* outgoing_var,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math1D::Vector<double>& forward);

  virtual double compute_forward(const AllInclusiveSepCumTRWSPairSeparator* incoming_sep, const AllInclusiveSepCumTRWSVariable* incoming_var, 
                                 const AllInclusiveSepCumTRWSPairSeparator* outgoing_sep,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math2D::Matrix<double>& forward);

protected:

  double eval_pair(uint pair_num, uint x, uint y, uint z) const;

  double eval_pair(uint pair_num, uint x, uint y, uint z, const Storage1D< Math2D::Matrix<double> >& pair_param) const;

  const Math3D::Tensor<float> cost_;
};

//4th order factor
class FourthOrderAllInclusiveSepCumTRWSFactor : public AllInclusiveSepCumTRWSFactor {
public:

  FourthOrderAllInclusiveSepCumTRWSFactor(const Storage1D<AllInclusiveSepCumTRWSVariable*>& vars, 
                                          const Storage1D<AllInclusiveSepCumTRWSPairSeparator*>& separators,
                                          const Storage1D<Math3D::Tensor<float> >& cost);

  virtual ~FourthOrderAllInclusiveSepCumTRWSFactor();
  
  virtual double compute_var_reparameterization(AllInclusiveSepCumTRWSVariable* var);

  virtual double compute_pair_reparameterization(AllInclusiveSepCumTRWSPairSeparator* pair);

  virtual double best_value();


  virtual double compute_forward(const AllInclusiveSepCumTRWSPairSeparator* incoming_sep, const AllInclusiveSepCumTRWSVariable* incoming_var, 
                                 const AllInclusiveSepCumTRWSVariable* outgoing_var,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math1D::Vector<double>& forward);

  virtual double compute_forward(const AllInclusiveSepCumTRWSPairSeparator* incoming_sep, const AllInclusiveSepCumTRWSVariable* incoming_var, 
                                 const AllInclusiveSepCumTRWSPairSeparator* outgoing_sep,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math2D::Matrix<double>& forward);

protected:

  double eval_pair(uint pair_num, uint x, uint y, uint z, uint w) const;

  double eval_pair(uint pair_num, uint x, uint y, uint z, uint w,
                   const Storage1D< Math2D::Matrix<double> >& pair_param) const;

  const Storage1D<Math3D::Tensor<float> > cost_;
};


/******************/

//this is the main class
class AllInclusiveSepCumTRWS {
public:

  AllInclusiveSepCumTRWS(uint nVars, uint nSeparators, uint nFactors);

  ~AllInclusiveSepCumTRWS();

  uint add_var(const Math1D::Vector<float>& cost);
  
  uint add_pair_separator(uint var1, uint var2);

  void add_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost);

  void add_ternary_factor(uint var1, uint var2, uint var3, const Storage1D<uint>& separators,
                          const Math3D::Tensor<float>& cost);

  void add_fourth_order_factor(uint var1, uint var2, uint var3, uint var4, const Storage1D<uint>& separators,
                               const Storage1D<Math3D::Tensor<float> >& cost);

  double optimize(uint nIter, bool quiet = false);

  const Math1D::Vector<uint>& labeling() const;

protected:

  void add_factor(AllInclusiveSepCumTRWSFactor* fac);

  double cur_bound(bool backward = false);

  Storage1D<AllInclusiveSepCumTRWSVariable*> var_;
  Storage1D<AllInclusiveSepCumTRWSPairSeparator*> separator_;
  Storage1D<AllInclusiveSepCumTRWSFactor*> factor_;

  Math1D::Vector<uint> labeling_;

  uint nUsedVars_;
  uint nUsedSeparators_;
  uint nUsedFactors_;

  bool optimize_called_;
};



#endif
