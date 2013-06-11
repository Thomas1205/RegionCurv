/******* written by Thomas Schoenemann as an employee of the University of Düsseldorf, Germany, 2012 *****/

/***** implements the solving of dual relaxations with pairwise separators *****/
/**** via the subgradient method for dual decomposition with factor (junction) chains *****/


#ifndef SEPARATOR_CHAIN_SUBGRADIENT_HH
#define SEPARATOR_CHAIN_SUBGRADIENT_HH

#include "vector.hh"
#include "matrix.hh"
#include "tensor.hh"

class SepChainDDFactor;
class SepChainDDPairSeparator;

// variable class
class SepChainDDVar {
public:

  SepChainDDVar(const Math1D::Vector<float>& cost);

  void add_cost(const Math1D::Vector<float>& cost);
  
  void add_factor(SepChainDDFactor* factor);
  
  void add_separator(SepChainDDPairSeparator* sep);

  const Math1D::Vector<float>& cost() const;
  
  uint nLabels() const;

  const Storage1D<SepChainDDFactor*>& neighboring_factor() const;

  uint nChains() const;

  void set_up_chains();

  //double dual_value(uint& arg_min);

protected:

  Math1D::Vector<float> cost_;
  
  Storage1D<SepChainDDFactor*> neighboring_factor_; 

  Storage1D<SepChainDDPairSeparator*> neighboring_separator_;
};

/***************************************/

// class for separators
class SepChainDDPairSeparator {
public:

  SepChainDDPairSeparator(SepChainDDVar* var1, SepChainDDVar* var2);

  SepChainDDVar* var1();

  SepChainDDVar* var2();

  const SepChainDDVar* var1() const;

  const SepChainDDVar* var2() const;

  void add_factor(SepChainDDFactor* factor);

  Storage1D<SepChainDDFactor*> neighboring_factor() const; 

protected:

  SepChainDDVar* var1_;
  
  SepChainDDVar* var2_;

  Storage1D<SepChainDDFactor*> neighboring_factor_; 
};

/***************************************/

// abstract base class for factors
/* abstract */ class SepChainDDFactor {
public:

  SepChainDDFactor(const Storage1D<SepChainDDVar*>& involved_vars, const Storage1D<SepChainDDPairSeparator*>& separators);

  virtual ~SepChainDDFactor();

  virtual double compute_forward(const SepChainDDPairSeparator* incoming_sep, const SepChainDDVar* incoming_var, 
                                 const SepChainDDVar* outgoing_var,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace) const = 0;


  virtual double compute_forward(const SepChainDDPairSeparator* incoming_sep, const SepChainDDVar* incoming_var, 
                                 const SepChainDDPairSeparator* outgoing_sep,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math2D::Matrix<double>& forward, 
                                 Math3D::Tensor<uint>& trace) const = 0;

  Math1D::Vector<double>& get_duals(const SepChainDDVar* var);

  Math2D::Matrix<double>& get_pair_duals(const SepChainDDPairSeparator* sep);

  SepChainDDVar* prev_var() const;

  SepChainDDVar* next_var() const;

  SepChainDDPairSeparator* prev_sep() const;

  SepChainDDPairSeparator* next_sep() const;

  SepChainDDFactor* prev_factor() const;

  SepChainDDFactor* next_factor() const;

  void set_prev_var(SepChainDDVar* var);

  void set_next_var(SepChainDDVar* var);

  void set_prev_sep(SepChainDDPairSeparator* sep);

  void set_next_sep(SepChainDDPairSeparator* sep);

  void set_prev_factor(SepChainDDFactor* factor);

  void set_next_factor(SepChainDDFactor* factor);

  const Storage1D<SepChainDDVar*>& involved_vars() const;

  const Storage1D<SepChainDDPairSeparator*>& involved_separators() const;

protected:

  SepChainDDVar* prev_var_;
  SepChainDDVar* next_var_;

  SepChainDDPairSeparator* prev_sep_;
  SepChainDDPairSeparator* next_sep_;

  SepChainDDFactor* prev_factor_;
  SepChainDDFactor* next_factor_;

  Storage1D<SepChainDDVar*> involved_var_;
  Storage1D<SepChainDDPairSeparator*> involved_separator_;

  Storage1D< Math1D::Vector<double> > dual_var_;
  Storage1D< Math2D::Matrix<double> > dual_pair_var_;
};

/***************************************/

//binary factor
class BinarySepChainDDFactor : public SepChainDDFactor {
public:
  BinarySepChainDDFactor(const Storage1D<SepChainDDVar*>& involved_vars,
                          const Math2D::Matrix<float>& cost);

  virtual ~BinarySepChainDDFactor();

  virtual double compute_forward(const SepChainDDPairSeparator* incoming_sep, const SepChainDDVar* incoming_var, 
                                 const SepChainDDVar* outgoing_var,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace) const;

  virtual double compute_forward(const SepChainDDPairSeparator* incoming_sep, const SepChainDDVar* incoming_var, 
                                 const SepChainDDPairSeparator* outgoing_sep,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math2D::Matrix<double>& forward, 
                                 Math3D::Tensor<uint>& trace) const;

protected:

  const Math2D::Matrix<float> cost_;
};

/***************************************/

//ternary factor
class TernarySepChainDDFactor : public SepChainDDFactor {
public:

  TernarySepChainDDFactor(const Storage1D<SepChainDDVar*>& involved_vars, Storage1D<SepChainDDPairSeparator*>& separators,
                          const Math3D::Tensor<float>& cost);

  virtual ~TernarySepChainDDFactor();

  virtual double compute_forward(const SepChainDDPairSeparator* incoming_sep, const SepChainDDVar* incoming_var, 
                                 const SepChainDDVar* outgoing_var,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace) const;

  virtual double compute_forward(const SepChainDDPairSeparator* incoming_sep, const SepChainDDVar* incoming_var, 
                                 const SepChainDDPairSeparator* outgoing_sep,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math2D::Matrix<double>& forward, 
                                 Math3D::Tensor<uint>& trace) const;
protected:

  const Math3D::Tensor<float> cost_;

  double eval_pair(uint s, uint l1, uint l2, uint l3, const Storage1D< Math2D::Matrix<double> >& pair_param) const;
};

/***************************************/

// 4th order factor
class FourthOrderSepChainDDFactor : public SepChainDDFactor {
public:

  FourthOrderSepChainDDFactor(const Storage1D<SepChainDDVar*>& involved_vars, Storage1D<SepChainDDPairSeparator*>& separators,
                              const Storage1D<Math3D::Tensor<float> >& cost);

  virtual ~FourthOrderSepChainDDFactor();

  virtual double compute_forward(const SepChainDDPairSeparator* incoming_sep, const SepChainDDVar* incoming_var, 
                                 const SepChainDDVar* outgoing_var,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math1D::Vector<double>& forward, 
                                 Math2D::Matrix<uint>& trace) const;

  virtual double compute_forward(const SepChainDDPairSeparator* incoming_sep, const SepChainDDVar* incoming_var, 
                                 const SepChainDDPairSeparator* outgoing_sep,
                                 const Math2D::Matrix<double>& prev_pair_forward, const Math1D::Vector<double>& prev_var_forward, 
                                 Math2D::Matrix<double>& forward, 
                                 Math3D::Tensor<uint>& trace) const;
protected:

  const Storage1D<Math3D::Tensor<float> > cost_;

  double eval_pair(uint s, uint l1, uint l2, uint l3, uint l4, const Storage1D< Math2D::Matrix<double> >& pair_param) const;
};

/***************************************/

//this is the main class
class SeparatorChainDualDecomposition {
public:
  
  //you need to provide upper bounds on the number of variables, pair-separators and factors you will add
  SeparatorChainDualDecomposition(uint nVars, uint nSeparators, uint nFactors);

  ~SeparatorChainDualDecomposition();

  uint add_var(const Math1D::Vector<float>& cost);

  uint add_pair_separator(uint var1, uint var2);

  void add_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost);

  void add_ternary_factor(uint var1, uint var2, uint var3, const Storage1D<uint>& separators, 
                          const Math3D::Tensor<float>& cost);

  void add_fourth_order_factor(uint var1, uint var2, uint var3, uint var4, const Storage1D<uint>& separators, 
                               const Storage1D<Math3D::Tensor<float> >& cost);

  double optimize(uint nIter, double start_step_size = 1.0);

  const Math1D::Vector<uint>& labeling();

protected:

  void add_factor(SepChainDDFactor* fac);

  void set_up_chains();

  Storage1D<SepChainDDVar*> var_;
  Storage1D<SepChainDDPairSeparator*> separator_;
  Storage1D<SepChainDDFactor*> factor_;

  Math1D::Vector<uint> labeling_;

  uint nUsedVars_;
  uint nUsedSeparators_;
  uint nUsedFactors_;

  bool optimize_called_;
};


#endif
