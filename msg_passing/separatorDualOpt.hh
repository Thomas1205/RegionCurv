/**** written by Thomas Schoenemann as an employee of the University of Pisa, Italy, Oct. 2011 ****/
/**** and continued at the University of Düsseldorf, Germany, 2012 ***/


/***** implements the solving of dual relaxations with pairwise separators *****/
/**** via the MSD with higher order factors (up to 4th order) *****/


#ifndef SEPARATOR_DUAL_OPT_HH
#define SEPARATOR_DUAL_OPT_HH

#include "vector.hh"
#include "matrix.hh"
#include "tensor.hh"
#include "factorDualOpt.hh"

#include <set>

//forward declarations
class sepDualOptFactor;
class sepDualOptPairSeparator;

//variable class
class sepDualOptVar {
public:
  sepDualOptVar(const Math1D::Vector<float>& cost);

  void add_cost(const Math1D::Vector<float>& cost);

  void add_factor(sepDualOptFactor* adjacent_fac);
  
  void add_pair_separator(sepDualOptPairSeparator* adjacent_sep);

  double dual_value(uint& arg_min) const;

  uint nLabels();

  void compute_message(const sepDualOptFactor* factor, Math1D::Vector<double>& msg);

  void compute_message(const sepDualOptPairSeparator* sep, Math1D::Vector<double>& msg);

  const Math1D::Vector<float>& cost() const;

protected:
  
  Math1D::Vector<float> cost_;

  Storage1D<sepDualOptFactor*> adjacent_factor_; //contains only those factors that do not subsume any separator with this var.
  Storage1D<std::pair<sepDualOptVar*,sepDualOptPairSeparator*> > adjacent_separator_;
};

//class for a pair separator
//Note: one could have a derived class that also handles its own cost matrix
class sepDualOptPairSeparator {
public:
  
  sepDualOptPairSeparator(sepDualOptVar* var1, sepDualOptVar* var2);

  virtual ~sepDualOptPairSeparator();

  const Math1D::Vector<double>& dual_var(const sepDualOptVar* var) const;

  virtual void update_duals(DualBCAMode mode);

  void add_factor(sepDualOptFactor* adjacent_fac);

  std::set<sepDualOptVar*> involved_vars();

  double dual_value() const;

  void compute_message(const sepDualOptFactor* factor, Math2D::Matrix<double>& msg);

  sepDualOptVar* var1();

  sepDualOptVar* var2();

protected:

  sepDualOptVar* var1_;
  sepDualOptVar* var2_;

  Storage1D<Math1D::Vector<double> > dual_var_;
  Storage1D<sepDualOptFactor*> adjacent_factor_;
};

//abstract base class for a factor
/* abstract */ class sepDualOptFactor {
public:
  
  sepDualOptFactor(const Storage1D<sepDualOptVar*>& vars, const Storage1D<sepDualOptPairSeparator*>& separators,
                   bool minimal_links);

  virtual ~sepDualOptFactor();

  const Math1D::Vector<double>& dual_var(const sepDualOptVar* var) const;

  const Math2D::Matrix<double>& pair_dual(const sepDualOptPairSeparator* pair_sep) const;

  virtual void update_duals(DualBCAMode mode) = 0;

  virtual double dual_value() const = 0;

  const Storage1D<sepDualOptVar*>& vars() const;

  const Storage1D<sepDualOptPairSeparator*>& separators() const;

  virtual void write_cost(std::ostream& out, double factor = 1.0) = 0;

protected:

  Storage1D<sepDualOptVar*> var_;
  Storage1D<sepDualOptPairSeparator*> separator_;

  //if a var is contained in a separator, the respective vector will have length 0
  Storage1D<Math1D::Vector<double> > dual_var_; 
  Storage1D<Math2D::Matrix<double> > pair_dual_;
};

//binary factor
class BinarySepDualOptFactor : public sepDualOptFactor {
public:

  BinarySepDualOptFactor(const Storage1D<sepDualOptVar*>& vars, const Math2D::Matrix<float>& cost, bool minimal_links);

  virtual ~BinarySepDualOptFactor();

  virtual void update_duals(DualBCAMode mode);

  virtual double dual_value() const;

  virtual void write_cost(std::ostream& out, double factor = 1.0);

protected:

  const Math2D::Matrix<float> cost_;
};

//ternary factor
class TernarySepDualOptFactor : public sepDualOptFactor {
public:

  TernarySepDualOptFactor(const Storage1D<sepDualOptVar*>& vars, const Storage1D<sepDualOptPairSeparator*>& separators,
                          const Math3D::Tensor<float>& cost, bool minimal_links);

  virtual ~TernarySepDualOptFactor();

  virtual void update_duals(DualBCAMode mode);

  virtual double dual_value() const;

  virtual void write_cost(std::ostream& out, double factor = 1.0);

protected:

  double eval_pair(uint pair_num, uint x, uint y, uint z) const;

  const Math3D::Tensor<float> cost_;
};

// 4th order factor
class FourthOrderSepDualOptFactor : public sepDualOptFactor {
public:

  FourthOrderSepDualOptFactor(const Storage1D<sepDualOptVar*>& vars, const Storage1D<sepDualOptPairSeparator*>& separators,
                              const Storage1D<Math3D::Tensor<float> >& cost, bool minimal_links);

  virtual ~FourthOrderSepDualOptFactor();

  virtual void update_duals(DualBCAMode mode);

  virtual double dual_value() const;

  virtual void write_cost(std::ostream& out, double factor = 1.0);

protected:

  double eval_pair(uint pair_num, uint x, uint y, uint z, uint w) const;

  const Storage1D<Math3D::Tensor<float> > cost_;
};

//this is the main class
class SeparatorDualOptimization {
public:
  
  //you have to provide upper bounds on the number of variables, pair separators and factors you will add
  SeparatorDualOptimization(uint nVars, uint nSeparators, uint nFactors, bool minimal_links = true);

  ~SeparatorDualOptimization();

  uint add_var(const Math1D::Vector<float>& cost);

  uint add_separator(uint v1, uint v2);

  void add_generic_binary_factor(uint v1, uint v2, const Math2D::Matrix<float>& cost);

  void add_generic_ternary_factor(uint v1, uint v2, uint v3, const Storage1D<uint>& separators,
                                  const Math3D::Tensor<float>& cost);

  void add_fourth_order_factor(uint v1, uint v2, uint v3, uint v4,
                               const Storage1D<uint>& separators,
                               const Storage1D<Math3D::Tensor<float> >& cost);

  sepDualOptVar* get_variable(uint v);

  double optimize(uint nIter, DualBCAMode mode = DUAL_BCA_MODE_MSD, bool quiet = true);

  const Math1D::Vector<uint>& labeling();

  //write files required for CMPLP by Globerson et al.
  void save_problem();

protected:

  void add_factor(sepDualOptFactor* fac);

  Storage1D<sepDualOptVar*> var_;
  Storage1D<sepDualOptPairSeparator*> separator_;
  Storage1D<sepDualOptFactor*> factor_;  

  Math1D::Vector<uint> labeling_;

  uint nUsedVars_;
  uint nUsedSeparators_;
  uint nUsedFactors_;

  bool minimal_links_;
};

#endif
