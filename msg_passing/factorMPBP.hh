/******* written by Thomas Schoenemann as a private person without employment, July 2011 *****/


/***** implements max-product belief propagation in negative log-space ****/

#ifndef FACTOR_MPBP_HH
#define FACTOR_MPBP_HH

#include "vector.hh"
#include "matrix.hh"
#include "tensor.hh"
#include "vardim_storage.hh"

class FactorNode;

class VariableNode {
public:
  
  VariableNode(const Math1D::Vector<float>& cost);
  
  void add_factor(FactorNode* node);
  
  double* get_message(FactorNode* node);
  
  void compute_messages();
  
  double cost(uint label) const;

  void compute_beliefs(Math1D::Vector<double>& beliefs);
  
  uint nLabels() const;

  void init_messages();

  const Storage1D<FactorNode*>& neighboring_factor() const;

protected:  

  Storage1D<FactorNode*> neighboring_factor_; 
  
  Math2D::Matrix<double> message_matrix_;

  const Math1D::Vector<float> cost_;
};


/*abstract*/ class FactorNode {
public:
  
  FactorNode(const Storage1D<VariableNode*>& participating_vars);

  virtual ~FactorNode();
  
  virtual void compute_messages() = 0;
  
  virtual double* get_message(VariableNode* node) = 0;
  
  virtual double cost(const Math1D::Vector<uint>& labels) const = 0;

  virtual void init_messages() = 0;

  //factors can modify the current labeling, e.g. to satisfy some constraints
  virtual bool process_labeling(Math1D::Vector<uint>& labels) const;

  const Storage1D<VariableNode*>& participating_nodes() const;
  
protected:
  
  Storage1D<VariableNode*> participating_var_;
};

class GenericFactorNode : public FactorNode {
public:

  GenericFactorNode(const Storage1D<VariableNode*>& participating_vars, const VarDimStorage<float>& cost);

  virtual ~GenericFactorNode();

  virtual void compute_messages();
  
  virtual double* get_message(VariableNode* node);
  
  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual void init_messages();

protected:
  
  Storage1D<Math1D::Vector<double> > message_;
  VarDimStorage<float> cost_;
};

/* abstract */ class BinaryFactorNodeBase : public FactorNode {
public:

  BinaryFactorNodeBase(const Storage1D<VariableNode*>& participating_vars);

  void compute_messages(const Math2D::Matrix<float>& cost);

  virtual double* get_message(VariableNode* node);

  virtual void init_messages();

protected:

  Math1D::Vector<double> message1_; 
  Math1D::Vector<double> message2_; 
};

class BinaryFactorNode : public BinaryFactorNodeBase {
public:

  BinaryFactorNode(const Storage1D<VariableNode*>& participating_vars, const Math2D::Matrix<float>& cost);

  virtual ~BinaryFactorNode();

  virtual void compute_messages();

  virtual double cost(const Math1D::Vector<uint>& labels) const;

protected:
  Math2D::Matrix<float> cost_;
};


//like above, but only a reference to the cost matrix is kept
class BinaryRefFactorNode : public BinaryFactorNodeBase {
public:

  BinaryRefFactorNode(const Storage1D<VariableNode*>& participating_vars, const Math2D::Matrix<float>& cost);

  virtual ~BinaryRefFactorNode();

  virtual void compute_messages();

  virtual double cost(const Math1D::Vector<uint>& labels) const;

protected:
  const Math2D::Matrix<float>& cost_;
};

//a Potts Factor is always binary
class PottsFactorNode: public BinaryFactorNodeBase {
public:

  PottsFactorNode(const Storage1D<VariableNode*>& participating_vars, float lambda);

  virtual ~PottsFactorNode();

  virtual void compute_messages();

  virtual double cost(const Math1D::Vector<uint>& labels) const;

protected:
  float lambda_;
};


/*abstract*/ class TernaryFactorNodeBase : public FactorNode {
public: 

  TernaryFactorNodeBase(const Storage1D<VariableNode*>& participating_vars);

  void compute_messages(const Math3D::Tensor<float>& cost);

  virtual void init_messages();

  virtual double* get_message(VariableNode* node);

protected:

  Math1D::Vector<double> message1_; 
  Math1D::Vector<double> message2_; 
  Math1D::Vector<double> message3_;
};

class TernaryFactorNode : public TernaryFactorNodeBase {
public: 

  TernaryFactorNode(const Storage1D<VariableNode*>& participating_vars, const Math3D::Tensor<float>& cost);

  virtual ~TernaryFactorNode();

  virtual void compute_messages();

  virtual double cost(const Math1D::Vector<uint>& labels) const;

protected:
  Math3D::Tensor<float> cost_;
};

class TernaryRefFactorNode : public TernaryFactorNodeBase {
public: 

  TernaryRefFactorNode(const Storage1D<VariableNode*>& participating_vars, const Math3D::Tensor<float>& cost);

  virtual ~TernaryRefFactorNode();

  virtual void compute_messages();

  virtual double cost(const Math1D::Vector<uint>& labels) const;

protected:
  const Math3D::Tensor<float>&  cost_;
};



/*abstract*/ class FourthOrderFactorNodeBase : public FactorNode {
public: 

  FourthOrderFactorNodeBase(const Storage1D<VariableNode*>& participating_vars);

  void compute_messages(const Storage1D<Math3D::Tensor<float> >& cost);

  virtual void init_messages();

  virtual double* get_message(VariableNode* node);

protected:

  Math1D::Vector<double> message1_; 
  Math1D::Vector<double> message2_; 
  Math1D::Vector<double> message3_;
  Math1D::Vector<double> message4_;
};

class FourthOrderFactorNode : public FourthOrderFactorNodeBase {
public: 

  FourthOrderFactorNode(const Storage1D<VariableNode*>& participating_vars, const Storage1D<Math3D::Tensor<float> >& cost);

  virtual ~FourthOrderFactorNode();

  virtual void compute_messages();

  virtual double cost(const Math1D::Vector<uint>& labels) const;

protected:
  Storage1D<Math3D::Tensor<float> > cost_;
};

class FourthOrderRefFactorNode : public FourthOrderFactorNodeBase {
public: 

  FourthOrderRefFactorNode(const Storage1D<VariableNode*>& participating_vars, const Storage1D<Math3D::Tensor<float> >& cost);

  virtual ~FourthOrderRefFactorNode();

  virtual void compute_messages();

  virtual double cost(const Math1D::Vector<uint>& labels) const;

protected:
  const Storage1D<Math3D::Tensor<float> >& cost_;
};


//TODO: Generalized Potts Factor

//1-of-N constraints (a special case of cardinality potentials)
class OneOfNFactorNode : public FactorNode {
public:

  //all variables must be binary
  OneOfNFactorNode(const Storage1D<VariableNode*>& participating_vars);

  virtual ~OneOfNFactorNode();

  virtual void init_messages();

  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual double* get_message(VariableNode* node);

  virtual void compute_messages();

  virtual bool process_labeling(Math1D::Vector<uint>& labels) const;
  
  Math2D::Matrix<double> message_matrix_;
};

class CardinalityFactorNode: public FactorNode {
public:

  //all variables must be binary
  CardinalityFactorNode(const Storage1D<VariableNode*>& participating_vars, const Math1D::Vector<float>& card_cost);

  virtual ~CardinalityFactorNode();

  virtual void compute_messages();

  virtual double* get_message(VariableNode* node);

  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual void init_messages();

protected:
  Math1D::Vector<float> card_cost_;

  Math2D::Matrix<double> message_matrix_;
};


//implements an integer linear programming constraint with 1/-1 - coefficients
class BILPConstraintFactorNode: public FactorNode {
public:

  //all variables must be binary
  BILPConstraintFactorNode(const Storage1D<VariableNode*>& participating_vars,
                           const Storage1D<bool>& positive, int rhs_lower = 0, int rhs_upper = 0);

  virtual ~BILPConstraintFactorNode();

  virtual void compute_messages();

  virtual double* get_message(VariableNode* node);

  virtual double cost(const Math1D::Vector<uint>& labels) const;

  virtual void init_messages();

  virtual bool process_labeling(Math1D::Vector<uint>& labels) const;

protected:
  
  int rhs_lower_;
  int rhs_upper_;

  Storage1D<bool> positive_;

  Math2D::Matrix<double> message_matrix_;
};

class FactorMPBP {
public:

  FactorMPBP(uint nVars = 0, uint nFactors=0);

  //delete all allocated owned var. and factor nodes
  ~FactorMPBP(); 

  /***** Setting up a factor graph ******/

  //return var. number
  uint add_var(const Math1D::Vector<float>& unary_cost);

  //NOTE: after calling this routine, owned vars can no longer be created
  VariableNode* get_var_node(uint v);

  uint add_generic_factor(const Math1D::Vector<uint> var, VarDimStorage<float>& cost);

  //if you set ref=true, make sure that the cost object exists (unmodified) for as long as this object exists
  uint add_generic_binary_factor(uint var1, uint var2, const Math2D::Matrix<float>& cost, bool ref=false);
  
  uint add_potts_factor(uint var1, uint var2, double lambda);

  //return factor number
  //if you set ref=true, make sure that the cost object exists (unmodified) for as long as this object exists
  uint add_generic_ternary_factor(uint var1, uint var2, uint var3, const Math3D::Tensor<float>& cost, bool ref=false);

  //return factor number
  //if you set ref=true, make sure that the cost object exists (unmodified) for as long as this object exists
  uint add_generic_fourth_order_factor(uint var1, uint var2, uint var3, uint var4,
                                       const Storage1D<Math3D::Tensor<float> >& cost, bool ref=false);


  //all participating variables must be binary
  uint add_one_of_N_factor(const Math1D::Vector<uint>& var);

  //all participating variables must be binary
  uint add_cardinality_factor(const Math1D::Vector<uint>& var, const Math1D::Vector<float>& card_cost);

  //all participating variables must be binary
  uint add_binary_ilp_factor(const Math1D::Vector<uint>& var, const Storage1D<bool>& positive, 
                             int rhs_lower = 0, int rhs_upper = 0);

  //CAUTION: the factor is deleted together with all truly owned factors
  uint pass_in_factor_node(FactorNode* factor);

  FactorNode* get_factor(uint f);

  /**** run inference ***/

  void mpbp(uint nIter);

  /**** get the state of the solver ****/

  double labeling_energy();

  const Math1D::Vector<uint>& labeling();

  void icm(uint iter);

protected:

  void process_labeling();
  
  uint add_factor(FactorNode* node);

  Storage1D<VariableNode*> var_;
  Storage1D<FactorNode*> factor_;

  Math1D::Vector<uint> labeling_;

  uint nUsedVars_;
  uint nUsedFactors_;
};


#endif
