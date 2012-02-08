#ifndef __eim_classes_h__
#define __eim_classes_h__

// local includes
#include "rb_eim_construction.h"
#include "assembly.h"

// A simple subclass of RBEIMEvaluation. Overload
// evaluate_parametrized_function to define the
// function that we "empirically" interpolate.
class SimpleEIMEvaluation : public RBEIMEvaluation
{
public:

  SimpleEIMEvaluation()
  {
    attach_paramerized_function(&sg);
  }

  /** 
   * Parametrized function that we approximate with EIM
   */
  ShiftedGaussian sg;

};

// A simple subclass of RBEIMConstruction.
class SimpleEIMConstruction : public RBEIMConstruction
{
public:

  /**
   * Constructor.
   */
  SimpleEIMConstruction (EquationSystems& es,
                         const std::string& name,
                         const unsigned int number)
  : Parent(es, name, number)
  {
  }
  
  /**
   * Destructor.
   */
  ~SimpleEIMConstruction()
  {
    for(unsigned int i=0; i<rb_eim_f_vector.size(); i++)
    {
      delete rb_eim_f_vector[i];
    }
    rb_eim_f_vector.clear();
  }
  
  /**
   * The type of the parent.
   */
  typedef RBEIMConstruction Parent;
  
  /**
   * Initialize data structures.
   */
  virtual void init_data()
  {
    u_var = this->add_variable ("f_EIM", FIRST);

    Parent::init_data();

    attach_inner_prod_assembly(&ip);
  }

  /**
   * Build and store a vector of EIM_F objects,
   * one for each EIM basis function.
   */
  void initialize_EIM_F_objects()
  {
    // Initialize the EIM_F objects
    rb_eim_f_vector.clear();
    for(unsigned int i=0; i<rb_eval->get_n_basis_functions(); i++)
    {
      rb_eim_f_vector.push_back(new EIM_F(*this, i));
    }
  }

  /**
   * Variable number for u.
   */
  unsigned int u_var;

  /**
   * Inner product assembly object
   */
  EIM_IP_assembly ip;
  
  /**
   * The vector of EIM_F objects that are created to point to
   * this RBEIMConstruction.
   */
  std::vector<ElemAssembly*> rb_eim_f_vector;
  
};

#endif