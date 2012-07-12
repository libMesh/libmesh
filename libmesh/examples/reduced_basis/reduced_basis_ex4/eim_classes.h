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
    attach_parametrized_function(&sg);
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
   * The type of the parent.
   */
  typedef RBEIMConstruction Parent;
  
  /**
   * Provide an implementation of build_eim_assembly
   */
  virtual AutoPtr<ElemAssembly> build_eim_assembly(unsigned int index)
  {
    return AutoPtr<ElemAssembly>(new EIM_F(*this, index));
  }
  
  /**
   * Initialize data structures.
   */
  virtual void init_data()
  {
    u_var = this->add_variable ("f_EIM", FIRST);

    Parent::init_data();

    set_inner_product_assembly(ip);
  }

  /**
   * Variable number for u.
   */
  unsigned int u_var;

  /**
   * Inner product assembly object
   */
  EIM_IP_assembly ip;
  
};

#endif