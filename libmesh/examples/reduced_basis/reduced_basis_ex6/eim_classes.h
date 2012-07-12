#ifndef __eim_classes_h__
#define __eim_classes_h__

// local includes
#include "rb_eim_construction.h"
#include "rb_eim_evaluation.h"
#include "assembly.h"

// A simple subclass of RBEIMEvaluation. Overload
// evaluate_parametrized_function to define the
// function that we "empirically" interpolate.
class SimpleEIMEvaluation : public RBEIMEvaluation
{
public:

  SimpleEIMEvaluation()
  {
    attach_parametrized_function(&g_0);
    attach_parametrized_function(&g_1);
  }

  /** 
   * Parametrized functions that we approximate with EIM
   */
  G_0 g_0;
  G_1 g_1;

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
    return AutoPtr<ElemAssembly>(new EIM_A(*this, index));
  }
  
  /**
   * Initialize data structures.
   */
  virtual void init_data()
  {
    g0_var = this->add_variable ("G_0", FIRST);
    g1_var = this->add_variable ("G_1", FIRST);

    Parent::init_data();

    set_inner_product_assembly(eim_ip);
  }

  /**
   * Variable numbers.
   */
  unsigned int g0_var;
  unsigned int g1_var;

  /**
   * Inner product assembly object
   */
  Ex6EIMInnerProduct eim_ip;
  
};

#endif