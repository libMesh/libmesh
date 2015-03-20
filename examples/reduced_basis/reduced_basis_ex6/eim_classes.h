#ifndef __eim_classes_h__
#define __eim_classes_h__

// local includes
#include "libmesh/rb_eim_construction.h"
#include "libmesh/rb_eim_evaluation.h"
#include "assembly.h"

// A simple subclass of RBEIMEvaluation. Overload
// evaluate_parametrized_function to define the
// function that we "empirically" interpolate.
class SimpleEIMEvaluation : public RBEIMEvaluation
{
public:

  SimpleEIMEvaluation(const libMesh::Parallel::Communicator& comm)
    : RBEIMEvaluation(comm)
  {
    attach_parametrized_function(&g_x);
    attach_parametrized_function(&g_y);
    attach_parametrized_function(&g_z);
  }

  /**
   * Build a ThetaEIM rather than an RBEIMTheta.
   */
  virtual UniquePtr<RBTheta> build_eim_theta(unsigned int index)
  {
    return UniquePtr<RBTheta>(new ThetaEIM(*this, index));
  }

  /**
   * Parametrized functions that we approximate with EIM
   */
  Gx g_x;
  Gy g_y;
  Gz g_z;

};

// A simple subclass of RBEIMConstruction.
class SimpleEIMConstruction : public RBEIMConstruction
{
public:

  /**
   * Constructor.
   */
  SimpleEIMConstruction (EquationSystems& es,
                         const std::string& name_in,
                         const unsigned int number_in)
    : Parent(es, name_in, number_in)
  {
  }

  /**
   * The type of the parent.
   */
  typedef RBEIMConstruction Parent;

  /**
   * Provide an implementation of build_eim_assembly
   */
  virtual UniquePtr<ElemAssembly> build_eim_assembly(unsigned int index)
  {
    return UniquePtr<ElemAssembly>(new AssemblyEIM(*this, index));
  }

  /**
   * Initialize data structures.
   */
  virtual void init_data()
  {
    Gx_var = this->add_variable ("x_comp_of_G", libMesh::FIRST);
    Gy_var = this->add_variable ("y_comp_of_G", libMesh::FIRST);
    Gz_var = this->add_variable ("z_comp_of_G", libMesh::FIRST);

    Parent::init_data();

    set_inner_product_assembly(eim_ip);
  }

  /**
   * Variable numbers.
   */
  unsigned int Gx_var;
  unsigned int Gy_var;
  unsigned int Gz_var;

  /**
   * Inner product assembly object
   */
  Ex6EIMInnerProduct eim_ip;

};

#endif
