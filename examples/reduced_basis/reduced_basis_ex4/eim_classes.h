#ifndef EIM_CLASSES_H
#define EIM_CLASSES_H

// libMesh includes
#include "libmesh/rb_eim_construction.h"
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// Example includes
#include "assembly.h"

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::EquationSystems;
using libMesh::RBEIMEvaluation;
#ifndef LIBMESH_HAVE_CXX14_MAKE_UNIQUE
using libMesh::make_unique;
#endif

// A simple subclass of RBEIMEvaluation. Overload
// evaluate_parametrized_function to define the
// function that we "empirically" interpolate.
class SimpleEIMEvaluation : public RBEIMEvaluation
{
public:

  SimpleEIMEvaluation(const libMesh::Parallel::Communicator & comm) :
    RBEIMEvaluation(comm)
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
  SimpleEIMConstruction (EquationSystems & es,
                         const std::string & name_in,
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
  virtual std::unique_ptr<ElemAssembly> build_eim_assembly(unsigned int index)
  {
    return libmesh_make_unique<EIM_F>(*this, index);
  }

  /**
   * Initialize data structures.
   */
  virtual void init_data()
  {
    Parent::init_data();

    set_inner_product_assembly(ip);
  }

  /**
   * Initialize the implicit system that is used to perform L2 projections.
   */
  virtual void init_implicit_system()
  {
    this->add_variable ("L2_proj_var", libMesh::FIRST);
  }

  /**
   * Initialize the explicit system that is used to store the basis functions.
   */
  virtual void init_explicit_system()
  {
    u_var = get_explicit_system().add_variable ("f_EIM", libMesh::FIRST);
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
