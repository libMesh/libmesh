#ifndef __eim_classes_h__
#define __eim_classes_h__

// local includes
#include "libmesh/rb_eim_construction.h"
#include "assembly.h"

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::EquationSystems;
using libMesh::RBEIMEvaluation;

// A simple subclass of RBEIMEvaluation. Overload
// evaluate_parametrized_function to define the
// function that we "empirically" interpolate.
class SimpleEIMEvaluation : public RBEIMEvaluation
{
public:

  SimpleEIMEvaluation(const libMesh::Parallel::Communicator& comm)
    : RBEIMEvaluation(comm)
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
    return UniquePtr<ElemAssembly>(new EIM_F(*this, index));
  }

  /**
   * Initialize data structures.
   */
  virtual void init_data()
  {
    u_var = this->add_variable ("f_EIM", libMesh::FIRST);

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
