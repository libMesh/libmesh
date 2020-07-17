#ifndef EIM_CLASSES_H
#define EIM_CLASSES_H

// libMesh includes
#include "libmesh/rb_eim_construction.h"
#include "libmesh/rb_eim_evaluation.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// Example includes
#include "assembly.h"

using libMesh::RBEIMConstruction;

#ifndef LIBMESH_HAVE_CXX14_MAKE_UNIQUE
using libMesh::make_unique;
#endif

// A simple subclass of RBEIMEvaluation. Overload
// evaluate_parametrized_function to define the
// function that we "empirically" interpolate.
class SimpleEIMEvaluation : public RBEIMEvaluation
{
public:

  SimpleEIMEvaluation(const libMesh::Parallel::Communicator & comm)
    : RBEIMEvaluation(comm)
  {
    set_parametrized_function(std::make_unique<Gxyz>());
  }

  /**
   * Build a ThetaEIM rather than an RBEIMTheta.
   */
  virtual std::unique_ptr<RBTheta> build_eim_theta(unsigned int index)
  {
    return libmesh_make_unique<ThetaEIM>(*this, index);
  }
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
    : RBEIMConstruction(es, name_in, number_in)
  {
  }

  /**
   * Provide an implementation of build_eim_assembly
   */
  virtual std::unique_ptr<ElemAssembly> build_eim_assembly(unsigned int index)
  {
    return libmesh_make_unique<AssemblyEIM>(get_rb_eim_evaluation(), index);
  }

  /**
   * Initialize data structures.
   */
  virtual void init_data()
  {
    this->add_variable ("eim_var", libMesh::FIRST);

    RBEIMConstruction::init_data();
  }

};

#endif
