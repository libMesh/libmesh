#ifndef EIM_CLASSES_H
#define EIM_CLASSES_H

// libMesh includes
#include "libmesh/rb_eim_construction.h"
#include "libmesh/rb_eim_evaluation.h"

// Example includes
#include "assembly.h"

// C++ includes
#include <memory>

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::EquationSystems;
using libMesh::RBEIMEvaluation;
using libMesh::RBEIMConstruction;

// A simple subclass of RBEIMEvaluation. Overload
// evaluate_parametrized_function to define the
// function that we "empirically" interpolate.
class SimpleEIMEvaluation : public RBEIMEvaluation
{
public:

  SimpleEIMEvaluation(const libMesh::Parallel::Communicator & comm) :
    RBEIMEvaluation(comm)
  {
    set_parametrized_function(std::make_unique<ShiftedGaussian>());
  }

  virtual bool use_eim_error_indicator() const override
  {
    // Indicate that we do use the EIM error indicator here.
    return true;
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
   * Initialize data structures.
   */
  virtual void init_data()
  {
    this->add_variable ("eim_var", libMesh::FIRST);

    RBEIMConstruction::init_data();
  }

  /**
   * Provide an implementation of build_eim_assembly
   */
  virtual std::unique_ptr<ElemAssembly> build_eim_assembly(unsigned int index)
  {
    return std::make_unique<EIM_F>(*this, index);
  }
};

#endif
