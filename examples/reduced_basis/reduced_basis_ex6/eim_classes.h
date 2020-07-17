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
    // FIXME: We can only attach a single RBParametrizedFunction by calling
    // set_parametrized_function(). Each RBParametrizedFunction can have multiple
    // components, though.
    // attach_parametrized_function(&g_x);
    // attach_parametrized_function(&g_y);
    // attach_parametrized_function(&g_z);
  }

  /**
   * Build a ThetaEIM rather than an RBEIMTheta.
   */
  virtual std::unique_ptr<RBTheta> build_eim_theta(unsigned int index)
  {
    return libmesh_make_unique<ThetaEIM>(*this, index);
  }

  /**
   * Parametrized functions that we approximate with EIM
   */
  // Gx g_x;
  // Gy g_y;
  // Gz g_z;
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
    return libmesh_make_unique<AssemblyEIM>(*this, index);
  }

  /**
   * Initialize data structures.
   */
  virtual void init_data()
  {
    Parent::init_data();

    set_inner_product_assembly(eim_ip);
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
    Gx_var = get_explicit_system().add_variable ("x_comp_of_G", libMesh::FIRST);
    Gy_var = get_explicit_system().add_variable ("y_comp_of_G", libMesh::FIRST);
    Gz_var = get_explicit_system().add_variable ("z_comp_of_G", libMesh::FIRST);
  }

  /**
   * Pre-request all relevant element data.
   */
  virtual void init_context(FEMContext & c)
  {
    // For efficiency, we should prerequest all
    // the data we will need to build the
    // linear system before doing an element loop.
    //
    // With equal variable types in x/y/z we only have the one elem_fe
    // to prerequest from.
    FEBase * elem_fe = nullptr;
    c.get_element_fe(Gx_var, elem_fe);

    elem_fe->get_JxW();
    elem_fe->get_phi();
    elem_fe->get_dphi();

    FEBase * side_fe = nullptr;
    c.get_side_fe(Gx_var, side_fe);
    side_fe->get_nothing();
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
