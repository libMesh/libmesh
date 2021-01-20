#ifndef RB_CLASSES_H
#define RB_CLASSES_H

// local includes
#include "assembly.h"

// rbOOmit includes
#include "libmesh/rb_construction.h"
#include "libmesh/rb_evaluation.h"

// libMesh includes
#include "libmesh/fe_base.h"
#include "libmesh/dof_map.h"

using namespace libMesh;

#ifdef LIBMESH_ENABLE_DIRICHLET

class ElasticityRBEvaluation : public RBEvaluation
{
public:

  /**
   * Constructor. Just set the theta expansion.
   */
  ElasticityRBEvaluation(const Parallel::Communicator & comm)
    : RBEvaluation(comm)
  {
    set_rb_theta_expansion(elasticity_theta_expansion);
  }

  /**
   * Return a "dummy" lower bound for the coercivity constant.
   * To do this rigorously we should use the SCM classes.
   */
  virtual Real get_stability_lower_bound() { return 1.; }

  /**
   * The object that stores the "theta" expansion of the parameter dependent PDE,
   * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
   */
  ElasticityThetaExpansion elasticity_theta_expansion;
};


class ElasticityRBConstruction : public RBConstruction
{
public:

  ElasticityRBConstruction (EquationSystems & es,
                            const std::string & name_in,
                            const unsigned int number_in) :
    Parent(es, name_in, number_in),
    u_var(0), v_var(0), w_var(0), var_order(FIRST),
    elasticity_assembly_expansion(*this),
    ip_assembly(*this)
  {}

  /**
   * Destructor.
   */
  virtual ~ElasticityRBConstruction () {}

  /**
   * The type of system.
   */
  typedef ElasticityRBConstruction sys_type;

  /**
   * The type of the parent.
   */
  typedef RBConstruction Parent;

  /**
   * Initialize data structures.
   */
  virtual void init_data()
  {
    unsigned int n_components = this->get_mesh().mesh_dimension();

    libmesh_assert_less_equal(n_components, 3);

    u_var = this->add_variable("u", var_order);
    if (n_components > 1)
      v_var = this->add_variable("v", var_order);
    if (n_components > 2)
      w_var = this->add_variable("w", var_order);

    // Generate a DirichletBoundary object
    dirichlet_bc = build_zero_dirichlet_boundary_object();

    // Set the Dirichlet boundary condition
    dirichlet_bc->b.insert(BOUNDARY_ID_MIN_X); // Dirichlet boundary at x=0
    dirichlet_bc->variables.push_back(u_var);
    if (n_components > 1)
      dirichlet_bc->variables.push_back(v_var);
    if (n_components > 2)
      dirichlet_bc->variables.push_back(w_var);

    // Attach dirichlet_bc (must do this _before_ Parent::init_data)
    get_dof_map().add_dirichlet_boundary(*dirichlet_bc);

    Parent::init_data();

    // Set the rb_assembly_expansion for this Construction object
    set_rb_assembly_expansion(elasticity_assembly_expansion);

    // We need to define an inner product matrix for this problem
    set_inner_product_assembly(ip_assembly);
  }

  /**
   * Pre-request all relevant element data.
   */
  virtual void init_context(FEMContext & c)
  {
    // For efficiency, we should prerequest all
    // the data we will need to build the
    // linear system before doing an element loop.
    FEBase * elem_fe = nullptr;
    c.get_element_fe(u_var, elem_fe);

    elem_fe->get_JxW();
    elem_fe->get_phi();
    elem_fe->get_dphi();

    FEBase * side_fe = nullptr;
    c.get_side_fe(u_var, side_fe);
    side_fe->get_JxW();
    side_fe->get_phi();
  }

  /**
   * Variable numbers and approximation order
   */
  unsigned int u_var;
  unsigned int v_var;
  unsigned int w_var;
  Order var_order;

  /**
   * The object that stores the "assembly" expansion of the parameter dependent PDE.
   */
  ElasticityAssemblyExpansion elasticity_assembly_expansion;

  /**
   * Object to assemble the inner product matrix
   */
  InnerProductAssembly ip_assembly;

  /**
   * The object that defines which degrees of freedom are on a Dirichlet boundary.
   */
  std::unique_ptr<DirichletBoundary> dirichlet_bc;
};

#endif // LIBMESH_ENABLE_DIRICHLET

#endif
