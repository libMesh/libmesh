// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic
//
//     This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef RB_CLASSES_H
#define RB_CLASSES_H

#if defined(LIBMESH_HAVE_SLEPC) && defined(LIBMESH_HAVE_GLPK)

// libMesh includes
#include "libmesh/rb_construction.h"
#include "libmesh/rb_scm_construction.h"
#include "libmesh/fe_base.h"
#include "libmesh/rb_evaluation.h"
#include "libmesh/rb_scm_evaluation.h"

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::EquationSystems;
using libMesh::FEMContext;
using libMesh::RBConstruction;
using libMesh::RBSCMConstruction;
using libMesh::RBEvaluation;
using libMesh::RBSCMEvaluation;
using libMesh::RBParameters;
using libMesh::RBThetaExpansion;
using libMesh::RBAssemblyExpansion;
using libMesh::DirichletBoundary;
using libMesh::Real;

// local include
#include "assembly.h"

// A simple subclass of RBEvaluation. We also store the theta expansion object
// for the affine expansion of the PDE as a member variable.
class SimpleRBEvaluation : public RBEvaluation
{
public:

  /**
   * Constructor. Just set the theta expansion.
   */
  SimpleRBEvaluation(const libMesh::Parallel::Communicator & comm)
    : RBEvaluation(comm)
  {
    set_rb_theta_expansion(ex02_rb_theta_expansion);
  }

  /**
   * We override get_stability_lower_bound so that it calls rb_scm_eval to return
   * a parameter-dependent lower bound for the coercivity constant.
   */
  virtual Real get_stability_lower_bound()
  {
    rb_scm_eval->set_parameters(get_parameters());
    return rb_scm_eval->get_SCM_LB();
  }

  /**
   * Pointer to the SCM object that will provide our coercivity constant lower bound.
   */
  RBSCMEvaluation * rb_scm_eval;

  /**
   * The object that stores the "theta" expansion of the parameter dependent PDE,
   * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
   */
  Ex02RBThetaExpansion ex02_rb_theta_expansion;

};

// A simple subclass of RBConstruction, which initializes libMesh-related data such
// as the number of variables and their finite element type. We also store the objects
// that define the affine expansion of the PDE as member variables.
class SimpleRBConstruction : public RBConstruction
{
public:

  SimpleRBConstruction (EquationSystems & es,
                        const std::string & name,
                        const unsigned int number)
    : Parent(es, name, number)
  {}

  /**
   * Destructor.
   */
  virtual ~SimpleRBConstruction () {}

  /**
   * The type of system.
   */
  typedef SimpleRBConstruction sys_type;

  /**
   * The type of the parent.
   */
  typedef RBConstruction Parent;

  /**
   * Initialize data structures.
   */
  virtual void init_data()
  {
    u_var = this->add_variable ("u", libMesh::FIRST);

    // Generate a DirichletBoundary object
    dirichlet_bc = build_zero_dirichlet_boundary_object();

    // Set the Dirichlet boundary IDs
    // and the Dirichlet boundary variable numbers
    dirichlet_bc->b.insert(3);
    dirichlet_bc->variables.push_back(u_var);

    // Attach dirichlet_bc (must do this _before_ Parent::init_data)
    get_dof_map().add_dirichlet_boundary(*dirichlet_bc);

    Parent::init_data();

    // Set the rb_assembly_expansion for this Construction object.
    set_rb_assembly_expansion(ex02_rb_assembly_expansion);

    // We need to define an inner product matrix for this problem
    set_inner_product_assembly(ex02_rb_assembly_expansion.B_assembly);
  }

  /**
   * Pre-request all relevant element data. (This is not essential, but it
   * allows libMesh to cache data and hence run faster.)
   */
  virtual void init_context(FEMContext & c)
  {
    // For efficiency, we should prerequest all
    // the data we will need to build the
    // linear system before doing an element loop.
    FEBase * elem_fe = libmesh_nullptr;
    c.get_element_fe(u_var, elem_fe);

    elem_fe->get_JxW();
    elem_fe->get_phi();
    elem_fe->get_dphi();
  }

  /**
   * Variable number for u.
   */
  unsigned int u_var;

  /**
   * The object that stores the "assembly" expansion of the parameter dependent PDE,
   * i.e. the objects that define how to assemble the set of parameter-independent
   * operators in the affine expansion of the PDE.
   */
  Ex02RBAssemblyExpansion ex02_rb_assembly_expansion;

  /**
   * The object that defines which degrees of freedom are on a Dirichlet boundary.
   */
  std::unique_ptr<DirichletBoundary> dirichlet_bc;

};

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

#endif
