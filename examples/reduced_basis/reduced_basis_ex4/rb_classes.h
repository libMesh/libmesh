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

#include "libmesh/rb_construction.h"
#include "libmesh/rb_evaluation.h"
#include "assembly.h"

#ifdef LIBMESH_ENABLE_DIRICHLET

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::DirichletBoundary;
using libMesh::RBConstruction;
using libMesh::RBEvaluation;

// A simple subclass of RBEvaluation.
class SimpleRBEvaluation : public RBEvaluation
{
public:

  /**
   * Constructor. Just set the theta expansion.
   */
  SimpleRBEvaluation(const libMesh::Parallel::Communicator & comm) :
    RBEvaluation(comm)
  {
    set_rb_theta_expansion(eim_test_rb_theta_expansion);
  }

  /**
   * The object that stores the "theta" expansion of the parameter dependent PDE,
   * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
   */
  EimTestRBThetaExpansion eim_test_rb_theta_expansion;
};

class SimpleRBConstruction : public RBConstruction
{
public:

  SimpleRBConstruction (EquationSystems & es,
                        const std::string & name_in,
                        const unsigned int number_in)
    : Parent(es, name_in, number_in),
      dirichlet_bc(std::unique_ptr<DirichletBoundary>())
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
    dirichlet_bc->b.insert(0);
    dirichlet_bc->b.insert(1);
    dirichlet_bc->b.insert(2);
    dirichlet_bc->b.insert(3);
    dirichlet_bc->variables.push_back(u_var);

    // Attach dirichlet_bc (must do this _before_ Parent::init_data)
    get_dof_map().add_dirichlet_boundary(*dirichlet_bc);

    Parent::init_data();

    // Set the rb_assembly_expansion for this Construction object.
    set_rb_assembly_expansion(eim_test_rb_assembly_expansion);

    // We need to define an inner product matrix for this problem
    set_inner_product_assembly(eim_test_rb_assembly_expansion.A0_assembly);
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
    side_fe->get_nothing();
  }

  /**
   * Variable number for u.
   */
  unsigned int u_var;

  /**
   * The object that stores the "theta" expansion of the parameter dependent PDE,
   * i.e. the set of parameter-dependent functions in the affine expansion of the PDE.
   */
  EimTestRBThetaExpansion eim_test_rb_theta_expansion;

  /**
   * The object that stores the "assembly" expansion of the parameter dependent PDE,
   * i.e. the objects that define how to assemble the set of parameter-independent
   * operators in the affine expansion of the PDE.
   */
  EimTestRBAssemblyExpansion eim_test_rb_assembly_expansion;

  /**
   * The object that defines which degrees of freedom are on a Dirichlet boundary.
   */
  std::unique_ptr<DirichletBoundary> dirichlet_bc;
};

#endif // LIBMESH_ENABLE_DIRICHLET

#endif
