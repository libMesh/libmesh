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

#ifndef __rb_classes_h__
#define __rb_classes_h__

#include "rb_construction.h"
#include "assembly.h"

// A simple subclass of RBEvaluation.
class SimpleRBEvaluation : public RBEvaluation
{
public:

  /**
   * Constructor. Just set the theta expansion.
   */
  SimpleRBEvaluation()
  {
    rb_theta_expansion = &eim_test_rb_theta_expansion;
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

  SimpleRBConstruction (EquationSystems& es,
                        const std::string& name,
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
    u_var = this->add_variable ("u", FIRST);

    Parent::init_data();

    // Attach rb_theta_expansion and rb_assembly_expansion
    // to this Construction object.
    // This also checks that the expansion objects are sized consistently
    attach_affine_expansion(eim_test_rb_theta_expansion,
                            eim_test_rb_assembly_expansion);


    // Attach the object that determines the Dirichlet boundary conditions for the PDE
    attach_dirichlet_dof_initialization(&eim_dirichlet_assembly);

    // We need to define an inner product matrix for this problem
    attach_inner_prod_assembly(&eim_test_rb_assembly_expansion.A0_assembly);
  }

  /**
   * Pre-request all relevant element data.
   */
  virtual void init_context(FEMContext &c)
  {
    // For efficiency, we should prerequest all
    // the data we will need to build the
    // linear system before doing an element loop.
    c.element_fe_var[u_var]->get_JxW();
    c.element_fe_var[u_var]->get_phi();
    c.element_fe_var[u_var]->get_dphi();
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
  EimDirichletDofAssembly eim_dirichlet_assembly;

};

#endif
