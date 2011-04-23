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

#ifndef __simple_rb_h__
#define __simple_rb_h__

#include "rb_system.h"
#include "fe_base.h"

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::EquationSystems;
using libMesh::FEMContext;
using libMesh::RBSystem;
using libMesh::Real;

// A simple subclass of RBSystem since we need to specialize
// get_SCM_lower_bound and get_SCM_upper_bound to return
// (a lower bound for) the coercivity constant for this problem,
// which is 0.05.

class SimpleRB : public RBSystem
{
public:

  SimpleRB (EquationSystems& es,
            const std::string& name,
            const unsigned int number)
  : Parent(es, name, number)
  {}

  /**
   * Destructor.
   */
  virtual ~SimpleRB () {}

  /**
   * The type of system.
   */
  typedef SimpleRB sys_type;

  /**
   * The type of the parent.
   */
  typedef RBSystem Parent;

  /**
   * The coercivity constant is bounded below by 0.05.
   */
  virtual Real get_SCM_lower_bound() { return 0.05; }

  /**
   * Use 0.05 here as well.
   */
  virtual Real get_SCM_upper_bound() { return 0.05; }

  /**
   * Initialize data structures.
   */
  virtual void init_data()
  {
    u_var = this->add_variable ("u", FIRST);

    Parent::init_data();
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

};

#endif
