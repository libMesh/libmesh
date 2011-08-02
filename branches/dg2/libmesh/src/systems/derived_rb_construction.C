// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic

// This file is part of rbOOmit.

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

#include "derived_rb_construction.h"
#include "derived_rb_evaluation.h"

#include "libmesh_logging.h"
#include "equation_systems.h"

namespace libMesh
{

template <class Base>
DerivedRBConstruction<Base>::DerivedRBConstruction (EquationSystems& es,
		                                    const std::string& name,
		                                    const unsigned int number)
  : Base(es, name, number)
{
  // We do not want to compute the output dual norms in
  // a derived system, we just copy them over from the
  // primary system
  Base::output_dual_norms_computed = true;
}

template <class Base>
std::string DerivedRBConstruction<Base>::system_type () const
{
  return "DerivedRBConstruction";
}

template <class Base>
Real DerivedRBConstruction<Base>::train_reduced_basis (const std::string& directory_name,
                                                       const bool resize_rb_eval_data)
{
  Real training_greedy_error = Base::train_reduced_basis(directory_name,
                                                         resize_rb_eval_data);
  
  generate_residual_terms_wrt_truth();
  
  return training_greedy_error;
}

template<class Base>
void DerivedRBConstruction<Base>::set_uber_current_parameters()
{
  EquationSystems& es = this->get_equation_systems();
  RBConstruction& uber_system = es.get_system<RBConstruction>(uber_system_name);

  uber_system.set_current_parameters( Base::get_current_parameters() );
}

template <class Base>
void DerivedRBConstruction<Base>::load_basis_function(unsigned int i)
{
  START_LOG("load_basis_function()", "DerivedRBConstruction");
  
  EquationSystems& es = Base::get_equation_systems();
  RBConstruction& uber_system = es.get_system<RBConstruction>(uber_system_name);

  DenseVector<Number> bf = get_derived_basis_function(i);

  for(unsigned int j=0; j<uber_system.rb_eval->get_n_basis_functions(); j++)
  {
    Base::solution->add(bf(j), uber_system.rb_eval->get_basis_function(j));
  }

  STOP_LOG("load_basis_function()", "DerivedRBConstruction");
}

// explicit instantiations
template class DerivedRBConstruction<RBConstruction>;

}
