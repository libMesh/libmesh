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

#include "derived_rb_system.h"
#include "libmesh_logging.h"
#include "equation_systems.h"
#include "derived_rb_evaluation.h"

namespace libMesh
{

template <class Base>
DerivedRBSystem<Base>::DerivedRBSystem (EquationSystems& es,
		    const std::string& name,
		    const unsigned int number)
  : Base(es, name, number),
    residual_type_flag(RESIDUAL_WRT_UBER)
  {
    // We do not want to compute the output dual norms in
    // a derived system, we just copy them over from the
    // primary system
    Base::output_dual_norms_computed = true;
  }

template <class Base>
std::string DerivedRBSystem<Base>::system_type () const
{
  return "DerivedRBSystem";
}


template<class Base>
void DerivedRBSystem<Base>::set_uber_current_parameters()
{
  EquationSystems& es = this->get_equation_systems();
  RBSystem& uber_system = es.get_system<RBSystem>(uber_system_name);

  uber_system.set_current_parameters( Base::get_current_parameters() );
}

template <class Base>
void DerivedRBSystem<Base>::generate_residual_terms_wrt_truth()
{
  START_LOG("generate_residual_terms_wrt_truth()", "DerivedRBSystem");

  if(residual_type_flag != RESIDUAL_WRT_TRUTH)
  {
    // Set flag to compute residual wrt truth space
    residual_type_flag = RESIDUAL_WRT_TRUTH;

    // Need to recompute _all_ residual terms
    Base::update_residual_terms_called = false;

    unsigned int saved_delta_N = Base::delta_N;
    Base::delta_N = Base::get_n_basis_functions();
  
    // Recompute all the residual terms
    update_residual_terms();
  
    Base::delta_N = saved_delta_N;
  }
  STOP_LOG("generate_residual_terms_wrt_truth()", "DerivedRBSystem");
}

template <class Base>
void DerivedRBSystem<Base>::load_basis_function(unsigned int i)
{
  START_LOG("load_basis_function()", "DerivedRBSystem");

  if(!Base::initialize_mesh_dependent_data)
  {
    libMesh::err << "Error: We must initialize the mesh dependent "
                 << "data structures in order to load basis function."
                 << std::endl;
    libmesh_error();
  }

  libmesh_assert(i < Base::get_n_basis_functions());
  
  EquationSystems& es = Base::get_equation_systems();
  RBSystem& uber_system = es.get_system<RBSystem>(uber_system_name);

  DenseVector<Number> bf = get_derived_basis_function(i);

  for(unsigned int j=0; j<uber_system.get_n_basis_functions(); j++)
  {
    Base::solution->add(bf(j), uber_system.get_basis_function(j));
  }

  STOP_LOG("load_basis_function()", "DerivedRBSystem");
}

template <class Base>
void DerivedRBSystem<Base>::write_offline_data_to_files(const std::string& directory_name)
{
  generate_residual_terms_wrt_truth();
  
  Base::write_offline_data_to_files(directory_name);
}

template <class Base>
void DerivedRBSystem<Base>::clear_basis_function_dependent_data()
{
  // Clear the basis functions
  Base::rb_eval->set_n_basis_functions(0);
  
  Base::clear_basis_function_dependent_data();
  
  // On clearing we restore to a residual wrt the uber system
  residual_type_flag = RESIDUAL_WRT_UBER;
}

// explicit instantiations
template class DerivedRBSystem<RBSystem>;

}
