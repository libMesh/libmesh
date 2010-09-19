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

#include <sys/stat.h>

namespace libMesh
{

template <class Base>
DerivedRBSystem<Base>::DerivedRBSystem (EquationSystems& es,
		    const std::string& name,
		    const unsigned int number)
  : Base(es, name, number),
    residual_type_flag(RESIDUAL_WRT_UBER)
  {}

template <class Base>
std::string DerivedRBSystem<Base>::system_type () const
{
  return "DerivedRBSystem";
}

template <class Base>
void DerivedRBSystem<Base>::generate_residual_terms_wrt_truth()
{
  if(residual_type_flag != RESIDUAL_WRT_TRUTH)
  {
    // Set flag to compute residual wrt truth space
    residual_type_flag = RESIDUAL_WRT_TRUTH;

    // Need to recompute _all_ residual terms
    Base::update_residual_terms_called = false;

    unsigned int saved_delta_N = Base::delta_N;
    Base::delta_N = get_n_basis_functions();
  
    // Recompute all the residual terms
    update_residual_terms();
  
    Base::delta_N = saved_delta_N;
  }
}

template <class Base>
void DerivedRBSystem<Base>::load_basis_function(unsigned int i)
{
  START_LOG("load_basis_function()", "DerivedRBSystem");

  if(!Base::initialize_calN_dependent_data)
  {
    std::cerr << "Error: We must initialize the calN dependent "
              << "data structures in order to load basis function."
              << std::endl;
    libmesh_error();
  }

  libmesh_assert(i < get_n_basis_functions());
  
  EquationSystems& es = Base::get_equation_systems();
  RBSystem& uber_system = es.get_system<RBSystem>(uber_system_name);

  for(unsigned int j=0; j<uber_system.get_n_basis_functions(); j++)
  {
    Base::solution->add(derived_basis_functions[i](j), uber_system.get_bf(j));
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
  set_n_basis_functions(0);
  
  Base::clear_basis_function_dependent_data();
  
  // On clearing we restore to a residual wrt the uber system
  residual_type_flag = RESIDUAL_WRT_UBER;
}

template <class Base>
void DerivedRBSystem<Base>::write_out_basis_functions(const std::string& directory_name,
                                                      const unsigned int precision_level)
{
  if( Base::store_basis_functions && (libMesh::processor_id() == 0) ) // Only write out on proc 0
  {
    std::cout << "Writing out the basis functions..." << std::endl;

    std::ostringstream file_name;

    for(unsigned int i=0; i<derived_basis_functions.size(); i++)
    {
      file_name.str(""); // reset the string
      file_name << directory_name << "/derived_bf" << i << ".dat";
      std::ofstream derived_bf_out(file_name.str().c_str());
      
      if ( !derived_bf_out.good() )
      {
        std::cerr << "Error opening derived_bf" << i << ".dat" << std::endl;
        libmesh_error();
      }
      
      derived_bf_out.precision(precision_level);
      for(unsigned int j=0; j<derived_basis_functions[i].size(); j++)
      {
        derived_bf_out << std::scientific << derived_basis_functions[i](j) << " ";
      }
      derived_bf_out.close();
    }
    
    // Also, need to write out the size of the derived basis functions
    {
      std::ofstream derived_bf_size_out;
      {
        OStringStream file_name;
        file_name << directory_name << "/derived_bf_size.dat";
        derived_bf_size_out.open(file_name.str().c_str());
      }
      if ( !derived_bf_size_out.good() )
      {
        std::cerr << "Error opening derived_bf_size.dat" << std::endl;
        libmesh_error();
      }
      derived_bf_size_out << derived_basis_functions[0].size();
      derived_bf_size_out.close();
    }
  }
}

template <class Base>
void DerivedRBSystem<Base>::read_in_basis_functions(const std::string& directory_name)
{
  if(Base::store_basis_functions)
  {
    std::cout << "Reading in the basis functions..." << std::endl;
    
    // First, get the number of size of the derived basis functions
    unsigned int derived_bf_size;
    {
      OStringStream file_name;
      file_name << directory_name << "/derived_bf_size.dat";
      std::ifstream derived_bf_size_in(file_name.str().c_str());

      if ( !derived_bf_size_in.good() )
      {
        std::cerr << "Error opening derived_bf_size.dat" << std::endl;
        libmesh_error();
      }

      derived_bf_size_in >> derived_bf_size;
      derived_bf_size_in.close();
    }

    std::ostringstream file_name;
    struct stat stat_info;

    for(unsigned int i=0; i<derived_basis_functions.size(); i++)
    {
      file_name.str(""); // reset the string
      file_name << directory_name << "/derived_bf" << i << ".dat";

      // On processor zero check to be sure the file exists
      if (libMesh::processor_id() == 0)
	{
	  int stat_result = stat(file_name.str().c_str(), &stat_info);

	  if (stat_result != 0)
	    {
	      std::cout << "File does not exist: " << file_name.str() << std::endl;
	      libmesh_error();
	    }
	}

      derived_basis_functions[i].resize(derived_bf_size);

      std::ifstream derived_bf_size_in(file_name.str().c_str());
      
      for(unsigned int j=0; j<derived_bf_size; j++)
      {
        Number  value;
        derived_bf_size_in >> value;
        derived_basis_functions[i](j) = value;
      }
      derived_bf_size_in.close();

    }
  }
}

// explicit instantiations
template class DerivedRBSystem<RBSystem>;

}