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

#include "derived_rb_evaluation.h"
#include "rb_system.h"

#include <sys/stat.h>

namespace libMesh
{

template <class Base>
DerivedRBEvaluation<Base>::DerivedRBEvaluation(RBSystem& rb_sys)
  :
  Base(rb_sys)
{}

template <class Base>
unsigned int DerivedRBEvaluation<Base>::get_n_basis_functions() const
{
  return derived_basis_functions.size();
}

template <class Base>
void DerivedRBEvaluation<Base>::set_n_basis_functions(unsigned int n_bfs)
{
  derived_basis_functions.resize(n_bfs);
}

template <class Base>
void DerivedRBEvaluation<Base>::write_out_basis_functions(const std::string& directory_name,
                                                          const unsigned int precision_level)
{
  if( Base::rb_sys.store_basis_functions && (libMesh::processor_id() == 0) ) // Only write out on proc 0
  {
    libMesh::out << "Writing out the basis functions..." << std::endl;

    std::ostringstream file_name;

    for(unsigned int i=0; i<derived_basis_functions.size(); i++)
    {
      file_name.str(""); // reset the string
      file_name << directory_name << "/derived_bf" << i << ".dat";
      std::ofstream derived_bf_out(file_name.str().c_str());
      
      if ( !derived_bf_out.good() )
      {
        libMesh::err << "Error opening derived_bf" << i << ".dat" << std::endl;
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
        libMesh::err << "Error opening derived_bf_size.dat" << std::endl;
        libmesh_error();
      }
      derived_bf_size_out << derived_basis_functions[0].size();
      derived_bf_size_out.close();
    }
  }
}

template <class Base>
void DerivedRBEvaluation<Base>::read_in_basis_functions(const std::string& directory_name)
{
  if(Base::rb_sys.store_basis_functions)
  {
    libMesh::out << "Reading in the basis functions..." << std::endl;
    
    // First, get the number of size of the derived basis functions
    unsigned int derived_bf_size;
    {
      OStringStream file_name;
      file_name << directory_name << "/derived_bf_size.dat";
      std::ifstream derived_bf_size_in(file_name.str().c_str());

      if ( !derived_bf_size_in.good() )
      {
        libMesh::err << "Error opening derived_bf_size.dat" << std::endl;
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
	      libMesh::out << "File does not exist: " << file_name.str() << std::endl;
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

template class DerivedRBEvaluation<RBEvaluation>;

}