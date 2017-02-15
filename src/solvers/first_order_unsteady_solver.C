// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "libmesh/first_order_unsteady_solver.h"
#include "libmesh/diff_system.h"

namespace libMesh
{

void FirstOrderUnsteadySolver::init()
{
  // First call parent init
  UnsteadySolver::init();

  // Now check for second order variables and add their velocities to the System.
  const std::set<unsigned int> & second_order_vars = this->get_second_order_vars();
  if( !second_order_vars.empty() )
    {
      for( std::set<unsigned int>::const_iterator var_it = second_order_vars.begin();
           var_it != second_order_vars.end(); ++var_it )
        {
          const Variable & var = _system.variable(*var_it);
          std::string new_var_name = std::string("dot_")+var.name();

          unsigned int v_var_idx;

          if( var.active_subdomains().empty() )
            v_var_idx = _system.add_variable( new_var_name, var.type() );
          else
            v_var_idx = _system.add_variable( new_var_name, var.type(), &var.active_subdomains() );

          _second_order_dot_vars.insert( std::pair<unsigned int,unsigned int>(*var_it,v_var_idx) );
        }
    }
}

} // namespace libMesh
