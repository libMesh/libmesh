// $Id: mesh_gmv_support.C,v 1.32 2004-03-24 05:49:12 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



// C/C++ includes
#include <fstream>

// Local includes
#include "equation_systems.h"


void MeshBase::write_discontinuous_gmv (const std::string& name,
					const EquationSystems& es,
					const bool write_partitioning) const
{
  std::vector<std::string> solution_names;
  std::vector<Number>      v;
  
  es.build_variable_names  (solution_names);
  es.build_discontinuous_solution_vector (v);
  
  if (this->processor_id() != 0) return;
  
  std::ofstream out(name.c_str());
  
  assert (out.good());

  // Begin interfacing with the GMV data file
  {

    // write the nodes    
    out << "gmvinput ascii" << std::endl << std::endl;
    
    // Compute the total weight
    {
      const_active_elem_iterator       it (this->elements_begin());
      const const_active_elem_iterator end(this->elements_end());

      unsigned int tw=0;
      
      for ( ; it != end; ++it)
	tw += (*it)->n_nodes();
      
      out << "nodes " << tw << std::endl;
    }
    


    // Write all the x values
    {
      const_active_elem_iterator       it (this->elements_begin());
      const const_active_elem_iterator end(this->elements_end());
      
      for ( ; it != end; ++it)   
	for (unsigned int n=0; n<(*it)->n_nodes(); n++)
	  out << (*it)->point(n)(0) << " ";

      out << std::endl;
    }
    
    
    // Write all the y values
    {
      const_active_elem_iterator       it (this->elements_begin());
      const const_active_elem_iterator end(this->elements_end());
      
      for ( ; it != end; ++it)   
	for (unsigned int n=0; n<(*it)->n_nodes(); n++)
	  out << (*it)->point(n)(1) << " ";

      out << std::endl;
    }
    
    
    // Write all the z values
    {
      const_active_elem_iterator       it (this->elements_begin());
      const const_active_elem_iterator end(this->elements_end());
      
      for ( ; it != end; ++it)   
	for (unsigned int n=0; n<(*it)->n_nodes(); n++)
	  out << (*it)->point(n)(2) << " ";

      out << std::endl << std::endl;
    }
  }


  
  {
    // write the connectivity
    
    out << "cells " << this->n_active_elem() << std::endl;

    const_active_elem_iterator       it (this->elements_begin());
    const const_active_elem_iterator end(this->elements_end());

    unsigned int nn=1;
    
    switch (this->mesh_dimension())
      {
      case 1:
	{
	  error();
	  break;
	}
	
      case 2:
	{
	  for ( ; it != end; ++it)
	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	      {
		if (((*it)->type() == QUAD4) ||
		    ((*it)->type() == QUAD8) ||
		    ((*it)->type() == QUAD9)
#ifdef ENABLE_INFINITE_ELEMENTS
		    || ((*it)->type() == INFQUAD4)
		    || ((*it)->type() == INFQUAD6)
#endif
		    )
		  {
		    out << "quad 4" << std::endl;
		    for (unsigned int i=0; i<(*it)->n_nodes(); i++)
		      out << nn++ << " ";
		    
		  }
		  else if (((*it)->type() == TRI3) ||
			   ((*it)->type() == TRI6))
		    {
		      out << "tri 3" << std::endl;
		    for (unsigned int i=0; i<(*it)->n_nodes(); i++)
		      out << nn++ << " ";
		    
		    }
		else
		  {
		    error();
		  }
		
		out << std::endl;
	      }
	  
	  break;
	}
	
	
      case 3:
	{
	  error();
	  break;
	}
	
      default:
	error();
      }
    
    out << std::endl;
  }
  

  
  // optionally write the partition information
  if (write_partitioning)
    {
      out << "material "
	  << this->n_processors()
	  << " 0"<< std::endl;

      for (unsigned int proc=0; proc<this->n_processors(); proc++)
	out << "proc_" << proc << std::endl;
      
      const_active_elem_iterator       it (this->elements_begin());
      const const_active_elem_iterator end(this->elements_end());

      for ( ; it != end; ++it)
	out << (*it)->processor_id()+1 << std::endl;
      
      out << std::endl;
    }


  
  // write the data
  {
    const unsigned int n_vars = solution_names.size();
      
    //    assert (v.size() == tw*n_vars);

    out << "variable" << std::endl;


      for (unsigned int c=0; c<n_vars; c++)
	{

#ifdef USE_COMPLEX_NUMBERS

	  // in case of complex data, write _tree_ data sets
	  // for each component

	  // this is the real part
	  out << "r_" << solution_names[c] << " 1" << std::endl;
	  {
	    const_active_elem_iterator       it (this->elements_begin());
	    const const_active_elem_iterator end(this->elements_end());
	    
	    for ( ; it != end; ++it)
	      for (unsigned int n=0; n<(*it)->n_nodes(); n++)
		out << std::setprecision(10) << v[(n++)*n_vars + c].real() << " ";	    
	  }	  	  
	  out << std::endl << std::endl;


	  // this is the imaginary part
	  out << "i_" << solution_names[c] << " 1" << std::endl;	  
	  {
	    const_active_elem_iterator       it (this->elements_begin());
	    const const_active_elem_iterator end(this->elements_end());
	    
	    for ( ; it != end; ++it)
	      for (unsigned int n=0; n<(*it)->n_nodes(); n++)
		out << std::setprecision(10) << v[(n++)*n_vars + c].imag() << " ";	    
	  }	  
	  out << std::endl << std::endl;

	  // this is the magnitude
	  out << "a_" << solution_names[c] << " 1" << std::endl;
	  {
	    const_active_elem_iterator       it (this->elements_begin());
	    const const_active_elem_iterator end(this->elements_end());
	    
	    for ( ; it != end; ++it)
	      for (unsigned int n=0; n<(*it)->n_nodes(); n++)
		out << std::setprecision(10)
		    << std::abs(v[(n++)*n_vars + c]) << " ";	    
	  }
	  out << std::endl << std::endl;

#else

	  out << solution_names[c] << " 1" << std::endl;
	  {
	    const_active_elem_iterator       it (this->elements_begin());
	    const const_active_elem_iterator end(this->elements_end());
	    
	    unsigned int nn=0;
	  
	    for ( ; it != end; ++it)
	      for (unsigned int n=0; n<(*it)->n_nodes(); n++)
		out << std::setprecision(10) << v[(nn++)*n_vars + c] << " ";	    
	  }	  
	  out << std::endl << std::endl;

#endif

	}
      
      out << "endvars" << std::endl;
    }

  
  // end of the file
  out << std::endl << "endgmv" << std::endl;
}




