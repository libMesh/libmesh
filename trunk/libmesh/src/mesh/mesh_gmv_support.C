// $Id: mesh_gmv_support.C,v 1.18 2003-05-23 00:30:38 ddreyer Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



// C++ includes
#include <iomanip>
#include <fstream>
#include <stdio.h>

// Local includes
#include "equation_systems.h"
#include "mesh_base.h"
#include "elem.h"



void MeshBase::write_gmv (const std::string& name,
			  EquationSystems& es,
			  const bool write_partitioning)
{
  std::ofstream out(name.c_str());

  write_gmv (out, es, write_partitioning);
}



void MeshBase::write_gmv (std::ostream& out,
			  EquationSystems& es,
			  const bool write_partitioning)
{
  std::vector<Number> soln;
  std::vector<std::string> names;
  
  es.build_variable_names  (names);
  es.build_solution_vector (soln);

  if (processor_id() == 0)
    write_gmv (out, &soln, &names, write_partitioning);
}



void MeshBase::write_gmv (const std::string& name,
			  const std::vector<Number>* v,
			  const std::vector<std::string>* solution_names,
			  const bool write_partitioning)
{
  std::ofstream out(name.c_str());

  write_gmv (out, v, solution_names, write_partitioning);
}




void MeshBase::write_gmv(std::ostream& out,
			 const std::vector<Number>*  v,
			 const std::vector<std::string>* solution_names,
			 const bool write_partitioning)
{
  assert (out);

  // Begin interfacing with the GMV data file
  {
    // write the nodes
    
    out << "gmvinput ascii" << std::endl << std::endl;
    out << "nodes " << n_nodes() << std::endl;
    for (unsigned int v=0; v<n_nodes(); v++)
      out << point(v)(0) << " ";
         
    out << std::endl;
    
    for (unsigned int v=0; v<n_nodes(); v++)
      out << point(v)(1) << " ";
    
    out << std::endl;
    
    for (unsigned int v=0; v<n_nodes(); v++)
      out << point(v)(2) << " ";
     
    out << std::endl << std::endl;
  }


  
  {
    // write the connectivity
    
    out << "cells " << n_active_sub_elem() << std::endl;

    const_active_elem_iterator       it (elements_begin());
    const const_active_elem_iterator end(elements_end());
    
    switch (_dim)
      {
      case 1:
	{
	  for ( ; it != end; ++it)
	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	      {
		out << "line 2" << std::endl;
		std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
		for (unsigned int i=0; i<conn.size(); i++)
		  out << conn[i] << " ";
		
		out << std::endl;
	      }
	  
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
		    std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
		    for (unsigned int i=0; i<conn.size(); i++)
		      out << conn[i] << " ";
		  }
		  else if (((*it)->type() == TRI3) ||
			   ((*it)->type() == TRI6))
		    {
		      out << "tri 3" << std::endl;
		      std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
		      for (unsigned int i=0; i<3; i++)
			out << conn[i] << " ";
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
	  for ( ; it != end; ++it)
	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	      {

#ifndef  ENABLE_INFINITE_ELEMENTS
		if (((*it)->type() == HEX8)   ||    
		    ((*it)->type() == HEX27))
		  {
		    out << "phex8 8" << std::endl;
		    std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
		    for (unsigned int i=0; i<conn.size(); i++)
		      out << conn[i] << " ";
		  }
		
		else if ((*it)->type() == HEX20)
		  {
		    out << "phex20 20" << std::endl;
		    out << (*it)->node(0)+1  << " "
			<< (*it)->node(1)+1  << " "
			<< (*it)->node(2)+1  << " "
			<< (*it)->node(3)+1  << " "
			<< (*it)->node(4)+1  << " "
			<< (*it)->node(5)+1  << " "
			<< (*it)->node(6)+1  << " "
			<< (*it)->node(7)+1  << " "
			<< (*it)->node(8)+1  << " "
			<< (*it)->node(9)+1  << " "
			<< (*it)->node(10)+1 << " "
			<< (*it)->node(11)+1 << " "
			<< (*it)->node(16)+1 << " "
			<< (*it)->node(17)+1 << " "
			<< (*it)->node(18)+1 << " "
			<< (*it)->node(19)+1 << " "
			<< (*it)->node(12)+1 << " "
			<< (*it)->node(13)+1 << " "
			<< (*it)->node(14)+1 << " "
			<< (*it)->node(15)+1 << " ";
		  }
#else
		/*
		 * In case of infinite elements, HEX20
		 * should be handled just like the
		 * INFHEX16, since these connect to each other
		 */
		if (((*it)->type() == HEX8)     ||
		    ((*it)->type() == HEX27)    ||
		    ((*it)->type() == INFHEX8)  ||
		    ((*it)->type() == INFHEX16) ||
		    ((*it)->type() == INFHEX18) ||
		    ((*it)->type() == HEX20))
		  {
		    out << "phex8 8" << std::endl;
		    std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
		    for (unsigned int i=0; i<conn.size(); i++)
		      out << conn[i] << " ";
		  }
#endif
		
		else if (((*it)->type() == TET4)  ||
			 ((*it)->type() == TET10))
		  {
		    out << "tet 4" << std::endl;
		    std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
		    out << conn[0] << " "
			<< conn[2] << " "
			<< conn[1] << " "
			<< conn[4] << " ";
		    }

		else if (((*it)->type() == PRISM6) ||
			 ((*it)->type() == INFPRISM6))
		  {
		    /**
		     * Note that here PRISM6 is treted as
		     * a degenerated phex8.
		     */
		    out << "phex8 8" << std::endl;
		    std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
		    for (unsigned int i=0; i<conn.size(); i++)
		      out << conn[i] << " ";
		  }
		
		else
		  {
		    error();
		  }
		
		out << std::endl;
	      }
	  
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
	  << n_processors()
	  << " 0"<< std::endl;

      for (unsigned int sbd=0; sbd<n_processors(); sbd++)
	out << "sbd_" << sbd << std::endl;
      
      const_active_elem_iterator       it (elements_begin());
      const const_active_elem_iterator end(elements_end());

      for ( ; it != end; ++it)
	for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	    out << (*it)->processor_id()+1 << std::endl;
      
      out << std::endl;
    }


  
  // optionally write the data
  if ((solution_names != NULL) &&
      (v != NULL))
    {      
      const unsigned int n_vars = solution_names->size();
      
      assert (v->size() == n_nodes()*n_vars);

      out << "variable" << std::endl;


      for (unsigned int c=0; c<n_vars; c++)
	{

#ifdef USE_COMPLEX_NUMBERS

	  // in case of complex data, write _tree_ data sets
	  // for each component

	  // this is the real part
	  out << "r_" << (*solution_names)[c] << " 1" << std::endl;
	  
	  for (unsigned int n=0; n<n_nodes(); n++)
	    out << std::setprecision(10) << (*v)[n*n_vars + c].real() << " ";
	  
	  out << std::endl << std::endl;


	  // this is the imaginary part
	  out << "i_" << (*solution_names)[c] << " 1" << std::endl;
	  
	  for (unsigned int n=0; n<n_nodes(); n++)
	    out << std::setprecision(10) << (*v)[n*n_vars + c].imag() << " ";
	  
	  out << std::endl << std::endl;

	  // this is the magnitude
	  out << "a_" << (*solution_names)[c] << " 1" << std::endl;
	  for (unsigned int n=0; n<n_nodes(); n++)
	    out << std::setprecision(10)
		<< std::abs((*v)[n*n_vars + c]) << " ";

	  out << std::endl << std::endl;

#else

	  out << (*solution_names)[c] << " 1" << std::endl;
	  
	  for (unsigned int n=0; n<n_nodes(); n++)
	    out << std::setprecision(10) << (*v)[n*n_vars + c] << " ";
	  
	  out << std::endl << std::endl;

#endif

	}
      
      out << "endvars" << std::endl;
    }

  
  // end of the file
  out << std::endl << "endgmv" << std::endl;
}




void MeshBase::write_gmv_binary (const std::string& name,
				 EquationSystems& es,
				 const bool write_partitioning)
{
  std::ofstream out(name.c_str());

  write_gmv_binary (out, es, write_partitioning);
}



void MeshBase::write_gmv_binary (std::ostream& out,
				 EquationSystems& es,
				 const bool write_partitioning)
{
  std::vector<Number> soln;
  std::vector<std::string> names;

  es.build_variable_names  (names);
  es.build_solution_vector (soln);
  
  if (processor_id() == 0)
    write_gmv_binary (out, &soln, &names, write_partitioning);
}



void MeshBase::write_gmv_binary (const std::string& name,
				 const std::vector<Number>* v,
				 const std::vector<std::string>* solution_names,
				 const bool write_partitioning)
{
  std::ofstream out (name.c_str(), std::ios::out|std::ios::binary);
  
  write_gmv_binary (out, v, solution_names, write_partitioning);
}




void MeshBase::write_gmv_binary(std::ostream& out,
				const std::vector<Number>* v,
				const std::vector<std::string>* solution_names,
				const bool write_partitioning)
{
  assert (out);

  char buf[80];

  // Begin interfacing with the GMV data file
  {
    // write the nodes
    strcpy(buf, "gmvinput");
    out.write(buf, strlen(buf));
	      
    strcpy(buf, "ieeei4r4");
    out.write(buf, strlen(buf));

  }


  
  // write the nodes
  {
    strcpy(buf, "nodes   ");
    out.write(buf, strlen(buf));

    unsigned int tempint = n_nodes();
    
    memcpy(buf, &tempint, sizeof(unsigned int));

    out.write(buf, sizeof(unsigned int));

    // write the x coordinate
    float *temp = new float[n_nodes()];
    for (unsigned int v=0; v<n_nodes(); v++)
      temp[v] = static_cast<float>(point(v)(0));
    out.write(reinterpret_cast<char *>(temp), sizeof(float)*n_nodes());

    for (unsigned int v=0; v<n_nodes(); v++)
      temp[v] = static_cast<float>(point(v)(1));
    out.write(reinterpret_cast<char *>(temp), sizeof(float)*n_nodes());

    for (unsigned int v=0; v<n_nodes(); v++)
      temp[v] = static_cast<float>(point(v)(2));
    out.write(reinterpret_cast<char *>(temp), sizeof(float)*n_nodes());

    delete [] temp;
  }


  
  // write the connectivity
  {
    strcpy(buf, "cells   ");
    out.write(buf, strlen(buf));

    unsigned int tempint = n_active_sub_elem();
    
    memcpy(buf, &tempint, sizeof(unsigned int));
    
    out.write(buf, sizeof(unsigned int));

    const_active_elem_iterator       it (elements_begin());
    const const_active_elem_iterator end(elements_end());

    switch (_dim)
      {

      case 1:
       
	for ( ; it != end; ++it)
	  for(unsigned se = 0; se < (*it)->n_sub_elem(); ++se)
	    {
	      strcpy(buf, "line    ");
	      out.write(buf, strlen(buf));
	      
	      tempint = 2;
	      memcpy(buf, &tempint, sizeof(unsigned int));
	      out.write(buf, sizeof(unsigned int));
	      
	      std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
	      
	      out.write(reinterpret_cast<char*>(&conn[0]), sizeof(unsigned int)*tempint);
	    }
	
	break;
	
      case 2:
       
      for ( ; it != end; ++it)
	for(unsigned se = 0; se < (*it)->n_sub_elem(); ++se)
	  {
	    strcpy(buf, "quad    ");
	    out.write(buf, strlen(buf));
	    
	    tempint = 4;
	    memcpy(buf, &tempint, sizeof(unsigned int));
	    out.write(buf, sizeof(unsigned int));
		
	    std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
	    
	    out.write(reinterpret_cast<char*>(&conn[0]), sizeof(unsigned int)*tempint);
	  }
	
	break;
	
      case 3:
       
      for ( ; it != end; ++it)
	for(unsigned se = 0; se < (*it)->n_sub_elem(); ++se)
	  {
	    strcpy(buf, "phex8   ");
	    out.write(buf, strlen(buf));
	    
	    tempint = 8;
	    memcpy(buf, &tempint, sizeof(unsigned int));
	    out.write(buf, sizeof(unsigned int));
	    
	    std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
	    
	    out.write(reinterpret_cast<char*>(&conn[0]), sizeof(unsigned int)*tempint);
	  }
      
	break;
	
      default:
	error();
	
      }
  }
  
  
  
  // optionally write the partition information
  if (write_partitioning)
    {
      strcpy(buf, "material");
      out.write(buf, strlen(buf));
      
      unsigned int tmpint = n_processors();
      memcpy(buf, &tmpint, sizeof(unsigned int));
      out.write(buf, sizeof(unsigned int));

      tmpint = 0; // IDs are cell based
      memcpy(buf, &tmpint, sizeof(unsigned int));
      out.write(buf, sizeof(unsigned int));


      for (unsigned int sbd=0; sbd<n_processors(); sbd++)
	{
	  sprintf(buf, "sbd_%d", sbd);
	  out.write(buf, 8);
	}

      unsigned int* sbd_id = new unsigned int[n_active_sub_elem()];
      
      unsigned int n=0;
      
      const_active_elem_iterator       it (elements_begin());
      const const_active_elem_iterator end(elements_end());

      for ( ; it != end; ++it)
	for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	  sbd_id[n++] = (*it)->processor_id()+1;
      
      
      out.write(reinterpret_cast<char *>(sbd_id), sizeof(unsigned int)*n_active_sub_elem());

      delete [] sbd_id;
    }


  
  // optionally write the data
  if ((solution_names != NULL) &&
      (v != NULL))
    {
      strcpy(buf, "variable");
      out.write(buf, strlen(buf));
      
      float *temp = new float[n_nodes()];

      const unsigned int n_vars = solution_names->size();
      
      for (unsigned int c=0; c<n_vars; c++)
	{

#ifdef USE_COMPLEX_NUMBERS
	  // for complex data, write two datasets


	  // Real part
	  strcpy(buf, "r_");
	  out.write(buf, 2);
	  strcpy(buf, (*solution_names)[c].c_str());
	  out.write(buf, 6);
	  
	  unsigned int tempint = 1; // always do nodal data
	  memcpy(buf, &tempint, sizeof(unsigned int));
	  out.write(buf, sizeof(unsigned int));
	  
	  for (unsigned int n=0; n<n_nodes(); n++)
	    temp[n] = static_cast<float>( (*v)[n*n_vars + c].real() );
	  
	  out.write(reinterpret_cast<char *>(temp), sizeof(float)*n_nodes());


	  // imaginary part
	  strcpy(buf, "i_");
	  out.write(buf, 2);
	  strcpy(buf, (*solution_names)[c].c_str());
	  out.write(buf, 6);
	  
	  memcpy(buf, &tempint, sizeof(unsigned int));
	  out.write(buf, sizeof(unsigned int));
	  
	  for (unsigned int n=0; n<n_nodes(); n++)
	    temp[n] = static_cast<float>( (*v)[n*n_vars + c].imag() );
	  
	  out.write(reinterpret_cast<char *>(temp), sizeof(float)*n_nodes());


#else


	  strcpy(buf, (*solution_names)[c].c_str());
	  out.write(buf, 8);
	  
	  unsigned int tempint = 1; // always do nodal data
	  memcpy(buf, &tempint, sizeof(unsigned int));
	  out.write(buf, sizeof(unsigned int));
	  
	  for (unsigned int n=0; n<n_nodes(); n++)
	    temp[n] = static_cast<float>((*v)[n*n_vars + c]);
	  
	  out.write(reinterpret_cast<char *>(temp), sizeof(float)*n_nodes());


#endif

	  
	}
    
      delete [] temp;
      
      strcpy(buf, "endvars ");
      out.write(buf, strlen(buf));

    }

  // end the file
  strcpy(buf, "endgmv  ");
  out.write(buf, strlen(buf));
}




