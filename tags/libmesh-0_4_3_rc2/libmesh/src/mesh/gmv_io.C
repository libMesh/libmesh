// $Id: gmv_io.C,v 1.12 2004-10-19 12:44:10 benkirk Exp $

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



// C++ includes
#include <fstream>
#include <stdio.h>  // for sprintf

// Local includes
#include "libmesh_config.h"
#include "gmv_io.h"



// ------------------------------------------------------------
// GMVIO  members
void GMVIO::write (const std::string& fname)
{
  if (libMesh::processor_id() == 0)
    if (this->binary())
      this->write_binary (fname);
    else
      this->write_ascii  (fname);
}



void GMVIO::write_nodal_data (const std::string& fname,
			      const std::vector<Number>& soln,
			      const std::vector<std::string>& names)
{
  if (libMesh::processor_id() == 0)
    if (this->binary())
      this->write_binary (fname, &soln, &names);
    else
      this->write_ascii  (fname, &soln, &names);
}



void GMVIO::write_ascii (const std::string& fname,
			 const std::vector<Number>* v,
			 const std::vector<std::string>* solution_names)
{
  // Open the output file stream
  std::ofstream out (fname.c_str());
  
  assert (out.good());

  // Get a reference to the mesh
  const MeshBase& mesh = this->cmesh();
  
  // Begin interfacing with the GMV data file
  {
    // write the nodes
    
    out << "gmvinput ascii\n\n";
    out << "nodes " << mesh.n_nodes() << '\n';
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      out << mesh.point(v)(0) << " ";
         
    out << '\n';
    
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      out << mesh.point(v)(1) << " ";
    
    out << '\n';
    
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      out << mesh.point(v)(2) << " ";
     
    out << '\n' << '\n';
  }


  
  {
    // write the connectivity
    
    out << "cells " << mesh.n_active_sub_elem() << '\n';

    const_active_elem_iterator       it (mesh.elements_begin());
    const const_active_elem_iterator end(mesh.elements_end());
    
    switch (mesh.mesh_dimension())
      {
      case 1:
	{
	  // The same temporary storage will be used for each element
	  std::vector<unsigned int> conn;

	  for ( ; it != end; ++it)
	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	      {
		out << "line 2\n";
		(*it)->connectivity(se, TECPLOT, conn);
		for (unsigned int i=0; i<conn.size(); i++)
		  out << conn[i] << " ";
		
		out << '\n';
	      }
	  break;
	}
	
      case 2:
	{
	  // The same temporary storage will be used for each element
	  std::vector<unsigned int> conn;
	  
	  for ( ; it != end; ++it)
	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	      {
		// Quad elements
		if (((*it)->type() == QUAD4) ||
		    ((*it)->type() == QUAD8) ||
		    ((*it)->type() == QUAD9)
#ifdef ENABLE_INFINITE_ELEMENTS
		    || ((*it)->type() == INFQUAD4)
		    || ((*it)->type() == INFQUAD6)
#endif
		    )
		  {
		    out << "quad 4\n";
		    (*it)->connectivity(se, TECPLOT, conn);
		    for (unsigned int i=0; i<conn.size(); i++)
		      out << conn[i] << " ";
		  }

		// Triangle elements
		else if (((*it)->type() == TRI3) ||
			 ((*it)->type() == TRI6))
		  {
		    out << "tri 3\n";
		    (*it)->connectivity(se, TECPLOT, conn);
		    for (unsigned int i=0; i<3; i++)
		      out << conn[i] << " ";
		  }
		
		else
		  {
		    error();
		  }
		
		out << '\n';
	      }
	  
	  break;
	}
	
	
      case 3:
	{
	  // The same temporary storage will be used for each element
	  std::vector<unsigned int> conn;
	  
	  for ( ; it != end; ++it)
	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	      {

#ifndef  ENABLE_INFINITE_ELEMENTS
		if (((*it)->type() == HEX8)   ||    
		    ((*it)->type() == HEX27))
		  {
		    out << "phex8 8\n";
		    (*it)->connectivity(se, TECPLOT, conn);
		    for (unsigned int i=0; i<conn.size(); i++)
		      out << conn[i] << " ";
		  }
		
		else if ((*it)->type() == HEX20)
		  {
		    out << "phex20 20\n";
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
		    out << "phex8 8\n";
		    (*it)->connectivity(se, TECPLOT, conn);
		    for (unsigned int i=0; i<conn.size(); i++)
		      out << conn[i] << " ";
		  }
#endif
		
		else if (((*it)->type() == TET4)  ||
			 ((*it)->type() == TET10))
		  {
		    out << "tet 4\n";
		    (*it)->connectivity(se, TECPLOT, conn);
		    out << conn[0] << " "
			<< conn[2] << " "
			<< conn[1] << " "
			<< conn[4] << " ";
		    }
#ifndef  ENABLE_INFINITE_ELEMENTS
		else if (((*it)->type() == PRISM6)  ||
			 ((*it)->type() == PRISM15) ||
			 ((*it)->type() == PRISM18))
#else
		else if (((*it)->type() == PRISM6)     ||
			 ((*it)->type() == PRISM15)    ||
			 ((*it)->type() == PRISM18)    ||
			 ((*it)->type() == INFPRISM6)  ||
			 ((*it)->type() == INFPRISM12))
#endif
		  {
		    /**
		     * Note that the prisms are treated as
		     * degenerated phex8's.
		     */
		    out << "phex8 8\n";
		    (*it)->connectivity(se, TECPLOT, conn);
		    for (unsigned int i=0; i<conn.size(); i++)
		      out << conn[i] << " ";
		  }
		
		else
		  {
		    error();
		  }
		
		out << '\n';
	      }
	  
	  break;
	}
	
      default:
	error();
      }
    
    out << '\n';
  }
  

  
  // optionally write the partition information
  if (this->partitioning())
    {
      out << "material "
	  << mesh.n_partitions()
	  << " 0"<< '\n';

      for (unsigned int proc=0; proc<mesh.n_partitions(); proc++)
	out << "proc_" << proc << '\n';
      
      const_active_elem_iterator       it (mesh.elements_begin());
      const const_active_elem_iterator end(mesh.elements_end());

      for ( ; it != end; ++it)
	for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	    out << (*it)->processor_id()+1 << '\n';
      
      out << '\n';
    }


  
  // optionally write the data
  if ((solution_names != NULL) &&
      (v != NULL))
    {      
      const unsigned int n_vars = solution_names->size();

      if (!(v->size() == mesh.n_nodes()*n_vars))
	std::cerr << "ERROR: v->size()=" << v->size()
		  << ", mesh.n_nodes()=" << mesh.n_nodes()
		  << ", n_vars=" << n_vars
		  << ", mesh.n_nodes()*n_vars=" << mesh.n_nodes()*n_vars
		  << std::endl;
      
      assert (v->size() == mesh.n_nodes()*n_vars);

      out << "variable\n";


      for (unsigned int c=0; c<n_vars; c++)
	{

#ifdef USE_COMPLEX_NUMBERS

	  // in case of complex data, write _tree_ data sets
	  // for each component

	  // this is the real part
	  out << "r_" << (*solution_names)[c] << " 1\n";
	  
	  for (unsigned int n=0; n<mesh.n_nodes(); n++)
	    out << std::setprecision(10) << (*v)[n*n_vars + c].real() << " ";
	  
	  out << '\n' << '\n';


	  // this is the imaginary part
	  out << "i_" << (*solution_names)[c] << " 1\n";
	  
	  for (unsigned int n=0; n<mesh.n_nodes(); n++)
	    out << std::setprecision(10) << (*v)[n*n_vars + c].imag() << " ";
	  
	  out << '\n' << '\n';

	  // this is the magnitude
	  out << "a_" << (*solution_names)[c] << " 1\n";
	  for (unsigned int n=0; n<mesh.n_nodes(); n++)
	    out << std::setprecision(10)
		<< std::abs((*v)[n*n_vars + c]) << " ";

	  out << '\n' << '\n';

#else

	  out << (*solution_names)[c] << " 1\n";
	  
	  for (unsigned int n=0; n<mesh.n_nodes(); n++)
	    out << std::setprecision(10) << (*v)[n*n_vars + c] << " ";
	  
	  out << '\n' << '\n';

#endif
	}
      
      out << "endvars\n";
    }

  
  // end of the file
  out << '\n' << "endgmv\n";
}



void GMVIO::write_binary (const std::string& fname,
			  const std::vector<Number>* vec,
			  const std::vector<std::string>* solution_names)
{
  std::ofstream out (fname.c_str(), std::ofstream::binary);
  
  assert (out.good());

  // get a reference to the mesh
  const MeshBase& mesh = this->cmesh();

  // Temporary buffer for storage
  char buf[80];


  // Begin interfacing with the GMV data file
  {
    std::string title = "gmvinput";
    out.write(title.c_str(), title.size());

    std::string format = "ieeei4r4";
    out.write(format.c_str(), format.size());
  }


  
  // write the nodes
  {
    std::string nodestr = "nodes   ";
    out.write(nodestr.c_str(), nodestr.size());

    // Write the number of nodes to the binary stream
    this->to_binary_stream(out, mesh.n_nodes());
    
    // write the x coordinates.  
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      this->to_binary_stream(out, static_cast<float>(mesh.point(v)(0)));
    
    // write the y coordinates.
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      this->to_binary_stream(out, static_cast<float>(mesh.point(v)(1)));

    // write the z coordinates.
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      this->to_binary_stream(out, static_cast<float>(mesh.point(v)(2)));
  }


  
  // write the connectivity
  {
    std::string cellstr = "cells   ";
    out.write(cellstr.c_str(), cellstr.size());

    // Write the number of active sub-elements
    this->to_binary_stream(out, mesh.n_active_sub_elem());

    
    const_active_elem_iterator       it (mesh.elements_begin());
    const const_active_elem_iterator end(mesh.elements_end());
    switch (mesh.mesh_dimension())
      {

      case 1:
	{
	  // The same temporary storage will be used for each element
	  std::vector<unsigned int> conn;
	  
	  // This string will be output for each element
	  std::string linestr = "line    ";
	  
	  for ( ; it != end; ++it)
	    for(unsigned se = 0; se < (*it)->n_sub_elem(); ++se)
	      {
		out.write(linestr.c_str(), linestr.size());

		// Write the number 2 to the file
		this->to_binary_stream(out, 2);
		
		(*it)->connectivity(se, TECPLOT, conn);

		for (unsigned int i=0; i<conn.size(); ++i)
		  this->to_binary_stream(out, conn[i]);
	      }
	
	  break;
	}
	
      case 2:
	{
	  // The same temporary storage will be used for each element
	  std::vector<unsigned int> conn;

	  // This string will be output for each element
	  std::string quadstr = "quad    ";
	  
	  for ( ; it != end; ++it)
	    for(unsigned se = 0; se < (*it)->n_sub_elem(); ++se)
	      {
		out.write(quadstr.c_str(), quadstr.size());

		// Write the number 4 to the file
		this->to_binary_stream(out, 4);
		
		(*it)->connectivity(se, TECPLOT, conn);

		for (unsigned int i=0; i<conn.size(); ++i)
		  this->to_binary_stream(out, conn[i]);
	      }
	
	  break;
	}
	
      case 3:
	{
	  // The same temporary storage will be used for each element
	  std::vector<unsigned int> conn;

	  // This string will be output for each element
	  std::string phexstr = "phex8   ";
	  
	  for ( ; it != end; ++it)
	    for(unsigned se = 0; se < (*it)->n_sub_elem(); ++se)
	      {
		out.write(phexstr.c_str(), phexstr.size());

		// Write the number 8 to the file
		this->to_binary_stream(out, 8);
		
		(*it)->connectivity(se, TECPLOT, conn);

		for (unsigned int i=0; i<conn.size(); ++i)
		  this->to_binary_stream(out, conn[i]);
	      }
      
	  break;
	}
	
      default:
	error();
	
      }
  }
  
  
  
  // optionally write the partition information
  if (this->partitioning())
    {
      std::string matstr = "material";
      out.write(matstr.c_str(), matstr.size());

      // Write the number of processors to the file
      this->to_binary_stream(out, mesh.n_processors());
      
      // IDs are cell based, so write the number 0 to the file
      this->to_binary_stream(out, 0);

      for (unsigned int proc=0; proc<mesh.n_processors(); proc++)
	{
	  sprintf(buf, "proc_%d", proc);
	  out.write(buf, 8);
	}

      std::vector<unsigned int> proc_id (mesh.n_active_sub_elem());
      
      unsigned int n=0;
      
      const_active_elem_iterator       it (mesh.elements_begin());
      const const_active_elem_iterator end(mesh.elements_end());

      for ( ; it != end; ++it)
	for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	  proc_id[n++] = (*it)->processor_id()+1;

      for (unsigned int i=0; i<proc_id.size(); ++i)
	this->to_binary_stream(out, proc_id[i]);
    }


  
  // optionally write the data
  if ((solution_names != NULL) &&
      (vec != NULL))
    {
      std::string varstr = "variable";
      out.write(varstr.c_str(), varstr.size());
      
      std::vector<float> temp (mesh.n_nodes());
      
      const unsigned int n_vars = solution_names->size();
      
      for (unsigned int c=0; c<n_vars; c++)
	{

#ifdef USE_COMPLEX_NUMBERS
	  // for complex data, write three datasets


	  // Real part
	  std::string rstr = "r_";
	  out.write(rstr.c_str(), rstr.size());
	  out.write((*solution_names)[c].c_str(), 6); // Write a max of 8 characters
	  
	  // Always do nodal data, so write the number 1 to file
	  this->to_binary_stream(out, 1);

	  // Get the real part into temporary storage
	  for (unsigned int n=0; n<mesh.n_nodes(); n++)
	    temp[n] = static_cast<float>( (*vec)[n*n_vars + c].real() );

	  for (unsigned int i=0; i<mesh.n_nodes(); ++i)
	    this->to_binary_stream(out, temp[i]);

	  
	  // Imaginary part
	  std::string istr = "i_";
	  out.write(istr.c_str(), istr.size());
	  out.write((*solution_names)[c].c_str(), 6); // Write a max of 8 characters

	  // Always do nodal data, so write the number 1 to file
	  this->to_binary_stream(out, 1);

	  // Get the imag part into temporary storage
	  for (unsigned int n=0; n<mesh.n_nodes(); n++)
	    temp[n] = static_cast<float>( (*vec)[n*n_vars + c].imag() );

	  for (unsigned int i=0; i<mesh.n_nodes(); ++i)
	    this->to_binary_stream(out, temp[i]);



	  
	  // Magnitude
	  std::string astr = "a_";
	  out.write(astr.c_str(), astr.size());
	  out.write((*solution_names)[c].c_str(), 6); // Write a max of 8 characters

	  // Always do nodal data, so write the number 1 to file
	  this->to_binary_stream(out, 1);

	  // Get the magnitude into temporary storage
	  for (unsigned int n=0; n<mesh.n_nodes(); n++)
	    temp[n] = static_cast<float>(std::abs((*vec)[n*n_vars + c]));

	  for (unsigned int i=0; i<mesh.n_nodes(); ++i)
	    this->to_binary_stream(out, temp[i]);

#else


	  out.write((*solution_names)[c].c_str(), 8); // Write a max of 8 characters
	  
	  // Always do nodal data, so write the number 1 to file
	  this->to_binary_stream(out, 1);

	  // Get the magnitudes into the temporary vector
	  for (unsigned int n=0; n<mesh.n_nodes(); n++)
	    temp[n] = static_cast<float>((*vec)[n*n_vars + c]);

	  for (unsigned int i=0; i<mesh.n_nodes(); ++i)
	    this->to_binary_stream(out, temp[i]);


#endif

	  
	}
      std::string endvarstr = "endvars ";
      out.write(endvarstr.c_str(), endvarstr.size());

    }

  // end the file
  std::string endstr = "endgmv  ";
  out.write(endstr.c_str(), endstr.size());
}





void GMVIO::write_discontinuous_gmv (const std::string& name,
				     const EquationSystems& es,
				     const bool write_partitioning) const
{
  std::vector<std::string> solution_names;
  std::vector<Number>      v;

  // Get a reference to the mesh
  const MeshBase& mesh = this->cmesh();
  
  es.build_variable_names  (solution_names);
  es.build_discontinuous_solution_vector (v);
  
  if (mesh.processor_id() != 0) return;
  
  std::ofstream out(name.c_str());
  
  assert (out.good());

  // Begin interfacing with the GMV data file
  {

    // write the nodes    
    out << "gmvinput ascii\n\n";
    
    // Compute the total weight
    {
      const_active_elem_iterator       it (mesh.elements_begin());
      const const_active_elem_iterator end(mesh.elements_end());

      unsigned int tw=0;
      
      for ( ; it != end; ++it)
	tw += (*it)->n_nodes();
      
      out << "nodes " << tw << '\n';
    }
    


    // Write all the x values
    {
      const_active_elem_iterator       it (mesh.elements_begin());
      const const_active_elem_iterator end(mesh.elements_end());
      
      for ( ; it != end; ++it)   
	for (unsigned int n=0; n<(*it)->n_nodes(); n++)
	  out << (*it)->point(n)(0) << " ";

      out << '\n';
    }
    
    
    // Write all the y values
    {
      const_active_elem_iterator       it (mesh.elements_begin());
      const const_active_elem_iterator end(mesh.elements_end());
      
      for ( ; it != end; ++it)   
	for (unsigned int n=0; n<(*it)->n_nodes(); n++)
	  out << (*it)->point(n)(1) << " ";

      out << '\n';
    }
    
    
    // Write all the z values
    {
      const_active_elem_iterator       it (mesh.elements_begin());
      const const_active_elem_iterator end(mesh.elements_end());
      
      for ( ; it != end; ++it)   
	for (unsigned int n=0; n<(*it)->n_nodes(); n++)
	  out << (*it)->point(n)(2) << " ";

      out << '\n' << '\n';
    }
  }


  
  {
    // write the connectivity
    
    out << "cells " << mesh.n_active_elem() << '\n';

    const_active_elem_iterator       it (mesh.elements_begin());
    const const_active_elem_iterator end(mesh.elements_end());

    unsigned int nn=1;
    
    switch (mesh.mesh_dimension())
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
		    out << "quad 4\n";
		    for (unsigned int i=0; i<(*it)->n_nodes(); i++)
		      out << nn++ << " ";
		    
		  }
		  else if (((*it)->type() == TRI3) ||
			   ((*it)->type() == TRI6))
		    {
		      out << "tri 3\n";
		    for (unsigned int i=0; i<(*it)->n_nodes(); i++)
		      out << nn++ << " ";
		    
		    }
		else
		  {
		    error();
		  }
		
		out << '\n';
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
    
    out << '\n';
  }
  

  
  // optionally write the partition information
  if (write_partitioning)
    {
      out << "material "
	  << mesh.n_processors()
	  << " 0"<< '\n';

      for (unsigned int proc=0; proc<mesh.n_processors(); proc++)
	out << "proc_" << proc << '\n';
      
      const_active_elem_iterator       it (mesh.elements_begin());
      const const_active_elem_iterator end(mesh.elements_end());

      for ( ; it != end; ++it)
	out << (*it)->processor_id()+1 << '\n';
      
      out << '\n';
    }


  
  // write the data
  {
    const unsigned int n_vars = solution_names.size();
      
    //    assert (v.size() == tw*n_vars);

    out << "variable\n";


      for (unsigned int c=0; c<n_vars; c++)
	{

#ifdef USE_COMPLEX_NUMBERS

	  // in case of complex data, write _tree_ data sets
	  // for each component

	  // this is the real part
	  out << "r_" << solution_names[c] << " 1\n";
	  {
	    const_active_elem_iterator       it (mesh.elements_begin());
	    const const_active_elem_iterator end(mesh.elements_end());
	    
	    for ( ; it != end; ++it)
	      for (unsigned int n=0; n<(*it)->n_nodes(); n++)
		out << std::setprecision(10) << v[(n++)*n_vars + c].real() << " ";	    
	  }	  	  
	  out << '\n' << '\n';


	  // this is the imaginary part
	  out << "i_" << solution_names[c] << " 1\n";	  
	  {
	    const_active_elem_iterator       it (mesh.elements_begin());
	    const const_active_elem_iterator end(mesh.elements_end());
	    
	    for ( ; it != end; ++it)
	      for (unsigned int n=0; n<(*it)->n_nodes(); n++)
		out << std::setprecision(10) << v[(n++)*n_vars + c].imag() << " ";	    
	  }	  
	  out << '\n' << '\n';

	  // this is the magnitude
	  out << "a_" << solution_names[c] << " 1\n";
	  {
	    const_active_elem_iterator       it (mesh.elements_begin());
	    const const_active_elem_iterator end(mesh.elements_end());
	    
	    for ( ; it != end; ++it)
	      for (unsigned int n=0; n<(*it)->n_nodes(); n++)
		out << std::setprecision(10)
		    << std::abs(v[(n++)*n_vars + c]) << " ";	    
	  }
	  out << '\n' << '\n';

#else

	  out << solution_names[c] << " 1\n";
	  {
	    const_active_elem_iterator       it (mesh.elements_begin());
	    const const_active_elem_iterator end(mesh.elements_end());
	    
	    unsigned int nn=0;
	  
	    for ( ; it != end; ++it)
	      for (unsigned int n=0; n<(*it)->n_nodes(); n++)
		out << std::setprecision(10) << v[(nn++)*n_vars + c] << " ";	    
	  }	  
	  out << '\n' << '\n';

#endif

	}
      
      out << "endvars\n";
    }

  
  // end of the file
  out << '\n' << "endgmv\n";
}

