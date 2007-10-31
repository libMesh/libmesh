// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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

// Changes: 
// o no more subelements, all elements are written down to GMV directly
// o Some nodes have to be left out, eg node 8 in QUAD9 
// o


// C++ includes
#include <iomanip>
#include <fstream>
#include <cstring> // for std::strcpy, std::memcpy
#include <cstdio>  // for std::sprintf
#include <vector>

// Local includes
#include "libmesh_config.h"
#include "gmv_io.h"
#include "mesh_base.h"
#include "elem.h"
#include "equation_systems.h"
#include "numeric_vector.h"

// Wrap everything in a GMV namespace and
// use extern "C" to avoid name mangling.
#ifdef HAVE_GMV
namespace GMV
{
  extern "C"
  {
#include "gmvread.h"
  }
}
#endif

// anonymous namespace to hold local data
namespace
{
  /**
   * Defines mapping from libMesh element types to GMV element types.
   */
  struct elementDefinition {
    std::string label;
    std::vector<unsigned int> nodes;
  };


  // maps from a libMesh element type to the proper
  // GMV elementDefinition.  Placing the data structure
  // here in this anonymous namespace gives us the
  // benefits of a global variable without the nasty
  // side-effects
  std::map<ElemType, elementDefinition> eletypes;



  // ------------------------------------------------------------
  // helper function to initialize the eletypes map
  void init_eletypes ()
  {
    if (eletypes.empty())
      {
	// This should happen only once.  The first time this method
	// is called the eletypes data struture will be empty, and
	// we will fill it.  Any subsequent calls will find an initialized
	// eletypes map and will do nothing.

	//==============================
	// setup the element definitions
	elementDefinition eledef;

	// use "swap trick" from Scott Meyer's "Effective STL" to initialize
	// eledef.nodes vector
	
	// EDGE2
	{
	  eledef.label = "line 2";
	  const unsigned int nodes[] = {0,1};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);
	
	  eletypes[EDGE2] = eledef;
	}
  
	// LINE3
	{
	  eledef.label = "3line 3";
	  const unsigned int nodes[] = {0,1,2};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);
	  
	  eletypes[EDGE3] = eledef;
	}
      
	// TRI3
	{
	  eledef.label = "tri3 3";
	  const unsigned int nodes[] = {0,1,2};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);
	  
	  eletypes[TRI3] = eledef;
	}
      
	// TRI6
	{
	  eledef.label = "6tri 6";
	  const unsigned int nodes[] = {0,1,2,3,4,5};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);

	  eletypes[TRI6] = eledef;
	}
      
	// QUAD4
	{
	  eledef.label = "quad 4";
	  const unsigned int nodes[] = {0,1,2,3};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);

	  eletypes[QUAD4] = eledef;
	}
      
	// QUAD8, QUAD9
	{
	  eledef.label = "8quad 8";
	  const unsigned int nodes[] = {0,1,2,3,4,5,6,7};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);

	  eletypes[QUAD8] = eledef;
	  eletypes[QUAD9] = eledef;
	}
      
	// HEX8
	{
	  eledef.label = "phex8 8";
	  const unsigned int nodes[] = {0,1,2,3,4,5,6,7};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);

	  eletypes[HEX8] = eledef;
	}
      
	// HEX20, HEX27
	{
	  eledef.label = "phex20 20";
	  const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);

	  eletypes[HEX20] = eledef;
	  eletypes[HEX27] = eledef;
	}
      
	// TET4
	{
	  eledef.label = "tet 4";
	  const unsigned int nodes[] = {0,1,2,3};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);

	  eletypes[TET4] = eledef;
	}
      
	// TET10
	{
	  eledef.label = "ptet10 10";
	  const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,9};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);
	  
	  eletypes[TET10] = eledef;
	}
      
	// PRISM6
	{
	  eledef.label = "pprism6 6";
	  const unsigned int nodes[] = {0,1,2,3,4,5};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);
	  
	  eletypes[PRISM6] = eledef;
	}
      
	// PRISM15, PRISM18
	{
	  eledef.label = "pprism15 15";
	  const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,12,13,14,9,10,11};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);
	  
	  eletypes[PRISM15] = eledef;
	  eletypes[PRISM18] = eledef;
	}
	//==============================      
      }
  }
  
} // end anonymous namespace


// ------------------------------------------------------------
// GMVIO  members
void GMVIO::write (const std::string& fname)
{
  if (libMesh::processor_id() == 0)
    if (this->binary())
      this->write_binary (fname);
    else
      this->write_ascii_old_impl  (fname);
}



void GMVIO::write_nodal_data (const std::string& fname,
                              const std::vector<Number>& soln,
                              const std::vector<std::string>& names)
{
  if (libMesh::processor_id() == 0)
    if (this->binary())
      this->write_binary (fname, &soln, &names);
    else
      this->write_ascii_old_impl  (fname, &soln, &names);
}



void GMVIO::write_ascii_new_impl (const std::string& fname,
				  const std::vector<Number>* v,
				  const std::vector<std::string>* solution_names)
{
#ifdef ENABLE_INFINITE_ELEMENTS

  std::cerr << "WARNING:  GMVIO::write_ascii_new_impl() not infinite-element aware!"
	    << std::endl;
  here();

  this->write_ascii_old_impl (fname, v, solution_names);

#else
  
  // Open the output file stream
  std::ofstream out (fname.c_str());
  
  assert (out.good());

  // Get a reference to the mesh
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

  unsigned int mesh_max_p_level = 0;

  // Begin interfacing with the GMV data file
  {
    out << "gmvinput ascii\n\n";

    // write the nodes
    out << "nodes " << mesh.n_nodes() << "\n";
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      out << mesh.point(v)(0) << " ";
    out << "\n";
    
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      out << mesh.point(v)(1) << " ";
    out << "\n";
    
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      out << mesh.point(v)(2) << " ";
    out << "\n\n";
  }

  {
    // write the connectivity
    out << "cells " << mesh.n_active_elem() << "\n";
    
    // initialize the eletypes map
    init_eletypes();

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end(); 
    
    for ( ; it != end; ++it)
      {
        const Elem* elem = *it;

        mesh_max_p_level = std::max(mesh_max_p_level,
                                    elem->p_level());

	// Make sure we have a valid entry for
	// the current element type.
	assert (eletypes.count(elem->type()));

        const elementDefinition& ele = eletypes[elem->type()];

	// The element mapper better not require any more nodes
	// than are present in the current element!
	assert (ele.nodes.size() <= elem->n_nodes());
	
        out << ele.label << "\n";
        for (unsigned int i=0; i < ele.nodes.size(); i++)
          out << elem->node(ele.nodes[i])+1 << " ";
        out << "\n";
      }
    out << "\n";
  }
  
  // optionally write the partition information
  if (this->partitioning())
    {
      out << "material "
	  << mesh.n_partitions()
	// Note: GMV may give you errors like
	// Error, material for cell 1 is greater than 1
        // Error, material for cell 2 is greater than 1
        // Error, material for cell 3 is greater than 1
	// ... because you put the wrong number of partitions here.
	// To ensure you write the correct number of materials, call
	// mesh.recalculate_n_partitions() before writing out the
	// mesh.
	// Note: we can't call it now because the Mesh is const here and
	// it is a non-const function.
          << " 0\n";

      for (unsigned int proc=0; proc<mesh.n_partitions(); proc++)
        out << "proc_" << proc << "\n";
      
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      // FIXME - don't we need to use an elementDefinition here? - RHS
      for ( ; it != end; ++it)
        out << (*it)->processor_id()+1 << "\n";
      out << "\n";
    }

  // If there are *any* variables at all in the system (including
  // p level, or arbitrary cell-based data)
  // to write, the gmv file needs to contain the word "variable"
  // on a line by itself.
  bool write_variable = false;

  // 1.) p-levels
  if (this->p_levels() && mesh_max_p_level)
    write_variable = true;

  // 2.) solution data
  if ((solution_names != NULL) && (v != NULL))
    write_variable = true;

  // 3.) cell-centered data
  if ( !(this->_cell_centered_data.empty()) )
    write_variable = true;
  
  if (write_variable)
    out << "variable\n";
  
//   if ((this->p_levels() && mesh_max_p_level) || 
//     ((solution_names != NULL) && (v != NULL)))
//     out << "variable\n";

  // optionally write the polynomial degree information
  if (this->p_levels() && mesh_max_p_level)
    {
      out << "p_level 0\n";

      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      for ( ; it != end; ++it)
        {
          const Elem* elem = *it;

          const elementDefinition& ele = eletypes[elem->type()];

	  // The element mapper better not require any more nodes
	  // than are present in the current element!
	  assert (ele.nodes.size() <= elem->n_nodes());

          for (unsigned int i=0; i < ele.nodes.size(); i++)
            out << elem->p_level() << " ";
        }
      out << "\n\n";
    }


  // optionally write cell-centered data
  if ( !(this->_cell_centered_data.empty()) )
    {
      std::map<std::string, const std::vector<Real>* >::iterator       it  = this->_cell_centered_data.begin();
      const std::map<std::string, const std::vector<Real>* >::iterator end = this->_cell_centered_data.end();      

      for (; it != end; ++it)
	{
	  // write out the variable name, followed by a zero.
	  out << (*it).first << " 0\n";

	  const std::vector<Real>* the_array = (*it).second;
	  
	  // Loop over active elements, write out cell data.  If second-order cells
	  // are split into sub-elements, the sub-elements inherit their parent's
	  // cell-centered data.
	  MeshBase::const_element_iterator       elem_it  = mesh.active_elements_begin();
	  const MeshBase::const_element_iterator elem_end = mesh.active_elements_end(); 
	  
	  for (; elem_it != elem_end; ++elem_it)
	    {
	      const Elem* e = *elem_it;
	      
	      // Use the element's ID to find the value.
	      assert (e->id() < the_array->size());
	      const Real the_value = the_array->operator[](e->id());
	      
	      if (this->subdivide_second_order())
		for (unsigned int se=0; se < e->n_sub_elem(); se++)
		  out << the_value << " ";
	      else
		out << the_value << " ";
	    }
	  
	  out << "\n\n";
	}
    }

  
  // optionally write the data
  if ((solution_names != NULL) && (v != NULL))
    {      
      const unsigned int n_vars = solution_names->size();

      if (!(v->size() == mesh.n_nodes()*n_vars))
        std::cerr << "ERROR: v->size()=" << v->size()
                  << ", mesh.n_nodes()=" << mesh.n_nodes()
                  << ", n_vars=" << n_vars
                  << ", mesh.n_nodes()*n_vars=" << mesh.n_nodes()*n_vars
                  << "\n";
      
      assert (v->size() == mesh.n_nodes()*n_vars);

      for (unsigned int c=0; c<n_vars; c++)
        {

#ifdef USE_COMPLEX_NUMBERS

          // in case of complex data, write _three_ data sets
          // for each component

          // this is the real part
          out << "r_" << (*solution_names)[c] << " 1\n";
	  
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            out << std::setprecision(10) << (*v)[n*n_vars + c].real() << " ";

          out << "\n\n";

          // this is the imaginary part
          out << "i_" << (*solution_names)[c] << " 1\n";
	  
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            out << std::setprecision(10) << (*v)[n*n_vars + c].imag() << " ";

          out << "\n\n";

          // this is the magnitude
          out << "a_" << (*solution_names)[c] << " 1\n";
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            out << std::setprecision(10)
                << std::abs((*v)[n*n_vars + c]) << " ";

          out << "\n\n";

#else

          out << (*solution_names)[c] << " 1\n";
	  
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            out << std::setprecision(10) << (*v)[n*n_vars + c] << " ";
	  
          out << "\n\n";

#endif
        }
      
    }

  // If we wrote any variables, we have to close the variable section now
  if (write_variable)
    out << "endvars\n";

  
  // end of the file
  out << "\nendgmv\n";

#endif
}






void GMVIO::write_ascii_old_impl (const std::string& fname,
				  const std::vector<Number>* v,
				  const std::vector<std::string>* solution_names)
{
  // Open the output file stream
  std::ofstream out (fname.c_str());
  
  assert (out.good());

  // Get a reference to the mesh
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

  unsigned int mesh_max_p_level = 0;
  
  // Begin interfacing with the GMV data file

  // FIXME - if subdivide_second_order() is off,
  // we probably should only be writing the
  // vertex nodes - RHS
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
    
    out << "cells ";
    if (this->subdivide_second_order())
      out << mesh.n_active_sub_elem();
    else
      out << mesh.n_active_elem();
    out << '\n';

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

    switch (mesh.mesh_dimension())
      {
      case 1:
	{
	  // The same temporary storage will be used for each element
	  std::vector<unsigned int> conn;

	  for ( ; it != end; ++it)
            {
              mesh_max_p_level = std::max(mesh_max_p_level,
                                          (*it)->p_level());

              if (this->subdivide_second_order())
	        for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	          {
		    out << "line 2\n";
		    (*it)->connectivity(se, TECPLOT, conn);
		    for (unsigned int i=0; i<conn.size(); i++)
		      out << conn[i] << " ";
		
		    out << '\n';
	          }
              else
                {
		  out << "line 2\n";
                  if ((*it)->default_order() == FIRST)
		    (*it)->connectivity(0, TECPLOT, conn);
                  else
                    {
                      AutoPtr<Elem> lo_elem = Elem::build(
                        Elem::first_order_equivalent_type((*it)->type()));
                      for (unsigned int i = 0; i != lo_elem->n_nodes(); ++i)
                        lo_elem->set_node(i) = (*it)->get_node(i);
		      lo_elem->connectivity(0, TECPLOT, conn);
                    }
		  for (unsigned int i=0; i<conn.size(); i++)
		    out << conn[i] << " ";
		
		  out << '\n';
                }
            }
	  break;
	}
	
      case 2:
	{
	  // The same temporary storage will be used for each element
	  std::vector<unsigned int> conn;
	  
	  for ( ; it != end; ++it)
            {
              mesh_max_p_level = std::max(mesh_max_p_level,
                                          (*it)->p_level());

              if (this->subdivide_second_order())
	        for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	          {
		    // Quad elements
		    if (((*it)->type() == QUAD4) ||
		        ((*it)->type() == QUAD8) || // Note: QUAD8 will be output as one central quad and
			                            // four surrounding triangles (though they will be written
			                            // to GMV as QUAD4s).  
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
		      error();
                  }
              else // !this->subdivide_second_order()
                {
		  // Quad elements
		  if (((*it)->type() == QUAD4)
#ifdef ENABLE_INFINITE_ELEMENTS
		      || ((*it)->type() == INFQUAD4)
#endif
		      )
                    {
		      (*it)->connectivity(0, TECPLOT, conn);
		      out << "quad 4\n";
		      for (unsigned int i=0; i<conn.size(); i++)
		        out << conn[i] << " ";
		    }
		  else if (((*it)->type() == QUAD8) ||
		           ((*it)->type() == QUAD9)
#ifdef ENABLE_INFINITE_ELEMENTS
		           || ((*it)->type() == INFQUAD6)
#endif
		          )
		    {
                      AutoPtr<Elem> lo_elem = Elem::build(
                        Elem::first_order_equivalent_type((*it)->type()));
                      for (unsigned int i = 0; i != lo_elem->n_nodes(); ++i)
                        lo_elem->set_node(i) = (*it)->get_node(i);
		      lo_elem->connectivity(0, TECPLOT, conn);
		      out << "quad 4\n";
		      for (unsigned int i=0; i<conn.size(); i++)
		        out << conn[i] << " ";
		    }
		  else if ((*it)->type() == TRI3)
		    {
		      (*it)->connectivity(0, TECPLOT, conn);
		      out << "tri 3\n";
		      for (unsigned int i=0; i<3; i++)
		        out << conn[i] << " ";
		    }
		  else if ((*it)->type() == TRI6)
		    {
                      AutoPtr<Elem> lo_elem = Elem::build(
                        Elem::first_order_equivalent_type((*it)->type()));
                      for (unsigned int i = 0; i != lo_elem->n_nodes(); ++i)
                        lo_elem->set_node(i) = (*it)->get_node(i);
		      lo_elem->connectivity(0, TECPLOT, conn);
		      out << "tri 3\n";
		      for (unsigned int i=0; i<3; i++)
		        out << conn[i] << " ";
		    }
		
		  out << '\n';
	        }
            }
	  
	  break;
	}
	
	
      case 3:
	{
	  // The same temporary storage will be used for each element
	  std::vector<unsigned int> conn;
	  
	  for ( ; it != end; ++it)
            {
              mesh_max_p_level = std::max(mesh_max_p_level,
                                          (*it)->p_level());

              if (this->subdivide_second_order())
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
		        std::cout << "Encountered an unrecognized element "
			          << "type: " << (*it)->type()
				  << "\nPossibly a dim-1 dimensional "
			          << "element?  Aborting..."
			          << std::endl;
		        error();
		      }
		
		    out << '\n';
	          }
              else // !this->subdivide_second_order()
                {
                  AutoPtr<Elem> lo_elem = Elem::build(
                    Elem::first_order_equivalent_type((*it)->type()));
                  for (unsigned int i = 0; i != lo_elem->n_nodes(); ++i)
                    lo_elem->set_node(i) = (*it)->get_node(i);
		  if ((lo_elem->type() == HEX8)
#ifdef  ENABLE_INFINITE_ELEMENTS
                      || (lo_elem->type() == HEX27)
#endif
                     )
		    {
		      out << "phex8 8\n";
		      lo_elem->connectivity(0, TECPLOT, conn);
		      for (unsigned int i=0; i<conn.size(); i++)
		        out << conn[i] << " ";
		    }
		
		  else if (lo_elem->type() == TET4)
		    {
		      out << "tet 4\n";
		      lo_elem->connectivity(0, TECPLOT, conn);
		      out << conn[0] << " "
			  << conn[2] << " "
			  << conn[1] << " "
			  << conn[4] << " ";
		    }
		  else if ((lo_elem->type() == PRISM6)
#ifdef  ENABLE_INFINITE_ELEMENTS
                           || (lo_elem->type() == INFPRISM6)
#endif
			   )
		    {
		      /**
		       * Note that the prisms are treated as
		       * degenerated phex8's.
		       */
		      out << "phex8 8\n";
		      lo_elem->connectivity(0, TECPLOT, conn);
		      for (unsigned int i=0; i<conn.size(); i++)
		        out << conn[i] << " ";
		    }
		
		  else
		    {
		      std::cout << "Encountered an unrecognized element "
			        << "type.  Possibly a dim-1 dimensional "
			        << "element?  Aborting..."
			        << std::endl;
		      error();
		    }
		
		  out << '\n';
	        }
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
      
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      for ( ; it != end; ++it)
        if (this->subdivide_second_order())
	  for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	    out << (*it)->processor_id()+1 << '\n';
        else
	  out << (*it)->processor_id()+1 << '\n';
      
      out << '\n';
    }


  // If there are *any* variables at all in the system (including
  // p level, or arbitrary cell-based data)
  // to write, the gmv file needs to contain the word "variable"
  // on a line by itself.
  bool write_variable = false;

  // 1.) p-levels
  if (this->p_levels() && mesh_max_p_level)
    write_variable = true;

  // 2.) solution data
  if ((solution_names != NULL) && (v != NULL))
    write_variable = true;

  // 3.) cell-centered data
  if ( !(this->_cell_centered_data.empty()) )
    write_variable = true;
  
  if (write_variable)
    out << "variable\n";


  // optionally write the p-level information
  if (this->p_levels() && mesh_max_p_level)
    {
      out << "p_level 0\n";

      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      for ( ; it != end; ++it)
        if (this->subdivide_second_order())
	  for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
            out << (*it)->p_level() << " ";
        else
          out << (*it)->p_level() << " ";
      out << "\n\n";
    }



  
  // optionally write cell-centered data
  if ( !(this->_cell_centered_data.empty()) )
    {
      std::map<std::string, const std::vector<Real>* >::iterator       it  = this->_cell_centered_data.begin();
      const std::map<std::string, const std::vector<Real>* >::iterator end = this->_cell_centered_data.end();      

      for (; it != end; ++it)
	{
	  // write out the variable name, followed by a zero.
	  out << (*it).first << " 0\n";

	  const std::vector<Real>* the_array = (*it).second;
	  
	  // Loop over active elements, write out cell data.  If second-order cells
	  // are split into sub-elements, the sub-elements inherit their parent's
	  // cell-centered data.
	  MeshBase::const_element_iterator       elem_it  = mesh.active_elements_begin();
	  const MeshBase::const_element_iterator elem_end = mesh.active_elements_end(); 
	  
	  for (; elem_it != elem_end; ++elem_it)
	    {
	      const Elem* e = *elem_it;
	      
	      // Use the element's ID to find the value...
	      assert (e->id() < the_array->size());
	      const Real the_value = the_array->operator[](e->id());
	      
	      if (this->subdivide_second_order())
		for (unsigned int se=0; se < e->n_sub_elem(); se++)
		  out << the_value << " ";
	      else
		out << the_value << " ";
	    }
	  
	  out << "\n\n";
	}
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
      
    }
  
  // If we wrote any variables, we have to close the variable section now
  if (write_variable)
    out << "endvars\n";

  
  // end of the file
  out << "\nendgmv\n";
}







void GMVIO::write_binary (const std::string& fname,
                          const std::vector<Number>* vec,
                          const std::vector<std::string>* solution_names)
{
  std::ofstream out (fname.c_str());
  
  assert (out.good());

  // get a reference to the mesh
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

  unsigned int mesh_max_p_level = 0;
  
  char buf[80];

  // Begin interfacing with the GMV data file
  {
    // write the nodes
    std::strcpy(buf, "gmvinput");
    out.write(buf, std::strlen(buf));
	      
    std::strcpy(buf, "ieeei4r4");
    out.write(buf, std::strlen(buf));
  }


  
  // write the nodes
  {
    std::strcpy(buf, "nodes   ");
    out.write(buf, std::strlen(buf));

    unsigned int tempint = mesh.n_nodes();
    
    std::memcpy(buf, &tempint, sizeof(unsigned int));

    out.write(buf, sizeof(unsigned int));

    // write the x coordinate
    float *temp = new float[mesh.n_nodes()];
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      temp[v] = static_cast<float>(mesh.point(v)(0));
    out.write(reinterpret_cast<char *>(temp), sizeof(float)*mesh.n_nodes());

    // write the y coordinate
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      temp[v] = static_cast<float>(mesh.point(v)(1));
    out.write(reinterpret_cast<char *>(temp), sizeof(float)*mesh.n_nodes());

    // write the z coordinate
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      temp[v] = static_cast<float>(mesh.point(v)(2));
    out.write(reinterpret_cast<char *>(temp), sizeof(float)*mesh.n_nodes());

    delete [] temp;
  }


  // write the connectivity
  {
    std::strcpy(buf, "cells   ");
    out.write(buf, std::strlen(buf));

    unsigned int tempint = mesh.n_active_elem();
    
    std::memcpy(buf, &tempint, sizeof(unsigned int));
    
    out.write(buf, sizeof(unsigned int));

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

    switch (mesh.mesh_dimension())
      {

      case 1:
        for ( ; it != end; ++it)
          {
            mesh_max_p_level = std::max(mesh_max_p_level,
                                        (*it)->p_level());

            for(unsigned se = 0; se < (*it)->n_sub_elem(); ++se)
              {
                std::strcpy(buf, "line    ");
                out.write(buf, std::strlen(buf));
	      
                tempint = 2;
                std::memcpy(buf, &tempint, sizeof(unsigned int));
                out.write(buf, sizeof(unsigned int));

                std::vector<unsigned int> conn;
                (*it)->connectivity(se,TECPLOT,conn);
	      
                out.write(reinterpret_cast<char*>(&conn[0]), sizeof(unsigned int)*tempint);
              }
          }
        break;
	
      case 2:
        for ( ; it != end; ++it)
          {
            mesh_max_p_level = std::max(mesh_max_p_level,
                                        (*it)->p_level());

            for(unsigned se = 0; se < (*it)->n_sub_elem(); ++se)
              {
                std::strcpy(buf, "quad    ");
                out.write(buf, std::strlen(buf));
                tempint = 4;
                std::memcpy(buf, &tempint, sizeof(unsigned int));
                out.write(buf, sizeof(unsigned int));
                std::vector<unsigned int> conn;
                (*it)->connectivity(se,TECPLOT,conn);
                out.write(reinterpret_cast<char*>(&conn[0]), sizeof(unsigned int)*tempint);
              }
          }
        break;
      case 3:
        for ( ; it != end; ++it)
          {
            mesh_max_p_level = std::max(mesh_max_p_level,
                                        (*it)->p_level());

            for(unsigned se = 0; se < (*it)->n_sub_elem(); ++se)
              {
                std::strcpy(buf, "phex8   ");
                out.write(buf, std::strlen(buf));
                tempint = 8;
                std::memcpy(buf, &tempint, sizeof(unsigned int));
                out.write(buf, sizeof(unsigned int));
                std::vector<unsigned int> conn;
                (*it)->connectivity(se,TECPLOT,conn);
                out.write(reinterpret_cast<char*>(&conn[0]), sizeof(unsigned int)*tempint);
              }
          }
        break;
      default:
        error();
	
      }
  }
  
  
  
  // optionally write the partition information
  if (this->partitioning())
    {
      std::strcpy(buf, "material");
      out.write(buf, std::strlen(buf));
      
      unsigned int tmpint = mesh.n_processors();
      std::memcpy(buf, &tmpint, sizeof(unsigned int));
      out.write(buf, sizeof(unsigned int));

      tmpint = 0; // IDs are cell based
      std::memcpy(buf, &tmpint, sizeof(unsigned int));
      out.write(buf, sizeof(unsigned int));


      for (unsigned int proc=0; proc<mesh.n_processors(); proc++)
        {
          std::sprintf(buf, "proc_%d", proc);
          out.write(buf, 8);
        }

      std::vector<unsigned int> proc_id (mesh.n_active_elem());
      
      unsigned int n=0;
      
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      for ( ; it != end; ++it)
        for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
          proc_id[n++] = (*it)->processor_id()+1;
      
      
      out.write(reinterpret_cast<char *>(&proc_id[0]),
                sizeof(unsigned int)*proc_id.size());
    }

  // If there are *any* variables at all in the system (including
  // p level, or arbitrary cell-based data)
  // to write, the gmv file needs to contain the word "variable"
  // on a line by itself.
  bool write_variable = false;

  // 1.) p-levels
  if (this->p_levels() && mesh_max_p_level)
    write_variable = true;

  // 2.) solution data
  if ((solution_names != NULL) && (vec != NULL))
    write_variable = true;

  //   // 3.) cell-centered data - unsupported
  //   if ( !(this->_cell_centered_data.empty()) )
  //     write_variable = true;

  if (write_variable)
    {
      std::strcpy(buf, "variable");
      out.write(buf, std::strlen(buf));
    }

  // optionally write the partition information
  if (this->p_levels() && mesh_max_p_level)
    {
      unsigned int n_floats = mesh.n_active_elem();
      for (unsigned int i=0; i != mesh.mesh_dimension(); ++i)
        n_floats *= 2;

      float *temp = new float[n_floats];

      std::strcpy(buf, "p_level");
      out.write(buf, std::strlen(buf));

      unsigned int tempint = 0; // p levels are cell data

      std::memcpy(buf, &tempint, sizeof(unsigned int));
      out.write(buf, sizeof(unsigned int));

      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 
      unsigned int n=0;

      for (; it != end; ++it)
        for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
          temp[n++] = static_cast<float>( (*it)->p_level() );

      out.write(reinterpret_cast<char *>(temp),
                sizeof(float)*n_floats);

      delete [] temp;
    }

  
   // optionally write cell-centered data
   if ( !(this->_cell_centered_data.empty()) )
     {
       std::cerr << "Cell-centered data not (yet) supported in binary I/O mode!" << std::endl;
       
//        std::map<std::string, const std::vector<Real>* >::iterator       it  = this->_cell_centered_data.begin();
//        const std::map<std::string, const std::vector<Real>* >::iterator end = this->_cell_centered_data.end();      

//        for (; it != end; ++it)
//  	{
//  	  // Write out the variable name ...
//  	  std::strcpy(buf, (*it).first.c_str());
//  	  out.write(buf, std::strlen(buf));
	  
//  	  // ... followed by a zero.
//  	  unsigned int tempint = 0; // 0 signifies cell data
//  	  std::memcpy(buf, &tempint, sizeof(unsigned int));
//  	  out.write(buf, sizeof(unsigned int));

//  	  // Get a pointer to the array of cell-centered data values
//  	  const std::vector<Real>* the_array = (*it).second;

// 	  // Since the_array might contain zeros (for inactive elements) we need to
// 	  // make a copy of it containing just values for active elements.
// 	  const unsigned int n_floats = mesh.n_active_elem() * (1<<mesh.mesh_dimension());
// 	  float *temp = new float[n_floats];

// 	  MeshBase::const_element_iterator       elem_it  = mesh.active_elements_begin();
// 	  const MeshBase::const_element_iterator elem_end = mesh.active_elements_end(); 
// 	  unsigned int n=0;

// 	  for (; elem_it != elem_end; ++elem_it)
// 	    {
// 	      // If there's a seg-fault, it will probably be here!
// 	      const float the_value = static_cast<float>(the_array->operator[]((*elem_it)->id()));
	      
// 	      for (unsigned int se=0; se<(*elem_it)->n_sub_elem(); se++)
// 		temp[n++] = the_value;
// 	    }

	  
//  	  // Write "the_array" directly to the file
//  	  out.write(reinterpret_cast<char *>(temp),
//  		    sizeof(float)*n_floats);

// 	  delete [] temp;
//  	}
     }

  
  
  
  // optionally write the data
  if ((solution_names != NULL) &&
      (vec != NULL))
    {
      float *temp = new float[mesh.n_nodes()];

      const unsigned int n_vars = solution_names->size();
      
      for (unsigned int c=0; c<n_vars; c++)
        {

#ifdef USE_COMPLEX_NUMBERS
          // for complex data, write three datasets


          // Real part
          std::strcpy(buf, "r_");
          out.write(buf, 2);
          std::strcpy(buf, (*solution_names)[c].c_str());
          out.write(buf, 6);
	  
          unsigned int tempint = 1; // always do nodal data
          std::memcpy(buf, &tempint, sizeof(unsigned int));
          out.write(buf, sizeof(unsigned int));
	  
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            temp[n] = static_cast<float>( (*vec)[n*n_vars + c].real() );
	  
          out.write(reinterpret_cast<char *>(temp), sizeof(float)*mesh.n_nodes());


          // imaginary part
          std::strcpy(buf, "i_");
          out.write(buf, 2);
          std::strcpy(buf, (*solution_names)[c].c_str());
          out.write(buf, 6);
	  
          std::memcpy(buf, &tempint, sizeof(unsigned int));
          out.write(buf, sizeof(unsigned int));
	  
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            temp[n] = static_cast<float>( (*vec)[n*n_vars + c].imag() );
	  
          out.write(reinterpret_cast<char *>(temp), sizeof(float)*mesh.n_nodes());

          // magnitude
          std::strcpy(buf, "a_");
          out.write(buf, 2);
          std::strcpy(buf, (*solution_names)[c].c_str());
          out.write(buf, 6);
	  
          std::memcpy(buf, &tempint, sizeof(unsigned int));
          out.write(buf, sizeof(unsigned int));
	  
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            temp[n] = static_cast<float>(std::abs((*vec)[n*n_vars + c]));
	  
          out.write(reinterpret_cast<char *>(temp), sizeof(float)*mesh.n_nodes());

#else


          std::strcpy(buf, (*solution_names)[c].c_str());
          out.write(buf, 8);
	  
          unsigned int tempint = 1; // always do nodal data
          std::memcpy(buf, &tempint, sizeof(unsigned int));
          out.write(buf, sizeof(unsigned int));
	  
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            temp[n] = static_cast<float>((*vec)[n*n_vars + c]);
	  
          out.write(reinterpret_cast<char *>(temp), sizeof(float)*mesh.n_nodes());


#endif

	  
        }
    
      delete [] temp;
      
    }

  // If we wrote any variables, we have to close the variable section now
  if (write_variable)
    {
      std::strcpy(buf, "endvars ");
      out.write(buf, std::strlen(buf));
    }

  // end the file
  std::strcpy(buf, "endgmv  ");
  out.write(buf, std::strlen(buf));
}









void GMVIO::write_discontinuous_gmv (const std::string& name,
				     const EquationSystems& es,
				     const bool write_partitioning) const
{
  std::vector<std::string> solution_names;
  std::vector<Number>      v;

  // Get a reference to the mesh
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();
  
  es.build_variable_names  (solution_names);
  es.build_discontinuous_solution_vector (v);
  
  if (mesh.processor_id() != 0) return;
  
  std::ofstream out(name.c_str());
  
  assert (out.good());

  // Begin interfacing with the GMV data file
  {

    // write the nodes    
    out << "gmvinput ascii" << std::endl << std::endl;
    
    // Compute the total weight
    {
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      unsigned int tw=0;
      
      for ( ; it != end; ++it)
	tw += (*it)->n_nodes();
      
      out << "nodes " << tw << std::endl;
    }
    


    // Write all the x values
    {
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 
      
      for ( ; it != end; ++it)   
	for (unsigned int n=0; n<(*it)->n_nodes(); n++)
	  out << (*it)->point(n)(0) << " ";

      out << std::endl;
    }
    
    
    // Write all the y values
    {
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      for ( ; it != end; ++it)   
	for (unsigned int n=0; n<(*it)->n_nodes(); n++)
	  out << (*it)->point(n)(1) << " ";

      out << std::endl;
    }
    
    
    // Write all the z values
    {
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 
      
      for ( ; it != end; ++it)   
	for (unsigned int n=0; n<(*it)->n_nodes(); n++)
	  out << (*it)->point(n)(2) << " ";

      out << std::endl << std::endl;
    }
  }


  
  {
    // write the connectivity
    
    out << "cells " << mesh.n_active_elem() << std::endl;

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

    unsigned int nn=1;
    
    switch (mesh.mesh_dimension())
      {
      case 1:
	{
	  for ( ; it != end; ++it)
	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	      {
		if (((*it)->type() == EDGE2) ||
		    ((*it)->type() == EDGE3) ||
		    ((*it)->type() == EDGE4)
#ifdef ENABLE_INFINITE_ELEMENTS
		    || ((*it)->type() == INFEDGE2)
#endif
		    )
		  {
		    out << "line 2" << std::endl;
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
	
      case 2:
	{
	  for ( ; it != end; ++it)
	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	      {
		if (((*it)->type() == QUAD4) ||
		    ((*it)->type() == QUAD8) || // Note: QUAD8 will be output as one central quad and
			                        // four surrounding triangles (though they will be written
			                        // to GMV as QUAD4s).  
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
	  for ( ; it != end; ++it)
	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	      {
		if (((*it)->type() == HEX8) ||
		    ((*it)->type() == HEX20) ||
		    ((*it)->type() == HEX27)
#ifdef ENABLE_INFINITE_ELEMENTS
		    || ((*it)->type() == INFHEX8)
		    || ((*it)->type() == INFHEX16)
		    || ((*it)->type() == INFHEX18)
#endif
		    )
		  {
		    out << "phex8 8" << std::endl;
		    for (unsigned int i=0; i<(*it)->n_nodes(); i++)
		      out << nn++ << " ";
		  }
                else if (((*it)->type() == PRISM6) ||
                         ((*it)->type() == PRISM15) ||
                         ((*it)->type() == PRISM18)
#ifdef ENABLE_INFINITE_ELEMENTS
		         || ((*it)->type() == INFPRISM6)
		         || ((*it)->type() == INFPRISM12)
#endif
			 )
                  {
                    out << "pprism6 6" << std::endl;
		    for (unsigned int i=0; i<(*it)->n_nodes(); i++)
		      out << nn++ << " ";
                  }
                else if (((*it)->type() == TET4) ||
                         ((*it)->type() == TET10))
                  {
                    out << "tet 4" << std::endl;
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
	
      default:
	error();
      }
    
    out << std::endl;
  }
  

  
  // optionally write the partition information
  if (write_partitioning)
    {
      out << "material "
	  << mesh.n_processors()
	  << " 0"<< std::endl;

      for (unsigned int proc=0; proc<mesh.n_processors(); proc++)
	out << "proc_" << proc << std::endl;
      
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      for ( ; it != end; ++it)
	out << (*it)->processor_id()+1 << std::endl;
      
      out << std::endl;
    }


  // Writing cell-centered data is not yet supported in discontinuous GMV files.
  if ( !(this->_cell_centered_data.empty()) )
    {
      std::cerr << "Cell-centered data not (yet) supported for discontinuous GMV files!" << std::endl;
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
          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

          for ( ; it != end; ++it)
            for (unsigned int n=0; n<(*it)->n_nodes(); n++)
              out << std::setprecision(10) << v[(n++)*n_vars + c].real() << " ";	    
        }	  	  
        out << std::endl << std::endl;


        // this is the imaginary part
        out << "i_" << solution_names[c] << " 1" << std::endl;	  
        {
          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end(); 
	    
          for ( ; it != end; ++it)
            for (unsigned int n=0; n<(*it)->n_nodes(); n++)
              out << std::setprecision(10) << v[(n++)*n_vars + c].imag() << " ";	    
        }	  
        out << std::endl << std::endl;

        // this is the magnitude
        out << "a_" << solution_names[c] << " 1" << std::endl;
        {
          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end(); 
	    
          for ( ; it != end; ++it)
            for (unsigned int n=0; n<(*it)->n_nodes(); n++)
              out << std::setprecision(10)
                  << std::abs(v[(n++)*n_vars + c]) << " ";	    
        }
        out << std::endl << std::endl;

#else

        out << solution_names[c] << " 1" << std::endl;
        {
          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end(); 
	    
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





void GMVIO::add_cell_centered_data (const std::string&       cell_centered_data_name,
				    const std::vector<Real>* cell_centered_data_vals)
{
  assert (cell_centered_data_vals != NULL);

  // Make sure there are *at least* enough entries for all the active elements.
  // There can also be entries for inactive elements, they will be ignored.
  assert (cell_centered_data_vals->size() >=
	  MeshOutput<MeshBase>::mesh().n_active_elem());
  this->_cell_centered_data[cell_centered_data_name] = cell_centered_data_vals;
}






void GMVIO::read (const std::string& name)
{

  untested();
  
#ifndef HAVE_GMV

  std::cerr << "Cannot read a GMV file without the GMV API." << std::endl;
  error();

#else
  // Clear the mesh so we are sure to start from a pristeen state.
  MeshInput<MeshBase>::mesh().clear();
  
  // It is apparently possible for gmv files to contain
  // a "fromfile" directive (?) But we currently don't make
  // any use of this feature in LibMesh.  Nonzero return val
  // from any function usually means an error has occurred.
  int ierr = GMV::gmvread_open_fromfileskip(const_cast<char*>(name.c_str()));
  if (ierr != 0)
    {
      std::cerr << "GMV::gmvread_open_fromfileskip failed!" << std::endl;
      error();
    }

  
  // Loop through file until GMVEND.  
  int iend = 0;
  while (iend == 0)
    {
      GMV::gmvread_data();

      /*  Check for GMVEND.  */
      if (GMV::gmv_data.keyword == GMVEND)
        {
	  iend = 1;
	  GMV::gmvread_close();
	  break;
        }

      /*  Check for GMVERROR.  */
      if (GMV::gmv_data.keyword == GMVERROR)
        {
	  std::cerr << "Encountered GMVERROR while reading!" << std::endl;
	  error();
        }

      /*  Process the data.  */
      switch (GMV::gmv_data.keyword)
        {
	case NODES:
	  {
	    //std::cout << "Reading nodes." << std::endl;

	    if (GMV::gmv_data.num2 == NODES)
	      this->_read_nodes();
	    
	    else if (GMV::gmv_data.num2 == NODE_V)
	      {
		std::cerr << "Unsupported GMV data type NODE_V!" << std::endl;
		error();
	      }
	    break;
	  }
	  
	case CELLS:
	  {
	    // Read 1 cell at a time
	    // std::cout << "\nReading one cell." << std::endl;
	    this->_read_one_cell();
	    break;
	  }

	case MATERIAL:
	  {
	    // keyword == 6
	    // These are the materials, which we use to specify the mesh
	    // partitioning.
	    this->_read_materials();
	    break;
	  }
	  
	case VARIABLE:
	  {
	    // keyword == 8
	    // This is a field variable.  

	    // Check to see if we're done reading variables and break out.
	    if (GMV::gmv_data.datatype == ENDKEYWORD)
	      {
		// std::cout << "Done reading GMV variables." << std::endl;
		break;
	      }

	    if (GMV::gmv_data.datatype == NODE)
	      {
		// std::cout << "Reading node field data for variable "
		// 	  << GMV::gmv_data.name1 << std::endl;
		this->_read_var();
		break;
	      }
	    
	    else
	      {
		std::cerr << "Warning: Skipping variable: "
			  << GMV::gmv_data.name1
			  << " which is of unsupported GMV datatype "
			  << GMV::gmv_data.datatype
			  << ".  Nodal field data is currently the only type currently supported."
			  << std::endl;
		break;
	      }

	  }
	  
	default:
	  {
	    std::cerr << "Encountered unknown GMV keyword "
		      << GMV::gmv_data.keyword
		      << std::endl;
	    error();
	  }
        } // end switch
    } // end while

  // Done reading in the mesh, now call find_neighbors, etc.
  // MeshInput<MeshBase>::mesh().find_neighbors();
  
  // Pass true flag to skip renumbering nodes and elements
  MeshInput<MeshBase>::mesh().prepare_for_use(true);
#endif
}




void GMVIO::_read_var()
{
#ifdef HAVE_GMV
  
  // Copy all the variable's values into a local storage vector.
  _nodal_data.insert ( std::make_pair(std::string(GMV::gmv_data.name1),
				      std::vector<Number>(GMV::gmv_data.doubledata1, GMV::gmv_data.doubledata1+GMV::gmv_data.num) ) );
#endif							 
}



void GMVIO::_read_materials()
{
#ifdef HAVE_GMV
  
  // LibMesh assigns materials on a per-cell basis
  assert (GMV::gmv_data.datatype == CELL);

  //   // Material names: LibMesh has no use for these currently...
  //   std::cout << "Number of material names="
  // 	    << GMV::gmv_data.num
  // 	    << std::endl;
  
  //   for (int i = 0; i < GMV::gmv_data.num; i++)
  //     {
  //       // Build a 32-char string from the appropriate entries
  //       std::string mat_string(&GMV::gmv_data.chardata1[i*33], 32);

  //       std::cout << "Material name " << i << ": " << mat_string << std::endl;
  //     }

  //   // Material labels: These correspond to (1-based) CPU IDs, and
  //   // there should be 1 of these for each element.
  //   std::cout << "Number of material labels = "
  // 	    << GMV::gmv_data.nlongdata1
  // 	    << std::endl;
  
  for (int i = 0; i < GMV::gmv_data.nlongdata1; i++)
    {
      // Debugging Info
      // std::cout << "Material ID " << i << ": "
      // << GMV::gmv_data.longdata1[i]
      // << std::endl;
      
      MeshInput<MeshBase>::mesh().elem(i)->processor_id() =
	GMV::gmv_data.longdata1[i]-1;
    }
  
#endif
}




void GMVIO::_read_nodes()
{
#ifdef HAVE_GMV
  //   // Debug Info 
  //   std::cout << "gmv_data.datatype="
  // 	    <<  GMV::gmv_data.datatype
  // 	    << std::endl;

  // LibMesh writes UNSTRUCT=100 node data 
  assert (GMV::gmv_data.datatype == UNSTRUCT);

  // The nodal data is stored in gmv_data.doubledata{1,2,3}
  // and is nnodes long
  for (int i = 0; i < GMV::gmv_data.num; i++)
    {
      //       std::cout << "(x,y,z)="
      // 		<< "("
      // 		<< GMV::gmv_data.doubledata1[i]
      // 		<< ","
      // 		<< GMV::gmv_data.doubledata2[i]
      // 		<< ","
      // 		<< GMV::gmv_data.doubledata3[i]
      // 		<< ")"
      // 		<< std::endl;
      
      // Add the point to the Mesh
      MeshInput<MeshBase>::mesh().add_point( Point(GMV::gmv_data.doubledata1[i],
						   GMV::gmv_data.doubledata2[i],
						   GMV::gmv_data.doubledata3[i]) );
    }
#endif  
}


void GMVIO::_read_one_cell()
{
#ifdef HAVE_GMV
  //   // Debug Info 
  //   std::cout << "gmv_data.datatype="
  // 	    <<  GMV::gmv_data.datatype
  // 	    << std::endl;

  // This is either a REGULAR=111 cell or
  // the ENDKEYWORD=207 of the cells
  {
    bool recognized =
      (GMV::gmv_data.datatype==REGULAR) ||
      (GMV::gmv_data.datatype==ENDKEYWORD);
    
    assert (recognized);
  }

  if (GMV::gmv_data.datatype == REGULAR)
    {
      //       std::cout << "Name of the cell is: "
      // 		<< GMV::gmv_data.name1
      // 		<< std::endl;

      //       std::cout << "Cell has "
      // 		<< GMV::gmv_data.num2
      // 		<< " vertices."
      // 		<< std::endl;

      // We need a mapping from GMV element types to LibMesh
      // ElemTypes.  Basically the reverse of the eletypes
      // std::map above.
      //
      // FIXME: Since Quad9's apparently don't exist for GMV, and since
      // In general we write linear sub-elements to GMV files, we need
      // to be careful to read back in exactly what we wrote out...
      ElemType type = this->_gmv_elem_to_libmesh_elem(GMV::gmv_data.name1);
      
      Elem* elem = Elem::build(type).release();

      // Print out the connectivity information for
      // this cell.
      for (int i = 0; i < GMV::gmv_data.num2; i++)
	{
	  // 	  // Debugging info
	  // 	  std::cout << "Vertex " << i << " is node "
	  // 		    << GMV::gmv_data.longdata1[i]
	  // 		    << std::endl;
	  
	  // Note: Node numbers are 1-based
	  elem->set_node(i) = MeshInput<MeshBase>::mesh().node_ptr(GMV::gmv_data.longdata1[i]-1);
	}

      // Add the newly-created element to the mesh
      MeshInput<MeshBase>::mesh().add_elem(elem);
    }


  if (GMV::gmv_data.datatype == ENDKEYWORD)
    {
      // There isn't a cell to read, so we just return
      return;
    }

#endif
}


ElemType GMVIO::_gmv_elem_to_libmesh_elem(const char* elemname)
{
  //
  // Linear Elements
  //
  if (!std::strncmp(elemname,"line",4))
    return EDGE2;
  
  if (!std::strncmp(elemname,"tri",3))
    return TRI3;

  if (!std::strncmp(elemname,"quad",4))
    return QUAD4;
  
  // FIXME: tet or ptet4?
  if ((!std::strncmp(elemname,"tet",3)) ||
      (!std::strncmp(elemname,"ptet4",5)))
    return TET4;

  // FIXME: hex or phex8?
  if ((!std::strncmp(elemname,"hex",3)) ||
      (!std::strncmp(elemname,"phex8",5)))
    return HEX8;

  // FIXME: prism or pprism6?
  if ((!std::strncmp(elemname,"prism",5)) ||
      (!std::strncmp(elemname,"pprism6",7)))
    return PRISM6;

  //
  // Quadratic Elements
  //
  if (!std::strncmp(elemname,"phex20",6))
    return HEX20;

  if (!std::strncmp(elemname,"phex27",6))
    return HEX27;

  if (!std::strncmp(elemname,"pprism15",8))
    return PRISM15;

  if (!std::strncmp(elemname,"ptet10",6))
    return TET10;

  if (!std::strncmp(elemname,"6tri",4))
    return TRI6;

  if (!std::strncmp(elemname,"8quad",5))
    return QUAD8;

  if (!std::strncmp(elemname,"3line",5))
    return EDGE3;
  
  // Unsupported/Unused types
  // if (!std::strncmp(elemname,"vface2d",7))
  // if (!std::strncmp(elemname,"vface3d",7))
  // if (!std::strncmp(elemname,"pyramid",7))
  // if (!std::strncmp(elemname,"ppyrmd5",7))
  // if (!std::strncmp(elemname,"ppyrmd13",8))

  // If we didn't return yet, then we didn't find the right cell!
  std::cerr << "Uknown/unsupported element: "
	    << elemname
	    << " was read."
	    << std::endl;
  error();
}




void GMVIO::copy_nodal_solution(EquationSystems& es)
{
  // Check for easy return if there isn't any nodal data
  if (_nodal_data.empty())
    {
      std::cerr << "Unable to copy nodal solution: No nodal "
		<< "solution has been read in from file." << std::endl;
      return;
    }

  // Be sure there is at least one system
  assert (es.n_systems());

  // Keep track of variable names which have been found and
  // copied already.  This could be used to prevent us from
  // e.g. copying the same var into 2 different systems ...
  // but this seems unlikely.  Also, it is used to tell if
  // any variables which were read in were not successfully
  // copied to the EquationSystems.
  std::set<std::string> vars_copied;
  
  // For each entry in the nodal data map, try to find a system
  // that has the same variable key name.
  for (unsigned int sys=0; sys<es.n_systems(); ++sys)
    {
      // Get a generic refernence to the current System
      System& system = es.get_system(sys);

      // And a reference to that system's dof_map
      // const DofMap & dof_map = system.get_dof_map();
	
      // For each var entry in the _nodal_data map, try to find
      // that var in the system
      std::map<std::string, std::vector<Number> >::iterator it = _nodal_data.begin();
      const std::map<std::string, std::vector<Number> >::iterator end = _nodal_data.end();
      for (; it != end; ++it)
	{
	  std::string var_name = (*it).first;
	  // std::cout << "Searching for var " << var_name << " in system " << sys << std::endl;

	  if (system.has_variable(var_name))
	    {
	      // Check if there are as many nodes in the mesh as there are entries
	      // in the stored nodal data vector
	      assert ( (*it).second.size() == MeshInput<MeshBase>::mesh().n_nodes() );
	      
	      const unsigned int var_num = system.variable_number(var_name);
	      
	      // std::cout << "Variable "
	      // 			<< var_name
	      // 			<< " is variable "
	      // 			<< var_num 
	      // 			<< " in system " << sys << std::endl;

	      // The only type of nodal data we can read in from GMV is for
	      // linear LAGRANGE type elements.
	      const FEType& fe_type = system.variable_type(var_num);
	      if ((fe_type.order != FIRST) || (fe_type.family != LAGRANGE))
		{
		  std::cerr << "Only FIRST-order LAGRANGE variables can be read from GMV files. "
			    << "Skipping variable " << var_name << std::endl;
		  //error();
		  break;
		}

	      
	      // Loop over the stored vector's entries, inserting them into
	      // the System's solution if appropriate.
	      for (unsigned int i=0; i<(*it).second.size(); ++i)
		{
		  // Since this var came from a GMV file, the index i corresponds to
		  // the (single) DOF value of the current variable for node i.
		  const unsigned int dof_index =
		    MeshInput<MeshBase>::mesh().node_ptr(i)->dof_number(sys,      /*system #*/
									var_num,  /*var # */
									0);       /*component #, always zero for LAGRANGE */

		  // std::cout << "Value " << i << ": "
		  // 			    << (*it).second [i]
		  // 			    << ", dof index="
		  // 			    << dof_index << std::endl;

		  // If the dof_index is local to this processor, set the value
		  if ((dof_index >= system.solution->first_local_index()) &&
		      (dof_index <  system.solution->last_local_index()))
		    system.solution->set (dof_index, (*it).second [i]);
		} // end loop over my GMVIO's copy of the solution

	      // Add the most recently copied var to the set of copied vars
	      vars_copied.insert (var_name);
	    } // end if (system.has_variable)
	} // end for loop over _nodal_data

      // Communicate parallel values before going to the next system.
      system.update();

    } // end loop over all systems


  
  // Warn the user if any GMV variables were not successfully copied over to the EquationSystems object
  {
    std::map<std::string, std::vector<Number> >::iterator it = _nodal_data.begin();
    const std::map<std::string, std::vector<Number> >::iterator end = _nodal_data.end();

    for (; it != end; ++it)
      {
	if (vars_copied.find( (*it).first ) == vars_copied.end())
	  {
	    std::cerr << "Warning: Variable "
		      << (*it).first
		      << " was not copied to the EquationSystems object."
		      << std::endl;
	  }
      }
  }
  
}
