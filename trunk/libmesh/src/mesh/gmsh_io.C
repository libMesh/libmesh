// $Id: gmsh_io.C,v 1.8 2005-02-22 22:17:39 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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

// Local includes
#include "libmesh_config.h"
#include "gmsh_io.h"
#include "elem.h"


// ------------------------------------------------------------
// GmshIO  members
void GmshIO::read (const std::string& name)
{
  std::ifstream in (name.c_str());

  this->read_stream (in);
}




void GmshIO::read_stream(std::istream& in)
{

  assert(in.good());
  // Not yet implemented
  error();
}



void GmshIO::write (const std::string& name)
{
  if (libMesh::processor_id() == 0)
    {
      std::ofstream out (name.c_str());
      this->write_stream (out);
    }
}




void GmshIO::write_stream (std::ostream& out)
{
  // Be sure that the stream is valid.
  assert (out.good());
  
  // Get a const reference to the mesh
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();
  
  // Note: we are using version 2.0 of the gmsh output format.
  
  {
    // Write the file header.
    out << "$MeshFormat\n";
    out << "2.0 0 " << sizeof(Real) << '\n';
    out << "$EndMeshFormat\n";
  }

  {
    // write the nodes in (n x y z) format
    out << "$Nodes\n";
    out << mesh.n_nodes() << '\n';
    
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      out << mesh.node(v).id()+1 << " "
	  << mesh.node(v)(0) << " "
	  << mesh.node(v)(1) << " "
	  << mesh.node(v)(2) << '\n';
    out << "$EndNodes\n";
  }


  
  {
    // write the connectivity
    
    out << "$Elements\n";
    out << mesh.n_active_sub_elem() << '\n';


//     const_active_elem_iterator       it (mesh.elements_begin());
//     const const_active_elem_iterator end(mesh.elements_end());

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end(); 
    
    switch (mesh.mesh_dimension())
      {

	// 1D Meshes
      case 1:
	{
// 	  for ( ; it != end; ++it)
// 	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
// 	      {
// 		out << "line 2\n";
// 		std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
// 		for (unsigned int i=0; i<conn.size(); i++)
// 		  out << conn[i] << " ";
		
// 		out << '\n';
// 	      }
	  error();  // Not yet implemented
	  break;
	}



	
	// 2D Meshes
      case 2:
	{
	  // The same temporary storage will be used for each element
	  std::vector<unsigned int> conn;

	  unsigned int ctr=1;
	  for ( ; it != end; ++it)
	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	      {
		// Write element number
		out << ctr++ << " ";

		// Quadrilateral elements
		if (((*it)->type() == QUAD4) ||
		    ((*it)->type() == QUAD8) ||
		    ((*it)->type() == QUAD9)
#ifdef ENABLE_INFINITE_ELEMENTS
		    || ((*it)->type() == INFQUAD4)
		    || ((*it)->type() == INFQUAD6)
#endif
		    )
		  {
		    out << "3 ";
		    
		    // Write the number of tags 
		    out << "0 ";

		    // Write connectivity
		    (*it)->connectivity(se, TECPLOT, conn);
		    for (unsigned int i=0; i<conn.size(); i++)
		      out << conn[i] << " ";

		  }

		// Triangular elements
		else if (((*it)->type() == TRI3) ||
			 ((*it)->type() == TRI6))
		  {
		    out << "2 ";
		    
		    // Write the number of tags 
		    out << "0 ";
		    
		    // Write connectivity (last entry of conn is skipped!)
		    (*it)->connectivity(se, TECPLOT, conn);
		    for (unsigned int i=0; i<3; i++)
		      out << conn[i] << " ";
		  }
		
		else
		  {
		    error();
		  }

		// Go to next line
		out << '\n';
	      }
	  
	  break;
	}
	
	
      case 3:
	{
// 	  for ( ; it != end; ++it)
// 	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
// 	      {

// #ifndef  ENABLE_INFINITE_ELEMENTS
// 		if (((*it)->type() == HEX8)   ||    
// 		    ((*it)->type() == HEX27))
// 		  {
// 		    out << "phex8 8\n";
// 		    std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
// 		    for (unsigned int i=0; i<conn.size(); i++)
// 		      out << conn[i] << " ";
// 		  }
		
// 		else if ((*it)->type() == HEX20)
// 		  {
// 		    out << "phex20 20\n";
// 		    out << (*it)->node(0)+1  << " "
// 			<< (*it)->node(1)+1  << " "
// 			<< (*it)->node(2)+1  << " "
// 			<< (*it)->node(3)+1  << " "
// 			<< (*it)->node(4)+1  << " "
// 			<< (*it)->node(5)+1  << " "
// 			<< (*it)->node(6)+1  << " "
// 			<< (*it)->node(7)+1  << " "
// 			<< (*it)->node(8)+1  << " "
// 			<< (*it)->node(9)+1  << " "
// 			<< (*it)->node(10)+1 << " "
// 			<< (*it)->node(11)+1 << " "
// 			<< (*it)->node(16)+1 << " "
// 			<< (*it)->node(17)+1 << " "
// 			<< (*it)->node(18)+1 << " "
// 			<< (*it)->node(19)+1 << " "
// 			<< (*it)->node(12)+1 << " "
// 			<< (*it)->node(13)+1 << " "
// 			<< (*it)->node(14)+1 << " "
// 			<< (*it)->node(15)+1 << " ";
// 		  }
// #else
// 		/*
// 		 * In case of infinite elements, HEX20
// 		 * should be handled just like the
// 		 * INFHEX16, since these connect to each other
// 		 */
// 		if (((*it)->type() == HEX8)     ||
// 		    ((*it)->type() == HEX27)    ||
// 		    ((*it)->type() == INFHEX8)  ||
// 		    ((*it)->type() == INFHEX16) ||
// 		    ((*it)->type() == INFHEX18) ||
// 		    ((*it)->type() == HEX20))
// 		  {
// 		    out << "phex8 8\n";
// 		    std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
// 		    for (unsigned int i=0; i<conn.size(); i++)
// 		      out << conn[i] << " ";
// 		  }
// #endif
		
// 		else if (((*it)->type() == TET4)  ||
// 			 ((*it)->type() == TET10))
// 		  {
// 		    out << "tet 4\n";
// 		    std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
// 		    out << conn[0] << " "
// 			<< conn[2] << " "
// 			<< conn[1] << " "
// 			<< conn[4] << " ";
// 		    }
// #ifndef  ENABLE_INFINITE_ELEMENTS
// 		else if (((*it)->type() == PRISM6)  ||
// 			 ((*it)->type() == PRISM15) ||
// 			 ((*it)->type() == PRISM18))
// #else
// 		else if (((*it)->type() == PRISM6)     ||
// 			 ((*it)->type() == PRISM15)    ||
// 			 ((*it)->type() == PRISM18)    ||
// 			 ((*it)->type() == INFPRISM6)  ||
// 			 ((*it)->type() == INFPRISM12))
// #endif
// 		  {
// 		    /**
// 		     * Note that the prisms are treated as
// 		     * degenerated phex8's.
// 		     */
// 		    out << "phex8 8\n";
// 		    std::vector<unsigned int> conn = (*it)->tecplot_connectivity(se);
// 		    for (unsigned int i=0; i<conn.size(); i++)
// 		      out << conn[i] << " ";
// 		  }
		
// 		else
// 		  {
// 		    error();
// 		  }
		
// 		out << '\n';
// 	      }
	  
	  error(); // Not yet implemented
	  break;
	}
	
      default:
	error();
      }
    
  }

  
  // end of the file
  out << '\n' << "$EndElements\n";
}




