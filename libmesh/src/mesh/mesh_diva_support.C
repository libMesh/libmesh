// $Id: mesh_diva_support.C,v 1.12 2004-03-18 15:43:37 jwpeterson Exp $

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


// Local includes
#include "libmesh_common.h"
#include "mesh.h"
#include "boundary_mesh.h"



void Mesh::write_diva (const std::string& name) const
{
  std::ofstream out(name.c_str());

  write_diva (out);
}



void Mesh::write_diva (std::ostream& out) const
{
  /*
    From Kelly: (kelly@tacc.utexas.edu)

    Ok, the following is the format:

    #points #triangles #quads #tets #prisms #pyramids #hexs
    loop over all points (written out x y z x y z ...)
    loop over all triangles (written out i1 i2 i3) (These are indices into
    the points going from
    1 to #points) 
    loop over all quads (written out i1 i2 i3 i4) (Same numbering scheme)
    loop over all triangles and quads (write out b1) (This is a boundary
    condition for each
    triangle and each
    hex. You can put
    anything you want
    here)
    loop over all tets (written out i1 i2 i3 i4) (Same)
    loop over all pyramids (written out i1 i2 i3 i4 i5) (Same)
    loop over all prisms (written out i1 i2 i3 i4 i5 i6) (Same)
    loop over all hexs (written out i1 i2 i3 i4 i5 i6 i7 i8) (Same)
    
  */
  
  assert (out);


  if (_dim < 3)
    {
      std::cerr << "WARNING: DIVA only supports 3D meshes.\n\n"
		<< "Exiting without producing output.\n";
      return;
    }
  


  BoundaryMesh boundary_mesh (this->mesh_dimension()-1);
  boundary_info.sync(boundary_mesh);
  

  /**
   * Write the header
   */ 
  {
    out << this->n_nodes() << " "
      
	<< (boundary_mesh.n_active_elem_of_type(TRI3) +
	    boundary_mesh.n_active_elem_of_type(TRI6)*4) << " "
      
	<< (boundary_mesh.n_active_elem_of_type(QUAD4)  +
	    boundary_mesh.n_active_elem_of_type(QUAD8) +
	    boundary_mesh.n_active_elem_of_type(QUAD9)*4) << " "
      
	<< (n_active_elem_of_type(TET4) +
	    n_active_elem_of_type(TET10)*8) << " "
      
	<< n_active_elem_of_type(PYRAMID5) << " "
      
	<< (n_active_elem_of_type(PRISM6) +
	    n_active_elem_of_type(PRISM18)*8) << " "
      
	<< (n_active_elem_of_type(HEX8)   +
	    n_active_elem_of_type(HEX20) +
	    n_active_elem_of_type(HEX27)*8) << " "
      
	<< std::endl;
  }
  

  boundary_mesh.clear();

  
  /**
   * Write the nodes
   */ 
  for (unsigned int v=0; v<n_nodes(); v++)
    point(v).write_unformatted(out);
  
  
  /**
   * Write the BC faces
   */
  {
    /**
     * Write the triangles
     */
    for(unsigned int e=0; e<n_elem(); e++)
      if (elem(e)->active())
	for (unsigned int s=0; s<elem(e)->n_sides(); s++)
	  if (elem(e)->neighbor(s) == NULL)
	    {
	      const AutoPtr<Elem> side(elem(e)->build_side(s));

	      if (side->type() == TRI3)
		{
		  out << side->node(0)+1 << " "
		      << side->node(1)+1 << " "
		      << side->node(2)+1 << std::endl;
		}
	      else if (side->type() == TRI6)
		{
		  out << side->node(0)+1 << " "
		      << side->node(3)+1 << " "
		      << side->node(5)+1 << std::endl

		      << side->node(3)+1 << " "
		      << side->node(1)+1 << " "
		      << side->node(4)+1 << std::endl

		      << side->node(5)+1 << " "
		      << side->node(4)+1 << " "
		      << side->node(2)+1 << std::endl

		      << side->node(3)+1 << " "
		      << side->node(4)+1 << " "
		      << side->node(5)+1 << std::endl;
		}
	    }

    
    /**
     * Write the quadrilaterals
     */
    for(unsigned int e=0; e<n_elem(); e++)
      if (elem(e)->active())
	for (unsigned int s=0; s<elem(e)->n_sides(); s++)
	  if (elem(e)->neighbor(s) == NULL)
	    {
	      const AutoPtr<Elem> side(elem(e)->build_side(s));

	      if ((side->type() == QUAD4) ||
		  (side->type() == QUAD8)  )		
		{
		  out << side->node(0)+1 << " "
		      << side->node(1)+1 << " "
		      << side->node(2)+1 << " "
		      << side->node(3)+1 << std::endl;
		}
	      else if (side->type() == QUAD9)
		{
		  out << side->node(0)+1 << " "
		      << side->node(4)+1 << " "
		      << side->node(8)+1 << " "
		      << side->node(7)+1 << std::endl

		      << side->node(4)+1 << " "
		      << side->node(1)+1 << " "
		      << side->node(5)+1 << " "
		      << side->node(8)+1 << std::endl

		      << side->node(7)+1 << " "
		      << side->node(8)+1 << " "
		      << side->node(6)+1 << " "
		      << side->node(3)+1 << std::endl

		      << side->node(8)+1 << " "
		      << side->node(5)+1 << " "
		      << side->node(2)+1 << " "
		      << side->node(6)+1 << std::endl;
		}
	    }
  }
  
	  

  /**
   * Write the BC IDs
   */
  {
    /**
     * Write the triangles
     */
    for(unsigned int e=0; e<n_elem(); e++)
      if (elem(e)->active())
	for (unsigned int s=0; s<elem(e)->n_sides(); s++)
	  if (elem(e)->neighbor(s) == NULL)
	    {
	      const AutoPtr<Elem> side(elem(e)->build_side(s));
	      
	      if ((side->type() == TRI3) ||
		  (side->type() == TRI6)  )

		out << boundary_info.boundary_id(elem(e), s)
		    << std::endl;
	    }

    
    /**
     * Write the quadrilaterals
     */
    for(unsigned int e=0; e<n_elem(); e++)
      if (elem(e)->active())
	for (unsigned int s=0; s<elem(e)->n_sides(); s++)
	  if (elem(e)->neighbor(s) == NULL)
	    {
	      const AutoPtr<Elem> side(elem(e)->build_side(s));
	      
	      if ((side->type() == QUAD4)  ||
		  (side->type() == QUAD8) ||
		  (side->type() == QUAD9)  )
		
		out << boundary_info.boundary_id(elem(e), s);
	    }
  }


  
  /**
   * Write all the Tets
   */
  for (unsigned int e=0; e<n_elem(); e++)
    if (elem(e)->active())
      if (elem(e)->type() == TET4)
	{
	  out << elem(e)->node(0)+1 << " "
	      << elem(e)->node(1)+1 << " "
	      << elem(e)->node(2)+1 << " "
	      << elem(e)->node(3)+1 << std::endl;
	}
      else if (elem(e)->type() == TET10)
	{
	  out << elem(e)->node(0)+1 << " "
	      << elem(e)->node(4)+1 << " "
	      << elem(e)->node(6)+1 << " "
	      << elem(e)->node(7)+1 << std::endl;
	    
	  out << elem(e)->node(4)+1 << " "
	      << elem(e)->node(1)+1 << " "
	      << elem(e)->node(5)+1 << " "
	      << elem(e)->node(8)+1 << std::endl;	
	    
	  out << elem(e)->node(6)+1 << " "
	      << elem(e)->node(5)+1 << " "
	      << elem(e)->node(2)+1 << " "
	      << elem(e)->node(9)+1 << std::endl;	
	    
	  out << elem(e)->node(7)+1 << " "
	      << elem(e)->node(8)+1 << " "
	      << elem(e)->node(9)+1 << " "
	      << elem(e)->node(3)+1 << std::endl;	
	    
	  out << elem(e)->node(4)+1 << " "
	      << elem(e)->node(8)+1 << " "
	      << elem(e)->node(6)+1 << " "
	      << elem(e)->node(7)+1 << std::endl;	
	    
	  out << elem(e)->node(4)+1 << " "
	      << elem(e)->node(5)+1 << " "
	      << elem(e)->node(6)+1 << " "
	      << elem(e)->node(8)+1 << std::endl;	
	    
	  out << elem(e)->node(6)+1 << " "
	      << elem(e)->node(5)+1 << " "
	      << elem(e)->node(9)+1 << " "
	      << elem(e)->node(8)+1 << std::endl;	
	  
	  out << elem(e)->node(6)+1 << " "
	      << elem(e)->node(8)+1 << " "
	      << elem(e)->node(9)+1 << " "
	      << elem(e)->node(7)+1 << std::endl;	
	}



  /**
   * Write all the Pyramids
   */
  for (unsigned int e=0; e<n_elem(); e++)
    if (elem(e)->active())
      if (elem(e)->type() == PYRAMID5)
	{
	  out << elem(e)->node(0)+1 << " "
	      << elem(e)->node(1)+1 << " "
	      << elem(e)->node(2)+1 << " "
	      << elem(e)->node(3)+1 << " "
	      << elem(e)->node(4)+1 << std::endl;
	}



  /**
   * Write all the Prisms
   */
  for (unsigned int e=0; e<n_elem(); e++)
    if (elem(e)->active())
      if (elem(e)->type() == PRISM6)
	{
	  out << elem(e)->node(0)+1 << " "
	      << elem(e)->node(1)+1 << " "
	      << elem(e)->node(2)+1 << " "
	      << elem(e)->node(3)+1 << " "
	      << elem(e)->node(4)+1 << " "
	      << elem(e)->node(5)+1 << std::endl;
	}
      else if (elem(e)->type() == PRISM18)
	{
	  error();
	}



  /**
   * Write all the Hexes
   */
  for (unsigned int e=0; e<n_elem(); e++)
    if (elem(e)->active())
      if ((elem(e)->type() == HEX8)   ||
	  (elem(e)->type() == HEX20) ||
	  (elem(e)->type() == HEX27)   )
	{
	  for (unsigned int se=0; se<elem(e)->n_sub_elem(); se++)
	    {
	      std::vector<unsigned int> conn =
		elem(e)->tecplot_connectivity(se);

	      out << conn[0] << " "
		  << conn[1] << " "
		  << conn[2] << " "
		  << conn[3] << " "
		  << conn[4] << " "
		  << conn[5] << " "
		  << conn[6] << " "
		  << conn[7] << std::endl;
	    }
	}
}
