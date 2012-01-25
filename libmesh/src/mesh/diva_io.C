// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "diva_io.h"
#include "boundary_mesh.h"
#include "mesh_tools.h"
#include "elem.h"
#include "boundary_info.h"

namespace libMesh
{

// ------------------------------------------------------------
// DivaIO class members
void DivaIO::write (const std::string& fname)
{
  // We may need to gather a ParallelMesh to output it, making that
  // const qualifier in our constructor a dirty lie
  MeshSerializer serialize(const_cast<MeshBase&>(this->mesh()), !_is_parallel_format);

  // Open the output file stream
  std::ofstream out(fname.c_str());

  // Make sure it opened correctly
  if (!out.good())
    libmesh_file_error(fname.c_str());

  this->write_stream (out);
}




void DivaIO::write_stream (std::ostream& out)
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

  // Be sure the stream has been created successfully.
  libmesh_assert (out.good());

  // Can't use a constant mesh reference since we have to
  // sync the boundary info.
  libmesh_here();
  libMesh::err << "WARNING...  Sure you want to do this?"
	        << std::endl;
  MeshBase& mesh = const_cast<MeshBase&>
    (MeshOutput<MeshBase>::mesh());

  if (mesh.mesh_dimension() < 3)
    {
      libMesh::err << "WARNING: DIVA only supports 3D meshes.\n\n"
		    << "Exiting without producing output.\n";
      return;
    }



  BoundaryMesh boundary_mesh (mesh.mesh_dimension()-1);
  mesh.boundary_info->sync(boundary_mesh);


  /**
   * Write the header
   */
  out << mesh.n_nodes() << " "

      << (MeshTools::n_active_elem_of_type(boundary_mesh,TRI3) +
	  MeshTools::n_active_elem_of_type(boundary_mesh,TRI6)*4) << " "

      << (MeshTools::n_active_elem_of_type(boundary_mesh, QUAD4) +
	  MeshTools::n_active_elem_of_type(boundary_mesh, QUAD8) +
	  MeshTools::n_active_elem_of_type(boundary_mesh, QUAD9)*4) << " "

      << (MeshTools::n_active_elem_of_type(mesh, TET4) +
	  MeshTools::n_active_elem_of_type(mesh, TET10)*8) << " "

      << MeshTools::n_active_elem_of_type(mesh, PYRAMID5) << " "

      << (MeshTools::n_active_elem_of_type(mesh, PRISM6) +
	  MeshTools::n_active_elem_of_type(mesh, PRISM18)*8) << " "

      << (MeshTools::n_active_elem_of_type(mesh, HEX8)   +
	  MeshTools::n_active_elem_of_type(mesh, HEX20) +
	  MeshTools::n_active_elem_of_type(mesh, HEX27)*8) << " "

      << '\n';


  boundary_mesh.clear();


  /**
   * Write the nodes
   */
  for (unsigned int v=0; v<mesh.n_nodes(); v++)
    mesh.point(v).write_unformatted(out);


  /**
   * Write the BC faces
   */
  {
    /**
     * Write the triangles
     */
    for(unsigned int e=0; e<mesh.n_elem(); e++)
      if (mesh.elem(e)->active())
	for (unsigned int s=0; s<mesh.elem(e)->n_sides(); s++)
	  if (mesh.elem(e)->neighbor(s) == NULL)
	    {
	      const AutoPtr<Elem> side(mesh.elem(e)->build_side(s));

	      if (side->type() == TRI3)
		{
		  out << side->node(0)+1 << " "
		      << side->node(1)+1 << " "
		      << side->node(2)+1 << '\n';
		}
	      else if (side->type() == TRI6)
		{
		  out << side->node(0)+1 << " "
		      << side->node(3)+1 << " "
		      << side->node(5)+1 << '\n'

		      << side->node(3)+1 << " "
		      << side->node(1)+1 << " "
		      << side->node(4)+1 << '\n'

		      << side->node(5)+1 << " "
		      << side->node(4)+1 << " "
		      << side->node(2)+1 << '\n'

		      << side->node(3)+1 << " "
		      << side->node(4)+1 << " "
		      << side->node(5)+1 << '\n';
		}
	    }


    /**
     * Write the quadrilaterals
     */
    for(unsigned int e=0; e<mesh.n_elem(); e++)
      if (mesh.elem(e)->active())
	for (unsigned int s=0; s<mesh.elem(e)->n_sides(); s++)
	  if (mesh.elem(e)->neighbor(s) == NULL)
	    {
	      const AutoPtr<Elem> side(mesh.elem(e)->build_side(s));

	      if ((side->type() == QUAD4) ||
		  (side->type() == QUAD8)  )
		{
		  out << side->node(0)+1 << " "
		      << side->node(1)+1 << " "
		      << side->node(2)+1 << " "
		      << side->node(3)+1 << '\n';
		}
	      else if (side->type() == QUAD9)
		{
		  out << side->node(0)+1 << " "
		      << side->node(4)+1 << " "
		      << side->node(8)+1 << " "
		      << side->node(7)+1 << '\n'

		      << side->node(4)+1 << " "
		      << side->node(1)+1 << " "
		      << side->node(5)+1 << " "
		      << side->node(8)+1 << '\n'

		      << side->node(7)+1 << " "
		      << side->node(8)+1 << " "
		      << side->node(6)+1 << " "
		      << side->node(3)+1 << '\n'

		      << side->node(8)+1 << " "
		      << side->node(5)+1 << " "
		      << side->node(2)+1 << " "
		      << side->node(6)+1 << '\n';
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
    for(unsigned int e=0; e<mesh.n_elem(); e++)
      if (mesh.elem(e)->active())
	for (unsigned int s=0; s<mesh.elem(e)->n_sides(); s++)
	  if (mesh.elem(e)->neighbor(s) == NULL)
	    {
	      const AutoPtr<Elem> side(mesh.elem(e)->build_side(s));

	      if ((side->type() == TRI3) ||
		  (side->type() == TRI6)  )

		out << mesh.boundary_info->boundary_id(mesh.elem(e), s)
		    << '\n';
	    }


    /**
     * Write the quadrilaterals
     */
    for(unsigned int e=0; e<mesh.n_elem(); e++)
      if (mesh.elem(e)->active())
	for (unsigned int s=0; s<mesh.elem(e)->n_sides(); s++)
	  if (mesh.elem(e)->neighbor(s) == NULL)
	    {
	      const AutoPtr<Elem> side(mesh.elem(e)->build_side(s));

	      if ((side->type() == QUAD4)  ||
		  (side->type() == QUAD8) ||
		  (side->type() == QUAD9)  )

		out << mesh.boundary_info->boundary_id(mesh.elem(e), s);
	    }
  }



  /**
   * Write all the Tets
   */
  for (unsigned int e=0; e<mesh.n_elem(); e++)
    if (mesh.elem(e)->active())
      {
	if (mesh.elem(e)->type() == TET4)
	  {
	    out << mesh.elem(e)->node(0)+1 << " "
		<< mesh.elem(e)->node(1)+1 << " "
		<< mesh.elem(e)->node(2)+1 << " "
		<< mesh.elem(e)->node(3)+1 << '\n';
	  }
	else if (mesh.elem(e)->type() == TET10)
	  {
	    out << mesh.elem(e)->node(0)+1 << " "
		<< mesh.elem(e)->node(4)+1 << " "
		<< mesh.elem(e)->node(6)+1 << " "
		<< mesh.elem(e)->node(7)+1 << '\n';

	    out << mesh.elem(e)->node(4)+1 << " "
		<< mesh.elem(e)->node(1)+1 << " "
		<< mesh.elem(e)->node(5)+1 << " "
		<< mesh.elem(e)->node(8)+1 << '\n';

	  out << mesh.elem(e)->node(6)+1 << " "
	      << mesh.elem(e)->node(5)+1 << " "
	      << mesh.elem(e)->node(2)+1 << " "
	      << mesh.elem(e)->node(9)+1 << '\n';

	  out << mesh.elem(e)->node(7)+1 << " "
	      << mesh.elem(e)->node(8)+1 << " "
	      << mesh.elem(e)->node(9)+1 << " "
	      << mesh.elem(e)->node(3)+1 << '\n';

	  out << mesh.elem(e)->node(4)+1 << " "
	      << mesh.elem(e)->node(8)+1 << " "
	      << mesh.elem(e)->node(6)+1 << " "
	      << mesh.elem(e)->node(7)+1 << '\n';

	  out << mesh.elem(e)->node(4)+1 << " "
	      << mesh.elem(e)->node(5)+1 << " "
	      << mesh.elem(e)->node(6)+1 << " "
	      << mesh.elem(e)->node(8)+1 << '\n';

	  out << mesh.elem(e)->node(6)+1 << " "
	      << mesh.elem(e)->node(5)+1 << " "
	      << mesh.elem(e)->node(9)+1 << " "
	      << mesh.elem(e)->node(8)+1 << '\n';

	  out << mesh.elem(e)->node(6)+1 << " "
	      << mesh.elem(e)->node(8)+1 << " "
	      << mesh.elem(e)->node(9)+1 << " "
	      << mesh.elem(e)->node(7)+1 << '\n';
	}
      }


  /**
   * Write all the Pyramids
   */
  for (unsigned int e=0; e<mesh.n_elem(); e++)
    if (mesh.elem(e)->active())
      if (mesh.elem(e)->type() == PYRAMID5)
	{
	  out << mesh.elem(e)->node(0)+1 << " "
	      << mesh.elem(e)->node(1)+1 << " "
	      << mesh.elem(e)->node(2)+1 << " "
	      << mesh.elem(e)->node(3)+1 << " "
	      << mesh.elem(e)->node(4)+1 << '\n';
	}



  /**
   * Write all the Prisms
   */
  for (unsigned int e=0; e<mesh.n_elem(); e++)
    if (mesh.elem(e)->active())
      {
	if (mesh.elem(e)->type() == PRISM6)
	  {
	    out << mesh.elem(e)->node(0)+1 << " "
		<< mesh.elem(e)->node(1)+1 << " "
		<< mesh.elem(e)->node(2)+1 << " "
		<< mesh.elem(e)->node(3)+1 << " "
		<< mesh.elem(e)->node(4)+1 << " "
		<< mesh.elem(e)->node(5)+1 << '\n';
	  }
	else if (mesh.elem(e)->type() == PRISM18)
	  {
	    libmesh_error();
	  }
      }


  /**
   * Write all the Hexes
   */
  for (unsigned int e=0; e<mesh.n_elem(); e++)
    if (mesh.elem(e)->active())
      if ((mesh.elem(e)->type() == HEX8)   ||
	  (mesh.elem(e)->type() == HEX20) ||
	  (mesh.elem(e)->type() == HEX27)   )
	{
	  std::vector<unsigned int> conn;
	  for (unsigned int se=0; se<mesh.elem(e)->n_sub_elem(); se++)
	    {
	      mesh.elem(e)->connectivity(se, TECPLOT, conn);

	      out << conn[0] << " "
		  << conn[1] << " "
		  << conn[2] << " "
		  << conn[3] << " "
		  << conn[4] << " "
		  << conn[5] << " "
		  << conn[6] << " "
		  << conn[7] << '\n';
	    }
	}
}

} // namespace libMesh
