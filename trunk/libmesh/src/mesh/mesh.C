// $Id: mesh.C,v 1.31 2004-01-03 15:37:43 benkirk Exp $

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

// Local includes
#include "mesh.h"
#include "libmesh.h"
#include "mesh_communication.h"
#include "libmesh_logging.h"


#ifdef HAVE_SFCURVES
// prototype for SFC code
namespace sfc {
  extern "C" {
#include "sfcurves.h"
  }
}
#endif



// ------------------------------------------------------------
// Mesh class member functions
Mesh::Mesh (unsigned int d) :
  MeshBase      (d)
{
  assert (libMesh::initialized());
}



Mesh::~Mesh ()
{
  this->clear ();
  
  assert (!libMesh::closed());
}



void Mesh::clear ()
{
  // Clear other data structures
  MeshBase::clear();
}



void Mesh::read (const std::string& name)
{
  START_LOG("read()", "Mesh");

  
  // Read the file based on extension.  Only processor 0
  // needs to read the mesh.  It will then broadcast it and
  // the other processors will pick it up
  if (libMesh::processor_id() == 0)
    {
      if (name.rfind(".mat") < name.size())
	this->read_matlab (name);
      
      else if (name.rfind(".ucd") < name.size())
	this->read_ucd (name);
      
      else if (name.rfind(".exd") < name.size())
	this->read_exd (name);
      
      else if (name.rfind(".xda") < name.size())
	this->read_xdr (name);
      
      else if ((name.rfind(".off")  < name.size()) ||
	       (name.rfind(".ogl")  < name.size()) ||
	       (name.rfind(".oogl") < name.size()))
	this->read_off(name);
      
      else if ((name.rfind(".xdr")  < name.size()) ||
	       (name.rfind(".0000") < name.size()))
	this->read_xdr_binary (name);
      
      else if (name.rfind(".mesh") < name.size())
	this->read_shanee (name);
      
      else if (name.rfind(".unv") < name.size())
	this->read_unv (name);
      
      else if ((name.rfind(".node")  < name.size()) ||
	       (name.rfind(".ele")   < name.size()))
	this->read_tetgen (name);
      
      else
	{
	  std::cerr << " ERROR: Unrecognized file extension: " << name
		    << "\n   I understand the following:\n\n"
		    << "     *.mat  -- Matlab triangular ASCII file\n"
		    << "     *.ucd  -- AVS's ASCII UCD format\n"
		    << "     *.mesh -- Ben's \"shanee\" format\n"
		    << "     *.off  -- OOGL OFF surface format\n"
		    << "     *.exd  -- Sandia's ExodusII format\n"
		    << "     *.xda  -- Internal ASCII format\n"
		    << "     *.xdr  -- Internal binary format,\n"
		    << "               compatible with XdrMGF\n"
		    << "     *.unv  -- I-deas Universal format\n"
		    << std::endl;
	  error();	  
	}    
    }

  // Send the mesh & bcs (which are now only on processor 0) to the other
  // processors
  {
    MeshCommunication mesh_communication;
  
    mesh_communication.distribute (*this);
  }
  
  STOP_LOG("read()", "Mesh");


  // Done reading the mesh.  Now prepare it for use.
  this->prepare_for_use();
}



void Mesh::write (const std::string& name)
{
  START_LOG("write()", "Mesh");
  
  // Write the file based on extension
  {
    if (name.rfind(".dat") < name.size())
      this->write_tecplot (name);
    
    else if (name.rfind(".plt") < name.size())
      this->write_tecplot_binary (name);

    else if (name.rfind(".ucd") < name.size())
      this->write_ucd (name);

    else if (name.rfind(".gmv") < name.size())
      {
	if (this->n_subdomains() > 1)
	  this->write_gmv_binary(name, NULL, NULL, true);
	else
	  this->write_gmv_binary(name);
      }


    else if (name.rfind(".ugrid") < name.size())
      this->write_diva (name);
    
    else if (name.rfind(".xda") < name.size())
      this->write_xdr (name);
    
    else if (name.rfind(".xdr") < name.size())
      this->write_xdr_binary (name);

    else if (name.rfind(".unv") < name.size())
      this->write_unv (name);

    else
      {
	std::cerr << " ERROR: Unrecognized file extension: " << name
		  << "\n   I understand the following:\n\n"
		  << "     *.dat   -- Tecplot ASCII file\n"
		  << "     *.plt   -- Tecplot binary file\n"
		  << "     *.ucd   -- AVS's ASCII UCD format\n"
		  << "     *.ugrid -- Kelly's DIVA ASCII format\n"
		  << "     *.gmv   -- LANL's GMV (General Mesh Viewer) format\n"
		  << "     *.xda   -- Internal ASCII format\n"
		  << "     *.xdr   -- Internal binary format,\n"
		  << "                compatible with XdrMGF\n"
		  << "     *.unv   -- I-deas Universal format\n"
		  << std::endl
		  << "\n Exiting without writing output\n";
      }    
  }
  
  STOP_LOG("write()", "Mesh");
}



void Mesh::write (const std::string& name,
		  std::vector<Number>& v,
		  std::vector<std::string>& vn)
{
  START_LOG("write()", "Mesh");

  // Write the file based on extension
  {
    if (name.rfind(".dat") < name.size())
      this->write_tecplot (name, &v, &vn);
    
    else if (name.rfind(".plt") < name.size())
      this->write_tecplot_binary (name, &v, &vn);

    else if (name.rfind(".gmv") < name.size())
      {
	if (n_subdomains() > 1)
	  this->write_gmv_binary(name, &v, &vn, true);
	else
	  this->write_gmv_binary(name, &v, &vn);
      }    
    else
      {
	std::cerr << " ERROR: Unrecognized file extension: " << name
		  << "\n   I understand the following:\n\n"
		  << "     *.dat  -- Tecplot ASCII file\n"
		  << "     *.plt  -- Tecplot binary file\n"
		  << "     *.gmv  -- LANL's GMV (General Mesh Viewer) format\n"
		  << "\n Exiting without writing output\n";
      }
  }

  STOP_LOG("write()", "Mesh");
}





void Mesh::create_pid_mesh(Mesh& pid_mesh,
			   const unsigned int pid) const
{

  // Issue a warning if the number the number of processors
  // currently available is less that that requested for
  // partitioning.  This is not necessarily an error since
  // you may run on one processor and still partition the
  // mesh into several partitions.
#ifdef DEBUG
  if (this->n_processors() < pid)
    {
      std::cout << "WARNING:  You are creating a "
		<< "mesh for a processor id (="
		<< pid
		<< ") greater than "
		<< "the number of processors available for "
		<< "the calculation. (="
		<< libMesh::n_processors()
		<< ")."
		<< std::endl;
    }
#endif
  
  // Create iterators to loop over the list of elements
  const_active_pid_elem_iterator       it(this->elements_begin(),   pid);
  const const_active_pid_elem_iterator it_end(this->elements_end(), pid);
    
  this->create_submesh (pid_mesh, it, it_end);
}







void Mesh::create_submesh (Mesh& new_mesh,
			   const_elem_iterator& it,
			   const const_elem_iterator& it_end) const
{
 
  // Just in case the subdomain_mesh already has some information
  // in it, get rid of it.
  new_mesh.clear();

  // How the nodes on this mesh will be renumbered to nodes
  // on the new_mesh.  
  std::vector<unsigned int> new_node_numbers (this->n_nodes());

  std::fill (new_node_numbers.begin(),
	     new_node_numbers.end(),
	     libMesh::invalid_uint);

  
  
  // the number of nodes on the new mesh, will be incremented
  unsigned int n_new_nodes = 0;
  unsigned int n_new_elem  = 0;
    
  for (; it != it_end; ++it)
    {
      // increment the new element counter
      n_new_elem++;
	
      const Elem* old_elem = *it;

      // Add an equivalent element type to the new_mesh
      Elem* new_elem = new_mesh.add_elem (Elem::build (old_elem->type()));

      assert (new_elem != NULL);
	
      // Loop over the nodes on this element.  
      for (unsigned int n=0; n<old_elem->n_nodes(); n++)
	{
	  assert (old_elem->node(n) < new_node_numbers.size());

	  if (new_node_numbers[old_elem->node(n)] == libMesh::invalid_uint)
	    {
	      new_node_numbers[old_elem->node(n)] = n_new_nodes;

	      // Add this node to the new mesh
	      new_mesh.add_point (old_elem->point(n));

	      // Increment the new node counter
	      n_new_nodes++;
	    }

	  // Define this element's connectivity on the new mesh
	  assert (new_node_numbers[old_elem->node(n)] < new_mesh.n_nodes());
	    
	  new_elem->set_node(n) = new_mesh.node_ptr (new_node_numbers[old_elem->node(n)]);
	}

      // Maybe add boundary conditions for this element
      for (unsigned int s=0; s<old_elem->n_sides(); s++)
	if (old_elem->neighbor(s) == NULL)
	  if (this->boundary_info.boundary_id (old_elem, s) !=
	      this->boundary_info.invalid_id)
	    new_mesh.boundary_info.add_side (new_elem,
					     s,
					     this->boundary_info.boundary_id (old_elem, s));
    } // end loop over elements
  

  // Prepare the new_mesh for use
  new_mesh.prepare_for_use();
  
}
