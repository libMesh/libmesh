// $Id: mesh.C,v 1.64 2005-08-16 13:35:30 benkirk Exp $

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

// Local includes
#include "mesh.h"
#include "mesh_communication.h"
#include "libmesh_logging.h"
#include "elem.h"
#include "boundary_info.h"

#include "diva_io.h"
#include "exodusII_io.h"
#include "gmv_io.h"
#include "tecplot_io.h"
#include "tetgen_io.h"
#include "ucd_io.h"
#include "unv_io.h"
#include "matlab_io.h"
#include "off_io.h"
#include "medit_io.h"
#include "gmsh_io.h"
#include "xdr_io.h"

#if   defined(HAVE_HASH_MAP)
# include <hash_map>
#elif defined(HAVE_EXT_HASH_MAP)
# include <ext/hash_map>
#endif


// ------------------------------------------------------------
// Mesh class member functions
Mesh::Mesh (unsigned int d) :
  MeshBase (d)
{
  assert (libMesh::initialized());
}


Mesh::Mesh (const Mesh& other_mesh) :
  MeshBase (other_mesh),
  _nodes             (other_mesh._nodes),
  _elements          (other_mesh._elements)
{

}


Mesh::~Mesh ()
{
  this->clear ();
  
  assert (!libMesh::closed());
}





const Point& Mesh::point (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);  

  return (*_nodes[i]);
}





const Node& Mesh::node (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);  
  
  return (*_nodes[i]);
}





Node& Mesh::node (const unsigned int i)
{
  if (i >= this->n_nodes())
    {
      std::cout << " i=" << i
		<< ", n_nodes()=" << this->n_nodes()
		<< std::endl;
      error();
    }
  
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);

  return (*_nodes[i]);
}



const Node* Mesh::node_ptr (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);  
  
  return _nodes[i];
}




Node* & Mesh::node_ptr (const unsigned int i)
{
  assert (i < this->n_nodes());

  return _nodes[i];
}




Elem* Mesh::elem (const unsigned int i) const
{
  assert (i < this->n_elem());
  assert (_elements[i] != NULL);
  
  return _elements[i];
}






void Mesh::clear ()
{
  // Call parent clear function
  MeshBase::clear();

  
  // Clear our elements and nodes
  {
    std::vector<Elem*>::iterator       it  = _elements.begin();
    const std::vector<Elem*>::iterator end = _elements.end();

    // There is no need to remove the elements from
    // the BoundaryInfo data structure since we
    // already cleared it.
    for (; it != end; ++it)
      delete *it;

    _elements.clear();
  }

  // clear the nodes data structure
  {
    std::vector<Node*>::iterator       it  = _nodes.begin();
    const std::vector<Node*>::iterator end = _nodes.end();

    // There is no need to remove the nodes from
    // the BoundaryInfo data structure since we
    // already cleared it.
    for (; it != end; ++it)
      delete *it;
    
    _nodes.clear();
  }
}





Node* Mesh::add_point (const Point& p)
{  
  _nodes.push_back (Node::build(p, this->n_nodes()).release());
  
  return _nodes.back();
}



Elem* Mesh::add_elem (Elem* e)
{
  if (e != NULL)
    e->set_id (_elements.size());
  
  _elements.push_back(e);

  return e;
}



void Mesh::delete_elem(Elem* e)
{
  assert (e != NULL);

  // Initialize an iterator to eventually point to the element we want to delete
  std::vector<Elem*>::iterator pos = _elements.end();
  
  // In many cases, e->id() gives us a clue as to where e
  // is located in the _elements vector.  Try that first
  // before trying the O(n_elem) search.
  assert (e->id() < _elements.size());

  if (_elements[e->id()] == e)
    {
      // We found it!
      pos = _elements.begin();
      std::advance(pos, e->id());
    }

  else
    {
      // This search is O(n_elem)
      pos = std::find (_elements.begin(),
		       _elements.end(),
		       e);
    }

  // Huh? Element not in the vector?
  assert (pos != _elements.end());

  // delete the element
  delete e;
  
  // explicitly NULL the pointer
  e    = NULL;
  *pos = NULL;
}



void Mesh::delete_node(Node* n)
{
  assert (n != NULL);
  
  std::vector<Node*>::iterator pos = std::find (_nodes.begin(),
						_nodes.end(),
						n);
  
  // Huh? Node not in the vector?
  assert (pos != _nodes.end());
  
  // delete the element
  delete n;
  
  // explicitly NULL the pointer
  n    = NULL;
  *pos = NULL;
}








void Mesh::find_neighbors()
{
  assert(this->n_nodes() != 0);
  assert(this->n_elem()  != 0);

  START_LOG("find_neighbors()", "Mesh");
  
  
  //TODO:[BSK] This should be removed later?!
  for (unsigned int e=0; e<this->n_elem(); e++)
    for (unsigned int s=0; s<this->elem(e)->n_neighbors(); s++)
      this->elem(e)->set_neighbor(s,NULL);
  
  // Find neighboring elements by first finding elements
  // with identical side keys and then check to see if they
  // are neighbors
  {
    // data structures -- Use the hash_multimap if available
    typedef unsigned int                    key_type;
    typedef std::pair<Elem*, unsigned char> val_type;
    typedef std::pair<key_type, val_type>   key_val_pair;
    
#if   defined(HAVE_HASH_MAP)    
    typedef std::hash_multimap<key_type, val_type> map_type;    
#elif defined(HAVE_EXT_HASH_MAP)
# if    (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0   
    typedef std::hash_multimap<key_type, val_type> map_type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
    typedef __gnu_cxx::hash_multimap<key_type, val_type> map_type;
# else
# error     DIE A HORRIBLE DEATH
# endif
#else
    typedef std::multimap<key_type, val_type>  map_type;
#endif
    
    // A map from side keys to corresponding elements & side numbers  
    map_type side_to_elem_map;
  

    element_iterator       el  = this->elements_begin();
    const element_iterator end = this->elements_end();


    for (; el != end; ++el)
      {
	Elem* element = *el;
	
	for (unsigned int ms=0; ms<element->n_neighbors(); ms++)
	  {
	  next_side:
	    
	    if (element->neighbor(ms) == NULL)
	      {
		// Get the key for the side of this element
		const unsigned int key = element->key(ms);
		
		// Look for elements that have an identical side key
		std::pair <map_type::iterator, map_type::iterator>
		  bounds = side_to_elem_map.equal_range(key);
		
		// May be multiple keys, check all the possible
		// elements which _might_ be neighbors.
		if (bounds.first != bounds.second)
		  {
		    // Get the side for this element
		    const AutoPtr<DofObject> my_side(element->side(ms));

		    // Look at all the entries with an equivalent key
		    while (bounds.first != bounds.second)
		      {
			// Get the potential element
			Elem* neighbor = bounds.first->second.first;
			
			// Get the side for the neighboring element
			const unsigned int ns = bounds.first->second.second;
			const AutoPtr<DofObject> their_side(neighbor->side(ns));
                        //assert (my_side.get() != NULL);
                        //assert (their_side.get() != NULL);			

			// If found a match wbounds.firsth my side
                        //
                        // We need a special case here for 1D, since parents
                        // and children have an equal side (i.e. a node), 
                        // so need to check ns != ms in 1D as well
			if( (*my_side == *their_side) && 
                            ((_dim != 1) || (ns != ms)) )
			  {
			    // So share a side.  Is this a mixed pair
			    // of subactive and active/ancestor
			    // elements?
                            // If not, then we're neighbors.
			    // If so, then the subactive's neighbor is 

                              if (element->subactive() ==
                                  neighbor->subactive())
                              {
                              // an element is only subactive if it has
                              // been coarsened but not deleted
                                element->set_neighbor (ms,neighbor);
                                neighbor->set_neighbor(ns,element);
                              }
                              else if (element->subactive())
                              {
                                element->set_neighbor(ms,neighbor);
                              }
                              else if (neighbor->subactive())
                              {
                                neighbor->set_neighbor(ns,element);
                              }
                              side_to_elem_map.erase (bounds.first);

                              // get out of this nested crap
                              goto next_side; 
			  }

			++bounds.first;
		      }
		  }
		    
		// didn't find a match...
		// Build the map entry for this element
		key_val_pair kvp;
		
		kvp.first         = key;
		kvp.second.first  = element;
		kvp.second.second = ms;
		
		// use the lower bound as a hint for
		// where to put it.
#if defined(HAVE_HASH_MAP) || defined(HAVE_EXT_HASH_MAP)
		side_to_elem_map.insert (kvp);
#else
		side_to_elem_map.insert (bounds.first,kvp);
#endif
	      }
	  }
      }
  }

  
  
#ifdef ENABLE_AMR

  /**
   * Here we look at all of the child elements.
   * If a child element has a NULL neighbor it is 
   * either because it is on the boundary or because
   * its neighbor is at a different level.  In the
   * latter case we must get the neighbor from the
   * parent.
   *
   * Furthermore, that neighbor better be active,
   * otherwise we missed a child somewhere.
   */
  element_iterator el  = this->not_level_elements_begin(0);
  element_iterator end = this->not_level_elements_end(0);

  for (; el != end; ++el)
    {
      Elem* elem = *el;

      assert (elem->parent() != NULL);
      for (unsigned int s=0; s < elem->n_neighbors(); s++)
        if (elem->neighbor(s) == NULL)
        {	    
          elem->set_neighbor(s, elem->parent()->neighbor(s));

#ifdef DEBUG	    
          Elem *neigh = elem->neighbor(s);
          if (neigh != NULL)
            // We ignore subactive elements here because
            // we don't care about neighbors of subactive element.
            if ((!neigh->active()) && (!elem->subactive()))
            {
              std::cerr << "Bad element ID = " << elem->id() 
                << ", Bad neighbor ID = " << neigh->id();
              std::cerr << "ERROR: " 
                << (elem->active()?"Active":"Ancestor")
                << " Element at level "
                << elem->level() << " found "
                << (neigh->subactive()?"subactive":"ancestor")
                << " neighbor at level " << neigh->level()
                << std::endl;
              GMVIO(*dynamic_cast<Mesh*>(this)).write ("bad_mesh.gmv");
              error();
            }
#endif // DEBUG
        }
    }
  
#endif // AMR

  STOP_LOG("find_neighbors()", "Mesh");
}






void Mesh::renumber_nodes_and_elements ()
{
  // Only renumber the nodes & trim the element vector
  // if we have performed some local coarsening.  This
  // will be the case if the element vector is not the
  // same size as the number of active elements.
  //
  // In this case simply number the elements (in case this
  // was not already done) and return.
  if (_elements.size() == this->n_active_elem())
    { 
      START_LOG("renumber_nodes_and_elem()", "Mesh");
      
      // Number the elements
      {
	element_iterator       it  = this->elements_begin();
	const element_iterator end = this->elements_end();

	for (unsigned int id=0; it != end; ++it)
	  (*it)->set_id() = id++;
      }

      // Number the nodes
      {
	node_iterator       it  = this->nodes_begin();
	const node_iterator end = this->nodes_end();

	for (unsigned int id=0; it != end; ++it)
	  (*it)->set_id() = id++;
      }      
      
      STOP_LOG("renumber_nodes_and_elem()", "Mesh");
      return;
    }
  
  START_LOG("renumber_nodes_and_elem()", "Mesh");
  
  // node and element id counters
  unsigned int next_free_elem = 0;
  unsigned int next_free_node = 0;

  // Loop over the elements.  Note that there may
  // be NULLs in the _elements vector from the coarsening
  // process.  Pack the elements in to a contiguous array
  // and then trim any excess.
  {      
    std::vector<Elem*>::iterator in        = _elements.begin();
    std::vector<Elem*>::iterator out       = _elements.begin();
    const std::vector<Elem*>::iterator end = _elements.end();

    for (; in != end; ++in)
      if (*in != NULL)
	{
	  Elem* elem = *in;
	  
	  *out = *in;
	  ++out;
	  
	  // Increment the element counter
	  elem->set_id (next_free_elem++);
	  
	  // Loop over this element's nodes.  Number them,
	  // if they have not been numbered already.  Also,
	  // position them in the _nodes vector so that they
	  // are packed contiguously from the beginning.
	  for (unsigned int n=0; n<elem->n_nodes(); n++)
	    if (elem->node(n) == next_free_node)     // don't need to process
	      next_free_node++;                      // [(src == dst) below]

	    else if (elem->node(n) > next_free_node) // need to process
	      {
		// The source and destination indices
		// for this node
		const unsigned int src_idx = elem->node(n);
		const unsigned int dst_idx = next_free_node++;

		// ensure we want to swap valid nodes
		assert (_nodes[src_idx] != NULL);
		assert (_nodes[dst_idx] != NULL);
		
		// Swap the source and destination nodes
		std::swap (_nodes[src_idx],
			   _nodes[dst_idx] );

		// Set proper indices
		_nodes[src_idx]->set_id (src_idx);
		_nodes[dst_idx]->set_id (dst_idx);
	      }
	}

    // Erase any additional storage. These elements have been
    // copied into NULL voids by the procedure above, and are
    // thus repeated and unnecessary.
    _elements.erase (out, end);
  }

  // Any nodes in the vector >= _nodes[next_free_node]
  // are not connected to any elements and may be deleted
  // if desired.

  // (This code block will erase the unused nodes)
  // Now, delete the unused nodes
  {
    std::vector<Node*>::iterator nd        = _nodes.begin();
    const std::vector<Node*>::iterator end = _nodes.end();

    std::advance (nd, next_free_node);
    
    for (std::vector<Node*>::iterator it=nd;
	 it != end; ++it)
      {
	assert (*it != NULL);

	// remove any boundary information associated with
	// this node
	this->boundary_info->remove (*it);
	
	// delete the node
	delete *it;
	*it = NULL;
      }
    
    _nodes.erase (nd, end);
  }
  

  assert (next_free_elem == _elements.size());
  assert (next_free_node == _nodes.size());
  
  STOP_LOG("renumber_nodes_and_elem()", "Mesh");
}








void Mesh::read (const std::string& name,
		 MeshData* mesh_data)
{
  START_LOG("read()", "Mesh");
  
  // Set the read_xda_file flag on all processors.
  // This ensures that renumber_nodes_and_elements is *not* called
  // during prepare_for_use().  This is required in cases 
  // where there is a associated solution file which expect
  // a certain ordering of the nodes.
  const bool read_xda_file =
    name.rfind(".xda") < name.size();
  
  // Read the file based on extension.  Only processor 0
  // needs to read the mesh.  It will then broadcast it and
  // the other processors will pick it up
  if (libMesh::processor_id() == 0)
    {
      if (name.rfind(".mat") < name.size())
	MatlabIO(*this).read(name);
      
      else if (name.rfind(".ucd") < name.size())
	UCDIO(*this).read (name);
      
      else if (name.rfind(".exd") < name.size())
	ExodusII_IO(*this).read (name);
      
      else if ((name.rfind(".off")  < name.size()) ||
	       (name.rfind(".ogl")  < name.size()) ||
	       (name.rfind(".oogl") < name.size()))
	OFFIO(*this).read (name);
     
      else if (name.rfind(".xda") < name.size())
	XdrIO(*this).read (name);
      
      else if (name.rfind(".xdr")  < name.size())
	XdrIO(*this,true).read (name);
      
      else if ((name.rfind(".mgf")  < name.size()) ||
	       (name.rfind(".0000") < name.size()))
	XdrIO(*this,true).read_mgf (name);
      
      else if (name.rfind(".unv") < name.size())
	{
	  if (mesh_data == NULL)
	    {
	      std::cerr << "Error! You must pass a "
			<< "valid MeshData pointer to "
			<< "read UNV files!" << std::endl;
	      error();
	    }
	  UNVIO(*this, *mesh_data).read (name);
	}
      
      else if ((name.rfind(".node")  < name.size()) ||
	       (name.rfind(".ele")   < name.size()))
	TetGenIO(*this,mesh_data).read (name);

      else if (name.rfind(".msh") < name.size())
	GmshIO(*this).read (name);
      
      else
	{
	  std::cerr << " ERROR: Unrecognized file extension: " << name
		    << "\n   I understand the following:\n\n"
		    << "     *.mat  -- Matlab triangular ASCII file\n"
		    << "     *.ucd  -- AVS's ASCII UCD format\n"
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
  
  STOP_LOG("read()", "Mesh");

  // Send the mesh & bcs (which are now only on processor 0) to the other
  // processors
  {
    MeshCommunication mesh_communication;
  
    mesh_communication.broadcast (*this);
  }

  // Done reading the mesh.  Now prepare it for use.
  this->prepare_for_use(read_xda_file);

}



void Mesh::write (const std::string& name,
		  MeshData* mesh_data)
{
  START_LOG("write()", "Mesh");
  
  // Write the file based on extension
  if (name.rfind(".dat") < name.size())
    TecplotIO(*this).write (name);
    
  else if (name.rfind(".plt") < name.size())
    TecplotIO(*this,true).write (name);
    
  else if (name.rfind(".ucd") < name.size())
    UCDIO (*this).write (name);
    
  else if (name.rfind(".gmv") < name.size())
    if (this->n_partitions() > 1)
      GMVIO(*this).write (name);
    else
      {
	GMVIO io(*this);
	io.partitioning() = false;
	io.write (name);
      }
    
  else if (name.rfind(".ugrid") < name.size())
    DivaIO(*this).write(name);
    
  else if (name.rfind(".xda") < name.size())
    XdrIO(*this).write(name);
    
  else if (name.rfind(".xdr") < name.size())
    XdrIO(*this,true).write(name);
    
  else if (name.rfind(".mgf")  < name.size())
    XdrIO(*this,true).write_mgf(name);
    
  else if (name.rfind(".unv") < name.size())
    {
      if (mesh_data == NULL)
	{
	  std::cerr << "Error! You must pass a "
		    << "valid MeshData pointer to "
		    << "write UNV files!" << std::endl;
	  error();
	}
      UNVIO(*this, *mesh_data).write (name);
    }

  else if (name.rfind(".mesh") < name.size())
    MEDITIO(*this).write (name);

  else if (name.rfind(".poly") < name.size())
    TetGenIO(*this).write (name);

  else if (name.rfind(".msh") < name.size())
    GmshIO(*this).write (name);
    
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
		<< "     *.mesh  -- MEdit mesh format\n"
		<< "     *.poly  -- TetGen ASCII file\n"
		<< std::endl
		<< "\n Exiting without writing output\n";
    }    
  
  STOP_LOG("write()", "Mesh");
}



void Mesh::write (const std::string& name,
		  const std::vector<Number>& v,
		  const std::vector<std::string>& vn)
{
  START_LOG("write()", "Mesh");

  // Write the file based on extension
  if (name.rfind(".dat") < name.size())
    TecplotIO(*this).write_nodal_data (name, v, vn);
  
  if (name.rfind(".plt") < name.size())
    TecplotIO(*this,true).write_nodal_data (name, v, vn);
  
  else if (name.rfind(".gmv") < name.size())
    {
      if (n_subdomains() > 1)
	GMVIO(*this).write_nodal_data (name, v, vn);
      else
	{
	  GMVIO io(*this);
	  io.partitioning() = false;
	  io.write_nodal_data (name, v, vn);
	}
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
//   const_active_pid_elem_iterator       it(this->elements_begin(),   pid);
//   const const_active_pid_elem_iterator it_end(this->elements_end(), pid);

  const_element_iterator       it     = this->active_pid_elements_begin(pid);
  const const_element_iterator it_end = this->active_pid_elements_end(pid);
    
  this->create_submesh (pid_mesh, it, it_end);
}







void Mesh::create_submesh (Mesh& new_mesh,
			   const_element_iterator& it,
			   const const_element_iterator& it_end) const
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
      Elem* new_elem = 
	new_mesh.add_elem (Elem::build(old_elem->type()).release());

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
	  if (this->boundary_info->boundary_id (old_elem, s) !=
	      this->boundary_info->invalid_id)
	    new_mesh.boundary_info->add_side (new_elem,
					     s,
					     this->boundary_info->boundary_id (old_elem, s));
    } // end loop over elements
  

  // Prepare the new_mesh for use
  new_mesh.prepare_for_use();
  
}




bool Mesh::contract ()
{
  START_LOG ("contract()", "Mesh");

  // Flag indicating if this call actually changes the mesh
  bool mesh_changed = false;

  std::vector<Elem*>::iterator in        = _elements.begin();
  std::vector<Elem*>::iterator out       = _elements.begin();
  const std::vector<Elem*>::iterator end = _elements.end();

#ifdef DEBUG
  for ( ; in != end; ++in)
    if (*in != NULL)
      {
	Elem* elem = *in;
	assert(elem->active() || elem->subactive() || elem->ancestor());
      }
  in = _elements.begin();
#endif
  
  unsigned int next_free_elem = 0;
  
  // Loop over the elements.   
  for ( ; in != end; ++in)
    if (*in != NULL)
      {
	Elem* elem = *in;

	// Delete all the subactive ones
	if (elem->subactive())
	  {
	    // Huh?  no level-0 element should be subactive
	    assert (elem->level() != 0);

	    // Make sure we dealt with parents first
	    if (elem->parent()->has_children())
	      {
		std::cerr << "Element being deleted is still a child." << std::endl;
	      }

	    // Delete the element
	    delete elem;
	    *in = NULL;
	    
	    // the mesh has certainly changed
	    mesh_changed = true;
	  }
	else
	  {
	    // Compress all the active ones
	    if (elem->active())
	      elem->contract();
	    else
	      assert (elem->ancestor());

	    // These elements are kept, pack them
	    // contiguously in the _elements vector
	    elem->set_id(next_free_elem++);
	    *out = elem;
	    ++out;
	  }
      }

  // Erase any additional storage. These elements have been
  // copied into NULL voids by the procedure above, and are
  // thus repeated and unnecessary.
  _elements.erase (out, end);

  STOP_LOG ("contract()", "Mesh");
  
  return mesh_changed;
}




