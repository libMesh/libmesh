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



// C++ includes
#include <fstream>

// C includes
#include <unistd.h>  // for unlink()

// Local includes
#include "boundary_info.h"
#include "unstructured_mesh.h"
#include "mesh_communication.h"
#include "libmesh_logging.h"
#include "elem.h"
#include "mesh_tools.h" // For n_levels
#include "parallel.h"
#include "remote_elem.h"

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
#include "fro_io.h"
#include "xdr_io.h"
#include "legacy_xdr_io.h"
#include "vtk_io.h"

#if   defined(HAVE_TR1_UNORDERED_MAP)
# include <tr1/unordered_map>
#elif defined(HAVE_UNORDERED_MAP)
# include <unordered_map>
#elif defined(HAVE_HASH_MAP)
# include <hash_map>
#elif defined(HAVE_EXT_HASH_MAP)
# include <ext/hash_map>
#else
# include <map>
#endif



// ------------------------------------------------------------
// Anonymous namespace for implementation details
namespace {
  bool is_parallel_file_format (const std::string &name)
  {
    // Certain mesh formats support parallel I/O, including the
    // "new" Xdr format.
    return ((name.rfind(".xda") < name.size()) ||
	    (name.rfind(".xdr") < name.size())
	    );
  }
}



// ------------------------------------------------------------
// UnstructuredMesh class member functions
UnstructuredMesh::UnstructuredMesh (unsigned int d) :
  MeshBase (d)
{
  libmesh_assert (libMesh::initialized());
}



void UnstructuredMesh::copy_nodes_and_elements 
  (const UnstructuredMesh& other_mesh)
{
  // We're assuming our subclass data needs no copy
  libmesh_assert(_n_sbd == other_mesh._n_sbd);
  libmesh_assert(_n_parts == other_mesh._n_parts);
  libmesh_assert(_dim == other_mesh._dim);
  libmesh_assert(_is_prepared == other_mesh._is_prepared);

  //Copy in Nodes
  {
    //Preallocate Memory if necessary
    this->reserve_nodes(other_mesh.n_nodes());
    
    const_node_iterator it = other_mesh.nodes_begin();
    const_node_iterator end = other_mesh.nodes_end();

    for (; it != end; ++it)
      this->add_point(*(*it)); //Add new nodes in old node Point locations
  }
  
  //Copy in Elements
  {
    //Preallocate Memory if necessary
    this->reserve_elem(other_mesh.n_elem());
    
    // Loop over the elements
    MeshBase::const_element_iterator it = other_mesh.elements_begin();
    const MeshBase::const_element_iterator end = other_mesh.elements_end();

    // FIXME: Where do we set element IDs??
    for (; it != end; ++it)
    {
      //Look at the old element
      Elem *old = *it;
      //Build a new element
      Elem *newparent = old->parent() ?
          this->elem(old->parent()->id()) : NULL;
      AutoPtr<Elem> ap = Elem::build(old->type(), newparent);
      Elem * elem = ap.release();

#ifdef ENABLE_AMR
      //Create the parent's child pointers if necessary
      if (newparent)
        {
          // Make sure we have space for those child pointers
          newparent->add_child(elem);

          // We'd better be adding these in the correct order
          libmesh_assert (newparent->which_child_am_i(elem) ==
                  old->parent()->which_child_am_i(old));
        }
      
      // Copy the refinement flags
      elem->set_refinement_flag(old->refinement_flag());
      elem->set_p_refinement_flag(old->p_refinement_flag());
#endif // #ifdef ENABLE_AMR

      //Assign all the nodes
      for(unsigned int i=0;i<elem->n_nodes();i++)
        elem->set_node(i) = &this->node(old->node(i));
      
      //Hold onto it
      this->add_elem(elem);
    }
  }
  
  //Finally prepare the Mesh for use
  this->prepare_for_use();
}
 
 

UnstructuredMesh::~UnstructuredMesh ()
{
//  this->clear ();  // Nothing to clear at this level
  
  libmesh_assert (!libMesh::closed());
}



void UnstructuredMesh::find_neighbors(bool reset_remote_elements)
{
  libmesh_assert(this->n_nodes() != 0);
  libmesh_assert(this->n_elem()  != 0);

  // This function must be run on all processors at once
  parallel_only();

  START_LOG("find_neighbors()", "Mesh");
  
  
  //TODO:[BSK] This should be removed later?!
  const element_iterator el_end = this->elements_end();
  for (element_iterator el = this->elements_begin(); el != el_end; ++el)
    {
      Elem* elem = *el;
      for (unsigned int s=0; s<elem->n_neighbors(); s++)
        if (elem->neighbor(s) != remote_elem ||
            reset_remote_elements)
          elem->set_neighbor(s,NULL);
    }
  
  // Find neighboring elements by first finding elements
  // with identical side keys and then check to see if they
  // are neighbors
  {
    // data structures -- Use the hash_multimap if available
    typedef unsigned int                    key_type;
    typedef std::pair<Elem*, unsigned char> val_type;
    typedef std::pair<key_type, val_type>   key_val_pair;
    
#if   defined(HAVE_UNORDERED_MAP)
    typedef std::unordered_multimap<key_type, val_type> map_type;    
#elif defined(HAVE_TR1_UNORDERED_MAP)
    typedef std::tr1::unordered_multimap<key_type, val_type> map_type;    
#elif defined(HAVE_HASH_MAP)    
    typedef std::hash_multimap<key_type, val_type> map_type;    
#elif defined(HAVE_EXT_HASH_MAP)
# if    (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0   
    typedef std::hash_multimap<key_type, val_type> map_type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
    typedef __gnu_cxx::hash_multimap<key_type, val_type> map_type;
# else
// XLC and who knows what other compilers get here.
// Try the most standard thing we can:
    typedef std::multimap<key_type, val_type>  map_type;
# endif
#else
    typedef std::multimap<key_type, val_type>  map_type;
#endif
    
    // A map from side keys to corresponding elements & side numbers  
    map_type side_to_elem_map;
  


    for (element_iterator el = this->elements_begin(); el != el_end; ++el)
      {
	Elem* element = *el;

	for (unsigned int ms=0; ms<element->n_neighbors(); ms++)
	  {
	  next_side:
	    // If we haven't yet found a neighbor on this side, try.
	    // Even if we think our neighbor is remote, that
	    // information may be out of date.
	    if (element->neighbor(ms) == NULL ||
		element->neighbor(ms) == remote_elem)
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
                        //libmesh_assert (my_side.get() != NULL);
                        //libmesh_assert (their_side.get() != NULL);			

			// If found a match with my side
                        //
			// We need special tests here for 1D:
			// since parents and children have an equal
			// side (i.e. a node), we need to check 
			// ns != ms, and we also check level() to
			// avoid setting our neighbor pointer to
			// any of our neighbor's descendants
			if( (*my_side == *their_side) && 
                            (element->level() == neighbor->level()) &&
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
#if defined(HAVE_UNORDERED_MAP) || defined(HAVE_TR1_UNORDERED_MAP) || defined(HAVE_HASH_MAP) || defined(HAVE_EXT_HASH_MAP)
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
   *
   * If a child element has a NULL neighbor it is 
   * either because it is on the boundary or because
   * its neighbor is at a different level.  In the
   * latter case we must get the neighbor from the
   * parent.
   *
   * If a child element has a remote_elem neighbor
   * on a boundary it shares with its parent, that
   * info may be out of date - if the parent's
   * neighbor is active then the child should share
   * it.
   *
   * Furthermore, that neighbor better be active,
   * otherwise we missed a child somewhere.
   */
  const unsigned int n_levels = MeshTools::n_levels(*this);
  for (unsigned int level = 1; level < n_levels; ++level)
    {
      element_iterator end = this->level_elements_end(level);
      for (element_iterator el = this->level_elements_begin(level);
           el != end; ++el)
        {
          Elem* elem = *el;
          libmesh_assert(elem);
          libmesh_assert(elem->parent());

          for (unsigned int s=0; s < elem->n_neighbors(); s++)
            if (elem->neighbor(s) == NULL)
// This currently leads to an infinite loop in ex10?
//            if (elem->neighbor(s) == NULL ||
//		(elem->neighbor(s) == remote_elem &&
//		 parent->is_child_on_side(parent->which_child_am_i(elem), s)))
            {	    
              Elem *neigh = elem->parent()->neighbor(s);

	      // If neigh was refined and had non-subactive children
	      // made remote earlier, then a non-subactive elem should
	      // actually have one of those remote children as a
	      // neighbor
              if (neigh && (neigh->ancestor()) && (!elem->subactive()))
                {
#ifdef DEBUG	    
                  // Let's make sure that "had children made remote"
	          // situation is actually the case
		  libmesh_assert(neigh->has_children());
		  bool neigh_has_remote_children = false;
		  for (unsigned int c = 0; c != neigh->n_children(); ++c)
                    {
                      if (neigh->child(c) == remote_elem)
                        neigh_has_remote_children = true;
                    }
                  libmesh_assert(neigh_has_remote_children);

	          // And let's double-check that we don't have
		  // a remote_elem neighboring a local element
                  libmesh_assert(elem->processor_id() !=
				 libMesh::processor_id());
#endif // DEBUG
                  neigh = const_cast<RemoteElem*>(remote_elem);
                }

              elem->set_neighbor(s, neigh);
#ifdef DEBUG	    
              if (neigh != NULL && neigh != remote_elem)
                // We ignore subactive elements here because
                // we don't care about neighbors of subactive element.
                if ((!neigh->active()) && (!elem->subactive()))
                {
                  std::cerr << "On processor " << libMesh::processor_id() 
                            << std::endl;
                  std::cerr << "Bad element ID = " << elem->id() 
                    << ", Side " << s << ", Bad neighbor ID = " << neigh->id() << std::endl;
                  std::cerr << "Bad element proc_ID = " << elem->processor_id() 
                    << ", Bad neighbor proc_ID = " << neigh->processor_id() << std::endl;
                  std::cerr << "Bad element size = " << elem->hmin() 
                    << ", Bad neighbor size = " << neigh->hmin() << std::endl;
                  std::cerr << "Bad element center = " << elem->centroid() 
                    << ", Bad neighbor center = " << neigh->centroid() << std::endl;
                  std::cerr << "ERROR: " 
                    << (elem->active()?"Active":"Ancestor")
                    << " Element at level "
                    << elem->level() << std::endl;
                  std::cerr << "with "
                    << (elem->parent()->active()?"active":
                        (elem->parent()->subactive()?"subactive":"ancestor"))
                    << " parent share "
                    << (neigh->subactive()?"subactive":"ancestor")
                    << " neighbor at level " << neigh->level()
                    << std::endl;
                  GMVIO(*dynamic_cast<UnstructuredMesh*>(this)).write ("bad_mesh.gmv");
                  libmesh_error();
                }
#endif // DEBUG
            }
        }
    }
  
#endif // AMR

#ifdef DEBUG
MeshTools::libmesh_assert_valid_neighbors(*this);
#endif

  STOP_LOG("find_neighbors()", "Mesh");
}



void UnstructuredMesh::read (const std::string& name,
			     MeshData* mesh_data)
{
  // See if the file exists.  Perform this check on all processors
  // so that the code is terminated properly in the case that the
  // file does not exist.
  {
    std::ifstream in (name.c_str());
    
    if (!in.good())
      {
	std::cerr << "ERROR: cannot locate specified file:\n\t"
		  << name
		  << std::endl;
	libmesh_error();
      }
  }

  // Set the skip_renumber_nodes_and_elements flag on all processors.
  // This ensures that renumber_nodes_and_elements is *not* called
  // during prepare_for_use() for certain types of mesh files.
  // This is required in cases where there is an associated solution
  // file which expects a certain ordering of the nodes.
  bool skip_renumber_nodes_and_elements = (name.rfind(".gmv") < name.size());
  
  // Look for parallel formats first
  if (is_parallel_file_format(name))
    {      
      // no need to handling bz2 files here -- the Xdr class does that.
      if ((name.rfind(".xda") < name.size()) ||
	  (name.rfind(".xdr") < name.size()))
	{
	  XdrIO xdr_io(*this);

	  // .xda* ==> bzip2/gzip/ASCII flavors
	  if (name.rfind(".xda") < name.size())
	    {
	      xdr_io.binary() = false;
	      xdr_io.read (name);
	    }
	  else // .xdr* ==> true binary XDR file
	    {
	      xdr_io.binary() = true;
	      xdr_io.read (name);
	    }

	  // The xdr_io object gets constructed with legacy() == false.
	  // if legacy() == true then it means that a legacy file was detected and
	  // thus processor 0 performed the read. We therefore need to broadcast the
	  // mesh.  Further, for this flavor of mesh solution data ordering is tied
	  // to the node ordering, so we better not reorder the nodes!
	  if (xdr_io.legacy())
	    {
	      skip_renumber_nodes_and_elements = true;
	      MeshCommunication().broadcast(*this);
	    }
	}      
    }

  // Serial mesh formats
  else
    {
      START_LOG("read()", "Mesh");  
  
      // Read the file based on extension.  Only processor 0
      // needs to read the mesh.  It will then broadcast it and
      // the other processors will pick it up
      if (libMesh::processor_id() == 0)
	{
	  // Nasty hack for reading/writing zipped files
	  std::string new_name = name;
	  if (name.size() - name.rfind(".bz2") == 4)
	    {
	      new_name.erase(new_name.end() - 4, new_name.end());
	      std::string system_string = "bunzip2 -f -k ";
	      system_string += name;
	      START_LOG("system(bunzip2)", "Mesh");
	      std::system(system_string.c_str());
	      STOP_LOG("system(bunzip2)", "Mesh");
	    }

	  if (new_name.rfind(".mat") < new_name.size())
	    MatlabIO(*this).read(new_name);
	  
	  else if (new_name.rfind(".ucd") < new_name.size())
	    UCDIO(*this).read (new_name);
	  
	  else if (new_name.rfind(".exd") < new_name.size() ||
		   new_name.rfind(".e") < new_name.size())
	    ExodusII_IO(*this).read (new_name);
	  
	  else if ((new_name.rfind(".off")  < new_name.size()) ||
		   (new_name.rfind(".ogl")  < new_name.size()) ||
		   (new_name.rfind(".oogl") < new_name.size()))
	    OFFIO(*this).read (new_name);
     
	  else if (new_name.rfind(".mgf") < new_name.size())
	    LegacyXdrIO(*this,true).read_mgf (new_name);
      
	  else if (new_name.rfind(".unv") < new_name.size())
	    {
	      if (mesh_data == NULL)
		{
		  std::cerr << "Error! You must pass a "
			    << "valid MeshData pointer to "
			    << "read UNV files!" << std::endl;
		  libmesh_error();
		}
	      UNVIO(*this, *mesh_data).read (new_name);
	    }
      
	  else if ((new_name.rfind(".node")  < new_name.size()) ||
		   (new_name.rfind(".ele")   < new_name.size()))
	    TetGenIO(*this,mesh_data).read (new_name);

	  else if (new_name.rfind(".msh") < new_name.size())
	    GmshIO(*this).read (new_name);

	  else if (new_name.rfind(".gmv") < new_name.size())
	    GMVIO(*this).read (new_name);

	  else if (new_name.rfind(".vtu") < new_name.size())
	    VTKIO(*this).read(new_name);
      
	  else
	    {
	      std::cerr << " ERROR: Unrecognized file extension: " << name
			<< "\n   I understand the following:\n\n"
			<< "     *.mat  -- Matlab triangular ASCII file\n"
			<< "     *.ucd  -- AVS's ASCII UCD format\n"
			<< "     *.gmv  -- LANL's General Mesh Viewer format\n"
			<< "     *.off  -- OOGL OFF surface format\n"
			<< "     *.exd  -- Sandia's ExodusII format\n"
			<< "     *.e    -- Sandia's ExodusII format\n"
			<< "     *.xda  -- Internal ASCII format\n"
			<< "     *.xdr  -- Internal binary format,\n"
			<< "               compatible with XdrMGF\n"
			<< "     *.unv  -- I-deas Universal format\n"
			<< std::endl;
	      libmesh_error();	  
	    }    
	  
	  // If we temporarily decompressed a .bz2 file, remove the
	  // uncompressed version
	  if (name.size() - name.rfind(".bz2") == 4)
	    std::remove(new_name.c_str());
	}
      
  
      STOP_LOG("read()", "Mesh");

      // Send the mesh & bcs (which are now only on processor 0) to the other
      // processors
      MeshCommunication().broadcast (*this);
    }

  
  // Done reading the mesh.  Now prepare it for use.
  this->prepare_for_use(skip_renumber_nodes_and_elements);

}



void UnstructuredMesh::write (const std::string& name,
			      MeshData* mesh_data)
{
  // parallel formats are special -- they may choose to write
  // separate files, let's not try to handle the zipping here.
  if (is_parallel_file_format(name))
    {	
      // no need to handling bz2 files here -- the Xdr class does that.
      if (name.rfind(".xda") < name.size())
	XdrIO(*this).write(name);
	
      else if (name.rfind(".xdr") < name.size())
	XdrIO(*this,true).write(name);
    }

  // serial file formats
  else
    {
      START_LOG("write()", "Mesh");

      // Nasty hack for reading/writing zipped files
      std::string new_name = name;
      if (name.size() - name.rfind(".bz2") == 4)
	new_name.erase(new_name.end() - 4, new_name.end());
  
      // New scope so that io will close before we try to zip the file
      {
	// Write the file based on extension
	if (new_name.rfind(".dat") < new_name.size())
	  TecplotIO(*this).write (new_name);
	
	else if (new_name.rfind(".plt") < new_name.size())
	  TecplotIO(*this,true).write (new_name);
	
	else if (new_name.rfind(".ucd") < new_name.size())
	  UCDIO (*this).write (new_name);
	
	else if (new_name.rfind(".gmv") < new_name.size())
	  if (this->n_partitions() > 1)
	    GMVIO(*this).write (new_name);
	  else
	    {
	      GMVIO io(*this);
	      io.partitioning() = false;
	      io.write (new_name);
	    }
	
	else if (new_name.rfind(".ugrid") < new_name.size())
	  DivaIO(*this).write(new_name);
    
	else if (new_name.rfind(".mgf")  < new_name.size())
	  LegacyXdrIO(*this,true).write_mgf(new_name);
	
	else if (new_name.rfind(".unv") < new_name.size())
	  {
	    if (mesh_data == NULL)
	      {
		std::cerr << "Error! You must pass a "
			  << "valid MeshData pointer to "
			  << "write UNV files!" << std::endl;
		libmesh_error();
	      }
	    UNVIO(*this, *mesh_data).write (new_name);
	  }
	
	else if (new_name.rfind(".mesh") < new_name.size())
	  MEDITIO(*this).write (new_name);
	
	else if (new_name.rfind(".poly") < new_name.size())
	  TetGenIO(*this).write (new_name);
	
	else if (new_name.rfind(".msh") < new_name.size())
	  GmshIO(*this).write (new_name);
	
	else if (new_name.rfind(".fro") < new_name.size())
	  FroIO(*this).write (new_name);
	
	else if (new_name.rfind(".vtu") < new_name.size())
	  VTKIO(*this).write (new_name);
	
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
		      << "     *.msh   -- GMSH ASCII file\n"
		      << "     *.fro   -- ACDL's surface triangulation file\n"
		      << std::endl
		      << "\n Exiting without writing output\n";
	  }    
      }
  
      // Nasty hack for reading/writing zipped files
      if (name.size() - name.rfind(".bz2") == 4)
	{
	  START_LOG("system(bzip2)", "Mesh");
	  if (libMesh::processor_id() == 0)
	    {
	      std::string system_string = "bzip2 -f ";
	      system_string += new_name;
	      std::system(system_string.c_str());
	    }
	  Parallel::barrier();
	  STOP_LOG("system(bzip2)", "Mesh");
	}
      
      STOP_LOG("write()", "Mesh");
    }  
}



void UnstructuredMesh::write (const std::string& name,
		  const std::vector<Number>& v,
		  const std::vector<std::string>& vn)
{
  START_LOG("write()", "Mesh");

  // Write the file based on extension
  if (name.rfind(".dat") < name.size())
    TecplotIO(*this).write_nodal_data (name, v, vn);
  
  else if (name.rfind(".plt") < name.size())
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
  else if (name.rfind(".pvtu") < name.size())
    {
      VTKIO(*this).write_nodal_data (name, v, vn);
    }
  else
    {
      std::cerr << " ERROR: Unrecognized file extension: " << name
		<< "\n   I understand the following:\n\n"
		<< "     *.dat  -- Tecplot ASCII file\n"
		<< "     *.plt  -- Tecplot binary file\n"
		<< "     *.pvtu -- Paraview VTK file\n"
		<< "     *.gmv  -- LANL's GMV (General Mesh Viewer) format\n"
		<< "\n Exiting without writing output\n";
    }

  STOP_LOG("write()", "Mesh");
}





void UnstructuredMesh::create_pid_mesh(UnstructuredMesh& pid_mesh,
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







void UnstructuredMesh::create_submesh (UnstructuredMesh& new_mesh,
			   const_element_iterator& it,
			   const const_element_iterator& it_end) const
{
  // Just in case the subdomain_mesh already has some information
  // in it, get rid of it.
  new_mesh.clear();

  // Fail if (*this == new_mesh), we cannot create a submesh inside ourself!
  // This may happen if the user accidently passes the original mesh into
  // this function!  We will check this by making sure we did not just
  // clear ourself.
  libmesh_assert (this->n_nodes() != 0);
  libmesh_assert (this->n_elem()  != 0); 

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

      libmesh_assert (new_elem != NULL);
	
      // Loop over the nodes on this element.  
      for (unsigned int n=0; n<old_elem->n_nodes(); n++)
	{
	  libmesh_assert (old_elem->node(n) < new_node_numbers.size());

	  if (new_node_numbers[old_elem->node(n)] == libMesh::invalid_uint)
	    {
	      new_node_numbers[old_elem->node(n)] = n_new_nodes;

	      // Add this node to the new mesh
	      new_mesh.add_point (old_elem->point(n));

	      // Increment the new node counter
	      n_new_nodes++;
	    }

	  // Define this element's connectivity on the new mesh
	  libmesh_assert (new_node_numbers[old_elem->node(n)] < new_mesh.n_nodes());
	    
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



#ifdef ENABLE_AMR
bool UnstructuredMesh::contract ()
{
  START_LOG ("contract()", "Mesh");

  // Flag indicating if this call actually changes the mesh
  bool mesh_changed = false;

  element_iterator in        = elements_begin();
  element_iterator out       = elements_begin();
  const element_iterator end = elements_end();

#ifdef DEBUG
  for ( ; in != end; ++in)
    if (*in != NULL)
      {
	Elem* elem = *in;
	libmesh_assert(elem->active() || elem->subactive() || elem->ancestor());
      }
  in = elements_begin();
#endif
  
  // Loop over the elements.   
  for ( ; in != end; ++in)
    if (*in != NULL)
      {
	Elem* elem = *in;

	// Delete all the subactive ones
	if (elem->subactive())
	  {
	    // Huh?  no level-0 element should be subactive
	    libmesh_assert (elem->level() != 0);

	    // Delete the element
	    // This just sets a pointer to NULL, and doesn't
	    // invalidate any iterators
	    this->delete_elem(elem);
	    
	    // the mesh has certainly changed
	    mesh_changed = true;
	  }
	else
	  {
	    // Compress all the active ones
	    if (elem->active())
	      elem->contract();
	    else
	      libmesh_assert (elem->ancestor());
	  }
      }

  // Strip any newly-created NULL voids out of the element array
  this->renumber_nodes_and_elements();

  STOP_LOG ("contract()", "Mesh");
  
  return mesh_changed;
}
#endif // #ifdef ENABLE_AMR



