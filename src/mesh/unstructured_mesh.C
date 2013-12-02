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
#include <sstream>
#include <iomanip>

// C includes
#include <unistd.h>  // for unlink()

// Local includes
#include "libmesh/boundary_info.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_tools.h" // For n_levels
#include "libmesh/parallel.h"
#include "libmesh/remote_elem.h"

#include "libmesh/diva_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmv_io.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/tetgen_io.h"
#include "libmesh/ucd_io.h"
#include "libmesh/unv_io.h"
#include "libmesh/matlab_io.h"
#include "libmesh/off_io.h"
#include "libmesh/medit_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/fro_io.h"
#include "libmesh/xdr_io.h"
#include "libmesh/legacy_xdr_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/abaqus_io.h"
#include "libmesh/checkpoint_io.h"

#include LIBMESH_INCLUDE_UNORDERED_MAP



// ------------------------------------------------------------
// Anonymous namespace for implementation details
namespace {
  bool is_parallel_file_format (const std::string &name)
  {
    // Certain mesh formats can support parallel I/O, including the
    // "new" Xdr format and the Nemesis format.
    return ((name.rfind(".xda") < name.size()) ||
	    (name.rfind(".xdr") < name.size()) ||
	    (name.rfind(".nem") < name.size()) ||
	    (name.rfind(".n") < name.size())   ||
            (name.rfind(".cp") < name.size())
	    );
  }
}


namespace libMesh
{


// ------------------------------------------------------------
// UnstructuredMesh class member functions
UnstructuredMesh::UnstructuredMesh (const Parallel::Communicator &comm,
				    unsigned int d) :
  MeshBase (comm,d)
{
  libmesh_assert (libMesh::initialized());
}



#ifndef LIBMESH_DISABLE_COMMWORLD
UnstructuredMesh::UnstructuredMesh (unsigned int d) :
  MeshBase (d)
{
  libmesh_assert (libMesh::initialized());
}
#endif



void UnstructuredMesh::copy_nodes_and_elements
  (const UnstructuredMesh& other_mesh)
{
  // We're assuming our subclass data needs no copy
  libmesh_assert_equal_to (_n_parts, other_mesh._n_parts);
  libmesh_assert_equal_to (_dim, other_mesh._dim);
  libmesh_assert_equal_to (_is_prepared, other_mesh._is_prepared);

  // We're assuming the other mesh has proper element number ordering,
  // so that we add parents before their children.
#ifdef DEBUG
  MeshTools::libmesh_assert_valid_amr_elem_ids(other_mesh);
#endif

  //Copy in Nodes
  {
    //Preallocate Memory if necessary
    this->reserve_nodes(other_mesh.n_nodes());

    const_node_iterator it = other_mesh.nodes_begin();
    const_node_iterator end = other_mesh.nodes_end();

    for (; it != end; ++it)
      {
        const Node *oldn = *it;

        // Add new nodes in old node Point locations
        /*Node *newn =*/ this->add_point(*oldn, oldn->id(), oldn->processor_id());

        // And start them off in the same subdomain
//        newn->processor_id() = oldn->processor_id();
      }
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
      const Elem *old = *it;
      //Build a new element
      Elem *newparent = old->parent() ?
          this->elem(old->parent()->id()) : NULL;
      AutoPtr<Elem> ap = Elem::build(old->type(), newparent);
      Elem * el = ap.release();

      el->subdomain_id() = old->subdomain_id();

      for (unsigned int s=0; s != old->n_sides(); ++s)
        if (old->neighbor(s) == remote_elem)
          el->set_neighbor(s, const_cast<RemoteElem*>(remote_elem));

#ifdef LIBMESH_ENABLE_AMR
      if (old->has_children())
        for (unsigned int c=0; c != old->n_children(); ++c)
          if (old->child(c) == remote_elem)
            el->add_child(const_cast<RemoteElem*>(remote_elem), c);

      //Create the parent's child pointers if necessary
      if (newparent)
        {
          unsigned int oldc = old->parent()->which_child_am_i(old);
          newparent->add_child(el, oldc);
        }

      // Copy the refinement flags
      el->set_refinement_flag(old->refinement_flag());
      el->set_p_refinement_flag(old->p_refinement_flag());
#endif // #ifdef LIBMESH_ENABLE_AMR

      //Assign all the nodes
      for(unsigned int i=0;i<el->n_nodes();i++)
        el->set_node(i) = &this->node(old->node(i));

      // And start it off in the same subdomain
      el->processor_id() = old->processor_id();

      // Give it the same id
      el->set_id(old->id());

      //Hold onto it
      this->add_elem(el);
    }
  }

  //Finally prepare the new Mesh for use.  Keep the same numbering and
  //partitioning but also the same renumbering and partitioning
  //policies as our source mesh.
  this->allow_renumbering(false);
  this->skip_partitioning(true);
  this->prepare_for_use();
  this->allow_renumbering(other_mesh.allow_renumbering());
  this->skip_partitioning(other_mesh.skip_partitioning());
}



UnstructuredMesh::~UnstructuredMesh ()
{
//  this->clear ();  // Nothing to clear at this level

  libmesh_assert (!libMesh::closed());
}





void UnstructuredMesh::find_neighbors (const bool reset_remote_elements,
				       const bool reset_current_list)
{
  // We might actually want to run this on an empty mesh
  // (e.g. the boundary mesh for a nonexistant bcid!)
  // libmesh_assert_not_equal_to (this->n_nodes(), 0);
  // libmesh_assert_not_equal_to (this->n_elem(), 0);

  // This function must be run on all processors at once
  parallel_object_only();

  START_LOG("find_neighbors()", "Mesh");

  const element_iterator el_end = this->elements_end();

  //TODO:[BSK] This should be removed later?!
  if (reset_current_list)
    for (element_iterator el = this->elements_begin(); el != el_end; ++el)
      {
	Elem* e = *el;
	for (unsigned int s=0; s<e->n_neighbors(); s++)
	  if (e->neighbor(s) != remote_elem ||
	      reset_remote_elements)
	    e->set_neighbor(s,NULL);
      }

  // Find neighboring elements by first finding elements
  // with identical side keys and then check to see if they
  // are neighbors
  {
    // data structures -- Use the hash_multimap if available
    typedef unsigned int                    key_type;
    typedef std::pair<Elem*, unsigned char> val_type;
    typedef std::pair<key_type, val_type>   key_val_pair;

    typedef LIBMESH_BEST_UNORDERED_MULTIMAP<key_type, val_type> map_type;

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
		    const AutoPtr<Elem> my_side(element->side(ms));

		    // Look at all the entries with an equivalent key
		    while (bounds.first != bounds.second)
		      {
			// Get the potential element
			Elem* neighbor = bounds.first->second.first;

			// Get the side for the neighboring element
			const unsigned int ns = bounds.first->second.second;
			const AutoPtr<Elem> their_side(neighbor->side(ns));
                        //libmesh_assert(my_side.get());
                        //libmesh_assert(their_side.get());

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
#if defined(LIBMESH_HAVE_UNORDERED_MAP) || defined(LIBMESH_HAVE_TR1_UNORDERED_MAP) || defined(LIBMESH_HAVE_HASH_MAP) || defined(LIBMESH_HAVE_EXT_HASH_MAP)
		side_to_elem_map.insert (kvp);
#else
		side_to_elem_map.insert (bounds.first,kvp);
#endif
	      }
	  }
      }
  }

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Here we look at all of the child elements which
   * don't already have valid neighbors.
   *
   * If a child element has a NULL neighbor it is
   * either because it is on the boundary or because
   * its neighbor is at a different level.  In the
   * latter case we must get the neighbor from the
   * parent.
   *
   * If a child element has a remote_elem neighbor
   * on a boundary it shares with its parent, that
   * info may have become out-dated through coarsening
   * of the neighbor's parent.  In this case, if the
   * parent's neighbor is active then the child should
   * share it.
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
          Elem* current_elem = *el;
          libmesh_assert(current_elem);
	  Elem* parent = current_elem->parent();
          libmesh_assert(parent);
	  const unsigned int my_child_num = parent->which_child_am_i(current_elem);

          for (unsigned int s=0; s < current_elem->n_neighbors(); s++)
            {
              if (current_elem->neighbor(s) == NULL ||
		  (current_elem->neighbor(s) == remote_elem &&
		   parent->is_child_on_side(my_child_num, s)))
                {
                  Elem *neigh = parent->neighbor(s);

	          // If neigh was refined and had non-subactive children
	          // made remote earlier, then a non-subactive elem should
	          // actually have one of those remote children as a
	          // neighbor
                  if (neigh && (neigh->ancestor()) && (!current_elem->subactive()))
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
                      libmesh_assert_not_equal_to (current_elem->processor_id(),
				                  this->processor_id());
#endif // DEBUG
                      neigh = const_cast<RemoteElem*>(remote_elem);
                    }

                  current_elem->set_neighbor(s, neigh);
#ifdef DEBUG
                  if (neigh != NULL && neigh != remote_elem)
                    // We ignore subactive elements here because
                    // we don't care about neighbors of subactive element.
                    if ((!neigh->active()) && (!current_elem->subactive()))
                      {
                        libMesh::err << "On processor " << this->processor_id()
                                      << std::endl;
                        libMesh::err << "Bad element ID = " << current_elem->id()
                          << ", Side " << s << ", Bad neighbor ID = " << neigh->id() << std::endl;
                        libMesh::err << "Bad element proc_ID = " << current_elem->processor_id()
                          << ", Bad neighbor proc_ID = " << neigh->processor_id() << std::endl;
                        libMesh::err << "Bad element size = " << current_elem->hmin()
                          << ", Bad neighbor size = " << neigh->hmin() << std::endl;
                        libMesh::err << "Bad element center = " << current_elem->centroid()
                          << ", Bad neighbor center = " << neigh->centroid() << std::endl;
                        libMesh::err << "ERROR: "
                          << (current_elem->active()?"Active":"Ancestor")
                          << " Element at level "
                          << current_elem->level() << std::endl;
                        libMesh::err << "with "
                          << (parent->active()?"active":
                              (parent->subactive()?"subactive":"ancestor"))
                          << " parent share "
                          << (neigh->subactive()?"subactive":"ancestor")
                          << " neighbor at level " << neigh->level()
                          << std::endl;
                        GMVIO(*this).write ("bad_mesh.gmv");
                        libmesh_error();
                      }
#endif // DEBUG
                }
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
			     MeshData* mesh_data,
			     bool skip_renumber_nodes_and_elements)
{
  // See if the file exists.  Perform this check on all processors
  // so that the code is terminated properly in the case that the
  // file does not exist.

  // For Nemesis files, the name we try to read will have suffixes
  // identifying processor rank
  if (name.rfind(".nem") + 4 == name.size() ||
      name.rfind(".n") + 2 == name.size())
    {
      std::ostringstream full_name;

      // Find the length of a string which represents the highest processor ID
      full_name << (this->n_processors());
      unsigned field_width = full_name.str().size();

      // reset the string stream
      full_name.str("");

      // And build up the full filename
      full_name << name
                << '.' << this->n_processors()
                << '.' << std::setfill('0') << std::setw(field_width) << this->processor_id();

      std::ifstream in (full_name.str().c_str());

      if (!in.good())
        {
	  libMesh::err << "ERROR: cannot locate specified file:\n\t"
		        << full_name.str()
		        << std::endl;
	  libmesh_error();
        }
    }
  else if(name.rfind(".cp")) {} // Do error checking in the reader
  else
    {
      std::ifstream in (name.c_str());

      if (!in.good())
        {
	  libMesh::err << "ERROR: cannot locate specified file:\n\t"
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
  if(name.rfind(".gmv")+4==name.size())
    {
      skip_renumber_nodes_and_elements =  true;
    }

  // Look for parallel formats first
  if (is_parallel_file_format(name))
    {
      // no need to handle bz2 files here -- the Xdr class does that.
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
	      this->allow_renumbering(false);
	      MeshCommunication().broadcast(*this);
	    }

          // libHilbert-enabled libMesh builds should construct files
          // with a canonical node ordering, which libHilbert-enabled
          // builds will be able to read in again regardless of any
	  // renumbering.  So in that case we're free to renumber.
	  // However, if either the writer or the reader of this file
	  // don't have libHilbert, then we'll have to skip
	  // renumbering because we need the numbering to remain
	  // consistent with any solution file we read in next.
#ifdef LIBMESH_HAVE_LIBHILBERT
	  // if (!xdr_io.libhilbert_ordering())
	  //   skip_renumber_nodes_and_elements = true;
#else
	  this->allow_renumbering(false);
#endif
	}
      else if (name.rfind(".nem") < name.size() ||
	       name.rfind(".n")   < name.size())
        Nemesis_IO(*this).read (name);
      else if (name.rfind(".cp") < name.size())
      {
        if(name.rfind(".cpa") < name.size())
          CheckpointIO(*this, false).read(name);
        else
          CheckpointIO(*this, true).read(name);
      }
    }

  // Serial mesh formats
  else
    {
      START_LOG("read()", "Mesh");

      // Read the file based on extension.  Only processor 0
      // needs to read the mesh.  It will then broadcast it and
      // the other processors will pick it up
      if (this->processor_id() == 0)
	{
          std::ostringstream pid_suffix;
          pid_suffix << '_' << getpid();
	  // Nasty hack for reading/writing zipped files
	  std::string new_name = name;
	  if (name.size() - name.rfind(".bz2") == 4)
	    {
#ifdef LIBMESH_HAVE_BZIP
	      new_name.erase(new_name.end() - 4, new_name.end());
              new_name += pid_suffix.str();
	      std::string system_string = "bunzip2 -f -k -c ";
	      system_string += name + " > " + new_name;
	      START_LOG("system(bunzip2)", "Mesh");
	      if (std::system(system_string.c_str()))
	        libmesh_file_error(system_string);
	      STOP_LOG("system(bunzip2)", "Mesh");
#else
              libMesh::err << "ERROR: need bzip2/bunzip2 to open .bz2 file "
		           << name << std::endl;
              libmesh_error();
#endif
	    }
          else if (name.size() - name.rfind(".xz") == 3)
            {
#ifdef LIBMESH_HAVE_XZ
	      new_name.erase(new_name.end() - 3, new_name.end());
              new_name += pid_suffix.str();
	      std::string system_string = "xz -f -d -k -c ";
	      system_string += name + " > " + new_name;
	      START_LOG("system(xz -d)", "XdrIO");
	      if (std::system(system_string.c_str()))
                libmesh_file_error(system_string);
	      STOP_LOG("system(xz -d)", "XdrIO");
#else
              libMesh::err << "ERROR: need xz to open .xz file "
		           << name << std::endl;
              libmesh_error();
#endif
            }

	  if (new_name.rfind(".mat") < new_name.size())
	    MatlabIO(*this).read(new_name);

	  else if (new_name.rfind(".ucd") < new_name.size())
	    UCDIO(*this).read (new_name);

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
		  libMesh::err << "Error! You must pass a "
			        << "valid MeshData pointer to "
			        << "read UNV files!" << std::endl;
		  libmesh_error();
		}
	      UNVIO(*this, *mesh_data).read (new_name);
	    }

	  else if ((new_name.rfind(".node")  < new_name.size()) ||
		   (new_name.rfind(".ele")   < new_name.size()))
	    TetGenIO(*this,mesh_data).read (new_name);

	  else if (new_name.rfind(".exd") < new_name.size() ||
		   new_name.rfind(".e") < new_name.size())
	    ExodusII_IO(*this).read (new_name);

	  else if (new_name.rfind(".msh") < new_name.size())
	    GmshIO(*this).read (new_name);

	  else if (new_name.rfind(".gmv") < new_name.size())
	    GMVIO(*this).read (new_name);

	  else if (new_name.rfind(".vtu") < new_name.size())
	    VTKIO(*this).read(new_name);

	  else if (new_name.rfind(".inp") < new_name.size())
	    AbaqusIO(*this).read(new_name);

	  else
	    {
	      libMesh::err << " ERROR: Unrecognized file extension: " << name
			    << "\n   I understand the following:\n\n"
			    << "     *.e    -- Sandia's ExodusII format\n"
			<< "     *.exd  -- Sandia's ExodusII format\n"
			<< "     *.gmv  -- LANL's General Mesh Viewer format\n"
			<< "     *.mat  -- Matlab triangular ASCII file\n"
			<< "     *.n    -- Sandia's Nemesis format\n"
			<< "     *.nem  -- Sandia's Nemesis format\n"
			<< "     *.off  -- OOGL OFF surface format\n"
			<< "     *.ucd  -- AVS's ASCII UCD format\n"
			<< "     *.unv  -- I-deas Universal format\n"
			<< "     *.vtu  -- Paraview VTK format\n"
			<< "     *.inp  -- Abaqus .inp format\n"
			<< "     *.xda  -- libMesh ASCII format\n"
			<< "     *.xdr  -- libMesh binary format\n"
			<< "     *.gz   -- any above format gzipped\n"
			<< "     *.bz2  -- any above format bzip2'ed\n"
			<< "     *.xz   -- any above format xzipped\n"
                        << "     *.cpa  -- libMesh Checkpoint ASCII format\n"
                        << "     *.cpr  -- libMesh Checkpoint binary format\n"

			<< std::endl;
	      libmesh_error();
	    }

	  // If we temporarily decompressed a file, remove the
	  // uncompressed version
	  if (name.size() - name.rfind(".bz2") == 4)
	    std::remove(new_name.c_str());
	  if (name.size() - name.rfind(".xz") == 3)
	    std::remove(new_name.c_str());
	}


      STOP_LOG("read()", "Mesh");

      // Send the mesh & bcs (which are now only on processor 0) to the other
      // processors
      MeshCommunication().broadcast (*this);
    }

  if (skip_renumber_nodes_and_elements)
    {
      // Use MeshBase::allow_renumbering() yourself instead.
      libmesh_deprecated();
      this->allow_renumbering(false);
    }

  // Done reading the mesh.  Now prepare it for use.
  this->prepare_for_use();
}



void UnstructuredMesh::write (const std::string& name,
			      MeshData* mesh_data)
{
  // parallel formats are special -- they may choose to write
  // separate files, let's not try to handle the zipping here.
  if (is_parallel_file_format(name))
    {
      // no need to handle bz2 files here -- the Xdr class does that.
      if (name.rfind(".xda") < name.size())
	XdrIO(*this).write(name);

      else if (name.rfind(".xdr") < name.size())
	XdrIO(*this,true).write(name);

      else if (name.rfind(".nem") < name.size() ||
               name.rfind(".n")   < name.size())
        Nemesis_IO(*this).write(name);
    }

  // serial file formats
  else
    {
      START_LOG("write()", "Mesh");

      // Nasty hack for reading/writing zipped files
      std::string new_name = name;
      processor_id_type pid_0 = 0;
      if (this->processor_id() == 0)
        pid_0 = getpid();
      this->comm().broadcast(pid_0);
      std::ostringstream pid_suffix;
      pid_suffix << '_' << pid_0;

      if (name.size() - name.rfind(".bz2") == 4)
        {
	  new_name.erase(new_name.end() - 4, new_name.end());
          new_name += pid_suffix.str();
        }
      else if (name.size() - name.rfind(".xz") == 3)
        {
	  new_name.erase(new_name.end() - 3, new_name.end());
          new_name += pid_suffix.str();
        }

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
        else if (new_name.rfind(".exd") < new_name.size() ||
                 new_name.rfind(".e") < new_name.size())
          ExodusII_IO(*this).write(new_name);
	else if (new_name.rfind(".mgf")  < new_name.size())
	  LegacyXdrIO(*this,true).write_mgf(new_name);

	else if (new_name.rfind(".unv") < new_name.size())
	  {
	    if (mesh_data == NULL)
	      {
		libMesh::err << "Error! You must pass a "
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
	    libMesh::err
              << " ERROR: Unrecognized file extension: " << name
              << "\n   I understand the following:\n\n"
              << "     *.dat   -- Tecplot ASCII file\n"
              << "     *.e     -- Sandia's ExodusII format\n"
              << "     *.exd   -- Sandia's ExodusII format\n"
              << "     *.fro   -- ACDL's surface triangulation file\n"
              << "     *.gmv   -- LANL's GMV (General Mesh Viewer) format\n"
              << "     *.mesh  -- MEdit mesh format\n"
              << "     *.mgf   -- MGF binary mesh format\n"
              << "     *.msh   -- GMSH ASCII file\n"
              << "     *.n     -- Sandia's Nemesis format\n"
              << "     *.nem   -- Sandia's Nemesis format\n"
              << "     *.plt   -- Tecplot binary file\n"
              << "     *.poly  -- TetGen ASCII file\n"
              << "     *.ucd   -- AVS's ASCII UCD format\n"
              << "     *.ugrid -- Kelly's DIVA ASCII format\n"
              << "     *.unv   -- I-deas Universal format\n"
              << "     *.vtu   -- VTK (paraview-readable) format\n"
              << "     *.xda   -- libMesh ASCII format\n"
              << "     *.xdr   -- libMesh binary format,\n"
              << std::endl
              << "\n Exiting without writing output\n";
	  }
      }

      // Nasty hack for reading/writing zipped files
      if (name.size() - name.rfind(".bz2") == 4)
	{
	  START_LOG("system(bzip2)", "Mesh");
	  if (this->processor_id() == 0)
	    {
	      std::string system_string = "bzip2 -f -c ";
	      system_string += new_name + " > " + name;
	      if (std::system(system_string.c_str()))
		libmesh_file_error(system_string);
	      std::remove(new_name.c_str());
	    }
	  this->comm().barrier();
	  STOP_LOG("system(bzip2)", "Mesh");
	}
      if (name.size() - name.rfind(".xz") == 3)
	{
	  START_LOG("system(xz)", "Mesh");
	  if (this->processor_id() == 0)
	    {
	      std::string system_string = "xz -f -c ";
	      system_string += new_name + " > " + name;
	      if (std::system(system_string.c_str()))
		libmesh_file_error(system_string);
	      std::remove(new_name.c_str());
	    }
	  this->comm().barrier();
	  STOP_LOG("system(xz)", "Mesh");
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
      libMesh::err
        << " ERROR: Unrecognized file extension: " << name
	<< "\n   I understand the following:\n\n"
	<< "     *.dat  -- Tecplot ASCII file\n"
	<< "     *.gmv  -- LANL's GMV (General Mesh Viewer) format\n"
	<< "     *.plt  -- Tecplot binary file\n"
	<< "     *.pvtu -- Paraview VTK file\n"
	<< "\n Exiting without writing output\n";
    }

  STOP_LOG("write()", "Mesh");
}





void UnstructuredMesh::create_pid_mesh(UnstructuredMesh& pid_mesh,
			               const processor_id_type pid) const
{

  // Issue a warning if the number the number of processors
  // currently available is less that that requested for
  // partitioning.  This is not necessarily an error since
  // you may run on one processor and still partition the
  // mesh into several partitions.
#ifdef DEBUG
  if (this->n_processors() < pid)
    {
      libMesh::out << "WARNING:  You are creating a "
		    << "mesh for a processor id (="
		    << pid
		    << ") greater than "
		    << "the number of processors available for "
		    << "the calculation. (="
		    << this->n_processors()
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
  libmesh_assert_not_equal_to (this->n_nodes(), 0);
  libmesh_assert_not_equal_to (this->n_elem(), 0);

  // How the nodes on this mesh will be renumbered to nodes
  // on the new_mesh.
  std::vector<dof_id_type> new_node_numbers (this->n_nodes());

  std::fill (new_node_numbers.begin(),
	     new_node_numbers.end(),
	     DofObject::invalid_id);



  // the number of nodes on the new mesh, will be incremented
  dof_id_type n_new_nodes = 0;
  dof_id_type n_new_elem  = 0;

  for (; it != it_end; ++it)
    {
      // increment the new element counter
      n_new_elem++;

      const Elem* old_elem = *it;

      // Add an equivalent element type to the new_mesh
      Elem* new_elem =
	new_mesh.add_elem (Elem::build(old_elem->type()).release());

      libmesh_assert(new_elem);

      // Loop over the nodes on this element.
      for (unsigned int n=0; n<old_elem->n_nodes(); n++)
	{
	  libmesh_assert_less (old_elem->node(n), new_node_numbers.size());

	  if (new_node_numbers[old_elem->node(n)] == DofObject::invalid_id)
	    {
	      new_node_numbers[old_elem->node(n)] = n_new_nodes;

	      // Add this node to the new mesh
	      new_mesh.add_point (old_elem->point(n));

	      // Increment the new node counter
	      n_new_nodes++;
	    }

	  // Define this element's connectivity on the new mesh
	  libmesh_assert_less (new_node_numbers[old_elem->node(n)], new_mesh.n_nodes());

	  new_elem->set_node(n) = new_mesh.node_ptr (new_node_numbers[old_elem->node(n)]);
	}

      // Copy ids for this element
      new_elem->subdomain_id() = old_elem->subdomain_id();
      new_elem->processor_id() = old_elem->processor_id();

      // Maybe add boundary conditions for this element
      for (unsigned int s=0; s<old_elem->n_sides(); s++)
// We're supporting boundary ids on internal sides now
//	if (old_elem->neighbor(s) == NULL)
          {
            const std::vector<boundary_id_type>& bc_ids = this->boundary_info->boundary_ids(old_elem, s);
            for (std::vector<boundary_id_type>::const_iterator id_it=bc_ids.begin(); id_it!=bc_ids.end(); ++id_it)
              {
                const boundary_id_type bc_id = *id_it;
	        if (bc_id != this->boundary_info->invalid_id)
	        new_mesh.boundary_info->add_side (new_elem,
					          s,
					          bc_id);
              }
          }
    } // end loop over elements


  // Prepare the new_mesh for use
  new_mesh.prepare_for_use(/*skip_renumber =*/false);

}



#ifdef LIBMESH_ENABLE_AMR
bool UnstructuredMesh::contract ()
{
  START_LOG ("contract()", "Mesh");

  // Flag indicating if this call actually changes the mesh
  bool mesh_changed = false;

  element_iterator in        = elements_begin();
  const element_iterator end = elements_end();

#ifdef DEBUG
  for ( ; in != end; ++in)
    if (*in != NULL)
      {
	Elem* el = *in;
	libmesh_assert(el->active() || el->subactive() || el->ancestor());
      }
  in = elements_begin();
#endif

  // Loop over the elements.
  for ( ; in != end; ++in)
    if (*in != NULL)
      {
	Elem* el = *in;

	// Delete all the subactive ones
	if (el->subactive())
	  {
	    // No level-0 element should be subactive.
	    // Note that we CAN'T test elem->level(), as that
	    // touches elem->parent()->dim(), and elem->parent()
	    // might have already been deleted!
	    libmesh_assert(el->parent());

	    // Delete the element
	    // This just sets a pointer to NULL, and doesn't
	    // invalidate any iterators
	    this->delete_elem(el);

	    // the mesh has certainly changed
	    mesh_changed = true;
	  }
	else
	  {
	    // Compress all the active ones
	    if (el->active())
	      el->contract();
	    else
	      libmesh_assert (el->ancestor());
	  }
      }

  // Strip any newly-created NULL voids out of the element array
  this->renumber_nodes_and_elements();

  // FIXME: Need to understand why deleting subactive children
  // invalidates the point locator.  For now we will clear it explicitly
  this->clear_point_locator();

  STOP_LOG ("contract()", "Mesh");

  return mesh_changed;
}
#endif // #ifdef LIBMESH_ENABLE_AMR

} // namespace libMesh
