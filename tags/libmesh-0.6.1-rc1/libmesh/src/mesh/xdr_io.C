// $Id: xdr_io.C,v 1.25 2007-05-23 23:36:11 roystgnr Exp $

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
#include <iostream>
#include <iomanip>

#include <vector>
#include <string>

// Local includes
#include "xdr_io.h"
#include "mesh_base.h"
#include "mesh_data.h"
#include "mesh_tools.h" // MeshTools::n_levels(mesh)
#include "cell_hex27.h" // Needed for MGF-style Hex27 meshes
#include "boundary_info.h"
#include "libmesh_logging.h"
#include "xdr_mgf.h"
#include "xdr_mesh.h"
#include "xdr_mhead.h"
#include "xdr_soln.h"
#include "xdr_shead.h"

#ifdef USE_COMPLEX_NUMBERS
#include "utility.h"
#endif







// ------------------------------------------------------------
// XdrIO members
XdrIO::XdrIO (MeshBase& mesh, const bool binary) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase>(mesh),
  _binary (binary)
{
}




XdrIO::XdrIO (const MeshBase& mesh, const bool binary) :
  MeshOutput<MeshBase>(mesh),
  _binary (binary)
{
}




XdrIO::~XdrIO ()
{
}




bool & XdrIO::binary ()
{
  return _binary;
}




bool XdrIO::binary () const
{
  return _binary;
}


void XdrIO::read (const std::string& name)
{
  if (this->binary())
    this->read_binary (name);
  else
    this->read_ascii  (name);
}



void XdrIO::read_mgf (const std::string& name)
{
  if (this->binary())
    this->read_binary (name, XdrIO::MGF);
  else
    this->read_ascii  (name, XdrIO::MGF);
}



void XdrIO::write (const std::string& name)
{
  if (this->binary())
    this->write_binary (name);
  else
    this->write_ascii  (name);
}



void XdrIO::write_mgf (const std::string& name)
{
  if (this->binary())
    this->write_binary (name, XdrIO::MGF);
  else
    this->write_ascii  (name, XdrIO::MGF);
}



void XdrIO::read_mgf_soln (const std::string& name,
			   std::vector<Number>& soln,
			   std::vector<std::string>& var_names) const
{
  here();
  std::cerr << "WARNING: this method is deprecated and will disappear soon!"
	    << std::endl;
  
#ifdef USE_COMPLEX_NUMBERS
  
  // buffer for writing separately
  std::vector<Real> real_soln;
  std::vector<Real> imag_soln;
  
  Utility::prepare_complex_data (soln, real_soln, imag_soln);

  this->read_soln (Utility::complex_filename(name, 0), 
		   real_soln, 
		   var_names);
  
  this->read_soln (Utility::complex_filename(name, 1), 
		   imag_soln, 
		   var_names);
  
#else
  
  this->read_soln (name, soln, var_names);
      
#endif
}



void XdrIO::write_mgf_soln (const std::string& name,
			    std::vector<Number>& soln,
			    std::vector<std::string>& var_names) const
{
  here();
  std::cerr << "WARNING: this method is deprecated and will disappear soon!"
	    << std::endl;
  
#ifdef USE_COMPLEX_NUMBERS
  
  // buffer for writing separately
  std::vector<Real> real_soln;
  std::vector<Real> imag_soln;
  
  Utility::prepare_complex_data (soln, real_soln, imag_soln);
  
  this->write_soln (Utility::complex_filename(name, 0), 
		    real_soln, 
		    var_names);
  
  this->write_soln (Utility::complex_filename(name, 1), 
		    imag_soln, 
		    var_names);
  
#else
  
  this->write_soln (name, soln, var_names);
      
#endif
}



void XdrIO::read_ascii (const std::string& name, const XdrIO::FileFormat originator)
{
  // get a writeable reference to the underlying mesh
  MeshBase& mesh = MeshInput<MeshBase>::mesh();
  
  // clear any existing mesh data
  mesh.clear();
    
  // read the mesh
  this->read_mesh (name, originator);
}



void XdrIO::read_binary (const std::string& name, const XdrIO::FileFormat originator)
{
#ifndef HAVE_XDR

  std::cerr << "WARNING: Compiled without XDR binary support.\n"
	    << "Will try ASCII instead" << std::endl << std::endl;

  this->read_ascii (name);
  
#else
  
  // get a writeable reference to the underlying mesh
  MeshBase& mesh = MeshInput<MeshBase>::mesh();
  
  // clear any existing mesh data
  mesh.clear();

  // read the mesh
  this->read_mesh (name, originator);
  
#endif
}



void XdrIO::write_ascii (const std::string& name, const XdrIO::FileFormat originator)
{
  this->write_mesh (name, originator);
}



void XdrIO::write_binary (const std::string& name, const XdrIO::FileFormat originator)
{
#ifndef HAVE_XDR

  std::cerr << "WARNING: Compiled without XDR binary support.\n"
	    << "Will try ASCII instead" << std::endl << std::endl;

  this->write_ascii (name);

#else
  
  this->write_mesh (name, originator);  
  
#endif
}



void XdrIO::read_mesh (const std::string& name,
		       const XdrIO::FileFormat originator,
		       MeshData* mesh_data)
{
  // get a writeable reference to the mesh
  MeshBase& mesh = MeshInput<MeshBase>::mesh();

  // clear any data in the mesh
  mesh.clear();
  
  // Create an XdrMESH object.
  XdrMESH m;

  // Create a pointer
  // to an XdrMESH file
  // header.
  XdrMHEAD mh;

  // Open the XDR file for reading.
  // Note 1: Provide an additional argument
  // to specify the dimension.
  //
  // Note 2: Has to do the right thing for
  // both binary and ASCII files.
  m.set_orig_flag(originator);
  m.init((this->binary() ? XdrMGF::DECODE : XdrMGF::R_ASCII), name.c_str(), 0); // mesh files are always number 0 ...

  // From here on, things depend
  // on whether we are reading or
  // writing!  First, we define
  // header variables that may
  // be read OR written.
  unsigned int              n_blocks = 0;
  unsigned int              n_levels = 0;
  
  if (m.get_orig_flag() == XdrIO::LIBM)
    n_levels = m.get_num_levels();
  
  
  std::vector<ElemType>     etypes;
  std::vector<unsigned int> neeb;
	
  // Get the information from
  // the header, and place it
  // in the header pointer.
  m.header(&mh);
	
  // Read information from the
  // file header.  This depends on
  // whether its a libMesh or MGF mesh.
  const int numElem     = mh.getNumEl();
  const int numNodes    = mh.getNumNodes();
  const int totalWeight = mh.getSumWghts();
  const int numBCs      = mh.getNumBCs();
  
  // If a libMesh-type mesh, read the augmented mesh information
  if ((m.get_orig_flag() == XdrIO::DEAL) || (m.get_orig_flag() == XdrIO::LIBM))
    {
      // Read augmented header
      n_blocks = mh.get_n_blocks();
      
      etypes.resize(n_blocks);
      mh.get_block_elt_types(etypes);
      
      mh.get_num_elem_each_block(neeb); 
    }

  
  
  // Read the connectivity  
  std::vector<int> conn;
  
  // Now that we know the
  // number of nodes and elements,
  // we can resize the
  // appropriate vectors if we are
  // reading information in.
  mesh.reserve_nodes (numNodes);
  mesh.reserve_elem  (numElem);
  
  // Each element stores two extra
  // locations: one which tells
  // what type of element it is,
  // and one which tells how many
  // nodes it has. Therefore,
  // the total number of nodes
  // (totalWeight) must be augmented
  // by 2 times the number of elements
  // in order to read in the entire
  // connectivity array.
  
  // Note: This section now depends on
  // whether we are reading an old-style libMesh, 
  // MGF, or a new-style libMesh mesh.  
  if (m.get_orig_flag() == XdrIO::DEAL) 
    {
      conn.resize(totalWeight);
      m.Icon(&conn[0], 1, totalWeight);
    }
  
  else if (m.get_orig_flag() == XdrIO::MGF) 
    {
      conn.resize(totalWeight+(2*numElem));
      m.Icon(&conn[0], 1, totalWeight+(2*numElem));
    }

  else if (m.get_orig_flag() == XdrIO::LIBM)
    {
      conn.resize(totalWeight);
      m.Icon(&conn[0], 1, totalWeight);
    }
  
  else
    {
      // I don't know what type of mesh it is.
      error();
    }


  // read in the nodal coordinates and form points.
  {
    std::vector<Real> coords(numNodes*mesh.spatial_dimension()); // Always use three coords per node
    m.coord(&coords[0], mesh.spatial_dimension(), numNodes);


  
    // Form Nodes out of
    // the coordinates.  If the    
    // MeshData object is active,
    // add the nodes and ids also          
    // to its map.
    for (int innd=0; innd<numNodes; ++innd)      
      {
	Node* node = mesh.add_point (Point(coords[0+innd*3],  
					   coords[1+innd*3],
					   coords[2+innd*3]));
				       
	if (mesh_data != NULL)
	  if (mesh_data->active())
	    {
	      // add the id to the MeshData, so that
	      // it knows the foreign id, even when 
	      // the underlying mesh got re-numbered,
	      // refined, elements/nodes added...   
	      mesh_data->add_foreign_node_id(node, innd);
	    }
      }  
  }

  
  
  // Build the elements.
  // Note: If the originator was MGF, we don't
  // have to do much checking ...
  // all the elements are Hex27.
  // If the originator was
  // this code, we have to loop over
  // et and neeb to read in all the
  // elements correctly.
  //
  // (This used to be before the coords block, but it
  // had to change now that elements store pointers to
  // nodes.  The nodes must exist before we assign them to
  // the elements. BSK, 1/13/2003)
  if ((m.get_orig_flag() == XdrIO::DEAL) || (m.get_orig_flag() == XdrIO::LIBM)) 
    {
      unsigned int lastConnIndex = 0;
      unsigned int lastFaceIndex = 0;

      // This map keeps track of elements we've previously
      // constructed, to avoid O(n) lookup times for parent pointers
      // and to enable elements to be added in ascending ID order
      std::map<unsigned int, Elem*> parents;

      for (unsigned int level=0; level<=n_levels; level++)
      {
        for (unsigned int idx=0; idx<n_blocks; idx++)  
        {
          for (unsigned int e=lastFaceIndex; e<lastFaceIndex+neeb[level*n_blocks+idx]; e++)  
          {
            // Build a temporary element of the right type, so we know how
            // connectivity entries will be on the line for this element.
            AutoPtr<Elem> temp_elem = Elem::build(etypes[idx]);

            // A pointer to the element which will eventually be added to the mesh.
            Elem* elem;

            // New-style libMesh mesh
            if (m.get_orig_flag() == XdrIO::LIBM)
            {
              unsigned int self_ID   = conn[lastConnIndex + temp_elem->n_nodes()];

#ifdef ENABLE_AMR
              unsigned int parent_ID = conn[lastConnIndex + temp_elem->n_nodes()+1];

              if (level > 0)
              {
                // Do a linear search for the parent
                Elem* my_parent;

                // Search for parent in the parents map (log(n))
                START_LOG("log(n) search for parent", "XdrIO::read_mesh");
                std::map<unsigned int, Elem*>::iterator it = parents.find(parent_ID);
                STOP_LOG("log(n) search for parent", "XdrIO::read_mesh");
                
                // If the parent was not previously added, we cannot continue.
                if (it == parents.end())
                {
                  std::cerr << "Parent element with ID " << parent_ID 
                            << " not found." << std::endl; 
                  error();
                }

                // Set the my_parent pointer
                my_parent = (*it).second;

                // my_parent is now INACTIVE, since he has children
                my_parent->set_refinement_flag(Elem::INACTIVE);
               
                // Now that we know the parent, build the child
                elem = Elem::build(etypes[idx],my_parent).release();

                // The new child is marked as JUST_REFINED
                elem->set_refinement_flag(Elem::JUST_REFINED); 
                
                // Tell the parent about his new child
                my_parent->add_child(elem);

                // sanity check
                assert (my_parent->type() == elem->type());
              }

              // Add level-0 elements to the mesh 
              else
#endif // #ifdef ENABLE_AMR
              {
                elem = Elem::build(etypes[idx]).release();
              }

              // Assign the newly-added element's ID so that future 
              // children which may be added can find it correctly.
              elem->set_id() = self_ID;
                
              // Add this element to the map, it may be a parent for a future element
              START_LOG("insert elem into map", "XdrIO::read_mesh");
              parents[self_ID] = elem;
              STOP_LOG("insert elem into map", "XdrIO::read_mesh");
            }

            // MGF-style meshes
            else
            {
              elem = mesh.add_elem(Elem::build(etypes[idx]).release());
            }
            
            // Add elements with the same id as in libMesh.  
            // Provided the data files that MeshData reads    
            // were only written with MeshData, then this      
            // should work properly.  This is an inline
            // function, so that for disabled MeshData, this
            // should not induce too much cost
            if (mesh_data != NULL)
              mesh_data->add_foreign_elem_id (elem, e);

            // Set the node pointers of the newly-created element
            for (unsigned int innd=0; innd < elem->n_nodes(); innd++)
            {
              elem->set_node(innd) = mesh.node_ptr(conn[innd+lastConnIndex]);
            }

            lastConnIndex += (m.get_orig_flag() == XdrIO::LIBM) ? (elem->n_nodes()+2) : elem->n_nodes();
          }
          lastFaceIndex += neeb[idx];
        }
        
      }

      if (m.get_orig_flag() == XdrIO::LIBM)
        // Iterate in ascending elem ID order
        for (std::map<unsigned int, Elem *>::iterator i =
             parents.begin();
             i != parents.end(); ++i)
          {
            Elem *elem = i->second;
            if (elem)
              mesh.add_elem(elem);
            else
              // We can probably handle this, but we don't expect it
              error();
          }

#ifdef ENABLE_AMR
      // All the elements at each level have been added, and their node pointers
      // have been set.  Now compute the node keys to put the mesh into a state consistent
      // with the state after being constructed through normal refinements. 
      MeshBase::element_iterator it = mesh.elements_begin();
      const MeshBase::element_iterator end = mesh.elements_end();
      for (; it!=end; ++it)
        (*it)->compute_children_node_keys();
#endif // #ifdef ENABLE_AMR
    }
 
  // MGF-style (1) Hex27 mesh
  else if (m.get_orig_flag() == XdrIO::MGF) 
    {
      
#ifdef DEBUG
      if (mesh_data != NULL)
	if (mesh_data->active())
	  {
	    std::cerr << "ERROR: MeshData not implemented for MGF-style mesh."
		      << std::endl;
	    error();
	  }
#endif
      
      for (int ielm=0; ielm < numElem; ++ielm)
	{
	  Elem* elem = mesh.add_elem(new Hex27);
	  
	  for (int innd=0; innd < 27; ++innd)
	    elem->set_node(innd) = mesh.node_ptr(conn[innd+2+(27+2)*ielm]);	
	}
    }

  
  // tell the MeshData object that we are finished 
  // reading data
  if (mesh_data != NULL)
    mesh_data->close_foreign_id_maps ();
  
  // Free memory used in
  // the connectivity
  // vector.
  conn.clear();


  // If we are reading,
  // read in the BCs
  // from the mesh file,
  // otherwise write the
  // boundary conditions
  // if the BoundaryInfo
  // object exists.
  if (numBCs > 0)
    {
      std::vector<int> bcs(numBCs*3);

      // Read the BCs from the XDR file
      m.BC(&bcs[0], numBCs);
  
      // Add to the boundary_info
      for (int ibc=0; ibc < numBCs; ibc++)
	mesh.boundary_info->add_side(bcs[0+ibc*3], bcs[1+ibc*3], bcs[2+ibc*3]);
    }
}



void XdrIO::write_mesh (const std::string& name,
			const XdrIO::FileFormat originator)
{
  // get a read-only reference to the mesh
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();


  
  // Create an XdrMESH object.
  XdrMESH m;

  // Create a pointer
  // to an XdrMESH file
  // header.
  XdrMHEAD mh;

  // Open the XDR file for writing.
  // Note 1: Provide an additional argument
  // to specify the dimension.
  //
  // Note 2: Has to do the right thing for
  // both binary and ASCII files.
  m.set_orig_flag(originator);

  // From here on, things depend
  // on whether we are reading or
  // writing!  First, we define
  // header variables that may
  // be read OR written.
  std::vector<unsigned int> neeb;
  std::vector<ElemType> etypes;
  

  int n_non_subactive = 0;
  int non_subactive_weight = 0;

  // This map will associate 
  // the distance from the beginning of the set
  // to each node ID with the node ID itself.
  std::map<unsigned int, unsigned int> node_map;

  {
    // For each non-subactive element:
    // 1.) Increment the number of non subactive elements
    // 2.) Accumulate the total weight
    // 3.) Add the node ids to a set of non subactive node ids 
    std::set<unsigned int> not_subactive_node_ids;
    MeshBase::const_element_iterator el = mesh.elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.elements_end();
    for( ; el != end_el; ++el)
    {
      Elem* elem = (*el);
      if(!elem->subactive())
      {
        n_non_subactive++;
        non_subactive_weight += elem->n_nodes();

        for (unsigned int n=0; n<elem->n_nodes(); ++n)
          not_subactive_node_ids.insert(elem->node(n));
      }
    }

    // Now that the set is built, most of the hard work is done.  We build
    // the map next and let the set go out of scope.
    std::set<unsigned int>::iterator it = not_subactive_node_ids.begin();
    const std::set<unsigned int>::iterator end = not_subactive_node_ids.end();
    unsigned int cnt=0;
    for (; it!=end; ++it)
      node_map[*it] = cnt++;
  }


  const int                   numElem  = n_non_subactive;       
  const int                   numBCs   = mesh.boundary_info->n_boundary_conds();
  const unsigned int          n_levels = MeshTools::n_levels(mesh);
  
  // Fill the etypes vector with all of the element types found in the mesh
  MeshTools::elem_types(mesh, etypes);
  
  // store number of elements in each block at each refinement level
  neeb.resize((n_levels+1)*etypes.size()); 

  // Store a variable for the number of element types                   
  const unsigned int n_el_types = etypes.size();
	
  m.set_num_levels(n_levels);

  // The last argument is zero because mesh files are always number 0 ...
  m.init((this->binary() ? XdrMGF::ENCODE : XdrMGF::W_ASCII), name.c_str(), 0); 

  // Loop over all levels and all element types to set the entries of neeb
  for(unsigned int level=0; level<=n_levels; level++)
    for (unsigned int el_type=0; el_type<n_el_types; el_type++)
      neeb[level*n_el_types + el_type] = 
        MeshTools::n_non_subactive_elem_of_type_at_level(mesh, etypes[el_type], level);
        // gotta change this function name!!!


  // Now we check to see if we're doing
  // MGF-style headers or libMesh-style
  // "augmented" headers.  An
  // augmented header contains 
  // information about mesh blocks,
  // allowing us to optimize storage
  // and minimize IO requirements
  // for these meshes.
  if ((m.get_orig_flag() == XdrIO::DEAL) || (m.get_orig_flag() == XdrIO::LIBM))
    {
      mh.set_n_blocks(etypes.size());
      mh.set_block_elt_types(etypes);
      mh.set_num_elem_each_block(neeb);
    }
  else
    assert(etypes.size() == 1);
  
  mh.setNumEl(numElem);
  mh.setNumNodes(node_map.size());
  mh.setStrSize(65536);
 
  // set a local variable for the total weight of the mesh
  int totalWeight =0;

  if (m.get_orig_flag() == XdrIO::DEAL)  // old-style LibMesh
    totalWeight=MeshTools::total_weight(mesh);

  else if (m.get_orig_flag() == XdrIO::MGF) // MGF-style
    totalWeight = MeshTools::total_weight(mesh)+2*numElem;

  else if (m.get_orig_flag() == XdrIO::LIBM) // new-style LibMesh
    totalWeight = non_subactive_weight+2*numElem;

  else
    error();
    
  // Set the total weight in the header
  mh.setSumWghts(totalWeight);
        
  mh.setNumBCs(numBCs);
  mh.setId("Id String");       // You can put whatever you want, it will be ignored 
  mh.setTitle("Title String"); // You can put whatever you want, it will be ignored 
  
  // Put the information
  // in the XDR file.
  m.header(&mh); 
  
  
  // Write the connectivity  
  {
    std::vector<int> conn;
    XdrIO::FileFormat orig_type = m.get_orig_flag();
   
    // Resize the connectivity vector to hold all the connectivity for the mesh
    conn.resize(totalWeight);
    
    unsigned int lastConnIndex = 0;
    unsigned int nn = 0;
   
    // Loop over levels and types again, write connectivity information to conn.
    for (unsigned int level=0; level<=n_levels; level++)
      for (unsigned int idx=0; idx<etypes.size(); idx++)
      {
        nn = lastConnIndex = 0;

        for (unsigned int e=0; e<mesh.n_elem(); e++)
          if ((mesh.elem(e)->type()  == etypes[idx]) && 
              (mesh.elem(e)->level() == level)       &&
              !mesh.elem(e)->subactive())
          {
            int nstart=0;
            
            if (orig_type == XdrIO::DEAL)
              nn = mesh.elem(e)->n_nodes();

            else if (orig_type == XdrIO::MGF)
            {
              nstart=2; // ignore the 27 and 0 entries
              nn = mesh.elem(e)->n_nodes()+2;
              conn[lastConnIndex + 0] = 27;
              conn[lastConnIndex + 1] = 0;
            }

            else if (orig_type == XdrIO::LIBM) // LIBMESH format
              nn = mesh.elem(e)->n_nodes() + 2;

            else
              error();

            // Loop over the connectivity entries for this element and write to conn.
            START_LOG("set connectivity", "XdrIO::write_mesh");
            const unsigned int loopmax = (orig_type==XdrIO::LIBM) ? nn-2 : nn;
            for (unsigned int n=nstart; n<loopmax; n++)
            {
              unsigned int connectivity_value=0;

              // old-style Libmesh and MGF meshes
              if (orig_type != XdrIO::LIBM)
                connectivity_value = mesh.elem(e)->node(n-nstart);

              // new-style libMesh meshes: compress the connectivity entries to account for
              // subactive nodes that will not be in the mesh we write out.
              else
              {
                std::map<unsigned int, unsigned int>::iterator pos = 
                  node_map.find(mesh.elem(e)->node(n-nstart));

                assert (pos != node_map.end());

                connectivity_value = (*pos).second;
              }
              conn[lastConnIndex + n] = connectivity_value;
            }
            STOP_LOG("set connectivity", "XdrIO::write_mesh");

            // In the case of an adaptive mesh, set last 2 entries to this ID and parent ID
            if (orig_type == XdrIO::LIBM)
            {
              int self_ID = mesh.elem(e)->id();
              int parent_ID = -1;
              if(level != 0)
                parent_ID = mesh.elem(e)->parent()->id();

              // Self ID is the second-to-last entry, Parent ID is the last
              // entry on each connectivity line
              conn[lastConnIndex+nn-2] = self_ID;
              conn[lastConnIndex+nn-1] = parent_ID;
            }

            lastConnIndex += nn;
          }

        // Send conn to the XDR file.  If there are no elements of this level and type,
        // then nn will be zero, and we there is no connectivity to write. 
        if (nn != 0)
          m.Icon(&conn[0], nn, lastConnIndex/nn);
      }
  }
    
  // create the vector of coords and send
  // it to the XDR file.
  {
    std::vector<Real> coords;
    
    coords.resize(mesh.spatial_dimension()*node_map.size());
    int lastIndex=0;

    std::map<unsigned int,unsigned int>::iterator it = node_map.begin();
    const std::map<unsigned int,unsigned int>::iterator end = node_map.end();
    for (; it != end; ++it)
      {
        const Point& p = mesh.node((*it).first);

        coords[lastIndex+0] = p(0);
        coords[lastIndex+1] = p(1);
        coords[lastIndex+2] = p(2);
        lastIndex += 3;
      }
   
    // Put the nodes in the XDR file
    m.coord(&coords[0], mesh.spatial_dimension(), node_map.size()); 
  }

  
  // write the
  // boundary conditions
  // if the BoundaryInfo
  // object exists.
  if (numBCs > 0)
    {
      std::vector<int> bcs(numBCs*3);
    
      //std::cout << "numBCs=" << numBCs << std::endl;
    
      //std::cout << "Preparing to write boundary conditions." << std::endl;
      std::vector<unsigned int> elem_list;
      std::vector<unsigned short int> side_list;
      std::vector<short int> elem_id_list;
      
      mesh.boundary_info->build_side_list (elem_list, side_list, elem_id_list);
    
      for (int ibc=0;  ibc<numBCs; ibc++)
	{
	  bcs[0+ibc*3] = elem_list[ibc];
	  bcs[1+ibc*3] = side_list[ibc];
	  bcs[2+ibc*3] = elem_id_list[ibc];
	}
    
      // Put the BCs in the XDR file
      m.BC(&bcs[0], numBCs);
    }
}



void XdrIO::read_soln (const std::string& name,
		       std::vector<Real>& soln,
		       std::vector<std::string>& var_names) const
{
  // Create an XdrSOLN object.
  XdrSOLN s;
  
  // Create an XdrSHEAD object.
  XdrSHEAD sh;
  
  // Open the XDR file for
  // reading or writing.
  // Note 1: Provide an additional argument
  // to specify the dimension.
  //
  // Note 2: Has to do the right thing for
  // both binary and ASCII files.
  s.init((this->binary() ? XdrMGF::DECODE : XdrMGF::R_ASCII), name.c_str(), 0); // mesh files are always number 0 ...
  
  // From here on, things depend
  // on whether we are reading or
  // writing!  First, we define
  // header variables that may
  // be read OR written.
  int         numVar      = 0;       
  int         numNodes    = 0;
  const char* varNames;
	
  // Get the information from
  // the header, and place it
  // in the header pointer.
  s.header(&sh);
	
  // Read information from the
  // file header.  This depends on
  // whether its a libMesh or MGF mesh.
  numVar   = sh.getWrtVar();
  numNodes = sh.getNumNodes();
  varNames = sh.getVarTitle();
	
  // Get the variable names
  {	  
    var_names.resize(numVar);
    
    const char* p = varNames;
    
    for (int i=0; i<numVar; i++)
      {
	var_names[i] = p;
	p += std::strlen(p) + 1;
      }
  }
  
  // Read the soln vector
  soln.resize(numVar*numNodes);
    
  s.values(&soln[0], numNodes);	
}  



void XdrIO::write_soln (const std::string& name,
			std::vector<Real>& soln,
			std::vector<std::string>& var_names) const
{
  // get a read-only reference to the mesh
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();
  
  // Create an XdrSOLN object.
  XdrSOLN s;
  
  // Create an XdrSHEAD object.
  XdrSHEAD sh;
  
  // Open the XDR file for
  // reading or writing.
  // Note 1: Provide an additional argument
  // to specify the dimension.
  //
  // Note 2: Has to do the right thing for
  // both binary and ASCII files.
  s.init((this->binary() ? XdrMGF::ENCODE : XdrMGF::W_ASCII), name.c_str(), 0); // mesh files are always number 0 ...

  // Build the header
  sh.setWrtVar(var_names.size());
  sh.setNumVar(var_names.size());
  sh.setNumNodes(mesh.n_nodes());
  sh.setNumBCs(mesh.boundary_info->n_boundary_conds());
  sh.setMeshCnt(0);
  sh.setKstep(0);
  sh.setTime(0.);
  sh.setStrSize(65536);
  sh.setId("Id String");	               // Ignored
  sh.setTitle("Title String");          // Ignored
  sh.setUserTitle("User Title String"); // Ignored
  
  // create the variable array
  {
    std::string var_title;
    
    for (unsigned int var=0; var<var_names.size(); var++)
      {
	for (unsigned int c=0; c<var_names[var].size(); c++)
	  var_title += var_names[var][c];
	
	var_title += '\0';
      }
    
    sh.setVarTitle(var_title.c_str(), var_title.size());
  }
  
  // Put the informationin the XDR file.
  s.header(&sh); // Needs to work for both types of file
  
  // Write the solution vector
  assert (soln.size() == var_names.size()*mesh.n_nodes());
  
  s.values(&soln[0], mesh.n_nodes());
}
