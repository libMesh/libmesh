// $Id: mesh_xdr_support.C,v 1.15 2003-09-02 19:50:21 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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



#include "mesh_common.h"

// C++ includes
#include <string.h>
#include <stdio.h>
#ifdef HAVE_RPC_RPC_H
#include <rpc/rpc.h>
#endif

// Local includes
#include "mesh.h"
#include "xdrIO.h"
#include "mesh_xdr_support.h"
#include "elem.h"
#include "cell_hex27.h"


void XdrInterface::mesh_interface(const std::string& name,
				  const XdrIO::XdrIO_TYPE access,
				  std::vector<Node*>& nodes,
				  std::vector<Elem*>& elements,
				  BoundaryInfo& boundary_info,
				  Mesh& mesh)
{
#ifdef DEBUG
  /**
   * Check for unknown access type
   */
  switch (access)
    {      
#ifdef HAVE_RPC_RPC_H
      
    case (XdrIO::ENCODE):
      break;
    case (XdrIO::DECODE):
      break;  
#endif
    case (XdrIO::W_ASCII):
      break;
    case (XdrIO::R_ASCII):
      break;
    default:
      error();
    }
#endif


  /**
   * Create an XdrMESH
   * object.  See XdrIO.h
   * for more information.
   */
  XdrMESH m;

  /**
   * Create a pointer
   * to an XdrMESH file
   * header.
   */
  XdrMHEAD* mh = new XdrMHEAD;

  /**
   * Open the XDR file for
   * reading or writing.
   * Note 1: Provide an additional argument
   * to specify the dimension.
   *
   * Note 2: Has to do the right thing for
   * both binary and ASCII files.
   */
  m.init(access, name.c_str(), 0); // mesh files are always number 0 ...

  /**
   * From here on, things depend
   * on whether we are reading or
   * writing!  First, we define
   * header variables that may
   * be read OR written.
   */
  int                       numElem     = 0;       
  int                       numNodes    = 0;      
  int                       totalWeight = 0;   
  int                       numBCs      = 0;
  unsigned int              n_blocks    = 0;
  std::vector<ElemType>     etypes;
  std::vector<unsigned int> neeb;
	
  switch (access)
    {
    case (XdrIO::R_ASCII):
    case (XdrIO::DECODE):
      {
	/**
	 * Get the information from
	 * the header, and place it
	 * in the header pointer.
	 */
	m.header(mh);
	
	/**
	 * Read information from the
	 * file header.  This depends on
	 * whether its a DEAL or MGF mesh.
	 */
	numElem           = mh->getNumEl();
	numNodes          = mh->getNumNodes();
	totalWeight       = mh->getSumWghts();
	numBCs            = mh->getNumBCs();

	/**
	 * If a deal type mesh, read the augmented mesh information
	 */
	int orig_type = m.get_orig_flag();
	if (orig_type == 0)
	  {
	    // Read augmented header
	    n_blocks = mh->get_n_blocks();

	    etypes.resize(n_blocks);
	    etypes   = mh->get_block_elt_types();

	    neeb.resize(n_blocks);
	    neeb     = mh->get_num_elem_each_block();
	  }
	  
	
	break;
      }
      
    case (XdrIO::W_ASCII):
    case (XdrIO::ENCODE):
      {
	/**
	 * Set the header information.
	 * We assume now that it is only
	 * necessary to write DEAL-style
	 * "augmented" headers.  An
	 * augmented header contains 
	 * information about mesh blocks,
	 * allowing us to optimize storage
	 * and minimize IO requirements
	 * for these meshes.
	 */
	numBCs = boundary_info.n_boundary_conds();
	etypes = mesh.elem_types();
	neeb.resize(etypes.size());
	
	for (unsigned int i=0; i<etypes.size(); i++)
	  neeb[i] = mesh.n_elem_of_type(etypes[i]);

	mh->set_n_blocks(etypes.size());
	mh->set_block_elt_types(etypes);
	mh->set_num_elem_each_block(neeb);	
	mh->setNumEl(mesh.n_elem());
	mh->setNumNodes(mesh.n_nodes());
	mh->setStrSize(65536);
	mh->setSumWghts(mesh.total_weight());
	mh->setNumBCs(numBCs);
	mh->setId("Id String");	   // Ignored
	mh->setTitle("Title String"); // Ignored
	
	/**
	 * Put the information
	 * in the XDR file.
	 */
	m.header(mh); // Needs to work for both types of file
	break;
      }
      
    default:
      {
	// Shouldn't have gotten here.
	error();
      }
    }

  
  /**
   * Deallocate storage for the
   * filename and mesh_base.header.
   */
  delete mh;

  /**
   * Now that we know the
   * number of nodes and elements,
   * we can resize the
   * appropriate vectors if we are
   * reading information in.
   */
  if ((access == XdrIO::DECODE) || (access == XdrIO::R_ASCII))
    {
      nodes.resize(numNodes);
      elements.resize(numElem);
    }

  
  
  /**
   * Read/Write the connectivity
   */

  std::vector<int> conn;
  switch (access)
    {
    case (XdrIO::DECODE):
    case (XdrIO::R_ASCII):
      {
	/**
	 * Each element stores two extra
	 * locations: one which tells
	 * what type of element it is,
	 * and one which tells how many
	 * nodes it has. Therefore,
	 * the total number of nodes
	 * (totalWeight) must be augmented
	 * by 2 times the number of elements
	 * in order to read in the entire
	 * connectivity array.
	 */

	/**
	 * Note: This section now depends on
	 * whether we are reading a DEAL or
	 * MGF style mesh.
	 */
	int orig_flag = m.get_orig_flag();

	if (orig_flag == 0) // DEAL 
	  {
	    conn.resize(totalWeight);
	    m.Icon(&conn[0], 1, totalWeight);
	  }
	
	else if (orig_flag == 1) // MGF
	  {
	    conn.resize(totalWeight+2*numElem);
	    m.Icon(&conn[0], 1, totalWeight+(2*numElem));
	  }
	
	else
	  {
	    // I don't know what type of mesh it is.
	    error();
	  }
	
	break;
      }

    case (XdrIO::ENCODE):
    case (XdrIO::W_ASCII):
      { 
	conn.resize(mesh.total_weight());
	
	unsigned int lastConnIndex = 0;
	unsigned int nn = 0;
	
	for (unsigned int idx=0; idx<etypes.size(); idx++)
	  {
	    nn = lastConnIndex = 0;
	    
	    for (unsigned int e=0; e<mesh.n_elem(); e++)
	      if (mesh.elem(e)->type() == etypes[idx])
		{
		  nn = mesh.elem(e)->n_nodes();
		  
		  for (unsigned int n=0; n<nn; n++)
		    conn[lastConnIndex + n] = mesh.elem(e)->node(n);
		  
		  lastConnIndex += nn;
		}
	    
	    // Send conn to the XDR file
	    m.Icon(&conn[0], nn, lastConnIndex/nn);
	  }
	
	break;
      }

    default:
      {
	// How'd we get here? We have to be either
	// reading or writing.
	error();
      }
    }
    
  /**
   * If we are reading,
   * read in the nodal
   * coordinates and form points.
   * If we are writing, create
   * the vector of coords and send
   * it to the XDR file.
   */
  {
    std::vector<Real> coords;
    
    switch (access)
      {
      case (XdrIO::R_ASCII):
      case (XdrIO::DECODE):
	{
	  coords.resize(numNodes*mesh.spatial_dimension()); // Always use three coords per node
	  m.coord(&coords[0], mesh.spatial_dimension(), numNodes);
	  
	  /**
	   * Form Nodes out of
	   * the coordinates.  If the    
           * MeshData object is active,
	   * add the nodes and ids also          
	   * to its map.
	   */	
	  if (mesh.data.active())   
	    {     
	      for (int innd=0; innd<numNodes; ++innd)      
	        {          
		  nodes[innd] = Node::build(coords[0+innd*3],  
					    coords[1+innd*3],
					    coords[2+innd*3],
					    innd);
                      
		  /*                  
		   * add the id to the MeshData, so that
		   * it knows the foreign id, even when 
		   * the underlying mesh got re-numbered,
		   * refined, elements/nodes added...   
                   */                 
		  mesh.data.add_foreign_node_id(nodes[innd],
						innd); 
		}  
	    }   
          else
	      for (int innd=0; innd<numNodes; ++innd)
		  nodes[innd] = Node::build(coords[0+innd*3],
					    coords[1+innd*3],
					    coords[2+innd*3],
					    innd);
	  
	  break;
	}
	
      case (XdrIO::W_ASCII):
      case (XdrIO::ENCODE):
	{
	  coords.resize(mesh.spatial_dimension()*mesh.n_nodes()); 
	  int lastIndex=0;
	  for (unsigned int i=0; i<mesh.n_nodes(); i++)
	    {
	      const Point& p = *nodes[i];
	      
	      coords[lastIndex+0] = p(0);
	      coords[lastIndex+1] = p(1);
	      coords[lastIndex+2] = p(2);
	      
	      lastIndex += 3;
	    }
	  
	  // Put the nodes in the XDR file
	  m.coord(&coords[0], mesh.spatial_dimension(), mesh.n_nodes()); 
	  break;
	}
	
      default:
	{
	  // How'd we get here? We have to be either
	  // reading or writing.
	  error();
	}
      }
    
    /**
     * Free memory used in
     * the coords vector.
     */
    coords.clear();
  }

  
  /**
   * If we are reading or encoding, build the elements.
   * Note: If the originator was MGF, we don't
   * have to do much checking ...
   * all the elements are Hex27.
   * If the originator was
   * this code, we have to loop over
   * et and neeb to read in all the
   * elements correctly.
   *
   * (This used to be before the coords block, but it
   * had to change now that elements store pointers to
   * nodes.  The nodes must exist before we assign them to
   * the elements. BSK, 1/13/2003)
   */
  {
    if ((access == XdrIO::DECODE) || (access == XdrIO::R_ASCII))
      {
	int orig_type = m.get_orig_flag();
	
	if (orig_type == 0) // DEAL-style (0) hybrid mesh possible
	  {
	    unsigned int lastConnIndex = 0;
	    unsigned int lastFaceIndex = 0;

	    for (unsigned int idx=0; idx<etypes.size(); idx++)  
	      {             
		for (unsigned int e=lastFaceIndex; e<lastFaceIndex+neeb[idx]; e++)  
		  {       
		    elements[e] = Elem::build(etypes[idx]);    
                   
		    /*             
		     * Add elements with the same id as in libMesh.  
		     * Provided the data files that MeshData reads    
		     * were only written with MeshData, then this      
		     * should work properly.  This is an inline
		     * function, so that for disabled MeshData, this
		     * should not induce too much cost
		     */                   
		    mesh.data.add_foreign_elem_id (elements[e],      
						   e);

		    for (unsigned int innd=0; innd < elements[e]->n_nodes(); innd++)
			elements[e]->set_node(innd) = nodes[conn[innd+lastConnIndex]];

		    lastConnIndex += mesh.elem(e)->n_nodes();
		  }
		lastFaceIndex += neeb[idx];
	      }

	  }
  
	else if (orig_type == 1) // MGF-style (1) Hex27 mesh
	  {

#ifdef DEBUG
	    if (mesh.data.active())
	      {
		  std::cerr << "ERROR: MeshData not implemented for MGF-style mesh."
			    << std::endl;
		  error();
	      }
#endif

	    for (int ielm=0; ielm < numElem; ++ielm)
	      {
		elements[ielm] = new Hex27;
		for (int innd=0; innd < 27; ++innd)
		  elements[ielm]->set_node(innd) = nodes[conn[innd+2+(27+2)*ielm]];	
	      }
	  }


	/*
	 * tell the MeshData object that we are finished 
	 * reading data
	 */
	mesh.data.close_foreign_id_maps ();
      }
  
    /**
     * Free memory used in
     * the connectivity
     * vector.
     */
    conn.clear();
  }


  /**
   * If we are reading,
   * read in the BCs
   * from the mesh file,
   * otherwise write the
   * boundary conditions
   * if the BoundaryInfo
   * object exists.
   */
  std::vector<int> bcs;
  bcs.resize(numBCs*3);   
  switch (access)
    {
    case (XdrIO::R_ASCII):
    case (XdrIO::DECODE):
      {
	m.BC(&bcs[0], numBCs);
	
	// Create the boundary_info !!
	for (int ibc=0; ibc < numBCs; ibc++)
	  boundary_info.add_side(bcs[0+ibc*3], bcs[1+ibc*3], bcs[2+ibc*3]);
	  
	break;
      }

    case (XdrIO::W_ASCII):
    case (XdrIO::ENCODE):
      {
	std::cout << "numBCs=" << numBCs << std::endl;
	
	//std::cout << "Preparing to write boundary conditions." << std::endl;
	std::vector<unsigned int> elem_list;
	std::vector<unsigned short int> side_list;
	std::vector<short int> elem_id_list;

	boundary_info.build_side_list (elem_list, side_list, elem_id_list);
	
	for (int ibc=0;  ibc<numBCs; ibc++)
	  {
	    bcs[0+ibc*3] = elem_list[ibc];
	    bcs[1+ibc*3] = side_list[ibc];
	    bcs[2+ibc*3] = elem_id_list[ibc];
	  }
	
	m.BC(&bcs[0], numBCs);
	
	break;
      }
      

    default:
      {
	// How'd we get here? We have to be either
	// reading or writing.
	error();
      }
    }
}



void XdrInterface::soln_interface(const std::string& name,
				  const XdrIO::XdrIO_TYPE access,
				  std::vector<Number>& soln,
				  std::vector<std::string>& var_names,
				  Mesh& mesh)
{

#ifdef USE_COMPLEX_NUMBERS
  // buffer for writing separately
  std::vector<Real> Real_soln;
  std::vector<Real> imag_soln;

  mesh.prepare_complex_data(&soln, &Real_soln, &imag_soln);

  soln_interface_impl(mesh.complex_filename(name, 0), 
		      access, 
		      Real_soln, 
		      var_names, 
		      mesh);

  soln_interface_impl(mesh.complex_filename(name, 1), 
		      access, 
		      imag_soln, 
		      var_names, 
		      mesh);


#else

  soln_interface_impl(name, access, soln, var_names, mesh);

#endif

}




void XdrInterface::soln_interface_impl(const std::string& name,
				       const XdrIO::XdrIO_TYPE access,
				       std::vector<Real>& soln,
				       std::vector<std::string>& var_names,
				       Mesh& mesh)
{
#ifdef DEBUG
  /**
   * Check for unknown access type
   */
  switch (access)
    {      
#ifdef HAVE_RPC_RPC_H      
    case (XdrIO::ENCODE):
      break;
    case (XdrIO::DECODE):
      break;      
#endif     
    case (XdrIO::W_ASCII):
      break;
    case (XdrIO::R_ASCII):
      break;
    default:
      error();
    }
#endif

  /**
   * Create an XdrSOLN
   * object.  See XdrIO.h
   * for more information.
   */
  XdrSOLN s;

  /**
   * Create a pointer
   * to an XdrSOLN file
   * header.
   */
  XdrSHEAD* sh = new XdrSHEAD;

  /**
   * Open the XDR file for
   * reading or writing.
   * Note 1: Provide an additional argument
   * to specify the dimension.
   *
   * Note 2: Has to do the right thing for
   * both binary and ASCII files.
   */
  s.init(access, name.c_str(), 0); // mesh files are always number 0 ...

  /**
   * From here on, things depend
   * on whether we are reading or
   * writing!  First, we define
   * header variables that may
   * be read OR written.
   */
  int         numVar      = 0;       
  int         numNodes    = 0;
  //int         strSize     = 0;
  const char* varNames;
	
  switch (access)
    {
    case (XdrIO::R_ASCII):
    case (XdrIO::DECODE):
      {
	/**
	 * Get the information from
	 * the header, and place it
	 * in the header pointer.
	 */
	s.header(sh);
	
	/**
	 * Read information from the
	 * file header.  This depends on
	 * whether its a DEAL or MGF mesh.
	 */
	numVar   = sh->getWrtVar();
	numNodes = sh->getNumNodes();
	//strSize  = sh->getStrSize();
	varNames = sh->getVarTitle();

	// Get the variable names
	{	  
	  var_names.resize(numVar);
	  
	  const char* p = varNames;
	  
	  for (int i=0; i<numVar; i++)
	    {
	      var_names[i] = p;
	      p += strlen(p) + 1;
	    }
	}
	
	// Read the soln vector
	{
	  soln.resize(numVar*numNodes);

	  s.values(&soln[0], numNodes);
	}	
	
	break;
      }
      
    case (XdrIO::W_ASCII):
    case (XdrIO::ENCODE):
      {
	sh->setWrtVar(var_names.size());
	sh->setNumVar(var_names.size());
	sh->setNumNodes(mesh.n_nodes());
	sh->setNumBCs(mesh.boundary_info.n_boundary_conds());
	sh->setMeshCnt(0);
	sh->setKstep(0);
	sh->setTime(0.);
	sh->setStrSize(65536);
	sh->setId("Id String");	               // Ignored
	sh->setTitle("Title String");          // Ignored
	sh->setUserTitle("User Title String"); // Ignored

	// create the variable array
	{
	  std::string var_title;

	  for (unsigned int var=0; var<var_names.size(); var++)
	    {
	      for (unsigned int c=0; c<var_names[var].size(); c++)
		var_title += var_names[var][c];

	      var_title += '\0';
	    }

	  sh->setVarTitle(var_title.c_str(), var_title.size());
	}
	
	/**
	 * Put the information
	 * in the XDR file.
	 */
	s.header(sh); // Needs to work for both types of file

	// Write the solution vector
	{
	  assert (soln.size() == var_names.size()*mesh.n_nodes());
	  
	  s.values(&soln[0], mesh.n_nodes());
	}
	
	break;
      }
      
    default:
      {
	// Shouldn't have gotten here.
	error();
      }
    }

  
  /**
   * Deallocate storage for the
   * filename and mesh_base.header.
   */
  delete sh;
}








// -------------------------------------------------- 
// Read Methods
// -------------------------------------------------- 

void Mesh::read_xdr(const std::string& name)
{
  /**
   * Clear any existing mesh data
   */
  clear();
  
  /**
   * Instantiate the proper
   * interface.
   */
  XdrInterface interface;
  
  /**
   * Call the appropriate function in the interface
   */
  interface.mesh_interface(name,
			   XdrIO::R_ASCII, // <-- ASCII Read flag
			   _nodes,
			   _elements,
			   boundary_info,
			   *this);
}





void Mesh::read_xdr_binary(const std::string& name)
{
#ifndef HAVE_RPC_RPC_H

  std::cerr << "WARNING: Compiled without XDR binary support." << std::endl
	    << "Will try ASCII instead" << std::endl << std::endl;

  read_xdr(name);

#else
  
  /**
   * Clear any existing mesh data
   */
  clear();
  
  /**
   * Instantiate the proper
   * interface.
   */
  XdrInterface interface;
  
  /**
   * Call the appropriate function in the interface
   */
  interface.mesh_interface(name,
			   XdrIO::DECODE, // <-- Binary Read flag
			   _nodes,
			   _elements,
			   boundary_info,
			   *this);
  
#endif
}
			       
  





void Mesh::read_xdr_soln(const std::string& name,
			 std::vector<Number>& soln,
			 std::vector<std::string>& var_names)
{
  /**
   * Instantiate the proper
   * interface.
   */
  XdrInterface interface;
  
  /**
   * Call the appropriate function in the interface
   */
  interface.soln_interface(name,
			   XdrIO::R_ASCII,
			   soln,
			   var_names,
			   *this);
}




void Mesh::read_xdr_soln_binary(const std::string& name,
				std::vector<Number>& soln,
				std::vector<std::string>& var_names)
{
#ifndef HAVE_RPC_RPC_H

  std::cerr << "WARNING: Compiled without XDR binary support." << std::endl
	    << "Will try ASCII instead" << std::endl << std::endl;

  read_xdr_soln(name, soln, var_names);

#else
  
  /**
   * Instantiate the proper
   * interface.
   */
  XdrInterface interface;
  
  /**
   * Call the appropriate function in the interface
   */
  interface.soln_interface(name,
			   XdrIO::DECODE, // <-- Binary Read flag
			   soln,
			   var_names,
			   *this);
  
#endif
}



// -------------------------------------------------- 
// Write methods
// -------------------------------------------------- 

void Mesh::write_xdr(const std::string& name)
{
  /**
   * Instantiate the proper
   * interface.
   */
  XdrInterface interface;
  
  /**
   * Call the appropriate function in the interface
   */
  interface.mesh_interface(name,
			   XdrIO::W_ASCII, // <-- ASCII Write flag
			   _nodes,
			   _elements,
			   boundary_info,
			   *this);
}





void Mesh::write_xdr_binary(const std::string& name)
{
#ifndef HAVE_RPC_RPC_H

  std::cerr << "WARNING: Compiled without XDR binary support." << std::endl
	    << "Will try ASCII instead" << std::endl << std::endl;

  write_xdr(name);

#else
  
  /**
   * Instantiate the proper
   * interface.
   */
  XdrInterface interface;
  
  /**
   * Call the appropriate function in the interface
   */
  interface.mesh_interface(name,
			   XdrIO::ENCODE, // <-- Binary Write flag
			   _nodes,
			   _elements,
			   boundary_info,
			   *this);

#endif
}




void Mesh::write_xdr_soln(const std::string& name,
			  std::vector<Number>& soln,
			  std::vector<std::string>& var_names)
{
  /**
   * Instantiate the proper
   * interface.
   */
  XdrInterface interface;
  
  /**
   * Call the appropriate function in the interface
   */
  interface.soln_interface(name,
			   XdrIO::W_ASCII,
			   soln,
			   var_names,
			   *this);
}




void Mesh::write_xdr_soln_binary(const std::string& name,
				 std::vector<Number>& soln,
				 std::vector<std::string>& var_names)

{
#ifndef HAVE_RPC_RPC_H

  std::cerr << "WARNING: Compiled without XDR binary support." << std::endl
	    << "Will try ASCII instead" << std::endl << std::endl;

  write_xdr_soln(name, soln, var_names);

#else
  
  /**
   * Instantiate the proper
   * interface.
   */
  XdrInterface interface;
  
  /**
   * Call the appropriate function in the interface
   */
  interface.soln_interface(name,
			   XdrIO::ENCODE, // <-- Binary Write flag
			   soln,
			   var_names,
			   *this);
 
#endif
}



