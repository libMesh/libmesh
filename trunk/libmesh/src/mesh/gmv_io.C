// $Id: gmv_io.C,v 1.16 2004-11-15 00:20:52 benkirk Exp $

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

// Changes: 
// o no more subelements, all elements are written down to GMV directly
// o Some nodes have to be left out, eg node 8 in QUAD9 
// o


// C++ includes
#include <fstream>
#include <string.h> // for strcpy, memcpy
#include <stdio.h>  // for sprintf
#include <vector>

// Local includes
#include "libmesh_config.h"
#include "gmv_io.h"



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
	here();

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
	  eledef.label = "line 3";
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
	  const unsigned int nodes[] = {0,1,2,3,4,5,6};
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
	  eledef.label = "tet10 10";
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
  // Open the output file stream
  std::ofstream out (fname.c_str());
  
  assert (out.good());

  // Get a reference to the mesh
  const MeshBase& mesh = this->cmesh();

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

	// Make sure we have a valid entry for
	// the current element type.
	assert (eletypes.count(elem->type()));

        const elementDefinition& ele = eletypes[elem->type()];
	
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
          << " 0\n";

      for (unsigned int proc=0; proc<mesh.n_partitions(); proc++)
        out << "proc_" << proc << "\n";
      
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      for ( ; it != end; ++it)
        out << (*it)->processor_id()+1 << "\n";
      out << "\n";
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

      out << "variable" << "\n";

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
      
      out << "endvars\n";
    }

  
  // end of the file
  out << "\nendgmv\n";
}



void GMVIO::write_ascii_old_impl (const std::string& fname,
				  const std::vector<Number>* v,
				  const std::vector<std::string>* solution_names)
{
  // Open the output file stream
  std::ofstream out (fname.c_str());
  
  assert (out.good());

  // Get a reference to the mesh
  const MeshBase& mesh = this->cmesh();
  
  // Begin interfacing with the GMV data file
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
    
    out << "cells " << mesh.n_active_sub_elem() << '\n';

//     const_active_elem_iterator       it (mesh.elements_begin());
//     const const_active_elem_iterator end(mesh.elements_end());

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end(); 
    
    switch (mesh.mesh_dimension())
      {
      case 1:
	{
	  // The same temporary storage will be used for each element
	  std::vector<unsigned int> conn;

	  for ( ; it != end; ++it)
	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	      {
		out << "line 2\n";
		(*it)->connectivity(se, TECPLOT, conn);
		for (unsigned int i=0; i<conn.size(); i++)
		  out << conn[i] << " ";
		
		out << '\n';
	      }
	  break;
	}
	
      case 2:
	{
	  // The same temporary storage will be used for each element
	  std::vector<unsigned int> conn;
	  
	  for ( ; it != end; ++it)
	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	      {
		// Quad elements
		if (((*it)->type() == QUAD4) ||
		    ((*it)->type() == QUAD8) ||
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
		  {
		    error();
		  }
		
		out << '\n';
	      }
	  
	  break;
	}
	
	
      case 3:
	{
	  // The same temporary storage will be used for each element
	  std::vector<unsigned int> conn;
	  
	  for ( ; it != end; ++it)
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
			      << "type.  Possibly a dim-1 dimensional "
			      << "element?  Aborting..."
			      << std::endl;
		    error();
		  }
		
		out << '\n';
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
      
//       const_active_elem_iterator       it (mesh.elements_begin());
//       const const_active_elem_iterator end(mesh.elements_end());

      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      for ( ; it != end; ++it)
	for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	    out << (*it)->processor_id()+1 << '\n';
      
      out << '\n';
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

      out << "variable\n";


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
      
      out << "endvars\n";
    }

  
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
  const MeshBase& mesh = this->cmesh();
  
  char buf[80];

  // Begin interfacing with the GMV data file
  {
    // write the nodes
    strcpy(buf, "gmvinput");
    out.write(buf, strlen(buf));
	      
    strcpy(buf, "ieeei4r4");
    out.write(buf, strlen(buf));
  }


  
  // write the nodes
  {
    strcpy(buf, "nodes   ");
    out.write(buf, strlen(buf));

    unsigned int tempint = mesh.n_nodes();
    
    memcpy(buf, &tempint, sizeof(unsigned int));

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
    strcpy(buf, "cells   ");
    out.write(buf, strlen(buf));

    unsigned int tempint = mesh.n_active_elem();
    
    memcpy(buf, &tempint, sizeof(unsigned int));
    
    out.write(buf, sizeof(unsigned int));

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

    switch (mesh.mesh_dimension())
      {

      case 1:
        for ( ; it != end; ++it)
          for(unsigned se = 0; se < (*it)->n_sub_elem(); ++se)
            {
              strcpy(buf, "line    ");
              out.write(buf, strlen(buf));
	      
              tempint = 2;
              memcpy(buf, &tempint, sizeof(unsigned int));
              out.write(buf, sizeof(unsigned int));

              std::vector<unsigned int> conn;
              (*it)->connectivity(se,TECPLOT,conn);
	      
              out.write(reinterpret_cast<char*>(&conn[0]), sizeof(unsigned int)*tempint);
            }
	

        break;
	
      case 2:
       
        for ( ; it != end; ++it)
          for(unsigned se = 0; se < (*it)->n_sub_elem(); ++se)
            {
              strcpy(buf, "quad    ");
              out.write(buf, strlen(buf));
              tempint = 4;
              memcpy(buf, &tempint, sizeof(unsigned int));
              out.write(buf, sizeof(unsigned int));
              std::vector<unsigned int> conn;
              (*it)->connectivity(se,TECPLOT,conn);
              out.write(reinterpret_cast<char*>(&conn[0]), sizeof(unsigned int)*tempint);
            }
        break;
      case 3:
        for ( ; it != end; ++it)
          for(unsigned se = 0; se < (*it)->n_sub_elem(); ++se)
            {
              strcpy(buf, "phex8   ");
              out.write(buf, strlen(buf));
              tempint = 8;
              memcpy(buf, &tempint, sizeof(unsigned int));
              out.write(buf, sizeof(unsigned int));
              std::vector<unsigned int> conn;
              (*it)->connectivity(se,TECPLOT,conn);
              out.write(reinterpret_cast<char*>(&conn[0]), sizeof(unsigned int)*tempint);
            }
        break;
      default:
        error();
	
      }
  }
  
  
  
  // optionally write the partition information
  if (this->partitioning())
    {
      strcpy(buf, "material");
      out.write(buf, strlen(buf));
      
      unsigned int tmpint = mesh.n_processors();
      memcpy(buf, &tmpint, sizeof(unsigned int));
      out.write(buf, sizeof(unsigned int));

      tmpint = 0; // IDs are cell based
      memcpy(buf, &tmpint, sizeof(unsigned int));
      out.write(buf, sizeof(unsigned int));


      for (unsigned int proc=0; proc<mesh.n_processors(); proc++)
        {
          sprintf(buf, "proc_%d", proc);
          out.write(buf, 8);
        }

      std::vector<unsigned int> proc_id (mesh.n_active_elem());
      
      unsigned int n=0;
      
      //       const_active_elem_iterator       it (mesh.elements_begin());
      //       const const_active_elem_iterator end(mesh.elements_end());

      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      for ( ; it != end; ++it)
        for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
          proc_id[n++] = (*it)->processor_id()+1;
      
      
      out.write(reinterpret_cast<char *>(&proc_id[0]),
                sizeof(unsigned int)*proc_id.size());
    }


  
  // optionally write the data
  if ((solution_names != NULL) &&
      (vec != NULL))
    {
      strcpy(buf, "variable");
      out.write(buf, strlen(buf));
      
      float *temp = new float[mesh.n_nodes()];

      const unsigned int n_vars = solution_names->size();
      
      for (unsigned int c=0; c<n_vars; c++)
        {

#ifdef USE_COMPLEX_NUMBERS
          // for complex data, write three datasets


          // Real part
          strcpy(buf, "r_");
          out.write(buf, 2);
          strcpy(buf, (*solution_names)[c].c_str());
          out.write(buf, 6);
	  
          unsigned int tempint = 1; // always do nodal data
          memcpy(buf, &tempint, sizeof(unsigned int));
          out.write(buf, sizeof(unsigned int));
	  
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            temp[n] = static_cast<float>( (*vec)[n*n_vars + c].real() );
	  
          out.write(reinterpret_cast<char *>(temp), sizeof(float)*mesh.n_nodes());


          // imaginary part
          strcpy(buf, "i_");
          out.write(buf, 2);
          strcpy(buf, (*solution_names)[c].c_str());
          out.write(buf, 6);
	  
          memcpy(buf, &tempint, sizeof(unsigned int));
          out.write(buf, sizeof(unsigned int));
	  
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            temp[n] = static_cast<float>( (*vec)[n*n_vars + c].imag() );
	  
          out.write(reinterpret_cast<char *>(temp), sizeof(float)*mesh.n_nodes());

          // magnitude
          strcpy(buf, "a_");
          out.write(buf, 2);
          strcpy(buf, (*solution_names)[c].c_str());
          out.write(buf, 6);
	  
          memcpy(buf, &tempint, sizeof(unsigned int));
          out.write(buf, sizeof(unsigned int));
	  
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            temp[n] = static_cast<float>(std::abs((*vec)[n*n_vars + c]));
	  
          out.write(reinterpret_cast<char *>(temp), sizeof(float)*mesh.n_nodes());

#else


          strcpy(buf, (*solution_names)[c].c_str());
          out.write(buf, 8);
	  
          unsigned int tempint = 1; // always do nodal data
          memcpy(buf, &tempint, sizeof(unsigned int));
          out.write(buf, sizeof(unsigned int));
	  
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            temp[n] = static_cast<float>((*vec)[n*n_vars + c]);
	  
          out.write(reinterpret_cast<char *>(temp), sizeof(float)*mesh.n_nodes());


#endif

	  
        }
    
      delete [] temp;
      
      strcpy(buf, "endvars ");
      out.write(buf, strlen(buf));

    }

  // end the file
  strcpy(buf, "endgmv  ");
  out.write(buf, strlen(buf));
}





void GMVIO::write_discontinuous_gmv (const std::string& name,
				     const EquationSystems& es,
				     const bool write_partitioning) const
{
  std::vector<std::string> solution_names;
  std::vector<Number>      v;

  // Get a reference to the mesh
  const MeshBase& mesh = this->cmesh();
  
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
      //       const_active_elem_iterator       it (mesh.elements_begin());
      //       const const_active_elem_iterator end(mesh.elements_end());

      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      unsigned int tw=0;
      
      for ( ; it != end; ++it)
	tw += (*it)->n_nodes();
      
      out << "nodes " << tw << std::endl;
    }
    


    // Write all the x values
    {
      //       const_active_elem_iterator       it (mesh.elements_begin());
      //       const const_active_elem_iterator end(mesh.elements_end());

      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 
      
      for ( ; it != end; ++it)   
	for (unsigned int n=0; n<(*it)->n_nodes(); n++)
	  out << (*it)->point(n)(0) << " ";

      out << std::endl;
    }
    
    
    // Write all the y values
    {
      //       const_active_elem_iterator       it (mesh.elements_begin());
      //       const const_active_elem_iterator end(mesh.elements_end());

      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      for ( ; it != end; ++it)   
	for (unsigned int n=0; n<(*it)->n_nodes(); n++)
	  out << (*it)->point(n)(1) << " ";

      out << std::endl;
    }
    
    
    // Write all the z values
    {
      //       const_active_elem_iterator       it (mesh.elements_begin());
      //       const const_active_elem_iterator end(mesh.elements_end());

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

    //     const_active_elem_iterator       it (mesh.elements_begin());
    //     const const_active_elem_iterator end(mesh.elements_end());

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

    unsigned int nn=1;
    
    switch (mesh.mesh_dimension())
      {
      case 1:
	{
	  error();
	  break;
	}
	
      case 2:
	{
	  for ( ; it != end; ++it)
	    for (unsigned int se=0; se<(*it)->n_sub_elem(); se++)
	      {
		if (((*it)->type() == QUAD4) ||
		    ((*it)->type() == QUAD8) ||
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
	  error();
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
      
      //       const_active_elem_iterator       it (mesh.elements_begin());
      //       const const_active_elem_iterator end(mesh.elements_end());

      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      for ( ; it != end; ++it)
	out << (*it)->processor_id()+1 << std::endl;
      
      out << std::endl;
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
          // 	    const_active_elem_iterator       it (mesh.elements_begin());
          // 	    const const_active_elem_iterator end(mesh.elements_end());

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
          // 	    const_active_elem_iterator       it (mesh.elements_begin());
          // 	    const const_active_elem_iterator end(mesh.elements_end());

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
          // 	    const_active_elem_iterator       it (mesh.elements_begin());
          // 	    const const_active_elem_iterator end(mesh.elements_end());

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
          // 	    const_active_elem_iterator       it (mesh.elements_begin());
          // 	    const const_active_elem_iterator end(mesh.elements_end());

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

