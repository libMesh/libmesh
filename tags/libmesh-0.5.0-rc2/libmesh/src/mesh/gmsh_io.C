// $Id: gmsh_io.C,v 1.13 2005-06-06 14:53:19 jwpeterson Exp $

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

// This file was massively overhauled and extended by Martin Lüthi, mluthi@tnoo.net

// C++ includes
#include <fstream>

// Local includes
#include "libmesh_config.h"
#include "gmsh_io.h"
#include "elem.h"
#include "mesh_base.h"

// anonymous namespace to hold local data
namespace
{
  /**
   * Defines mapping from libMesh element types to Gmsh element types.
   */
  struct elementDefinition {
    std::string label;
    std::vector<unsigned int> nodes;
    ElemType type;
    unsigned int exptype;
  };


  // maps from a libMesh element type to the proper
  // Gmsh elementDefinition.  Placing the data structure
  // here in this anonymous namespace gives us the
  // benefits of a global variable without the nasty
  // side-effects
  std::map<ElemType, elementDefinition> eletypes_exp;
  std::map<unsigned int, elementDefinition> eletypes_imp;



  // ------------------------------------------------------------
  // helper function to initialize the eletypes map
  void init_eletypes ()
  {
    if (eletypes_exp.empty() && eletypes_imp.empty())
      {
	// This should happen only once.  The first time this method
	// is called the eletypes data struture will be empty, and
	// we will fill it.  Any subsequent calls will find an initialized
	// eletypes map and will do nothing.

	//==============================
	// setup the element definitions
	elementDefinition eledef;

	// use "swap trick" from Scott Meyer's "Effective STL" to initialize
	// eledef.nodes vector
	
	// POINT (only Gmsh)
	{
	  eledef.exptype = 15;
	  eledef.nodes.clear();
	
          // import only
          eletypes_imp[15] = eledef;
	}

	// EDGE2
	{
	  eledef.type    = EDGE2;
	  eledef.exptype = 1;
	  eledef.nodes.clear();
	
	  eletypes_exp[EDGE2] = eledef;
          eletypes_imp[1]     = eledef;
	}
  
	// EDGE3
	{
	  eledef.type    = EDGE3;
	  eledef.exptype = 8;
	  eledef.nodes.clear();
	  
	  eletypes_exp[EDGE3] = eledef;
          eletypes_imp[8]     = eledef;
	}
      
	// TRI3
	{
          eledef.type    = TRI3;
	  eledef.exptype = 2;
	  eledef.nodes.clear();
	  
	  eletypes_exp[TRI3] = eledef;
          eletypes_imp[2] = eledef;
	}
      
	// TRI6
	{
          eledef.type    = TRI6;
	  eledef.exptype = 9;
	  eledef.nodes.clear();

	  eletypes_exp[TRI6] = eledef;
          eletypes_imp[9]    = eledef;
	}
      
	// QUAD4
	{
          eledef.type    = QUAD4;
	  eledef.exptype = 3;
	  eledef.nodes.clear();

	  eletypes_exp[QUAD4] = eledef;
          eletypes_imp[3]     = eledef;
	}
      
	// QUAD8
        // TODO: what should be done with this on writing?
	{
          eledef.type    = QUAD8;
	  eledef.exptype = 100;
	  const unsigned int nodes[] = {1,2,3,4,5,6,7,8};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);

	  eletypes_exp[QUAD8] = eledef;
          eletypes_imp[10]    = eledef;
	}
      
	// QUAD9
	{
          eledef.type    = QUAD9;
	  eledef.exptype = 10;
	  eledef.nodes.clear();

	  eletypes_exp[QUAD9] = eledef;
          eletypes_imp[10]    = eledef;
	}
      
	// HEX8
	{
          eledef.type    = HEX8;
	  eledef.exptype = 5;
	  eledef.nodes.clear();

	  eletypes_exp[HEX8] = eledef;
          eletypes_imp[5]    = eledef;
	}
      
	// HEX20
        // TODO: what should be done with this on writing?
	{
          eledef.type    = HEX20;
	  eledef.exptype = 101;
	  const unsigned int nodes[] = {1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,16};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);

	  eletypes_exp[HEX20] = eledef;
          eletypes_imp[12]    = eledef;
	}
      
	// HEX27
	{
          eledef.type    = HEX27;
	  eledef.exptype = 12;
          // export
          const unsigned int nodes_exp[] = {0,1,2,3,4,5,6,7,8,11,12,9,13,10,14,15,16,19,17,18,20,21,24,22,23,25,26};
          unsigned int nnodes = sizeof(nodes_exp)/sizeof(nodes_exp[0]);
	  std::vector<unsigned int>(nodes_exp, nodes_exp+nnodes).swap(eledef.nodes);
	  eletypes_exp[HEX27] = eledef;

          // import 
          const unsigned int nodes_imp[] = {1,2,3,4,5,6,7,8,9,12,14,10,11,13,15,16,17,19,20,18,21,22,24,25,23,26,27};
          std::vector<unsigned int>(nodes_imp, nodes_imp+nnodes).swap(eledef.nodes);
          eletypes_imp[12]    = eledef;
	}
      
	// TET4
	{
          eledef.type    = TET4;
	  eledef.exptype = 4;
	  eledef.nodes.clear();

          eletypes_exp[TET4] = eledef;
          eletypes_imp[4]    = eledef;
	}
      
	// TET10
	{
          eledef.type    = TET10;
	  eledef.exptype = 11;
          // export
          const unsigned int nodes_exp[] = {0,1,2,3,4,5,6,7,9,8};
          unsigned int nnodes = sizeof(nodes_exp)/sizeof(nodes_exp[0]);
	  std::vector<unsigned int>(nodes_exp, nodes_exp+nnodes).swap(eledef.nodes);
	  eletypes_exp[TET10] = eledef;

          // import 
          const unsigned int nodes_imp[] = {1,2,3,4,5,6,7,8,10,9};
          std::vector<unsigned int>(nodes_imp, nodes_imp+nnodes).swap(eledef.nodes);
          eletypes_imp[11]    = eledef;
	}
      
	// PRISM6
	{
          eledef.type    = PRISM6;
	  eledef.exptype = 6;
	  eledef.nodes.clear();
	  
	  eletypes_exp[PRISM6] = eledef;
          eletypes_imp[6]      = eledef;
	}
      
	// PRISM15
        // TODO: what should be done with this on writing?
	{
          eledef.type    = PRISM15;
	  eledef.exptype = 103;
	  eledef.nodes.clear();
	  
	  eletypes_exp[PRISM15] = eledef;
          eletypes_imp[13] = eledef;
	}

	// PRISM18
	{
          eledef.type    = PRISM18;
	  eledef.exptype = 13;
          // export
          const unsigned int nodes_exp[] = {0,1,2,3,4,5,6,8,9,7,10,11,12,14,13,15,17,16};
          unsigned int nnodes = sizeof(nodes_exp)/sizeof(nodes_exp[0]);
	  std::vector<unsigned int>(nodes_exp, nodes_exp+nnodes).swap(eledef.nodes);
	  eletypes_exp[PRISM18] = eledef;

          // import 
          const unsigned int nodes_imp[] = {1,2,3,4,5,6,7,10,8,9,11,12,13,15,14,16,18,17};
          std::vector<unsigned int>(nodes_imp, nodes_imp+nnodes).swap(eledef.nodes);
          eletypes_imp[13]      = eledef;
	}

	// PYRAMID5
	{
          eledef.type    = PYRAMID5;
	  eledef.exptype = 7;
	  eledef.nodes.clear();
	  
	  eletypes_exp[PYRAMID5] = eledef;
          eletypes_imp[7]        = eledef;
	}

	//==============================      
      }
  }
  
} // end anonymous namespace


// ------------------------------------------------------------
// GmshIO  members
void GmshIO::read (const std::string& name)
{
  std::ifstream in (name.c_str());
  this->read_mesh (in);
}


void GmshIO::read_mesh(std::istream& in)
{

  assert(in.good());

  // initialize the map with element types
  init_eletypes();

  // clear any data in the mesh
  MeshBase& mesh = MeshInput<MeshBase>::mesh();
  mesh.clear();

  // some variables
  const int  bufLen = 256;
  char       buf[bufLen+1];
  int        format=0, size=0;
  Real       version = 1.0;
  
  // map to hold the node numbers for translation
  // note the the nodes can be non-consecutive
  std::map<unsigned int, unsigned int> nodetrans;
  
  {
    while (!in.eof()) {
      in >> buf;

      if (!strncmp(buf,"$MeshFormat",11))
        {
          in >> version >> format >> size;
          if(version != 2.0){
            std::cerr << "Error: Wrong msh file version " << version << "\n";
            error();
          }
          if(format){
            std::cerr << "Error: Unknown data format for mesh\n";
            error();
          }
        }

      // read the node block
      else if (!strncmp(buf,"$NOD",4) ||
          !strncmp(buf,"$NOE",4) ||
          !strncmp(buf,"$Nodes",6) 
          )
        {
          unsigned int numNodes = 0;
          in >> numNodes;
          mesh.reserve_nodes (numNodes);

          // read in the nodal coordinates and form points.
          Real x, y, z;
          unsigned int id;
        
          // add the nodal coordinates to the mesh
          for (unsigned int i=0; i<numNodes; ++i)      
            {
              in >> id >> x >> y >> z;
              mesh.add_point (Point(x, y, z));
              nodetrans[id] = i;
            }
          // read the $ENDNOD delimiter
          in >> buf;
        }

      // read the element block
      else if (!strncmp(buf,"$ELM",4) ||
               !strncmp(buf,"$Elements",9)
               )
        {
          unsigned int numElem = 0;
          in >> numElem;
          mesh.reserve_elem (numElem);

          // read the elements
          for (unsigned int iel=0; iel<numElem; ++iel)      
            {
              unsigned int id, type, physical, elementary, partition = 1, numNodes, numTags;
              if(version <= 1.0)
                {
                  in >> id >> type >> physical >> elementary >> numNodes;
                }
              else
                {
                  in >> id >> type >> numTags;
                  elementary = physical = partition = 1;
                  for(unsigned int j = 0; j < numTags; j++)
                    {
                      int tag;
                      in >> tag;
                      if(j == 0)
                        physical = tag;
                      else if(j == 1)
                        elementary = tag;
                      else if(j == 2)
                        partition = tag;
                      // ignore any other tags for now
                    }
                }

              // consult the import element table which element to build
              const elementDefinition& eletype = eletypes_imp[type];
    
              // add the elements to the mesh
              Elem* elem = mesh.add_elem(Elem::build(eletype.type).release());

              // check number of nodes. We cannot do that for version 2.0
              if (version <= 1.0) 
                {
                  if (elem->n_nodes() != numNodes)
                    {
                      std::cerr << "Number of nodes for element " << id
                                << " of type " << eletypes_imp[type].type
                                << " (Gmsh type " << type  
                                << ") does not match Libmesh definition. "
                                << "I expected " << elem->n_nodes()
                                << " nodes, but got " << numNodes << "\n";
                      error();
                    }
                }

              // add node pointers to the elements
              int nod = 0;
              // if there is a node translation table, use it
              if (eletype.nodes.size() > 0)
                for (unsigned int i=0; i<numNodes; i++)
                  {
                    in >> nod;
                    elem->set_node(i) = mesh.node_ptr(nodetrans[eletype.nodes[nod]]);
                  }
              else
                {
                  for (unsigned int i=0; i<numNodes; i++)
                    {
                      in >> nod;
                      elem->set_node(i) = mesh.node_ptr(nodetrans[nod]);
                    }
                }
            }
          // read the $ENDELM delimiter
          in >> buf;
        }
    
    }

  }
  //mesh.reserve_elem  (numElem);

}



void GmshIO::write (const std::string& name)
{
  if (libMesh::processor_id() == 0)
    {
      std::ofstream out (name.c_str());
      this->write_mesh (out);
    }
}


void GmshIO::write_nodal_data (const std::string& fname,
                               const std::vector<Number>& soln,
                               const std::vector<std::string>& names)
{
  //this->_binary = true;
  if (libMesh::processor_id() == 0)
    this->write_post  (fname, &soln, &names);
}


void GmshIO::write_mesh (std::ostream& out)
{
  // Be sure that the stream is valid.
  assert (out.good());
  
  // initialize the map with element types
  init_eletypes();

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
    out << mesh.n_active_elem() << '\n';

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end(); 
    
    // the element number
    for ( ; it != end; ++it)
      {
        const Elem* elem = *it;

	// Make sure we have a valid entry for
	// the current element type.
	assert (eletypes_exp.count(elem->type()));

        // consult the import element table 
        const elementDefinition& eletype = eletypes_exp[elem->type()];

	// The element mapper better not require any more nodes
	// than are present in the current element!
	assert (eletype.nodes.size() <= elem->n_nodes());
	
        // elements ids are 1 based in Gmsh
        out << elem->id()+1 << " ";

        // element type
        out << eletype.exptype;

        // write the number of tags and
        // tag1 (physical entity), and tag2 (geometric entity)
        out << " 3 1 1 ";

        // write the partition the element belongs to
        out << elem->processor_id()+1 << " ";

        // if there is a node translation table, use it
        if (eletype.nodes.size() > 0)
          for (unsigned int i=0; i < elem->n_nodes(); i++)
            out << elem->node(eletype.nodes[i])+1 << " ";   // gmsh is 1-based 
        // otherwise keep the same node order
        else
          for (unsigned int i=0; i < elem->n_nodes(); i++)
            out << elem->node(i)+1 << " ";                  // gmsh is 1-based 
        out << "\n";
      } // element loop
    out << "$EndElements\n";
  }
}


void GmshIO::write_post (const std::string& fname,
                         const std::vector<Number>* v,
                         const std::vector<std::string>* solution_names)
{

  // Should only do this on processor 0!
  assert (libMesh::processor_id() == 0);
  
  // Create an output stream
  std::ofstream out(fname.c_str());

  // initialize the map with element types
  init_eletypes();

  if (!out.good())
    {
      std::cerr << "ERROR: opening output file " << fname
		<< std::endl;
      error();
    }

  // create a character buffer
  char buf[80];

  // Get a constant reference to the mesh.
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

  //  write the data
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

      // write the header
      out << "$PostFormat\n";
      if (this->binary())
        out << "1.2 1 " << sizeof(double) << "\n";
      else
        out << "1.2 0 " << sizeof(double) << "\n";
      out << "$EndPostFormat\n";

      // Loop over the elements to see how much of each type there are
      unsigned int n_points=0, n_lines=0, n_triangles=0, n_quadrangles=0,
        n_tetrahedra=0, n_hexahedra=0, n_prisms=0, n_pyramids=0;
      unsigned int n_scalar=0, n_vector=0, n_tensor=0;
      unsigned int nb_text2d=0, nb_text2d_chars=0, nb_text3d=0, nb_text3d_chars=0;

      {
        MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
        const MeshBase::const_element_iterator end = mesh.active_elements_end(); 


        for ( ; it != end; ++it)
          {
            const ElemType elemtype = (*it)->type();
            
            switch (elemtype)
              {
              case EDGE2:
              case EDGE3:
              case EDGE4:
                {
                  n_lines += 1;
                  break;
                }
              case TRI3:
              case TRI6:
                {
                  n_triangles += 1;
                  break;
                }
              case QUAD4:
              case QUAD8:
              case QUAD9:
                {
                  n_quadrangles += 1;
                  break;
                }
              case TET4:
              case TET10:
                {
                  n_tetrahedra += 1;
                  break;
                }
              case HEX8:
              case HEX20:
              case HEX27:
                {
                  n_hexahedra += 1;
                  break;
                }
              case PRISM6:
              case PRISM15:
              case PRISM18:
                {
                  n_prisms += 1;
                  break;
                }
              case PYRAMID5:
                {
                  n_pyramids += 1;
                  break;
                }
              default:
                {
                  std::cerr << "ERROR: Not existant element type "
                            << (*it)->type() << std::endl;
                  error();
                }
              }
          }
      }

      // create a view for each variable
      for (unsigned int ivar=0; ivar < n_vars; ivar++)
        {
          std::string varname = (*solution_names)[ivar];

          // at the moment, we just write out scalar quantities
          // later this should be made configurable through
          // options to the writer class
          n_scalar = 1;
      
          // write the variable as a view, and the number of time steps
          out << "$View\n" << varname << " " << 1 << "\n";
      
          // write how many of each geometry type are written
          out << n_points * n_scalar << " "
              << n_points * n_vector << " "
              << n_points * n_tensor << " "
              << n_lines * n_scalar << " "
              << n_lines * n_vector << " "
              << n_lines * n_tensor << " "
              << n_triangles * n_scalar << " "
              << n_triangles * n_vector << " "
              << n_triangles * n_tensor << " "
              << n_quadrangles * n_scalar << " "
              << n_quadrangles * n_vector << " "
              << n_quadrangles * n_tensor << " "
              << n_tetrahedra * n_scalar << " "
              << n_tetrahedra * n_vector << " "
              << n_tetrahedra * n_tensor << " "
              << n_hexahedra * n_scalar << " "
              << n_hexahedra * n_vector << " "
              << n_hexahedra * n_tensor << " "
              << n_prisms * n_scalar << " "
              << n_prisms * n_vector << " "
              << n_prisms * n_tensor << " "
              << n_pyramids * n_scalar << " "
              << n_pyramids * n_vector << " "
              << n_pyramids * n_tensor << " "
              << nb_text2d << " "
              << nb_text2d_chars << " "
              << nb_text3d << " "
              << nb_text3d_chars << "\n";
      
          // if binary, write a marker to identify the endianness of the file
          if (this->binary())
            {
              const int one = 1;
              memcpy(buf, &one, sizeof(int));
              out.write(buf, sizeof(int));
            }

          // the time steps (there is just 1 at the moment)
          if (this->binary())
            {
              double one = 1;
              memcpy(buf, &one, sizeof(double));
              out.write(buf, sizeof(double));
            }
          else
            out << "1\n";

          // Loop over the elements and write out the data
          MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = mesh.active_elements_end(); 
    
          for ( ; it != end; ++it)
            {
              const Elem* elem = *it;
	
              // this is quite crappy, but I did not invent that file format!
              for (unsigned int d=0; d<3; d++)  // loop over the dimensions 
                {
                  for (unsigned int n=0; n < elem->n_vertices(); n++)   // loop over vertices
                    {
                      const Point vertex = elem->point(n);
                      if (this->binary())
                        {
                          double tmp = vertex(d);
                          memcpy(buf, &tmp, sizeof(double));
                          out.write(reinterpret_cast<char *>(buf), sizeof(double));
                        }
                      else
                        out << vertex(d) << " ";
                    }
                  if (!this->binary())
                    out << "\n";
                }

              // now finally write out the data
              for (unsigned int i=0; i < elem->n_vertices(); i++)   // loop over vertices
                if (this->binary())
                  {
#ifndef USE_COMPLEX_NUMBERS
                    double tmp = (*v)[elem->node(i)*n_vars + ivar];
#else
		    std::cout << "WARNING: Gmsh::write_post does not fully support "
			      << "complex numbers. Will only write the real part of "
			      << "variable " << varname << std::endl;

		    double tmp = (*v)[elem->node(i)*n_vars + ivar].real();
#endif
                    memcpy(buf, &tmp, sizeof(double));
                    out.write(reinterpret_cast<char *>(buf), sizeof(double));
                  }
                else
		  {
#ifndef USE_COMPLEX_NUMBERS
		    out << (*v)[elem->node(i)*n_vars + ivar] << "\n";
#else
		    std::cout << "WARNING: Gmsh::write_post does not fully support "
			      << "complex numbers. Will only write the real part of "
			      << "variable " << varname << std::endl;

		    out << (*v)[elem->node(i)*n_vars + ivar].real() << "\n";
#endif
		  }
            }
          if (this->binary())
            out << "\n";
          out << "$EndView\n";

        } // end variable loop (writing the views)
    }

}

