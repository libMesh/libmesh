/* $Id: ex12.C,v 1.4 2003-11-10 22:14:35 benkirk Exp $ */

/* The Next Great Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

 // Example 12 -- The \p MeshData class
 //
 // The previous examples covered the certainly involved
 // aspects of simulating multiple equation systems, and prior 
 // to this even using adaptive refinement.  I.e.: cool stuff.
 //
 // This example now reduces brain effort significantly,
 // and presents some supplements concerning data I/O,
 // especially how to handle the \p MeshData object.
 // The \p MeshData class may be used for handling input
 // data for a simulation, like actual material properties, 
 // (not indicators, like in the \p BoundaryInfo class),
 // or mode shapes that serve as input for acoustic radiation
 // simulations.  The use of the \p MeshData during simulation
 // is straightforward:  the 
 //
 // - \p Number               \p MeshData::operator()(const Node*, int),
 //   get the i-th floating-point value associated with the node
 // - \p bool                 \p MeshData::has_data(const Node*),
 //   verify whether a certain node has data associated
 // - \p std::vector<Number>& \p MeshData::get_data (const Node* node)
 //   to get read-access to the data associated with this node (better
 //   make sure first that this node @e has data, see \p has_data() )
 // - iterator for nodal data \p MeshData::const_node_data_iterator
 //   to directly iterate over the set of nodes that hold data, 
 //   instead of asking through \p has_data() each time again.
 //
 // (and corresponding methods for const Elem*) provide access to
 // the floating-point data associated with nodes/elements.
 // This example does @e not show how to use these aforementioned
 // methods, this is straightforward.  Simply check if the current 
 // Elem* has a node with data. If so, add this contribution to the 
 // RHS, or whatever. Or ask the \p MeshData for each Elem* for its 
 // material properties...
 //
 // Here, only the actual file I/O necessary to handle such 
 // nodal/element data is presented.  Lean back, relax...

// C++ include files that we need.
#include <math.h>

// Functions to initialize the library
// and provide some further features (e.g. 
// our own \p pi)
#include "libmesh.h"

// Basic include files needed for the mesh and 
// \p MeshData functionality.
#include "mesh.h"

// Function prototype for creating artificial nodal data
// that can be inserted into a \p MeshData object.
void create_artificial_data (const Mesh& mesh,
			     std::map<const Node*, std::vector<Number> >& art_data);

// The main program.
int main (int argc, char** argv)
{
  // Initialize the library.
  libMesh::init (argc, argv);

  // Force all our objects to have local scope.  
  {    
    // Check for proper usage. The program is designed to be run
    //  as follows:
    //
    // ./ex12 -d 3 in_mesh.unv
    //
    // where in_mesh.unv should better be a Universal file.
    if (argc < 4)
      {
	std::cerr << "Usage: " << argv[0] << " -d <dim> in_mesh.unv"
		  << std::endl;
	
	error();
      }

    // Get the dimensionality of the mesh from argv[2]
    const unsigned int dim = atoi(argv[2]);
    
    // The filename of the mesh
    const std::string mesh_file = argv[3];

    // The following example makes currently sense
    // only with a Universal file
    if (mesh_file.rfind(".unv") >= mesh_file.size())
      {

	std::cerr << "ERROR:  This example works only properly with a Universal mesh file!"
		  << std::endl;
	
	error();
      }
  
     // Some notes on \p MeshData:
     //
     // - The \p MeshData is @e not a mesh!  Consult the \p Mesh,
     //   \p MeshBase, and \p BoundaryMesh classes for details on
     //   handling meshes!
     //
     // - The \p MeshData is an object handling arbitrary floating-point 
     //   data associated with mesh entities (nodes, elements), currently
     //   most features only available for nodal data.  However,
     //   extending to element-data (does @e anybody need this?)
     //   is straightforward.
     //
     // - Currently std::vector<Number> is the data (associated
     //   with nodes/elements) that can be handled by \p MeshData,
     //
     // - In order to provide @e full functionality, the \p MeshData
     //   @e has to be activated prior to reading in a mesh.  However,
     //   there is a so-called compatibility mode available when you 
     //   (intentionally) forgot to "turn on" the \p MeshData.
     //
     // - It is possible to provide user-defined nodal data that
     //   may subsequently be used or written to file.
     //
     // - Translate the nodal-/element-associated data to formats
     //   suitable for visual inspection, e.g. GMV files.
     //   
     // - Two I/O file formats are currently supported: the Universal 
     //   file format with dataset 2414, and an extension of 
     //   the libMesh-own xda/xdr format, named xtr, xta.
     //   Some details on this:
     //
     //   - The xtr/xta format is simply an extension of the
     //     xdr/xda format to read/write also nodal/element-associated
     //     floating-point data.  The xtr interface of the \p MeshData
     //     creates stand-alone files that may be binary or ASCII.
     //     You cannot append files created by the \p MeshData I/O methods
     //     to a mesh file (xdr/xda).
     //     The xtr files may be seen as a "backdrop", especially when
     //     binary files are preferred.  Note that unv files are @e always
     //     ASCII and may become pretty big!
     //
     //   - The unv format is an extremely versatile text-based file format
     //     for arbitrary data.  Its functionality is @e large, and \p libMesh
     //     supports only a small share: namely the I/O for nodes, elements and
     //     arbitrary node-/element-associated data, stored in the 
     //     so-called datasets "2411", "2412", and "2414", respectively.
     //     Here, only the last dataset is covered.  The two first datasets are
     //     implemented in the Universal support for I/O of meshes.  A single
     //     unv file may hold @e multiple datasets of type "2414".  To
     //     distinguish data, each dataset has a header.  The \p libMesh pendant
     //     to this header is the \p MeshDataUnvHeader class.  When you
     //     read/write unv files using the \p MeshData, you @always
     //     automatically also read/write such headers.
     //
     // Enough babble for now.  Examples. 
    {
      // Create a mesh with the requested dimension.
      // Below we require a 3-dim mesh, therefore assert
      // it.
      assert (dim == 3);
      Mesh mesh(dim);

      // Activate the \p MeshData of the mesh, so that
      // we can remember the node and element ids used
      // in the file.  When we do not activate the \p MeshData,
      // then there is no chance to import node- or element-
      // associated data.
      mesh.data.activate();

      // Now we can safely read the input mesh.  Note that
      // this should better be a .unv or .xdr/xda file, otherwise
      // we cannot load/save associated data.
      mesh.read (mesh_file);
    
      // Print information about the mesh and the data
      // to the screen.  Obviously, there is no data
      // (apart from node & element ids) in it, yet.
      std::cout << std::endl 
		<< "Finished reading the mesh.  MeshData is active but empty:" << std::endl
		<< "---------------------------------------------------------" << std::endl;
      
      mesh.print_info();
      mesh.data.print_info();
    
      // Create some artificial node-associated data and store
      // it in the mesh's \p MeshData.  Use a \p std::map for this. 
      // Let the function create_artificial_data do the work.
      {
	std::map<const Node*, std::vector<Number> > artificial_data;
	
	create_artificial_data (mesh, artificial_data);
	
	// Before we let the map go out of scope, insert it into
	// the \p MeshData.  Note that by default the element-associated
	// data containers are closed, so that the \p MeshData is
	// ready for use.
	mesh.data.insert_node_data(artificial_data);

	// Let \p artificial_data go out of scope
      }
      
       // Print information about the data to the screen.
      std::cout << std::endl 
		<< "After inserting artificial data into the MeshData:" << std::endl
		<< "--------------------------------------------------" << std::endl;
      
      mesh.data.print_info();

      // Write the artificial data into a universal
      // file.  Use a default header for this.
      std::string first_out_data="data_first_out_with_default_header.unv";
      std::cout << "Writing MeshData to: " << first_out_data << std::endl;
      mesh.data.write(first_out_data);

      // Alternatively, create your own header.
      std::cout << std::endl 
		<< "Attach our own MeshDataUnvHeader to the MeshData:" << std::endl
		<< "-------------------------------------------------" << std::endl
		<< " (note the warning: the number of values per node in" << std::endl
		<< "  my_header is not correct)" << std::endl << std::endl;
      MeshDataUnvHeader my_header;

      // Specify an int for this dataset.
      // This is particularly helpful
      // when there are multiple datasets in
      // a file and you want to find a specific one.
      // For this, use MeshDataUnvHeader::which_dataset(int),
      // which is not covered here.
      my_header.dataset_label = 3;

      // Specify some text that helps the @e user 
      // identify the data.  This text is @e not used for
      // finding a specific dataset.  Leave default for
      // the remaining 2 lines.
      my_header.id_lines_1_to_5[0] = "Artificial data";
      my_header.id_lines_1_to_5[1] = "sin curve in z-direction";
      my_header.id_lines_1_to_5[2] = "line in x-direction";
      
      // Also some float data, not associated with nodes,
      // can be stored in the header.
      my_header.record_12[0] = libMesh::pi;
      
      // Now attach this header to the \p MeshData, and write
      // the same file again, but with the personalized header.
      mesh.data.set_unv_header(&my_header);

      // Write again to file.
      std::string second_out_data="data_second_with_header_out.unv";
      std::cout << "Writing MeshData to: " << second_out_data << std::endl;
      mesh.data.write(second_out_data);
      
      // Print information about the data to the screen.
      std::cout << std::endl 
		<< "Before clearing the MeshData:" << std::endl
		<< "-----------------------------" << std::endl;
      mesh.data.print_info();

      // Now clear only the data associated with nodes/elements,
      // but keep the node/element ids used in the mesh.  To
      // clear also these ids, use MeshData::slim() (not used here).
      mesh.data.clear();

      // Print information about the data to the screen.
      std::cout << std::endl 
		<< "After clearing the MeshData:" << std::endl
		<< "----------------------------" << std::endl;
      mesh.data.print_info();

      // Now the \p MeshData is open again to read data from
      // file.  Read the file that we created first.
      mesh.data.read(first_out_data);
      
      std::cout << std::endl 
		<< "After re-reading the first file:" << std::endl
		<< "--------------------------------" << std::endl;
      mesh.data.print_info();

      // Let the mesh, the unv_header etc go out of scope, and
      // do another example.
    }

    std::cout << std::endl 
	      << "----------------------------------------------" << std::endl
	      << "---------- next example with MeshData --------" << std::endl
	      << "----------------------------------------------" << std::endl;

    // Create a new mesh, read it again -- but this time
    // with de-activated \p MeshData.  Then we are
    // able to use the compatibility mode not only for
    // handling \p MeshData dat, but also to @e write 
    // a mesh in unv format.
    // The libMesh-internal node and element ids are used.
    {
      Mesh mesh(dim);

      // Read the input mesh, but with deactivated \p MeshData.
      mesh.read (mesh_file);
    
      // Print information about the mesh and the data
      // to the screen.
      std::cout << std::endl 
		<< "De-activated MeshData:" << std::endl
		<< "----------------------" << std::endl;
      mesh.print_info();
      mesh.data.print_info();
 
      // Write the @e mesh (not the MeshData!) as .unv file.
      // In general, the \p MeshBase interface for .unv I/O
      // needs an active \p MeshData.  However, use compatibility
      // mode to at least write a .unv file, but with the 
      // ids from libMesh.
      const std::string out_mesh = "mesh_with_libmesh_ids.unv";
      std::cout << "Writing _Mesh_ to: " << out_mesh << std::endl
		<< "Try 'diff " << out_mesh << " " << mesh_file << "'" << std::endl
		<< "to see the differences in node numbers." << std::endl
		<< "---------------------------------------" << std::endl
		<< std::endl;
      mesh.write(out_mesh);

      // Again create some artificial node-associated data,
      // as before.
      {
	std::map<const Node*, std::vector<Number> > artificial_data;
	create_artificial_data (mesh, artificial_data);
	mesh.data.insert_node_data(artificial_data);
      }

      // Note that even with (only) compatibility mode MeshData
      // data can be written.  But again, the ids from libMesh
      // are used.  Consult the warning messages issued in
      // DEBUG mode.  And the user @e has to specify that the
      // \p MeshData should change to compatibility mode.
      mesh.data.enable_compatibility_mode();

      // Now that compatibility mode is used, data can be written.
      // _Without_ explicitly enabling compatibility mode, we
      // would get an error message and (with the following write()
      // statement).
      std::string mesh_data_file = "data_third_with_libmesh_ids_out.unv";
      std::cout << std::endl 
		<< "Writing MeshData to: " << mesh_data_file << std::endl
		<< "----------------------------------------------------------" 
		<< std::endl << std::endl;
      mesh.data.write (mesh_data_file);

#ifdef HAVE_ZLIB_H

      // As may already seen, UNV files are text-based, so the may
      // become really big.  When \p ./configure found \p zlib.h,
      // then we may also @e read or @e write \p .unv files in gzip'ed 
      // format! -- Pretty cool, and also pretty fast, due to zlib.h.
      //
      // Note that this works also for mesh files, not only for 
      // meshdata files.
      //
      // In order to write a ".unv.gz" file instead of a ".unv" file,
      // simply provide the full name with ".gz" appended to the
      // write method; it will then figure out whether this file should
      // be gzip'ed or not.
      std::string packed_mesh_data_file =  "packed_" + mesh_data_file + ".gz";
      std::cout << std::endl 
		<< "Writing gzip'ed MeshData to: " << packed_mesh_data_file << std::endl
		<< "---------------------------------------------------------------------------" << std::endl
		<< " To verify the integrity of the packed version, type:" << std::endl << std::endl
		<< "   gunzip " << packed_mesh_data_file << "; " << std::endl
		<< "   diff packed_" << mesh_data_file << " " 
		<< mesh_data_file << std::endl << std::endl;
      
      mesh.data.write (packed_mesh_data_file);
      
#endif

      // And now a last gimmick: The \p MeshData::translate()
      // conveniently converts the nodal- or element-associated
      // data (currently only nodal) to vectors that may be used 
      // for writing a mesh with "solution" vectors, where the solution 
      // vector contains the data from the \p MeshData.  Particularly
      // useful for @e inspecting the data contained in \p MeshData.
      //
      // And even better: the user can choose the mesh for which
      // to export the data.  E.g. not only use the \p mesh
      // itself, but alternatively use the \p BoundaryMesh 
      // (any mesh that uses the @e same nodes, i.e. the \p Node* 
      // have to be the same.  Only exception that will not work:
      // A mesh created using Mesh::create_submesh()
      // actually will @e not work with \p MeshData::translate() ).
      //
      // All in all not bad, hm?
      {
	// have a vector for the actual values and a vector
	// for the names of the data available.
	std::vector<Number> translated_data;
	std::vector<std::string> data_names;
	
	// Use the \p mesh itself.  Alternatively, use the 
	// \p BoundaryMesh of \p mesh.
	mesh.data.translate (mesh,
			     translated_data,
			     data_names);


	// And write the data to a GMV file
	const std::string gmv_file = "data_and_mesh_out.gmv";
	std::cout << std::endl 
		  << "Writing the data from the MeshData to the GMV file " 
		  << gmv_file << std::endl
		  << "------------------------------------------------------------------------" 
		  << std::endl;
	
	mesh.write_gmv_binary (gmv_file,
			       &translated_data,
			       &data_names);
	
	 // Let the vectors with translated data
	 // go out of scope.
      }
      
      // Let the second mesh go out of scope.
    }
  }

   // All done.
  return libMesh::close();
}

// This function creates the data to populate the \p MeshData object
void create_artificial_data (const Mesh& mesh,
			     std::map<const Node*, std::vector<Number> >& art_data)
{
  // get the bounding box to have some sensible data
  std::pair<Point, Point> b_box = mesh.bounding_box();

  const Real z_min = b_box.first (2);
  const Real z_max = b_box.second(2);
  assert (fabs(z_max-z_min) > TOLERANCE);

  const Real x_min = b_box.first (0);
  const Real x_max = b_box.second(0);
  assert (fabs(x_max-x_min) > TOLERANCE);


  const_node_iterator node_it = mesh.nodes_begin();
  const const_node_iterator node_end = mesh.nodes_end();

  for (; node_it != node_end; ++node_it)
    {
      // All the vectors in \p artificial_data @e have to have the
      // same size.  Here we use only two entries per node,
      // but theoretically arbitrary size is possible.
      std::vector<Number> node_data;
      node_data.resize(2);

      // Use a sin curve in z-direction and a linear
      // increase in x-direction
      const Point& p = **node_it;
	
      const Real z_normalized = (p(2)-z_min)/(z_max-z_min);
      const Real x_normalized = (p(0)-x_min)/(x_max-x_min);

      node_data[0] = sin(2*libMesh::pi*z_normalized);
      node_data[1] = x_normalized;
      
      // Insert the current data together with a pointer to
      // the current node in the map.
      art_data.insert (std::make_pair(*node_it,node_data));
    }
}
