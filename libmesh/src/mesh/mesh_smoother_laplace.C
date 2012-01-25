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
#include <algorithm> // for std::copy, std::sort


// Local includes
#include "mesh_smoother_laplace.h"
#include "mesh_tools.h"
#include "elem.h"
#include "unstructured_mesh.h"
#include "parallel.h"
#include "parallel_ghost_sync.h" // sync_dofobject_data_by_id()
#include "parallel_algebra.h" // StandardType<Point>

namespace libMesh
{

  // Anonymous namespace to hold functor that works with Roy's
  // Parallel::sync_dofobject_data_by_id() utility.
  namespace {

    // This struct must be created and passed to the
    // Parallel::sync_dofobject_data_by_id() function.  The member
    // functions declared here are defined in an anonymous namespace
    // at the end of this file to make it easier to read.
    struct SyncNodalPositions
    {
      // The constructor.  You need a reference to the mesh where you will
      // be setting/getting nodal positions.
      SyncNodalPositions(MeshBase& m);

      // The datum typedef is required of this functor, so that the
      // Parallel::sync_dofobject_data_by_id() function can create e.g.
      // std::vector<datum>.
      typedef Point datum;

      // First required interface.  This function must fill up the data vector for the
      // ids specified in the ids vector.
      void gather_data (const std::vector<unsigned int>& ids, std::vector<datum>& data);

      // Second required interface.  This function must do something with the data in
      // the data vector for the ids in the ids vector.
      void act_on_data (const std::vector<unsigned int>& ids, std::vector<datum>& data);

      // Data to be used by the SyncNodalPositions struct
      MeshBase &mesh;
    };

  } // anonymous namespace





  // LaplaceMeshSmoother member functions
  LaplaceMeshSmoother::LaplaceMeshSmoother(UnstructuredMesh& mesh)
    : MeshSmoother(mesh),
      _initialized(false)
  {
  }




void LaplaceMeshSmoother::smooth(unsigned int n_iterations)
{
  if (!_initialized)
    this->init();

  // Don't smooth the nodes on the boundary...
  // this would change the mesh geometry which
  // is probably not something we want!
  std::vector<bool> on_boundary;
  MeshTools::find_boundary_nodes(_mesh, on_boundary);

  // Ensure that the find_boundary_nodes() function returned a properly-sized vector
  if (on_boundary.size() != _mesh.n_nodes())
    {
      libMesh::err << "MeshTools::find_boundary_nodes() returned incorrect length vector!" << std::endl;
      libmesh_error();
    }

  // We can only update the nodes after all new positions were
  // determined. We store the new positions here
  std::vector<Point> new_positions;

  for (unsigned int n=0; n<n_iterations; n++)
    {
      new_positions.resize(_mesh.n_nodes());

      {
	MeshBase::node_iterator       it     = _mesh.local_nodes_begin();
	const MeshBase::node_iterator it_end = _mesh.local_nodes_end();
	for (; it != it_end; ++it)
	  {
	    Node* node = *it;

	    if (node == NULL)
	      {
		libMesh::err << "[" << libMesh::processor_id() << "]: Node iterator returned NULL pointer." << std::endl;
		libmesh_error();
	      }

	    // leave the boundary intact
	    // Only relocate the nodes which are vertices of an element
	    // All other entries of _graph (the secondary nodes) are empty
	    if (!on_boundary[node->id()] && (_graph[node->id()].size() > 0))
	      {
		Point avg_position(0.,0.,0.);

		for (unsigned j=0; j<_graph[node->id()].size(); ++j)
		  {
		    // Will these nodal positions always be available
		    // or will they refer to remote nodes?  To be
		    // careful, we grab a pointer and test it against
		    // NULL.
		    Node* connected_node = _mesh.node_ptr(_graph[node->id()][j]);

		    if (connected_node == NULL)
		      {
			libMesh::err << "Error! Libmesh returned NULL pointer for node " << _graph[connected_node->id()][j] << std::endl;
			libmesh_error();
		      }

		    avg_position.add( *connected_node );
		  } // end for(j)

		// Compute the average, store in the new_positions vector
		new_positions[node->id()] = avg_position / static_cast<Real>(_graph[node->id()].size());
	      } // end if
	  } // end for
      } // end scope


      // now update the node positions (local node positions only)
      {
	MeshBase::node_iterator it           = _mesh.local_nodes_begin();
	const MeshBase::node_iterator it_end = _mesh.local_nodes_end();
	for (; it != it_end; ++it)
	  {
	    Node* node = *it;

	    if (!on_boundary[node->id()] && (_graph[node->id()].size() > 0))
	      {
		// Should call Point::op=
		// libMesh::out << "Setting node id " << node->id() << " to position " << new_positions[node->id()];
		_mesh.node(node->id()) = new_positions[node->id()];
	      }
	  } // end for
      } // end scope

      // Now the nodes which are ghosts on this processor may have been moved on
      // the processors which own them.  So we need to synchronize with our neighbors
      // and get the most up-to-date positions for the ghosts.
      SyncNodalPositions sync_object(_mesh);
      Parallel::sync_dofobject_data_by_id(_mesh.nodes_begin(), _mesh.nodes_end(), sync_object);

    } // end for n_iterations

  // finally adjust the second order nodes (those located between vertices)
  // these nodes will be located between their adjacent nodes
  // do this element-wise
  MeshBase::element_iterator       el  = _mesh.active_elements_begin();
  const MeshBase::element_iterator end = _mesh.active_elements_end();

  for (; el != end; ++el)
    {
      // Constant handle for the element
      const Elem* elem = *el;

      // get the second order nodes (son)
      // their element indices start at n_vertices and go to n_nodes
      const unsigned int son_begin = elem->n_vertices();
      const unsigned int son_end   = elem->n_nodes();

      // loop over all second order nodes (son)
      for (unsigned int son=son_begin; son<son_end; son++)
        {
	  // Don't smooth second-order nodes which are on the boundary
	  if (!on_boundary[elem->node(son)])
	    {
	      const unsigned int n_adjacent_vertices =
		elem->n_second_order_adjacent_vertices(son);

	      // calculate the new position which is the average of the
	      // position of the adjacent vertices
	      Point avg_position(0,0,0);
	      for (unsigned int v=0; v<n_adjacent_vertices; v++)
		avg_position +=
		  _mesh.point( elem->node( elem->second_order_adjacent_vertex(son,v) ) );

	      _mesh.node(elem->node(son)) = avg_position / n_adjacent_vertices;
	    }
        }
    }
}





void LaplaceMeshSmoother::init()
{
  switch (_mesh.mesh_dimension())
    {

      // TODO:[BSK] Fix this to work for refined meshes...  I think
      // the implementation was done quickly for Damien, who did not have
      // refined grids.  Fix it here and in the original Mesh member.

    case 2: // Stolen directly from build_L_graph in mesh_base.C
      {
	// Initialize space in the graph.  It is n_nodes long.  Each
	// node may be connected to an arbitrary number of other nodes
	// via edges.
	_graph.resize(_mesh.n_nodes());

	MeshBase::element_iterator       el  = _mesh.active_local_elements_begin();
	const MeshBase::element_iterator end = _mesh.active_local_elements_end();

	for (; el != end; ++el)
	  {
	    // Constant handle for the element
	    const Elem* elem = *el;

	    for (unsigned int s=0; s<elem->n_neighbors(); s++)
	      {
		// Only operate on sides which are on the
		// boundary or for which the current element's
		// id is greater than its neighbor's.
                // Sides get only built once.
		if ((elem->neighbor(s) == NULL) ||
		    (elem->id() > elem->neighbor(s)->id()))
		  {
		    AutoPtr<Elem> side(elem->build_side(s));
		    _graph[side->node(0)].push_back(side->node(1));
		    _graph[side->node(1)].push_back(side->node(0));
		  }
	      }
	  }
	_initialized = true;
	break;
      } // case 2

    case 3: // Stolen blatantly from build_L_graph in mesh_base.C
      {
	// Initialize space in the graph.
	_graph.resize(_mesh.n_nodes());

	MeshBase::element_iterator       el  = _mesh.active_local_elements_begin();
	const MeshBase::element_iterator end = _mesh.active_local_elements_end();

	for (; el != end; ++el)
	  {
	    // Shortcut notation for simplicity
	    const Elem* elem = *el;

	    for (unsigned int f=0; f<elem->n_neighbors(); f++) // Loop over faces
	      if ((elem->neighbor(f) == NULL) ||
		  (elem->id() > elem->neighbor(f)->id()))
		{
		  // We need a full (i.e. non-proxy) element for the face, since we will
		  // be looking at its sides as well!
		  AutoPtr<Elem> face = elem->build_side(f, /*proxy=*/false);

		  for (unsigned int s=0; s<face->n_neighbors(); s++) // Loop over face's edges
		    {
		      // Here we can use a proxy
		      AutoPtr<Elem> side = face->build_side(s);

		      // At this point, we just insert the node numbers
		      // again.  At the end we'll call sort and unique
		      // to make sure there are no duplicates
		      _graph[side->node(0)].push_back(side->node(1));
		      _graph[side->node(1)].push_back(side->node(0));
		    }
		}
	  }

	_initialized = true;
	break;
      } // case 3

    default:
      {
	libMesh::err << "At this time it is not possible "
		      << "to smooth a dimension "
		      << _mesh.mesh_dimension()
		      << "mesh.  Aborting..."
		      << std::endl;
	libmesh_error();
      }
    }

  // Done building graph from local elements.  Let's now allgather the
  // graph so that it is available on all processors for the actual
  // smoothing operation?
  this->allgather_graph();

  // In 3D, it's possible for > 2 processor partitions to meet
  // at a single edge, while in 2D only 2 processor partitions
  // share an edge.  Therefore the allgather'd graph in 3D may
  // now have duplicate entries and we need to remove them so
  // they don't foul up the averaging algorithm employed by the
  // Laplace smoother.
  for (unsigned i=0; i<_graph.size(); ++i)
    {
      // The std::unique algorithm removes duplicate *consecutive* elements from a range,
      // so it only makes sense to call it on a sorted range...
      std::sort(_graph[i].begin(), _graph[i].end());
      _graph[i].erase(std::unique(_graph[i].begin(), _graph[i].end()), _graph[i].end());
    }

} // init()




void LaplaceMeshSmoother::print_graph(std::ostream& out) const
{
  for (unsigned int i=0; i<_graph.size(); ++i)
    {
      out << i << ": ";
      std::copy(_graph[i].begin(),
		_graph[i].end(),
		std::ostream_iterator<unsigned>(out, " "));
      out << std::endl;
    }
}



  void LaplaceMeshSmoother::allgather_graph()
  {
    // The graph data structure is not well-suited for parallel communication,
    // so copy the graph into a single vector defined by:
    // NA A_0 A_1 ... A_{NA} | NB B_0 B_1 ... B_{NB} | NC C_0 C_1 ... C_{NC}
    // where:
    // * NA is the number of graph connections for node A
    // * A_0, A_1, etc. are the IDs connected to node A
    std::vector<unsigned> flat_graph;

    // Reserve at least enough space for each node to have zero entries
    flat_graph.reserve(_graph.size());

    for (unsigned i=0; i<_graph.size(); ++i)
      {
	// First push back the number of entries for this node
	flat_graph.push_back (_graph[i].size());

	// Then push back all the IDs
	for (unsigned j=0; j<_graph[i].size(); ++j)
	  flat_graph.push_back(_graph[i][j]);
      }

    // // A copy of the flat graph (for printing only, delete me later)
    // std::vector<unsigned> copy_of_flat_graph(flat_graph);

    // Use the allgather routine to combine all the flat graphs on all processors
    Parallel::allgather(flat_graph);

    // Now reconstruct _graph from the allgathered flat_graph.

    // // (Delete me later, the copy is just for printing purposes.)
    // std::vector<std::vector<unsigned > > copy_of_graph(_graph);

    // Make sure the old graph is cleared out
    _graph.clear();
    _graph.resize(_mesh.n_nodes());

    // Our current position in the allgather'd flat_graph
    unsigned cursor=0;

    // There are n_nodes * n_processors vectors to read in total
    for (unsigned p=0; p<libMesh::n_processors(); ++p)
      for (unsigned node_ctr=0; node_ctr<_mesh.n_nodes(); ++node_ctr)
	{
	  // Read the number of entries for this node, move cursor
	  unsigned n_entries = flat_graph[cursor++];

	  // Reserve space for that many more entries, then push back
	  _graph[node_ctr].reserve(_graph[node_ctr].size() + n_entries);

	  // Read all graph connections for this node, move the cursor each time
	  // Note: there might be zero entries but that's fine
	  for (unsigned i=0; i<n_entries; ++i)
	    _graph[node_ctr].push_back(flat_graph[cursor++]);
	}


//    // Print local graph to uniquely named file (debugging)
//    {
//      // Generate unique filename for this processor
//      std::ostringstream oss;
//      oss << "graph_filename_" << libMesh::processor_id() << ".txt";
//      std::ofstream graph_stream(oss.str().c_str());
//
//      // Print the local non-flat graph
//      std::swap(_graph, copy_of_graph);
//      print_graph(graph_stream);
//
//      // Print the (local) flat graph for verification
//      for (unsigned i=0; i<copy_of_flat_graph.size(); ++i)
//	graph_stream << copy_of_flat_graph[i] << " ";
//      graph_stream << "\n";
//
//      // Print the allgather'd grap for verification
//      for (unsigned i=0; i<flat_graph.size(); ++i)
//	graph_stream << flat_graph[i] << " ";
//      graph_stream << "\n";
//
//      // Print the global non-flat graph
//      std::swap(_graph, copy_of_graph);
//      print_graph(graph_stream);
//    }
  } // allgather_graph()




  // Re-open anonymous namepsace, define functor functions
  namespace {

    SyncNodalPositions::SyncNodalPositions(MeshBase& m)
      : mesh(m)
    {}



    void SyncNodalPositions::gather_data (const std::vector<unsigned int>& ids, std::vector<datum>& data)
    {
      data.resize(ids.size());

      // Gather (x,y,z) data for all node IDs in the ids vector
      for (unsigned i=0; i<ids.size(); ++i)
	{
	  // Look for this node in the mesh
	  Node *node = mesh.node_ptr(ids[i]);

	  if (node == NULL)
	    {
	      libMesh::err << "Error! Mesh returned a NULL node pointer in SyncNodalPosition::gather_data()." << std::endl;
	      libmesh_error();
	    }

	  // Store this node's position in the data array.
	  // This should call Point::op=
	  data[i] = *node;
	} // end for
    } // gather_data()



    void SyncNodalPositions::act_on_data (const std::vector<unsigned int>& ids, std::vector<datum>& data)
    {
      for (unsigned i=0; i<ids.size(); ++i)
	{

	  // Get a pointer to the node whose position is to be updated.
	  Node* node = mesh.node_ptr(ids[i]);

	  if (node == NULL)
	    {
	      libMesh::err << "Error! Mesh returned a NULL node pointer in SyncNodalPosition::act_on_data()." << std::endl;
	      libmesh_error();
	    }

	  // Update this node's position.  Should call Point::op=
	  *node = data[i];
	} // end for
    } // act_on_data()

  } // anonymous namespace


} // namespace libMesh
