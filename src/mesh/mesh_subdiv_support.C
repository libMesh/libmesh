// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/mesh_tools.h"
#include "libmesh/mesh_subdiv_support.h"
#include "libmesh/boundary_info.h"

namespace libMesh
{


void MeshTools::Subdiv::find_one_ring(const Tri3SD* elem, std::vector<Node*>& nodes)
{
	libmesh_assert(elem->is_subdiv_updated());
	libmesh_assert(elem->get_ordered_node(0));

	unsigned int valence = elem->get_ordered_valence(0);
	nodes.resize(valence + 6);

	// The first three vertices in the patch are the ones from the element triangle
	nodes[0]       = elem->get_ordered_node(0);
	nodes[1]       = elem->get_ordered_node(1);
	nodes[valence] = elem->get_ordered_node(2);

	const unsigned int nn0 = elem->local_node_number(nodes[0]->id());

	Tri3SD* nb = dynamic_cast<Tri3SD*>(elem->neighbor(nn0)); 
	libmesh_assert(nb);

	unsigned int j, i = 1;

	do
	{
		++i;
		j = nb->local_node_number(nodes[0]->id());
		nodes[i] = nb->get_node(next[j]);
		nb = static_cast<Tri3SD*>(nb->neighbor(j));
	} while (nb != elem);

	/* for nodes connected with N (= valence[0]) */
	nb = static_cast<Tri3SD*>(elem->neighbor(next[nn0]));
	j = nb->local_node_number(nodes[1]->id());
	nodes[valence+1] = nb->get_node(next[j]);

	nb = static_cast<Tri3SD*>(nb->neighbor(next[j]));
	j = nb->local_node_number(nodes[valence+1]->id());
	nodes[valence+4] = nb->get_node(next[j]);

	nb = static_cast<Tri3SD*>(nb->neighbor(next[j]));
	j = nb->local_node_number(nodes[valence+4]->id());
	nodes[valence+5] = nb->get_node(next[j]);

	/* for nodes connected with 1 */
	nb = static_cast<Tri3SD*>(elem->neighbor(next[nn0]));
	j = nb->local_node_number(nodes[1]->id());
	// nodes[valence+1] has been determined already

	nb = static_cast<Tri3SD*>(nb->neighbor(j));
	j = nb->local_node_number(nodes[1]->id());
	nodes[valence+2] = nb->get_node(next[j]);

	nb = static_cast<Tri3SD*>(nb->neighbor(j));
	j = nb->local_node_number(nodes[1]->id());
	nodes[valence+3] = nb->get_node(next[j]);

	return;
}


void MeshTools::Subdiv::all_subdiv(MeshBase& mesh)
{
	std::vector<Elem*> new_elements;
	new_elements.reserve(mesh.n_elem());
	const bool mesh_has_boundary_data = (mesh.boundary_info->n_boundary_ids() > 0);

	std::vector<Elem*> new_boundary_elements;
	std::vector<unsigned int> new_boundary_sides;
	std::vector<short int> new_boundary_ids;

	MeshBase::const_element_iterator       el     = mesh.elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.elements_end(); 
	for (; el != end_el; ++el)
	{
		const Elem* elem = *el;
		libmesh_assert_equal_to(elem->type(), TRI3);

		Elem* tri = new Tri3SD;
		tri->set_id(elem->id());
		tri->set_node(0) = (*el)->get_node(0);
		tri->set_node(1) = (*el)->get_node(1);
		tri->set_node(2) = (*el)->get_node(2);

		if (mesh_has_boundary_data)
		{
			for (unsigned int side = 0; side < elem->n_sides(); ++side)
			{
				const short int boundary_id = mesh.boundary_info->boundary_id(elem, side);
				if (boundary_id != BoundaryInfo::invalid_id)
				{
					// add the boundary id to the list of new boundary ids
					new_boundary_ids.push_back(boundary_id);
					new_boundary_elements.push_back(tri);
					new_boundary_sides.push_back(side);
				}
			}
				
			// remove the original element from the BoundaryInfo structure
			mesh.boundary_info->remove(elem);
		}

		new_elements.push_back(tri);
		mesh.insert_elem(tri);
	}
	mesh.prepare_for_use();

	if (mesh_has_boundary_data)
	{
		// If the old mesh had boundary data, the new mesh better have some too.
		libmesh_assert_greater(new_boundary_elements.size(), 0);

		// We should also be sure that the lengths of the new boundary data vectors
		// are all the same.
		libmesh_assert_equal_to(new_boundary_sides.size(), new_boundary_elements.size());
		libmesh_assert_equal_to(new_boundary_sides.size(), new_boundary_ids.size());

		// Add the new boundary info to the mesh.
		for (unsigned int s = 0; s < new_boundary_elements.size(); ++s)
			mesh.boundary_info->add_side(new_boundary_elements[s],
						     new_boundary_sides[s],
						     new_boundary_ids[s]);
	}

	mesh.prepare_for_use();
}


void MeshTools::Subdiv::prepare_subdiv_mesh(MeshBase& mesh, bool ghosted)
{
	mesh.prepare_for_use();

	// convert all mesh elements to subdivision elements
	all_subdiv(mesh);

	if (!ghosted)
	{
		// add the ghost elements for the boundaries
		add_boundary_ghosts(mesh);
	}
	else
	{
		// This assumes that the mesh already has the ghosts. Only tagging them is required here.
		tag_boundary_ghosts(mesh);
	}

	mesh.prepare_for_use();
   
	std::vector<std::vector<const Elem*> > nodes_to_elem_map;
	MeshTools::build_nodes_to_elem_map(mesh, nodes_to_elem_map);

	// compute the node valences
	MeshBase::const_node_iterator       nd     = mesh.nodes_begin();
	const MeshBase::const_node_iterator end_nd = mesh.nodes_end(); 
	for (; nd != end_nd; ++nd)
	{  
		Node* node = *nd;
		std::vector<const Node*> neighbors;
		MeshTools::find_nodal_neighbors(mesh, *node, nodes_to_elem_map, neighbors);
		const unsigned int valence = neighbors.size();
		libmesh_assert_greater(valence, 1);
		node->set_valence(valence);
	}

	MeshBase::const_element_iterator       el     = mesh.elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.elements_end(); 
	for (; el != end_el; ++el)
	{  
		Tri3SD* elem = dynamic_cast<Tri3SD*>(*el);
		libmesh_assert(elem);
		if (!elem->is_ghost())
			elem->prepare_subdiv_properties();
	}
}


void MeshTools::Subdiv::tag_boundary_ghosts(MeshBase& mesh)
{
	MeshBase::element_iterator       el     = mesh.elements_begin();
	const MeshBase::element_iterator end_el = mesh.elements_end(); 
	for (; el != end_el; ++el)
	{ 
		Elem* elem = *el;
		libmesh_assert_equal_to(elem->type(), TRI3SD);

		Tri3SD* sd_elem = static_cast<Tri3SD*>(elem);
		for (unsigned int i = 0; i < elem->n_sides(); ++i)
		{
			if (elem->neighbor(i) == NULL) 
			{
				sd_elem->set_ghost(true);
				// set all other neighbors to ghosts as well
				if (elem->neighbor(next[i]))
				{
					Tri3SD* nb = static_cast<Tri3SD*>(elem->neighbor(next[i]));
					nb->set_ghost(true);
				}
				if (elem->neighbor(prev[i]))
				{
					Tri3SD* nb = static_cast<Tri3SD*>(elem->neighbor(prev[i]));
					nb->set_ghost(true);
				}
			}
		}
	}
}


void MeshTools::Subdiv::add_boundary_ghosts(MeshBase& mesh)
{
	static const Real tol = 1e-5;

	// add the mirrored ghost elements (without using iterators, because the mesh is modified in the course)
	std::vector<Tri3SD*> ghost_elems;
	std::vector<Node*> ghost_nodes;
	const unsigned int n_elem = mesh.n_elem();
	for (unsigned int eid = 0; eid < n_elem; ++eid)
	{   
		Elem* elem = mesh.elem(eid);
		libmesh_assert_equal_to(elem->type(), TRI3SD);

		for (unsigned int i = 0; i < elem->n_sides(); ++i)
		{
			libmesh_assert_not_equal_to(elem->neighbor(i), elem);
			if (elem->neighbor(i) == NULL)
			{
				// this is the vertex to be mirrored
				Point point = elem->point(i) + elem->point(next[i]) - elem->point(prev[i]);

				// Check if the proposed vertex doesn't coincide with one of the existing vertices.
				// This is necessary because for some triangulations, it can happen that two mirrored
				// ghost vertices coincide, which would then lead to a zero size ghost element below.
				Node* node = NULL;
				for (unsigned int j = 0; j < ghost_nodes.size(); ++j)
				{
					if ((*ghost_nodes[j] - point).size() < tol * (elem->point(i) - point).size())
					{
						node = ghost_nodes[j];
						break;
					}
				}

				// add the new vertex only if no other is nearby
				if (node == NULL)
				{
					node = mesh.add_point(point);
					ghost_nodes.push_back(node);
				}

				Tri3SD* newelem = new Tri3SD();
				ghost_elems.push_back(newelem);
				
				newelem->set_node(0) = elem->get_node(next[i]);
				newelem->set_node(1) = elem->get_node(i);
				newelem->set_node(2) = node;
				newelem->set_neighbor(0,elem);
				newelem->set_ghost(true);
				elem->set_neighbor(i,newelem);

				mesh.add_elem(newelem);
				mesh.boundary_info->add_node(elem->get_node(i), 1);
				mesh.boundary_info->add_node(elem->get_node(next[i]), 1);
				mesh.boundary_info->add_node(elem->get_node(prev[i]), 1);
				mesh.boundary_info->add_node(node, 1);
			}
		}
	}
      
	// add the missing ghost elements (connecting new ghost nodes)
	std::vector<Tri3SD*> missing_ghost_elems;
	std::vector<Tri3SD*>::iterator       ghost_el     = ghost_elems.begin();
	const std::vector<Tri3SD*>::iterator end_ghost_el = ghost_elems.end();
	for (; ghost_el != end_ghost_el; ++ghost_el)
	{
		Tri3SD *elem = *ghost_el;
		libmesh_assert(elem->is_ghost());

		for (unsigned int i = 0; i < elem->n_sides(); ++i)
		{
			if (elem->neighbor(i) == NULL && elem->neighbor(prev[i]) != NULL) 
			{
				// go around counter-clockwise
				Tri3SD *nb1 = static_cast<Tri3SD *>(elem->neighbor(prev[i]));
				Tri3SD *nb2 = nb1;
				unsigned int j = i;
				while (nb1 != NULL && nb1->id() != elem->id())
				{
					j = nb1->local_node_number(elem->node(i));
					nb2 = nb1;
					nb1 = static_cast<Tri3SD *>(nb1->neighbor(prev[j]));
					libmesh_assert(nb1 == NULL || nb1->id() != nb2->id());
				}

				libmesh_assert_not_equal_to(nb2->id(), elem->id());

				// Above, we merged coinciding ghost vertices. Therefore, we need
				// to exclude the case where there is no ghost element to add between
				// these two (identical) ghost nodes.
				if (elem->get_node(next[i])->id() == nb2->get_node(prev[j])->id())
					break;

				Tri3SD *newelem = new Tri3SD();
				newelem->set_node(0) = elem->get_node(next[i]);
				newelem->set_node(1) = elem->get_node(i);
				newelem->set_node(2) = nb2->get_node(prev[j]);
				newelem->set_neighbor(0,elem);
				newelem->set_neighbor(1,nb2);
				newelem->set_neighbor(2,NULL);
				newelem->set_ghost(true);

				elem->set_neighbor(i,newelem);
				nb2->set_neighbor(prev[j],newelem);

				missing_ghost_elems.push_back(newelem);
				break;             
			}
		} // end side loop
	} // end ghost element loop

	// add the missing ghost elements to the mesh
	std::vector<Tri3SD*>::iterator       missing_el     = missing_ghost_elems.begin();
	const std::vector<Tri3SD*>::iterator end_missing_el = missing_ghost_elems.end();
	for (; missing_el != end_missing_el; ++missing_el)
		mesh.add_elem(*missing_el);
}

} // namespace libMesh
