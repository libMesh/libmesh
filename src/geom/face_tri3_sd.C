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
#include "libmesh/face_tri3_sd.h"
#include "libmesh/mesh_subdiv_support.h"

namespace libMesh
{



// ------------------------------------------------------------
// Tri3 subdivision class member functions

Tri3SD::Tri3SD(Elem *p) : Tri3(p), _subdiv_updated(true)
{
	if (p)
	{
		libmesh_assert_equal_to(p->type(), TRI3SD);
		Tri3SD* sd_elem = static_cast<Tri3SD*>(p);
		_is_ghost = sd_elem->is_ghost();

		if (!_is_ghost)
		{
			_ordered_nodes[0] = sd_elem->get_ordered_node(0);
			_ordered_nodes[1] = sd_elem->get_ordered_node(1);
			_ordered_nodes[2] = sd_elem->get_ordered_node(2);
		}
	}
}


void Tri3SD::prepare_subdiv_properties()
{
	/*
	 * Find the index of the irregular vertex, if any.
	 * The current implementation can only handle triangles with
	 * no more than one irregular vertex. That is, a vertex with
	 * valence != 6.
	 */
	unsigned int irregular_idx = 0;
	for (unsigned int i = 0; i < 3; ++i)
	{
		if (this->get_node(i)->valence() != 6)
		{
			irregular_idx = i;
			if (this->get_node(MeshTools::Subdiv::next[i])->valence() != 6 || this->get_node(MeshTools::Subdiv::prev[i])->valence() != 6)
			{
				std::cout << "Error: The mesh contains elements with more than one irregular vertex!" << std::endl;
				libmesh_error();
			}
		}
	}

	/*
	 * Rotate ordered vertices such that ordered_nodes[0] is the
	 * irregular vertex. Doing this once in advance lets the evaluation
	 * of subdivision interpolation be much more efficient afterwards.
	 */
	switch (irregular_idx)
	{
		case 0:
			_ordered_nodes[0] = this->get_node(0);
			_ordered_nodes[1] = this->get_node(1);
			_ordered_nodes[2] = this->get_node(2);
			break;
		case 1:
			_ordered_nodes[0] = this->get_node(1);
			_ordered_nodes[1] = this->get_node(2);
			_ordered_nodes[2] = this->get_node(0);
			break;
		case 2:
			_ordered_nodes[0] = this->get_node(2);
			_ordered_nodes[1] = this->get_node(0);
			_ordered_nodes[2] = this->get_node(1);
			break;
		default:
		libmesh_error();
	}
	
	_subdiv_updated = true;
}


unsigned int Tri3SD::local_node_number(unsigned int node_id) const
{
	return (_nodes[0]->id() == node_id) ? 0 : ( (_nodes[1]->id() == node_id) ? 1 : ( (_nodes[2]->id() == node_id) ? 2 : 3 ) );
}


unsigned int Tri3SD::get_ordered_valence(unsigned int node_id) const
{
	libmesh_assert_less(node_id, n_neighbors());
	libmesh_assert(_subdiv_updated);
	return get_ordered_node(node_id)->valence();
}


Node* Tri3SD::get_ordered_node(unsigned int node_id) const
{
	libmesh_assert_less(node_id, 3);
	libmesh_assert(_subdiv_updated);
	return _ordered_nodes[node_id];
}

} // namespace libMesh
