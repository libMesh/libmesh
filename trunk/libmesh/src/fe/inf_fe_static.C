// $Id: inf_fe_static.C,v 1.15 2003-04-02 21:58:44 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



// Local includes
#include "mesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS
#include "inf_fe.h"
#include "fe.h"
#include "fe_interface.h"
#include "elem.h"


// ------------------------------------------------------------
// InfFE class static member initialization
template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
ElemType InfFE<Dim,T_radial,T_map>::_compute_node_indices_fast_current_elem_type = INVALID_ELEM;




// ------------------------------------------------------------
// InfFE static class members
template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::n_shape_functions (const FEType& fet,
							   const ElemType t)
{
  return InfFE<Dim,T_radial,T_map>::n_dofs(fet, t);
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::n_dofs (const FEType& fet,
						const ElemType inf_elem_type)
{
  const ElemType base_et (Base::get_elem_type(inf_elem_type));
    
  if (Dim > 1)
    return FEInterface::n_dofs(Dim-1, fet, base_et) * Radial::n_dofs(fet.radial_order);
  else
    return Radial::n_dofs(fet.radial_order);
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::n_dofs_at_node (const FEType& fet,
							const ElemType inf_elem_type,
							const unsigned int n)
{
  const ElemType base_et (Base::get_elem_type(inf_elem_type));

  unsigned int n_base, n_radial;
  compute_node_indices(inf_elem_type, n, n_base, n_radial);
  
//   std::cout << "elem_type=" << inf_elem_type 
// 	    << ",  fet.radial_order=" << fet.radial_order
// 	    << ",  n=" << n 
// 	    << ",  n_radial=" << n_radial 
// 	    << ",  n_base=" << n_base
// 	    << std::endl;

  if (Dim > 1)
    return FEInterface::n_dofs_at_node(Dim-1, fet, base_et, n_base) 
        * Radial::n_dofs_at_node(fet.radial_order, n_radial);
  else
    return Radial::n_dofs_at_node(fet.radial_order, n_radial);
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::n_dofs_per_elem (const FEType& fet,
							 const ElemType inf_elem_type)
{
  const ElemType base_et (Base::get_elem_type(inf_elem_type));

  if (Dim > 1)
    return FEInterface::n_dofs_per_elem(Dim-1, fet, base_et) 
        * Radial::n_dofs_per_elem(fet.radial_order);
  else
    return Radial::n_dofs_per_elem(fet.radial_order);
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::nodal_soln(const FEType& /* fet */,
					   const Elem* elem,
					   const std::vector<Number>& /* elem_soln */,
					   std::vector<Number>&       nodal_soln)
{
  const unsigned int n_nodes = elem->n_nodes();
  
//  const ElemType inf_elem_type = elem->type();

  nodal_soln.resize(n_nodes);

//TODO:[DD] CONTINUE CONTINUE CONTINUE CONTINUE CONTINUE CONTINUE CONTINUE CONTINUE CONTINUE CONTINUE 

  std::cerr << "ERROR: The concept of a nodal solution is not "
	    << "applicable to infinite elements!" << std::endl;
  error();  
}





template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Real InfFE<Dim,T_radial,T_map>::shape(const FEType& fet,
				      const ElemType inf_elem_type,
				      const unsigned int i,
				      const Point& p)
{
  assert (Dim != 0);

  const ElemType     base_et  (Base::get_elem_type(inf_elem_type));
  const Order        o_radial (fet.radial_order);
  const Real         v        (p(Dim-1));

  unsigned int i_base, i_radial;
  compute_shape_indices(fet, inf_elem_type, i, i_base, i_radial);

  //TODO:[SP/DD]  exp(ikr) is still missing here!
  if (Dim > 1)
    return FEInterface::shape(Dim-1, fet, base_et, i_base, p)
        * InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial)
        * InfFE<Dim,T_radial,T_map>::Radial::decay(v);
  else
    return InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial)
        * InfFE<Dim,T_radial,T_map>::Radial::decay(v);
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Real InfFE<Dim,T_radial,T_map>::shape(const FEType& fet,
				      const Elem* inf_elem,
				      const unsigned int i,
				      const Point& p)
{
  assert (inf_elem != NULL);
  assert (Dim != 0);

  const Order        o_radial (fet.radial_order);
  const Real         v        (p(Dim-1));
  AutoPtr<Elem>      base_el  (inf_elem->build_side(0));

  unsigned int i_base, i_radial;
  compute_shape_indices(fet, inf_elem->type(), i, i_base, i_radial);

  //TODO:[SP/DD]  exp(ikr) is still missing here!
  if (Dim > 1)
    return FEInterface::shape(Dim-1, fet, base_el.get(), i_base, p)
        * InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial)
        * InfFE<Dim,T_radial,T_map>::Radial::decay(v);
  else
    return InfFE<Dim,T_radial,T_map>::eval(v, o_radial, i_radial)
        * InfFE<Dim,T_radial,T_map>::Radial::decay(v);
}

















template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::compute_node_indices (const ElemType inf_elem_type,
						      const unsigned int outer_node_index,
						      unsigned int& base_node,
						      unsigned int& radial_node)
{
  switch (inf_elem_type)
    {
    case INFEDGE2:
      {
	assert (outer_node_index < 2);
	base_node   = 0;
	radial_node = outer_node_index;
	return;
      }


    // linear base approximation, easy to determine
    case INFQUAD4:
      {
	assert (outer_node_index < 4);
	base_node   = outer_node_index % 2;
	radial_node = outer_node_index / 2;
	return;
      }

    case INFPRISM6:
      {
	assert (outer_node_index < 6);
	base_node   = outer_node_index % 3;
	radial_node = outer_node_index / 3;
	return;
      }

    case INFHEX8:
      {
	assert (outer_node_index < 8);
	base_node   = outer_node_index % 4;
	radial_node = outer_node_index / 4;
	return;
      }


    // higher order base approximation, more work necessary
    case INFQUAD6:
      {
	switch (outer_node_index)
	  {
	  case 0:
	  case 1:
	    {
	      radial_node = 0;
	      base_node   = outer_node_index;
	      return;
	    }

	  case 2:
	  case 3:
	    {
	      radial_node = 1;
	      base_node   = outer_node_index-2;
	      return;
	    }

	  case 4:
	    {
	      radial_node = 0;
	      base_node   = 2;
	      return;
	    }

	  case 5:
	    {
	      radial_node = 1;
	      base_node   = 2;
	      return;
	    }

	  default:
	    {
	      error();
	      return;
	    }
	  }
      }


    case INFHEX16:
    case INFHEX18:
      {
	switch (outer_node_index)
	  {
	  case 0:
	  case 1:
	  case 2:
	  case 3:
	    {
	      radial_node = 0;
	      base_node   = outer_node_index;
	      return;
	    }

	  case 4:
	  case 5:
	  case 6:
	  case 7:
	    {
	      radial_node = 1;
	      base_node   = outer_node_index-4;
	      return;
	    }

	  case 8:
	  case 9:
	  case 10:
	  case 11:
	    {
	      radial_node = 0;
	      base_node   = outer_node_index-4;
	      return;
	    }

	  case 12:
	  case 13:
	  case 14:
	  case 15:
	    {
	      radial_node = 1;
	      base_node   = outer_node_index-8;
	      return;
	    }

	  case 16:
	    {
	      assert (inf_elem_type == INFHEX18);
	      radial_node = 0;
	      base_node   = 8;
	      return;
	    }

	  case 17:
	    {
	      assert (inf_elem_type == INFHEX18);
	      radial_node = 1;
	      base_node   = 8;
	      return;
	    }

	  default:
	    {
	      error();
	      return;
	    }
	  }
      }


    case INFPRISM12:
      {
	switch (outer_node_index)
	  {
	  case 0:
	  case 1:
	  case 2:
	    {
	      radial_node = 0;
	      base_node   = outer_node_index;
	      return;
	    }

	  case 3:
	  case 4:
	  case 5:
	    {
	      radial_node = 1;
	      base_node   = outer_node_index-3;
	      return;
	    }

	  case 6:
	  case 7:
	  case 8:
	    {
	      radial_node = 0;
	      base_node   = outer_node_index-3;
	      return;
	    }

	  case 9:
	  case 10:
	  case 11:
	    {
	      radial_node = 1;
	      base_node   = outer_node_index-6;
	      return;
	    }

	  default:
	    {
	      error();
	      return;
	    }
	  }
      }


    default:
      { 
        std::cerr << "ERROR: Bad infinite element type=" << inf_elem_type 
		  << ", node=" << outer_node_index << std::endl;
	error();
	return;
      }
    }
}






template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::compute_node_indices_fast (const ElemType inf_elem_type,
							   const unsigned int outer_node_index,
							   unsigned int& base_node,
							   unsigned int& radial_node)
{
  assert (inf_elem_type != INVALID_ELEM);

  static std::vector<unsigned int> _static_base_node_index;
  static std::vector<unsigned int> _static_radial_node_index;

  /* 
   * fast counterpart to compute_node_indices(), uses local static buffers
   * to store the index maps.  The class member
   * \p _compute_node_indices_fast_current_elem_type remembers
   * the current element type.
   *
   * Note that there exist non-static members storing the
   * same data.  However, you never know what element type
   * is currently used by the \p InfFE object, and what
   * request is currently directed to the static \p InfFE
   * members (which use \p compute_node_indices_fast()).
   * So separate these.
   *
   * check whether the work for this elemtype has already 
   * been done.  If so, use this index.  Otherwise, refresh
   * the buffer to this element type.
   */
  if (inf_elem_type==_compute_node_indices_fast_current_elem_type)
    {
      base_node   = _static_base_node_index  [outer_node_index];
      radial_node = _static_radial_node_index[outer_node_index];
      return;
    }
  else
    {      
      // store the map for _all_ nodes for this element type
      _compute_node_indices_fast_current_elem_type = inf_elem_type;

      unsigned int n_nodes;

      switch (inf_elem_type)
        {
        case INFEDGE2:
	  {
	    n_nodes = 2;
	    break;
	  }
	case INFQUAD4:
	  {
	    n_nodes = 4;
	    break;
	  }
        case INFQUAD6:
	  {
	    n_nodes = 6;
	    break;
	  }
        case INFHEX8:
	  {
	    n_nodes = 8;
	    break;
	  }
        case INFHEX16:
	  {
	    n_nodes = 16;
	    break;
	  }
        case INFHEX18:
	  {
	    n_nodes = 18;
	    break;
	  }
        case INFPRISM6:
	  {
	    n_nodes = 6;
	    break;
	  }
        case INFPRISM12:
	  {
	    n_nodes = 12;
	    break;
	  }
	default:
	  { 
	    std::cerr << "ERROR: Bad infinite element type=" << inf_elem_type 
		      << ", node=" << outer_node_index << std::endl;
	    error();
	    break;
	  }
	}


      _static_base_node_index.resize  (n_nodes);
      _static_radial_node_index.resize(n_nodes);

      for (unsigned int n=0; n<n_nodes; n++)
	  compute_node_indices (inf_elem_type, 
				n, 
				_static_base_node_index  [outer_node_index], 
				_static_radial_node_index[outer_node_index]);

      // and return for the specified node
      base_node   = _static_base_node_index  [outer_node_index];
      radial_node = _static_radial_node_index[outer_node_index];
      return;
    }
}






template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
void InfFE<Dim,T_radial,T_map>::compute_shape_indices (const FEType& fet,
						       const ElemType inf_elem_type,
						       const unsigned int i,
						       unsigned int& base_shape,
						       unsigned int& radial_shape)
{

  /*
   * An example is provided:  the numbers in comments refer to
   * a fictitious InfHex18.  The numbers are chosen as exemplary
   * values.  There is currently no base approximation that
   * requires this many dof's at nodes, sides, faces and in the element.
   *
   * the order of the shape functions is heavily related with the
   * order the dofs are assigned in \p DofMap::distributed_dofs().
   * Due to the infinite elements with higher-order base approximation,
   * some more effort is necessary.
   *
   * numbering scheme:
   * 1. all vertices in the base, assign node->n_comp() dofs to each vertex
   * 2. all vertices further out: innermost loop: radial shapes,
   *    then the base approximation shapes
   * 3. all side nodes in the base, assign node->n_comp() dofs to each side node
   * 4. all side nodes further out: innermost loop: radial shapes,
   *    then the base approximation shapes
   * 5. (all) face nodes in the base, assign node->n_comp() dofs to each face node
   * 6. (all) face nodes further out: innermost loop: radial shapes,
   *    then the base approximation shapes
   * 7. element-associated dof in the base
   * 8. element-associated dof further out
   */

  const unsigned int radial_order       = static_cast<unsigned int>(fet.radial_order);             // 4
  const unsigned int radial_order_p_one = radial_order+1;                                          // 5

  const ElemType base_elem_type           (Base::get_elem_type(inf_elem_type));                    // QUAD9

  // assume that the number of dof is the same for all vertices
  unsigned int n_base_vertices;                                                                    // 4
  const unsigned int n_base_vertex_dof = FEInterface::n_dofs_at_node  (Dim-1, fet, base_elem_type, 0);// 2

  unsigned int n_base_side_nodes;                                                                  // 4
  unsigned int n_base_side_dof;                                                                    // 3

  unsigned int n_base_face_nodes;                                                                  // 1
  unsigned int n_base_face_dof;                                                                    // 5

  const unsigned int n_base_elem_dof   = FEInterface::n_dofs_per_elem (Dim-1, fet, base_elem_type);// 9


  switch (inf_elem_type)
    {
    case INFEDGE2:
      {
	n_base_vertices   = 1;
	n_base_side_nodes = 0;
	n_base_face_nodes = 0;
	n_base_side_dof   = 0;
	n_base_face_dof   = 0;
	break;
      }

    case INFQUAD4:
      {
	n_base_vertices   = 2;
	n_base_side_nodes = 0;
	n_base_face_nodes = 0;
	n_base_side_dof   = 0;
	n_base_face_dof   = 0;
	break;
      }

    case INFQUAD6:
      {
	n_base_vertices   = 2;
	n_base_side_nodes = 1;
	n_base_face_nodes = 0;
	n_base_side_dof   = FEInterface::n_dofs_at_node (Dim-1, fet,base_elem_type, n_base_vertices);
	n_base_face_dof   = 0;
	break;
      }
	
    case INFHEX8:
      {
	n_base_vertices   = 4;
	n_base_side_nodes = 0;
	n_base_face_nodes = 0;
	n_base_side_dof   = 0;
	n_base_face_dof   = 0;
	break;
      }

    case INFHEX16:
      {
	n_base_vertices   = 4;
	n_base_side_nodes = 4;
	n_base_face_nodes = 0;
	n_base_side_dof   = FEInterface::n_dofs_at_node (Dim-1, fet,base_elem_type, n_base_vertices);
	n_base_face_dof   = 0;
	break;
      }

    case INFHEX18:
      {
	n_base_vertices   = 4;
	n_base_side_nodes = 4;
	n_base_face_nodes = 1;
	n_base_side_dof   = FEInterface::n_dofs_at_node (Dim-1, fet,base_elem_type, n_base_vertices);
	n_base_face_dof   = FEInterface::n_dofs_at_node (Dim-1, fet,base_elem_type, 8);
	break;
      }

		 
    case INFPRISM6:
      {
	n_base_vertices   = 3;
	n_base_side_nodes = 0;
	n_base_face_nodes = 0;
	n_base_side_dof   = 0;
	n_base_face_dof   = 0;
	break;
      }

    case INFPRISM12:
      {
	n_base_vertices   = 3;
	n_base_side_nodes = 3;
	n_base_face_nodes = 0;
	n_base_side_dof   = FEInterface::n_dofs_at_node (Dim-1, fet,base_elem_type, n_base_vertices);
	n_base_face_dof   = 0;
	break;
      }

    default:
	error();
    }


  {
    // these are the limits describing the intervals where the shape function lies
    const unsigned int n_dof_at_base_vertices = n_base_vertices*n_base_vertex_dof;                 // 8
    const unsigned int n_dof_at_all_vertices  = n_dof_at_base_vertices*radial_order_p_one;         // 40

    const unsigned int n_dof_at_base_sides    = n_base_side_nodes*n_base_side_dof;                 // 12
    const unsigned int n_dof_at_all_sides     = n_dof_at_base_sides*radial_order_p_one;            // 60

    const unsigned int n_dof_at_base_face     = n_base_face_nodes*n_base_face_dof;                 // 5
    const unsigned int n_dof_at_all_faces     = n_dof_at_base_face*radial_order_p_one;             // 25


    // start locating the shape function
    if (i < n_dof_at_base_vertices)                                              // range of i: 0..7
      {
	// belongs to vertex in the base
	radial_shape = 0;
	base_shape   = i;
      }

    else if (i < n_dof_at_all_vertices)                                          // range of i: 8..39
      {
	/* belongs to vertex in the outer shell
	 *
	 * subtract the number of dof already counted,
	 * so that i_offset contains only the offset for the base
	 */
	const unsigned int i_offset = i - n_dof_at_base_vertices;                // 0..31
	
	// first the radial dof are counted, then the base dof
	radial_shape = (i_offset % radial_order) + 1;
	base_shape   = i_offset / radial_order;
      }

    else if (i < n_dof_at_all_vertices+n_dof_at_base_sides)                      // range of i: 40..51
      {
	// belongs to base, is a side node
	radial_shape = 0;
	base_shape = i - radial_order * n_dof_at_base_vertices;                  //  8..19
      }

    else if (i < n_dof_at_all_vertices+n_dof_at_all_sides)                       // range of i: 52..99
      {
	// belongs to side node in the outer shell
	const unsigned int i_offset = i - (n_dof_at_all_vertices
					   + n_dof_at_base_sides);               // 0..47
	radial_shape = (i_offset % radial_order) + 1;
	base_shape   = (i_offset / radial_order) + n_dof_at_base_vertices;
      }

    else if (i < n_dof_at_all_vertices+n_dof_at_all_sides+n_dof_at_base_face)    // range of i: 100..104
      {
	// belongs to the node in the base face
	radial_shape = 0;
	base_shape = i - radial_order*(n_dof_at_base_vertices
				       + n_dof_at_base_sides);                   //  20..24
      }

    else if (i < n_dof_at_all_vertices+n_dof_at_all_sides+n_dof_at_all_faces)    // range of i: 105..124
      {
	// belongs to the node in the outer face
	const unsigned int i_offset = i - (n_dof_at_all_vertices 
					   + n_dof_at_all_sides 
					   + n_dof_at_base_face);                // 0..19
	radial_shape = (i_offset % radial_order) + 1;
	base_shape   = (i_offset / radial_order) + n_dof_at_base_vertices + n_dof_at_base_sides;
      }

    else if (i < n_dof_at_all_vertices+n_dof_at_all_sides+n_dof_at_all_faces+n_base_elem_dof)      // range of i: 125..133
      {
	// belongs to the base and is an element associated shape
	radial_shape = 0;
	base_shape = i - (n_dof_at_all_vertices 
			  + n_dof_at_all_sides 
			  + n_dof_at_all_faces);                                 // 0..8 
      }

    else                                                                         // range of i: 134..169
      {
	assert (i < n_dofs(fet, inf_elem_type));
	// belongs to the outer shell and is an element associated shape
	const unsigned int i_offset = i - (n_dof_at_all_vertices 
					   + n_dof_at_all_sides 
					   + n_dof_at_all_faces 
					   + n_base_elem_dof);                   // 0..19
	radial_shape = (i_offset % radial_order) + 1;
	base_shape   = (i_offset / radial_order) + n_dof_at_base_vertices + n_dof_at_base_sides + n_dof_at_base_face;
      }
  }

  return;
}



//--------------------------------------------------------------
// Explicit instantiations
#include "inf_fe_instantiate_1D.h"
#include "inf_fe_instantiate_2D.h"
#include "inf_fe_instantiate_3D.h"



#endif //ifdef ENABLE_INFINITE_ELEMENTS

