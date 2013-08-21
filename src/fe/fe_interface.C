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



// Local includes
#include "libmesh/fe_interface.h"
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/fe_compute_data.h"
#include "libmesh/dof_map.h"

namespace libMesh
{

//------------------------------------------------------------
//FEInterface class members
FEInterface::FEInterface()
{
  libMesh::err << "ERROR: Do not define an object of this type."
	        << std::endl;
  libmesh_error();
}


#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
#define fe_family_switch(dim, func_and_args, prefix, suffix) \
      do { \
	switch (fe_t.family) \
	  { \
	  case CLOUGH: \
	    prefix FE<dim,CLOUGH>::func_and_args suffix \
	  case HERMITE: \
	    prefix FE<dim,HERMITE>::func_and_args suffix \
	  case HIERARCHIC: \
	    prefix FE<dim,HIERARCHIC>::func_and_args suffix \
	  case L2_HIERARCHIC: \
	    prefix FE<dim,L2_HIERARCHIC>::func_and_args suffix \
	  case LAGRANGE: \
	    prefix FE<dim,LAGRANGE>::func_and_args suffix \
	  case L2_LAGRANGE: \
	    prefix FE<dim,L2_LAGRANGE>::func_and_args suffix \
	  case MONOMIAL: \
	    prefix FE<dim,MONOMIAL>::func_and_args suffix \
          case SCALAR: \
            prefix FE<dim,SCALAR>::func_and_args suffix \
	  case BERNSTEIN: \
	    prefix FE<dim,BERNSTEIN>::func_and_args suffix \
	  case SZABAB: \
	    prefix FE<dim,SZABAB>::func_and_args suffix \
	  case XYZ: \
	    prefix FEXYZ<dim>::func_and_args suffix \
	  default: \
	    libmesh_error(); \
	  } \
      } while (0)

#define fe_family_with_vec_switch(dim, func_and_args, prefix, suffix) \
      do { \
	switch (fe_t.family) \
	  { \
	  case CLOUGH: \
	    prefix FE<dim,CLOUGH>::func_and_args suffix \
	  case HERMITE: \
	    prefix FE<dim,HERMITE>::func_and_args suffix \
	  case HIERARCHIC: \
	    prefix FE<dim,HIERARCHIC>::func_and_args suffix \
	  case L2_HIERARCHIC: \
	    prefix FE<dim,L2_HIERARCHIC>::func_and_args suffix \
	  case LAGRANGE: \
	    prefix FE<dim,LAGRANGE>::func_and_args suffix \
	  case LAGRANGE_VEC: \
	    prefix FELagrangeVec<dim>::func_and_args suffix \
	  case L2_LAGRANGE: \
	    prefix FE<dim,L2_LAGRANGE>::func_and_args suffix \
	  case MONOMIAL: \
	    prefix FE<dim,MONOMIAL>::func_and_args suffix \
          case SCALAR: \
            prefix FE<dim,SCALAR>::func_and_args suffix \
	  case BERNSTEIN: \
	    prefix FE<dim,BERNSTEIN>::func_and_args suffix \
	  case SZABAB: \
	    prefix FE<dim,SZABAB>::func_and_args suffix \
	  case XYZ: \
	    prefix FEXYZ<dim>::func_and_args suffix \
	  case NEDELEC_ONE: \
            prefix FENedelecOne<dim>::func_and_args suffix \
	  default: \
	    libmesh_error(); \
	  } \
      } while (0)

#define fe_scalar_vec_error_switch(dim, func_and_args, prefix, suffix) \
  do { \
  switch (fe_t.family)  \
    {  \
    case CLOUGH: \
      prefix FE<dim,CLOUGH>::func_and_args suffix \
    case HERMITE: \
      prefix FE<dim,HERMITE>::func_and_args suffix \
    case HIERARCHIC: \
      prefix FE<dim,HIERARCHIC>::func_and_args suffix \
    case L2_HIERARCHIC: \
      prefix FE<dim,L2_HIERARCHIC>::func_and_args suffix\
    case LAGRANGE: \
      prefix FE<dim,LAGRANGE>::func_and_args suffix\
    case L2_LAGRANGE: \
      prefix FE<dim,L2_LAGRANGE>::func_and_args suffix\
    case MONOMIAL: \
      prefix FE<dim,MONOMIAL>::func_and_args suffix\
    case SCALAR: \
      prefix FE<dim,SCALAR>::func_and_args suffix\
    case BERNSTEIN: \
      prefix FE<dim,BERNSTEIN>::func_and_args suffix\
    case SZABAB: \
      prefix FE<dim,SZABAB>::func_and_args suffix\
    case XYZ: \
      prefix FEXYZ<dim>::func_and_args suffix \
    case LAGRANGE_VEC: \
    case NEDELEC_ONE: \
      libMesh::err << "Error: Can only request scalar valued elements for Real FEInterface::func_and_args"\
		   << std::endl;\
      libmesh_error();\
    default: \
      libmesh_error(); \
    }\
    } while(0)


#define fe_vector_scalar_error_switch(dim, func_and_args, prefix, suffix) \
	do { \
	switch (fe_t.family) \
    { \
    case LAGRANGE_VEC: \
      prefix FELagrangeVec<dim>::func_and_args suffix \
    case NEDELEC_ONE:					    \
      prefix FENedelecOne<dim>::func_and_args suffix  \
    case HERMITE: \
    case HIERARCHIC: \
    case L2_HIERARCHIC: \
    case LAGRANGE: \
    case L2_LAGRANGE: \
    case MONOMIAL: \
    case SCALAR: \
    case BERNSTEIN: \
    case SZABAB: \
    case XYZ: \
      libMesh::err << "Error: Can only request vector valued elements for RealGradient FEInterface::shape" \
		   << std::endl; \
      libmesh_error();\
    default: \
      libmesh_error(); \
    } \
    } while(0)

#else
#define fe_family_switch(dim, func_and_args, prefix, suffix) \
      do { \
	switch (fe_t.family) \
	  { \
	  case CLOUGH: \
	    prefix FE<dim,CLOUGH>::func_and_args suffix \
	  case HERMITE: \
	    prefix FE<dim,HERMITE>::func_and_args suffix \
	  case HIERARCHIC: \
	    prefix FE<dim,HIERARCHIC>::func_and_args suffix \
	  case L2_HIERARCHIC: \
	    prefix FE<dim,L2_HIERARCHIC>::func_and_args suffix \
	  case LAGRANGE: \
	    prefix FE<dim,LAGRANGE>::func_and_args suffix \
	  case L2_LAGRANGE: \
	    prefix FE<dim,L2_LAGRANGE>::func_and_args suffix \
	  case MONOMIAL: \
	    prefix FE<dim,MONOMIAL>::func_and_args suffix \
          case SCALAR: \
            prefix FE<dim,SCALAR>::func_and_args suffix \
	  case XYZ: \
	    prefix FEXYZ<dim>::func_and_args suffix \
	  default: \
	    libmesh_error(); \
	  } \
      } while (0)

#define fe_family_with_vec_switch(dim, func_and_args, prefix, suffix) \
      do { \
	switch (fe_t.family) \
	  { \
	  case CLOUGH: \
	    prefix FE<dim,CLOUGH>::func_and_args suffix \
	  case HERMITE: \
	    prefix FE<dim,HERMITE>::func_and_args suffix \
	  case HIERARCHIC: \
	    prefix FE<dim,HIERARCHIC>::func_and_args suffix \
	  case L2_HIERARCHIC: \
	    prefix FE<dim,L2_HIERARCHIC>::func_and_args suffix \
	  case LAGRANGE: \
	    prefix FE<dim,LAGRANGE>::func_and_args suffix \
	  case LAGRANGE_VEC: \
	    prefix FELagrangeVec<dim>::func_and_args suffix \
	  case L2_LAGRANGE: \
	    prefix FE<dim,L2_LAGRANGE>::func_and_args suffix \
	  case MONOMIAL: \
	    prefix FE<dim,MONOMIAL>::func_and_args suffix \
          case SCALAR: \
            prefix FE<dim,SCALAR>::func_and_args suffix \
	  case XYZ: \
	    prefix FEXYZ<dim>::func_and_args suffix \
	  case NEDELEC_ONE: \
            prefix FENedelecOne<dim>::func_and_args suffix \
	  default: \
	    libmesh_error(); \
	  } \
      } while (0)

#define fe_scalar_vec_error_switch(dim, func_and_args, prefix, suffix) \
  do { \
  switch (fe_t.family)  \
    {  \
    case CLOUGH: \
      prefix  FE<dim,CLOUGH>::func_and_args suffix \
    case HERMITE: \
      prefix  FE<dim,HERMITE>::func_and_args suffix \
    case HIERARCHIC: \
      prefix  FE<dim,HIERARCHIC>::func_and_args suffix \
    case L2_HIERARCHIC: \
      prefix  FE<dim,L2_HIERARCHIC>::func_and_args suffix\
    case LAGRANGE: \
      prefix  FE<dim,LAGRANGE>::func_and_args suffix\
    case L2_LAGRANGE: \
      prefix  FE<dim,L2_LAGRANGE>::func_and_args suffix\
    case MONOMIAL: \
      prefix  FE<dim,MONOMIAL>::func_and_args suffix\
    case SCALAR: \
      prefix  FE<dim,SCALAR>::func_and_args suffix \
    case XYZ: \
      prefix  FEXYZ<dim>::func_and_args suffix \
    case LAGRANGE_VEC: \
    case NEDELEC_ONE: \
      libMesh::err << "Error: Can only request scalar valued elements for Real FEInterface::func_and_args"\
		   << std::endl;\
      libmesh_error();\
    default: \
      libmesh_error(); \
    }\
    } while(0)


#define fe_vector_scalar_error_switch(dim, func_and_args, prefix, suffix) \
	do { \
	switch (fe_t.family) \
    { \
    case LAGRANGE_VEC: \
      prefix FELagrangeVec<dim>::func_and_args suffix \
    case NEDELEC_ONE:				      \
      prefix FENedelecOne<dim>::func_and_args suffix \
    case HERMITE: \
    case HIERARCHIC: \
    case L2_HIERARCHIC: \
    case LAGRANGE: \
    case L2_LAGRANGE: \
    case MONOMIAL: \
    case SCALAR: \
    case XYZ: \
      libMesh::err << "Error: Can only request vector valued elements for RealGradient FEInterface::func_and_args" \
		   << std::endl; \
      libmesh_error();\
    default: \
      libmesh_error(); \
    } \
    } while(0)
#endif


#define fe_switch(func_and_args) \
  do { \
    switch (dim) \
      { \
        /* 0D */ \
      case 0: \
        fe_family_switch (0, func_and_args, return, ;); \
        /* 1D */ \
      case 1: \
        fe_family_switch (1, func_and_args, return, ;); \
        /* 2D */ \
      case 2: \
        fe_family_switch (2, func_and_args, return, ;); \
        /* 3D */ \
      case 3: \
        fe_family_switch (3, func_and_args, return, ;); \
      default: \
        libmesh_error(); \
      } \
  } while (0)

#define fe_with_vec_switch(func_and_args) \
  do { \
    switch (dim) \
      { \
        /* 0D */ \
      case 0: \
        fe_family_with_vec_switch (0, func_and_args, return, ;); \
        /* 1D */ \
      case 1: \
        fe_family_with_vec_switch (1, func_and_args, return, ;); \
        /* 2D */ \
      case 2: \
        fe_family_with_vec_switch (2, func_and_args, return, ;); \
        /* 3D */ \
      case 3: \
        fe_family_with_vec_switch (3, func_and_args, return, ;); \
      default: \
        libmesh_error(); \
      } \
  } while (0)


#define void_fe_switch(func_and_args) \
  do { \
    switch (dim) \
      { \
        /* 0D */ \
      case 0: \
        fe_family_switch (0, func_and_args, ;, ; return;); \
        /* 1D */ \
      case 1: \
        fe_family_switch (1, func_and_args, ;, ; return;); \
        /* 2D */ \
      case 2: \
        fe_family_switch (2, func_and_args, ;, ; return;); \
        /* 3D */ \
      case 3: \
        fe_family_switch (3, func_and_args, ;, ; return;); \
      default: \
        libmesh_error(); \
      } \
  } while (0)

#define void_fe_with_vec_switch(func_and_args) \
  do { \
    switch (dim) \
      { \
        /* 0D */ \
      case 0: \
        fe_family_with_vec_switch (0, func_and_args, ;, ; return;); \
        /* 1D */ \
      case 1: \
        fe_family_with_vec_switch (1, func_and_args, ;, ; return;); \
        /* 2D */ \
      case 2: \
        fe_family_with_vec_switch (2, func_and_args, ;, ; return;); \
        /* 3D */ \
      case 3: \
        fe_family_with_vec_switch (3, func_and_args, ;, ; return;); \
      default: \
        libmesh_error(); \
      } \
  } while (0)



unsigned int FEInterface::n_shape_functions(const unsigned int dim,
					    const FEType& fe_t,
					    const ElemType t)
{

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  /*
   * Since the FEType, stored in DofMap/(some System child), has to
   * be the _same_ for InfFE and FE, we have to catch calls
   * to infinite elements through the element type.
   */

  if ( is_InfFE_elem(t) )
    return ifem_n_shape_functions(dim, fe_t, t);

#endif

  const Order o = fe_t.order;

  fe_with_vec_switch(n_shape_functions(t, o));

  libmesh_error();
  return 0;
}





unsigned int FEInterface::n_dofs(const unsigned int dim,
				 const FEType& fe_t,
				 const ElemType t)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(t) )
    return ifem_n_dofs(dim, fe_t, t);

#endif

  const Order o = fe_t.order;

  fe_with_vec_switch(n_dofs(t, o));

  libmesh_error();
  return 0;
}




unsigned int FEInterface::n_dofs_at_node(const unsigned int dim,
					 const FEType& fe_t,
					 const ElemType t,
					 const unsigned int n)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(t) )
    return ifem_n_dofs_at_node(dim, fe_t, t, n);

#endif

  const Order o = fe_t.order;

  fe_with_vec_switch(n_dofs_at_node(t, o, n));

  libmesh_error();
  return 0;
}





unsigned int FEInterface::n_dofs_per_elem(const unsigned int dim,
					  const FEType& fe_t,
					  const ElemType t)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(t) )
    return ifem_n_dofs_per_elem(dim, fe_t, t);

#endif

  const Order o = fe_t.order;

  fe_with_vec_switch(n_dofs_per_elem(t, o));

  libmesh_error();
  return 0;
}




void FEInterface::dofs_on_side(const Elem* const elem,
			       const unsigned int dim,
			       const FEType& fe_t,
			       unsigned int s,
			       std::vector<unsigned int>& di)
{
  const Order o = fe_t.order;

  void_fe_with_vec_switch(dofs_on_side(elem, o, s, di));

  libmesh_error();
}



void FEInterface::dofs_on_edge(const Elem* const elem,
			       const unsigned int dim,
			       const FEType& fe_t,
			       unsigned int e,
			       std::vector<unsigned int>& di)
{
  const Order o = fe_t.order;

  void_fe_with_vec_switch(dofs_on_edge(elem, o, e, di));

  libmesh_error();
}




void FEInterface::nodal_soln(const unsigned int dim,
			     const FEType& fe_t,
			     const Elem* elem,
			     const std::vector<Number>& elem_soln,
			     std::vector<Number>&       nodal_soln)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(elem->type()) )
  {
    ifem_nodal_soln(dim, fe_t, elem, elem_soln, nodal_soln);
    return;
  }

#endif

  const Order order = fe_t.order;

  void_fe_with_vec_switch(nodal_soln(elem, order, elem_soln, nodal_soln));
}




Point FEInterface::map(unsigned int dim,
                       const FEType& fe_t,
                       const Elem* elem,
                       const Point& p)
{
  fe_with_vec_switch(map(elem, p));

  libmesh_error();
  return Point();
}





Point FEInterface::inverse_map (const unsigned int dim,
				const FEType& fe_t,
				const Elem* elem,
				const Point& p,
				const Real tolerance,
				const bool secure)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(elem->type()) )
    return ifem_inverse_map(dim, fe_t, elem, p,tolerance, secure);

#endif

  fe_with_vec_switch(inverse_map(elem, p, tolerance, secure));

  libmesh_error();
  return Point();
}




void FEInterface::inverse_map (const unsigned int dim,
			       const FEType& fe_t,
			       const Elem* elem,
			       const std::vector<Point>& physical_points,
			       std::vector<Point>&       reference_points,
			       const Real tolerance,
			       const bool secure)
{
  const std::size_t n_pts = physical_points.size();

  // Resize the vector
  reference_points.resize(n_pts);

  if (n_pts == 0)
    {
      libMesh::err << "WARNING: empty vector physical_points!"
		    << std::endl;
      libmesh_here();
      return;
    }

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(elem->type()) )
    {
      ifem_inverse_map(dim, fe_t, elem, physical_points, reference_points, tolerance, secure);
      return;

//       libMesh::err << "ERROR: Not implemented!"
// 		<< std::endl;
//       libmesh_error();
    }

#endif

  void_fe_with_vec_switch(inverse_map(elem, physical_points, reference_points, tolerance, secure));

  libmesh_error();
  return;
}



bool FEInterface::on_reference_element(const Point& p,
				       const ElemType t,
				       const Real eps)
{
  return FEBase::on_reference_element(p,t,eps);
}




Real FEInterface::shape(const unsigned int dim,
			const FEType& fe_t,
			const ElemType t,
			const unsigned int i,
			const Point& p)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(t) )
    return ifem_shape(dim, fe_t, t, i, p);

#endif

  const Order o = fe_t.order;

  fe_switch(shape(t,o,i,p));

  libmesh_error();
  return 0.;
}

Real FEInterface::shape(const unsigned int dim,
			const FEType& fe_t,
			const Elem* elem,
			const unsigned int i,
			const Point& p)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(elem->type()) )
    return ifem_shape(dim, fe_t, elem, i, p);

#endif

  const Order o = fe_t.order;

  fe_switch(shape(elem,o,i,p));

  libmesh_error();
  return 0.;
}

template<>
void FEInterface::shape<Real>(const unsigned int dim,
			      const FEType& fe_t,
			      const ElemType t,
			      const unsigned int i,
			      const Point& p,
			      Real& phi)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(t) )
    phi = ifem_shape(dim, fe_t, t, i, p);

#endif

  const Order o = fe_t.order;

  switch(dim)
    {
    case 0:
      fe_scalar_vec_error_switch(0, shape(t,o,i,p), phi = , ; break;);
      break;
    case 1:
      fe_scalar_vec_error_switch(1, shape(t,o,i,p), phi = , ; break;);
      break;
    case 2:
      fe_scalar_vec_error_switch(2, shape(t,o,i,p), phi = , ; break;);
      break;
    case 3:
      fe_scalar_vec_error_switch(3, shape(t,o,i,p), phi = , ; break;);
      break;
    }

  return;
}

template<>
void FEInterface::shape<Real>(const unsigned int dim,
			      const FEType& fe_t,
			      const Elem* elem,
			      const unsigned int i,
			      const Point& p,
			      Real& phi)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(elem->type()) )
    phi = ifem_shape(dim, fe_t, elem, i, p);

#endif

  const Order o = fe_t.order;

  switch(dim)
    {
    case 0:
      fe_scalar_vec_error_switch(0, shape(elem,o,i,p), phi = , ; break;);
      break;
    case 1:
      fe_scalar_vec_error_switch(1, shape(elem,o,i,p), phi = , ; break;);
      break;
    case 2:
      fe_scalar_vec_error_switch(2, shape(elem,o,i,p), phi = , ; break;);
      break;
    case 3:
      fe_scalar_vec_error_switch(3, shape(elem,o,i,p), phi = , ; break;);
      break;
    }

  return;
}

template<>
void FEInterface::shape<RealGradient>(const unsigned int dim,
				      const FEType& fe_t,
				      const ElemType t,
				      const unsigned int i,
				      const Point& p,
				      RealGradient& phi)
{
  const Order o = fe_t.order;

  switch(dim)
    {
    case 0:
      fe_vector_scalar_error_switch(0, shape(t,o,i,p), phi = , ; break;);
      break;
    case 1:
      fe_vector_scalar_error_switch(1, shape(t,o,i,p), phi = , ; break;);
      break;
    case 2:
      fe_vector_scalar_error_switch(2, shape(t,o,i,p), phi = , ; break;);
      break;
    case 3:
      fe_vector_scalar_error_switch(3, shape(t,o,i,p), phi = , ; break;);
      break;
    }

  return;
}

template<>
void FEInterface::shape<RealGradient>(const unsigned int dim,
				      const FEType& fe_t,
				      const Elem* elem,
				      const unsigned int i,
				      const Point& p,
				      RealGradient& phi)
{
  const Order o = fe_t.order;

  switch(dim)
    {
    case 0:
      fe_vector_scalar_error_switch(0, shape(elem,o,i,p), phi = , ; break;);
      break;
    case 1:
      fe_vector_scalar_error_switch(1, shape(elem,o,i,p), phi = , ; break;);
      break;
    case 2:
      fe_vector_scalar_error_switch(2, shape(elem,o,i,p), phi = , ; break;);
      break;
    case 3:
      fe_vector_scalar_error_switch(3, shape(elem,o,i,p), phi = , ; break;);
      break;
    }

  return;
}

void FEInterface::compute_data(const unsigned int dim,
			       const FEType& fe_t,
			       const Elem* elem,
			       FEComputeData& data)
{
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  if ( is_InfFE_elem(elem->type()) )
    {
      data.init();
      ifem_compute_data(dim, fe_t, elem, data);
      return;
    }

#endif

  FEType p_refined = fe_t;
  p_refined.order = static_cast<Order>(p_refined.order + elem->p_level());

  const unsigned int n_dof = n_dofs (dim, p_refined, elem->type());
  const Point&       p     = data.p;
  data.shape.resize(n_dof);

  // set default values for all the output fields
  data.init();

  for (unsigned int n=0; n<n_dof; n++)
      data.shape[n] = shape(dim, p_refined, elem, n, p);

   return;
}



#ifdef LIBMESH_ENABLE_AMR

void FEInterface::compute_constraints (DofConstraints &constraints,
				       DofMap &dof_map,
				       const unsigned int variable_number,
				       const Elem* elem)
{
  libmesh_assert(elem);

  const FEType& fe_t = dof_map.variable_type(variable_number);

  switch (elem->dim())
    {
    case 0:
    case 1:
      {
	// No constraints in 0D/1D.
	return;
      }


    case 2:
      {
	switch (fe_t.family)
	  {
	  case CLOUGH:
	    FE<2,CLOUGH>::compute_constraints (constraints,
					       dof_map,
					       variable_number,
					       elem); return;

	  case HERMITE:
	    FE<2,HERMITE>::compute_constraints (constraints,
					        dof_map,
					        variable_number,
					        elem); return;

	  case LAGRANGE:
	    FE<2,LAGRANGE>::compute_constraints (constraints,
						 dof_map,
						 variable_number,
						 elem); return;

	  case HIERARCHIC:
	    FE<2,HIERARCHIC>::compute_constraints (constraints,
						   dof_map,
						   variable_number,
						   elem); return;

	  case L2_HIERARCHIC:
	    FE<2,L2_HIERARCHIC>::compute_constraints (constraints,
						      dof_map,
						      variable_number,
						      elem); return;

	  case LAGRANGE_VEC:
	    FE<2,LAGRANGE_VEC>::compute_constraints (constraints,
						     dof_map,
						     variable_number,
						     elem); return;


	  default:
	    return;
	  }
      }


    case 3:
      {
	switch (fe_t.family)
	  {
	  case HERMITE:
	    FE<3,HERMITE>::compute_constraints (constraints,
					        dof_map,
					        variable_number,
					        elem); return;

	  case LAGRANGE:
	    FE<3,LAGRANGE>::compute_constraints (constraints,
					         dof_map,
						 variable_number,
						 elem); return;

	  case HIERARCHIC:
	    FE<3,HIERARCHIC>::compute_constraints (constraints,
						   dof_map,
						   variable_number,
						   elem); return;

	  case L2_HIERARCHIC:
	    FE<3,L2_HIERARCHIC>::compute_constraints (constraints,
						      dof_map,
						      variable_number,
						      elem); return;

	  case LAGRANGE_VEC:
	    FE<3,LAGRANGE_VEC>::compute_constraints (constraints,
						     dof_map,
						     variable_number,
						     elem); return;
	  default:
	    return;
	  }
      }


    default:
      libmesh_error();
    }
}

#endif // #ifdef LIBMESH_ENABLE_AMR



#ifdef LIBMESH_ENABLE_PERIODIC

void FEInterface::compute_periodic_constraints (DofConstraints &constraints,
				                DofMap &dof_map,
                                                const PeriodicBoundaries &boundaries,
						const MeshBase &mesh,
                                                const PointLocatorBase* point_locator,
				                const unsigned int variable_number,
				                const Elem* elem)
{
  // No element-specific optimizations currently exist
  FEBase::compute_periodic_constraints (constraints,
                                        dof_map,
                                        boundaries,
                                        mesh,
                                        point_locator,
				        variable_number,
				        elem);
}

#endif // #ifdef LIBMESH_ENABLE_PERIODIC



unsigned int FEInterface::max_order(const FEType& fe_t,
			            const ElemType& el_t)
{
  // Yeah, I know, infinity is much larger than 11, but our
  // solvers don't seem to like high degree polynomials, and our
  // quadrature rules and number_lookups tables
  // need to go up higher.
  const unsigned int unlimited = 11;

  // If we used 0 as a default, then elements missing from this
  // table (e.g. infinite elements) would be considered broken.
  const unsigned int unknown = unlimited;

  switch (fe_t.family)
    {
      case LAGRANGE:
      case L2_LAGRANGE: // TODO: L2_LAGRANGE can have higher "max_order" than LAGRANGE
      case LAGRANGE_VEC:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	    case EDGE4:
	      return 3;
	    case TRI3:
	      return 1;
	    case TRI6:
	      return 2;
	    case QUAD4:
	      return 1;
	    case QUAD8:
	    case QUAD9:
	      return 2;
	    case TET4:
	      return 1;
	    case TET10:
	      return 2;
	    case HEX8:
	      return 1;
	    case HEX20:
	    case HEX27:
	      return 2;
	    case PRISM6:
	    case PRISM15:
	      return 1;
	    case PRISM18:
	      return 2;
	    case PYRAMID5:
	      return 1;
	    default:
	      return unknown;
	  }
	break;
      case MONOMIAL:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	    case EDGE4:
	    case TRI3:
	    case TRI6:
	    case QUAD4:
	    case QUAD8:
	    case QUAD9:
	    case TET4:
	    case TET10:
	    case HEX8:
	    case HEX20:
	    case HEX27:
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return unlimited;
	    default:
	      return unknown;
	  }
	break;
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
      case BERNSTEIN:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	    case EDGE4:
	      return unlimited;
	    case TRI3:
	      return 0;
	    case TRI6:
	      return 6;
	    case QUAD4:
	      return 0;
	    case QUAD8:
	    case QUAD9:
	      return unlimited;
	    case TET4:
	      return 1;
	    case TET10:
	      return 2;
	    case HEX8:
	      return 0;
	    case HEX20:
	      return 2;
	    case HEX27:
	      return 4;
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return 0;
	    default:
	      return unknown;
	  }
	break;
      case SZABAB:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	    case EDGE4:
	      return 7;
	    case TRI3:
	      return 0;
	    case TRI6:
	      return 7;
	    case QUAD4:
	      return 0;
	    case QUAD8:
	    case QUAD9:
	      return 7;
	    case TET4:
	    case TET10:
	    case HEX8:
	    case HEX20:
	    case HEX27:
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return 0;
	    default:
	      return unknown;
	  }
	break;
#endif
      case XYZ:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	    case EDGE4:
	    case TRI3:
	    case TRI6:
	    case QUAD4:
	    case QUAD8:
	    case QUAD9:
	    case TET4:
	    case TET10:
	    case HEX8:
	    case HEX20:
	    case HEX27:
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return unlimited;
	    default:
	      return unknown;
	  }
	break;
      case CLOUGH:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	      return 3;
	    case EDGE4:
	    case TRI3:
	      return 0;
	    case TRI6:
	      return 3;
	    case QUAD4:
	    case QUAD8:
	    case QUAD9:
	    case TET4:
	    case TET10:
	    case HEX8:
	    case HEX20:
	    case HEX27:
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return 0;
	    default:
	      return unknown;
	  }
	break;
      case HERMITE:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	      return unlimited;
	    case EDGE4:
	    case TRI3:
	    case TRI6:
	      return 0;
	    case QUAD4:
	      return 3;
	    case QUAD8:
	    case QUAD9:
	      return unlimited;
	    case TET4:
	    case TET10:
	      return 0;
	    case HEX8:
	      return 3;
	    case HEX20:
	    case HEX27:
	      return unlimited;
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return 0;
	    default:
	      return unknown;
	  }
	break;
      case HIERARCHIC:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	    case EDGE4:
	      return unlimited;
	    case TRI3:
	      return 1;
	    case TRI6:
	      return unlimited;
	    case QUAD4:
	      return 1;
	    case QUAD8:
	    case QUAD9:
	      return unlimited;
	    case TET4:
	    case TET10:
	      return 0;
	    case HEX8:
	    case HEX20:
	      return 1;
	    case HEX27:
	      return unlimited;
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return 0;
	    default:
	      return unknown;
	  }
	break;
      case L2_HIERARCHIC:
	switch (el_t)
	  {
	    case EDGE2:
	    case EDGE3:
	    case EDGE4:
	      return unlimited;
	    case TRI3:
	      return 1;
	    case TRI6:
	      return unlimited;
	    case QUAD4:
	      return 1;
	    case QUAD8:
	    case QUAD9:
	      return unlimited;
	    case TET4:
	    case TET10:
	      return 0;
	    case HEX8:
	    case HEX20:
	      return 1;
	    case HEX27:
	      return unlimited;
	    case PRISM6:
	    case PRISM15:
	    case PRISM18:
	    case PYRAMID5:
	      return 0;
	    default:
	      return unknown;
	  }
	break;
    case NEDELEC_ONE:
      switch (el_t)
	{
	case TRI6:
	case QUAD8:
	case QUAD9:
        case HEX20:
        case HEX27:
	  return 1;
	default:
	  return 0;
	}
      break;
      default:
	return 0;
	break;
    }
}



bool FEInterface::extra_hanging_dofs(const FEType& fe_t)
{
  switch (fe_t.family)
    {
    case LAGRANGE:
    case L2_LAGRANGE:
    case MONOMIAL:
    case L2_HIERARCHIC:
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
    case BERNSTEIN:
    case SZABAB:
#endif
    case XYZ:
    case LAGRANGE_VEC:
    case NEDELEC_ONE:
      return false;
    case CLOUGH:
    case HERMITE:
    case HIERARCHIC:
    default:
      return true;
    }
}

FEFieldType FEInterface::field_type( const FEType& fe_type )
{
  return FEInterface::field_type( fe_type.family );
}

FEFieldType FEInterface::field_type( const FEFamily& fe_family )
{
  switch (fe_family)
    {
    case LAGRANGE_VEC:
    case NEDELEC_ONE:
      return TYPE_VECTOR;
    default:
      return TYPE_SCALAR;
    }
}

unsigned int FEInterface::n_vec_dim( const MeshBase& mesh, const FEType& fe_type )
{
  switch (fe_type.family)
    {
      //FIXME: We currently assume that the number of vector components is tied
      //       to the mesh dimension. This will break for mixed-dimension meshes.
    case LAGRANGE_VEC:
    case NEDELEC_ONE:
      return mesh.mesh_dimension();
    default:
      return 1;
    }
}

} // namespace libMesh
