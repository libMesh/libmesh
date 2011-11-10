// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include <cmath> // for std::sqrt


// Local includes
#include "libmesh_common.h"
#include "fe.h"
#include "quadrature.h"
#include "elem.h"
#include "libmesh_logging.h"

namespace libMesh
{

//-------------------------------------------------------
// Full specializations for useless methods in 0D, 1D
#define REINIT_ERROR(_dim, _type, _func)       \
template <>                                    \
void FE<_dim,_type>::_func(const Elem*,        \
			   const unsigned int, \
                           const Real,         \
                           const std::vector<Point>* const,     \
                           const std::vector<Real>* const)      \
{                                              \
  libMesh::err << "ERROR: This method makes no sense for low-D elements!" \
	        << std::endl;                      \
  libmesh_error();                                     \
}

#define SIDEMAP_ERROR(_dim, _type, _func)       \
template <>                                    \
void FE<_dim,_type>::_func(const Elem*,        \
                           const Elem*,        \
			   const unsigned int, \
                           const std::vector<Point>&,     \
                           std::vector<Point>&)      \
{                                              \
  libMesh::err << "ERROR: This method makes no sense for low-D elements!" \
	        << std::endl;                      \
  libmesh_error();                                     \
}

REINIT_ERROR(0, CLOUGH, reinit)
REINIT_ERROR(0, CLOUGH, edge_reinit)
SIDEMAP_ERROR(0, CLOUGH, side_map)
REINIT_ERROR(0, HERMITE, reinit)
REINIT_ERROR(0, HERMITE, edge_reinit)
SIDEMAP_ERROR(0, HERMITE, side_map)
REINIT_ERROR(0, HIERARCHIC, reinit)
REINIT_ERROR(0, HIERARCHIC, edge_reinit)
SIDEMAP_ERROR(0, HIERARCHIC, side_map)
REINIT_ERROR(0, L2_HIERARCHIC, reinit)
REINIT_ERROR(0, L2_HIERARCHIC, edge_reinit)
SIDEMAP_ERROR(0, L2_HIERARCHIC, side_map)
REINIT_ERROR(0, LAGRANGE, reinit)
REINIT_ERROR(0, LAGRANGE, edge_reinit)
SIDEMAP_ERROR(0, LAGRANGE, side_map)
REINIT_ERROR(0, MONOMIAL, reinit)
REINIT_ERROR(0, MONOMIAL, edge_reinit)
SIDEMAP_ERROR(0, MONOMIAL, side_map)
REINIT_ERROR(0, SCALAR, reinit)
REINIT_ERROR(0, SCALAR, edge_reinit)
SIDEMAP_ERROR(0, SCALAR, side_map)
REINIT_ERROR(0, XYZ, reinit)
REINIT_ERROR(0, XYZ, edge_reinit)
SIDEMAP_ERROR(0, XYZ, side_map)
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
REINIT_ERROR(0, BERNSTEIN, reinit)
REINIT_ERROR(0, BERNSTEIN, edge_reinit)
SIDEMAP_ERROR(0, BERNSTEIN, side_map)
REINIT_ERROR(0, SZABAB, reinit)
REINIT_ERROR(0, SZABAB, edge_reinit)
SIDEMAP_ERROR(0, SZABAB, side_map)
#endif

REINIT_ERROR(1, CLOUGH, edge_reinit)
REINIT_ERROR(1, HERMITE, edge_reinit)
REINIT_ERROR(1, HIERARCHIC, edge_reinit)
REINIT_ERROR(1, L2_HIERARCHIC, edge_reinit)
REINIT_ERROR(1, LAGRANGE, edge_reinit)
REINIT_ERROR(1, XYZ, edge_reinit)
REINIT_ERROR(1, MONOMIAL, edge_reinit)
REINIT_ERROR(1, SCALAR, edge_reinit)
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
REINIT_ERROR(1, BERNSTEIN, edge_reinit)
REINIT_ERROR(1, SZABAB, edge_reinit)
#endif


//-------------------------------------------------------
// Methods for 2D, 3D
template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::reinit(const Elem* elem,
		       const unsigned int s,
		       const Real /* tolerance */,
                       const std::vector<Point>* const pts,
                       const std::vector<Real>* const weights)
{
  libmesh_assert (elem  != NULL);
  libmesh_assert (qrule != NULL || pts != NULL);
  // We now do this for 1D elements!
  // libmesh_assert (Dim != 1);

  // Build the side of interest 
  const AutoPtr<Elem> side(elem->build_side(s));

  // Find the max p_level to select 
  // the right quadrature rule for side integration
  unsigned int side_p_level = elem->p_level();
  if (elem->neighbor(s) != NULL)
    side_p_level = std::max(side_p_level, elem->neighbor(s)->p_level());

  // Initialize the shape functions at the user-specified
  // points
  if (pts != NULL)
    {
      // The shape functions do not correspond to the qrule
      shapes_on_quadrature = false;

      // Initialize the face shape functions
      this->init_face_shape_functions (*pts, side.get());

      // Compute the Jacobian*Weight on the face for integration
      if (weights != NULL)
        {
          this->compute_face_map (*weights, side.get());
        }
      else
        {
          std::vector<Real> dummy_weights (pts->size(), 1.);
          this->compute_face_map (dummy_weights, side.get());
        }
    }
  // If there are no user specified points, we use the
  // quadrature rule
  else
    {
      // initialize quadrature rule
      qrule->init(side->type(), side_p_level);

      if(qrule->shapes_need_reinit())
        shapes_on_quadrature = false;

      // FIXME - could this break if the same FE object was used
      // for both volume and face integrals? - RHS
      // We might not need to reinitialize the shape functions
      if ((this->get_type() != elem->type())    ||
          (side->type() != last_side)           ||
          (this->get_p_level() != side_p_level) ||
          this->shapes_need_reinit()            ||
          !shapes_on_quadrature)
        {
          // Set the element type and p_level
          elem_type = elem->type();

          // Set the last_side
          last_side = side->type();

          // Set the last p level
          _p_level = side_p_level;

          // Initialize the face shape functions
          this->init_face_shape_functions (qrule->get_points(),  side.get());
        }

      // Compute the Jacobian*Weight on the face for integration
      this->compute_face_map (qrule->get_weights(), side.get());

      // The shape functions correspond to the qrule
      shapes_on_quadrature = true;
    }

  // make a copy of the Jacobian for integration
  const std::vector<Real> JxW_int(JxW);

  // make a copy of shape on quadrature info
  bool shapes_on_quadrature_side = shapes_on_quadrature;

  // Find where the integration points are located on the
  // full element.
  const std::vector<Point>* ref_qp;
  if (pts != NULL)
    ref_qp = pts;
  else
    ref_qp = &qrule->get_points();

  std::vector<Point> qp;
  this->side_map(elem, side.get(), s, *ref_qp, qp);

  // compute the shape function and derivative values
  // at the points qp
  this->reinit  (elem, &qp);

  shapes_on_quadrature = shapes_on_quadrature_side;

  // copy back old data
  JxW = JxW_int;
}



template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::edge_reinit(const Elem* elem,
		            const unsigned int e,
			    const Real tolerance,
                            const std::vector<Point>* const pts,
                            const std::vector<Real>* const weights)
{
  libmesh_assert (elem  != NULL);
  libmesh_assert (qrule != NULL || pts != NULL);
  // We don't do this for 1D elements!
  libmesh_assert (Dim != 1);

  // Build the side of interest 
  const AutoPtr<Elem> edge(elem->build_edge(e));

  // Initialize the shape functions at the user-specified
  // points
  if (pts != NULL)
    {
      // The shape functions do not correspond to the qrule
      shapes_on_quadrature = false;

      // Initialize the edge shape functions
      this->init_edge_shape_functions (*pts, edge.get());
  
      // Compute the Jacobian*Weight on the face for integration
      if (weights != NULL)
        {
          this->compute_edge_map (*weights, edge.get());
        }
      else
        {
          std::vector<Real> dummy_weights (pts->size(), 1.);
          this->compute_edge_map (dummy_weights, edge.get());
        }
    }
  // If there are no user specified points, we use the
  // quadrature rule
  else
    {  
      // initialize quadrature rule
      qrule->init(edge->type(), elem->p_level());

      if(qrule->shapes_need_reinit())
        shapes_on_quadrature = false;

      // We might not need to reinitialize the shape functions
      if ((this->get_type() != elem->type())                   ||
          (edge->type() != static_cast<int>(last_edge))        || // Comparison between enum and unsigned, cast the unsigned to int
          this->shapes_need_reinit()                           ||
          !shapes_on_quadrature)
        {
          // Set the element type
          elem_type = elem->type();

          // Set the last_edge
          last_edge = edge->type();
      
          // Initialize the edge shape functions
          this->init_edge_shape_functions (qrule->get_points(), edge.get());
        }
  
      // Compute the Jacobian*Weight on the face for integration
      this->compute_edge_map (qrule->get_weights(), edge.get());

      // The shape functions correspond to the qrule
      shapes_on_quadrature = true;
    }

  // make a copy of the Jacobian for integration
  const std::vector<Real> JxW_int(JxW);

  // Find where the integration points are located on the
  // full element.
  std::vector<Point> qp; this->inverse_map (elem, xyz, qp, tolerance);
  
  // compute the shape function and derivative values
  // at the points qp
  this->reinit  (elem, &qp);

  // copy back old data
  JxW = JxW_int;
}

template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::side_map (const Elem* elem,
	                  const Elem* side,
                          const unsigned int s,
                          const std::vector<Point>& reference_side_points,
	                  std::vector<Point>&       reference_points)
{
  unsigned int side_p_level = elem->p_level();
  if (elem->neighbor(s) != NULL)
    side_p_level = std::max(side_p_level, elem->neighbor(s)->p_level());

  if (side->type() != last_side ||
      side_p_level != _p_level ||
      !shapes_on_quadrature)
    {
      // Set the element type
      elem_type = elem->type();
      _p_level = side_p_level;

      // Set the last_side
      last_side = side->type();

      // Initialize the face shape functions
      this->init_face_shape_functions(reference_side_points, side);
    }
  
  const unsigned int n_points = reference_side_points.size();
  reference_points.resize(n_points);
  for (unsigned int i = 0; i < n_points; i++)
    reference_points[i].zero();
  
  std::vector<unsigned int> elem_nodes_map;
  elem_nodes_map.resize(side->n_nodes());
  for (unsigned int j = 0; j < side->n_nodes(); j++)
    for (unsigned int i = 0; i < elem->n_nodes(); i++)
      if (side->node(j) == elem->node(i))
         elem_nodes_map[j] = i;
  std::vector<Point> refspace_nodes;
  this->get_refspace_nodes(elem->type(), refspace_nodes);

  for (unsigned int i=0; i<psi_map.size(); i++) // sum over the nodes
  {
    const Point& side_node = refspace_nodes[elem_nodes_map[i]]; 
    for (unsigned int p=0; p<n_points; p++)
      reference_points[p].add_scaled (side_node, psi_map[i][p]);
  }
}


template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::init_face_shape_functions(const std::vector<Point>& qp,
					  const Elem* side)
{
  libmesh_assert (side  != NULL);
  
  /**
   * Start logging the shape function initialization
   */
  START_LOG("init_face_shape_functions()", "FE");

  // The element type and order to use in
  // the map
  const Order    mapping_order     (side->default_order()); 
  const ElemType mapping_elem_type (side->type());

  // The number of quadrature points.
  const unsigned int n_qp = qp.size();
	
  const unsigned int n_mapping_shape_functions =
    FE<Dim,LAGRANGE>::n_shape_functions (mapping_elem_type,
					 mapping_order);
  
  // resize the vectors to hold current data
  // Psi are the shape functions used for the FE mapping
  psi_map.resize        (n_mapping_shape_functions);

  if (Dim > 1)
    {
      dpsidxi_map.resize    (n_mapping_shape_functions);
      d2psidxi2_map.resize  (n_mapping_shape_functions);
    }
  
  if (Dim == 3)
    {
      dpsideta_map.resize     (n_mapping_shape_functions);
      d2psidxideta_map.resize (n_mapping_shape_functions);
      d2psideta2_map.resize   (n_mapping_shape_functions);
    }
  
  for (unsigned int i=0; i<n_mapping_shape_functions; i++)
    {
      // Allocate space to store the values of the shape functions
      // and their first and second derivatives at the quadrature points.
      psi_map[i].resize        (n_qp);
      if (Dim > 1)
        {
          dpsidxi_map[i].resize    (n_qp);
          d2psidxi2_map[i].resize  (n_qp);
        }
      if (Dim == 3)
	{
	  dpsideta_map[i].resize     (n_qp);
	  d2psidxideta_map[i].resize (n_qp);
	  d2psideta2_map[i].resize   (n_qp);
	}
  
      // Compute the value of shape function i, and its first and
      // second derivatives at quadrature point p
      // (Lagrange shape functions are used for the mapping)
      for (unsigned int p=0; p<n_qp; p++)
	{
	  psi_map[i][p]        = FE<Dim-1,LAGRANGE>::shape             (mapping_elem_type, mapping_order, i,    qp[p]);
          if (Dim > 1)
	    {
	      dpsidxi_map[i][p]    = FE<Dim-1,LAGRANGE>::shape_deriv       (mapping_elem_type, mapping_order, i, 0, qp[p]);
	      d2psidxi2_map[i][p]  = FE<Dim-1,LAGRANGE>::shape_second_deriv(mapping_elem_type, mapping_order, i, 0, qp[p]);
	    }
	  // libMesh::out << "d2psidxi2_map["<<i<<"][p]=" << d2psidxi2_map[i][p] << std::endl;

	  // If we are in 3D, then our sides are 2D faces.
	  // For the second derivatives, we must also compute the cross
	  // derivative d^2() / dxi deta
	  if (Dim == 3)
	    {
	      dpsideta_map[i][p]     = FE<Dim-1,LAGRANGE>::shape_deriv       (mapping_elem_type, mapping_order, i, 1, qp[p]);
	      d2psidxideta_map[i][p] = FE<Dim-1,LAGRANGE>::shape_second_deriv(mapping_elem_type, mapping_order, i, 1, qp[p]); 
	      d2psideta2_map[i][p]   = FE<Dim-1,LAGRANGE>::shape_second_deriv(mapping_elem_type, mapping_order, i, 2, qp[p]);
	    }
	}
    }

  
  /**
   * Stop logging the shape function initialization
   */
  STOP_LOG("init_face_shape_functions()", "FE");
}


  
template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::init_edge_shape_functions(const std::vector<Point>& qp,
					  const Elem* edge)
{
  libmesh_assert (edge != NULL);
  
  /**
   * Start logging the shape function initialization
   */
  START_LOG("init_edge_shape_functions()", "FE");

  // The element type and order to use in
  // the map
  const Order    mapping_order     (edge->default_order()); 
  const ElemType mapping_elem_type (edge->type());

  // The number of quadrature points.
  const unsigned int n_qp = qp.size();
	
  const unsigned int n_mapping_shape_functions =
    FE<Dim,LAGRANGE>::n_shape_functions (mapping_elem_type,
					 mapping_order);
  
  // resize the vectors to hold current data
  // Psi are the shape functions used for the FE mapping
  psi_map.resize        (n_mapping_shape_functions);
  dpsidxi_map.resize    (n_mapping_shape_functions);
  d2psidxi2_map.resize  (n_mapping_shape_functions);
  
  for (unsigned int i=0; i<n_mapping_shape_functions; i++)
    {
      // Allocate space to store the values of the shape functions
      // and their first and second derivatives at the quadrature points.
      psi_map[i].resize        (n_qp);
      dpsidxi_map[i].resize    (n_qp);
      d2psidxi2_map[i].resize  (n_qp);
  
      // Compute the value of shape function i, and its first and
      // second derivatives at quadrature point p
      // (Lagrange shape functions are used for the mapping)
      for (unsigned int p=0; p<n_qp; p++)
	{
	  psi_map[i][p]        = FE<1,LAGRANGE>::shape             (mapping_elem_type, mapping_order, i,    qp[p]);
	  dpsidxi_map[i][p]    = FE<1,LAGRANGE>::shape_deriv       (mapping_elem_type, mapping_order, i, 0, qp[p]);
	  d2psidxi2_map[i][p]  = FE<1,LAGRANGE>::shape_second_deriv(mapping_elem_type, mapping_order, i, 0, qp[p]);
	}
    }
  
  /**
   * Stop logging the shape function initialization
   */
  STOP_LOG("init_edge_shape_functions()", "FE");
}

  

void FEBase::compute_face_map(const std::vector<Real>& qw,
			      const Elem* side)
{
  libmesh_assert (side  != NULL);

  START_LOG("compute_face_map()", "FE");

  // The number of quadrature points.
  const unsigned int n_qp = qw.size();
  
  
  switch (dim)
    {
    case 1:
      {
	// A 1D finite element, currently assumed to be in 1D space
	// This means the boundary is a "0D finite element", a
        // NODEELEM.

	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize(n_qp);
	  normals.resize(n_qp);

	  JxW.resize(n_qp);
        }

        // If we have no quadrature points, there's nothing else to do
        if (!n_qp)
          break;

        // We need to look back at the full edge to figure out the normal
        // vector
        const Elem *elem = side->parent();
        libmesh_assert (elem);
        if (side->node(0) == elem->node(0))
          normals[0] = Point(-1.);
        else
          {
            libmesh_assert (side->node(0) == elem->node(1));
            normals[0] = Point(1.);
          }

        // Calculate x at the point
	libmesh_assert (psi_map.size() == 1);
        // In the unlikely event we have multiple quadrature
        // points, they'll be in the same place
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    xyz[p].zero();
	    xyz[p].add_scaled          (side->point(0), psi_map[0][p]);
            normals[p] = normals[0];
	    JxW[p] = 1.0*qw[p];
          }

	// done computing the map
	break;
      }
      
    case 2:
      {
	// A 2D finite element living in either 2D or 3D space.
	// This means the boundary is a 1D finite element, i.e.
	// and EDGE2 or EDGE3.
	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize(n_qp);
	  dxyzdxi_map.resize(n_qp);
	  d2xyzdxi2_map.resize(n_qp);
	  tangents.resize(n_qp);
	  normals.resize(n_qp);
	  curvatures.resize(n_qp);
	  
	  JxW.resize(n_qp);
	}
	
	// Clear the entities that will be summed
	// Compute the tangent & normal at the quadrature point
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    tangents[p].resize(LIBMESH_DIM-1); // 1 Tangent in 2D, 2 in 3D
	    xyz[p].zero();
	    dxyzdxi_map[p].zero();
	    d2xyzdxi2_map[p].zero();
	  }
	
	// compute x, dxdxi at the quadrature points    
	for (unsigned int i=0; i<psi_map.size(); i++) // sum over the nodes
	  {
	    const Point& side_point = side->point(i);
	    
	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
	      {	  
		xyz[p].add_scaled          (side_point, psi_map[i][p]);
		dxyzdxi_map[p].add_scaled  (side_point, dpsidxi_map[i][p]);
		d2xyzdxi2_map[p].add_scaled(side_point, d2psidxi2_map[i][p]);
	      }
	  }

	// Compute the tangent & normal at the quadrature point
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    // The first tangent comes from just the edge's Jacobian
	    tangents[p][0] = dxyzdxi_map[p].unit();
	    
#if LIBMESH_DIM == 2
	    // For a 2D element living in 2D, the normal is given directly
	    // from the entries in the edge Jacobian.
	    normals[p] = (Point(dxyzdxi_map[p](1), -dxyzdxi_map[p](0), 0.)).unit();
	    
#elif LIBMESH_DIM == 3
	    // For a 2D element living in 3D, there is a second tangent.
	    // For the second tangent, we need to refer to the full
	    // element's (not just the edge's) Jacobian.
	    const Elem *elem = side->parent();
	    libmesh_assert (elem != NULL);

	    // Inverse map xyz[p] to a reference point on the parent...
	    Point reference_point = FE<2,LAGRANGE>::inverse_map(elem, xyz[p]);
	    
	    // Get dxyz/dxi and dxyz/deta from the parent map.
	    Point dx_dxi  = FE<2,LAGRANGE>::map_xi (elem, reference_point);
	    Point dx_deta = FE<2,LAGRANGE>::map_eta(elem, reference_point);

	    // The second tangent vector is formed by crossing these vectors.
	    tangents[p][1] = dx_dxi.cross(dx_deta).unit();

	    // Finally, the normal in this case is given by crossing these
	    // two tangents.
	    normals[p] = tangents[p][0].cross(tangents[p][1]).unit();
#endif 
	    

	    // The curvature is computed via the familiar Frenet formula:
	    // curvature = [d^2(x) / d (xi)^2] dot [normal]
	    // For a reference, see:
	    // F.S. Merritt, Mathematics Manual, 1962, McGraw-Hill, p. 310
	    //
	    // Note: The sign convention here is different from the
	    // 3D case.  Concave-upward curves (smiles) have a positive
	    // curvature.  Concave-downward curves (frowns) have a
	    // negative curvature.  Be sure to take that into account!
	    const Real numerator   = d2xyzdxi2_map[p] * normals[p];
	    const Real denominator = dxyzdxi_map[p].size_sq();
	    libmesh_assert (denominator != 0);
	    curvatures[p] = numerator / denominator;
	  }
	
	// compute the jacobian at the quadrature points
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    const Real jac = dxyzdxi_map[p].size();
	    
	    libmesh_assert (jac > 0.);
	    
	    JxW[p] = jac*qw[p];
	  }
	
	// done computing the map
	break;
      }


      
    case 3:
      {
	// A 3D finite element living in 3D space.
	// Resize the vectors to hold data at the quadrature points
	{  
	  xyz.resize(n_qp);
	  dxyzdxi_map.resize(n_qp);
	  dxyzdeta_map.resize(n_qp);
	  d2xyzdxi2_map.resize(n_qp);
	  d2xyzdxideta_map.resize(n_qp);
	  d2xyzdeta2_map.resize(n_qp);
	  tangents.resize(n_qp);
	  normals.resize(n_qp);
	  curvatures.resize(n_qp);

	  JxW.resize(n_qp);
	}
    
	// Clear the entities that will be summed
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    tangents[p].resize(LIBMESH_DIM-1); // 1 Tangent in 2D, 2 in 3D
	    xyz[p].zero();
	    dxyzdxi_map[p].zero();
	    dxyzdeta_map[p].zero();
	    d2xyzdxi2_map[p].zero();
	    d2xyzdxideta_map[p].zero();
	    d2xyzdeta2_map[p].zero();
	  }
	
	// compute x, dxdxi at the quadrature points    
	for (unsigned int i=0; i<psi_map.size(); i++) // sum over the nodes
	  {
	    const Point& side_point = side->point(i);
	    
	    for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
	      {
		xyz[p].add_scaled         (side_point, psi_map[i][p]);
		dxyzdxi_map[p].add_scaled (side_point, dpsidxi_map[i][p]);
		dxyzdeta_map[p].add_scaled(side_point, dpsideta_map[i][p]);
		d2xyzdxi2_map[p].add_scaled   (side_point, d2psidxi2_map[i][p]);
		d2xyzdxideta_map[p].add_scaled(side_point, d2psidxideta_map[i][p]);
		d2xyzdeta2_map[p].add_scaled  (side_point, d2psideta2_map[i][p]);
	      }
	  }

	// Compute the tangents, normal, and curvature at the quadrature point
	for (unsigned int p=0; p<n_qp; p++)
	  {	    
	    const Point n  = dxyzdxi_map[p].cross(dxyzdeta_map[p]);
	    normals[p]     = n.unit();
	    tangents[p][0] = dxyzdxi_map[p].unit();
	    tangents[p][1] = n.cross(dxyzdxi_map[p]).unit();
	    
	    // Compute curvature using the typical nomenclature
	    // of the first and second fundamental forms.
	    // For reference, see:
	    // 1) http://mathworld.wolfram.com/MeanCurvature.html
	    //    (note -- they are using inward normal)
	    // 2) F.S. Merritt, Mathematics Manual, 1962, McGraw-Hill
	    const Real L  = -d2xyzdxi2_map[p]    * normals[p];
	    const Real M  = -d2xyzdxideta_map[p] * normals[p];
	    const Real N  = -d2xyzdeta2_map[p]   * normals[p];
	    const Real E  =  dxyzdxi_map[p].size_sq();
	    const Real F  =  dxyzdxi_map[p]      * dxyzdeta_map[p];
	    const Real G  =  dxyzdeta_map[p].size_sq();
	    
	    const Real numerator   = E*N -2.*F*M + G*L;
	    const Real denominator = E*G - F*F;
	    libmesh_assert (denominator != 0.);
	    curvatures[p] = 0.5*numerator/denominator;
	  }  
    	
	// compute the jacobian at the quadrature points, see
	// http://sp81.msi.umn.edu:999/fluent/fidap/help/theory/thtoc.htm
	for (unsigned int p=0; p<n_qp; p++)
	  {
	    const Real g11 = (dxdxi_map(p)*dxdxi_map(p) +
			      dydxi_map(p)*dydxi_map(p) +
			      dzdxi_map(p)*dzdxi_map(p));
	    
	    const Real g12 = (dxdxi_map(p)*dxdeta_map(p) +
			      dydxi_map(p)*dydeta_map(p) +
			      dzdxi_map(p)*dzdeta_map(p));
	    
	    const Real g21 = g12;
	    
	    const Real g22 = (dxdeta_map(p)*dxdeta_map(p) +
			      dydeta_map(p)*dydeta_map(p) +
			      dzdeta_map(p)*dzdeta_map(p));
	    
	    
	    const Real jac = std::sqrt(g11*g22 - g12*g21);
	    
	    libmesh_assert (jac > 0.);

	    JxW[p] = jac*qw[p];
	  }
	
	// done computing the map
	break;
      }


    default:
      libmesh_error();
      
    }
  STOP_LOG("compute_face_map()", "FE");
}




void FEBase::compute_edge_map(const std::vector<Real>& qw,
			      const Elem* edge)
{
  libmesh_assert (edge != NULL);

  if (dim == 2)
    {
      // A 2D finite element living in either 2D or 3D space.
      // The edges here are the sides of the element, so the
      // (misnamed) compute_face_map function does what we want
      FEBase::compute_face_map(qw, edge);
      return;
    }

  libmesh_assert (dim == 3);  // 1D is unnecessary and currently unsupported

  START_LOG("compute_edge_map()", "FE");

  // The number of quadrature points.
  const unsigned int n_qp = qw.size();
  
  // Resize the vectors to hold data at the quadrature points
  xyz.resize(n_qp);
  dxyzdxi_map.resize(n_qp);
  dxyzdeta_map.resize(n_qp);
  d2xyzdxi2_map.resize(n_qp);
  d2xyzdxideta_map.resize(n_qp);
  d2xyzdeta2_map.resize(n_qp);
  tangents.resize(n_qp);
  normals.resize(n_qp);
  curvatures.resize(n_qp);

  JxW.resize(n_qp);
    
  // Clear the entities that will be summed
  for (unsigned int p=0; p<n_qp; p++)
    {
      tangents[p].resize(1);
      xyz[p].zero();
      dxyzdxi_map[p].zero();
      dxyzdeta_map[p].zero();
      d2xyzdxi2_map[p].zero();
      d2xyzdxideta_map[p].zero();
      d2xyzdeta2_map[p].zero();
    }

  // compute x, dxdxi at the quadrature points    
  for (unsigned int i=0; i<psi_map.size(); i++) // sum over the nodes
    {
      const Point& edge_point = edge->point(i);
      
      for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
        {
	  xyz[p].add_scaled             (edge_point, psi_map[i][p]);
	  dxyzdxi_map[p].add_scaled     (edge_point, dpsidxi_map[i][p]);
	  d2xyzdxi2_map[p].add_scaled   (edge_point, d2psidxi2_map[i][p]);
        }
    }

  // Compute the tangents at the quadrature point
  // FIXME: normals (plural!) and curvatures are uncalculated
  for (unsigned int p=0; p<n_qp; p++)
    {    
      const Point n  = dxyzdxi_map[p].cross(dxyzdeta_map[p]);
      tangents[p][0] = dxyzdxi_map[p].unit();

      // compute the jacobian at the quadrature points
      const Real jac = std::sqrt(dxdxi_map(p)*dxdxi_map(p) +
				 dydxi_map(p)*dydxi_map(p) +
				 dzdxi_map(p)*dzdxi_map(p));
	    
      libmesh_assert (jac > 0.);

      JxW[p] = jac*qw[p];
    }

  STOP_LOG("compute_edge_map()", "FE");
}




//--------------------------------------------------------------
// Explicit instantiations
template void FE<1,LAGRANGE>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<1,LAGRANGE>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<1,HIERARCHIC>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<1,HIERARCHIC>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<1,L2_HIERARCHIC>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<1,L2_HIERARCHIC>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<1,CLOUGH>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<1,CLOUGH>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<1,HERMITE>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<1,HERMITE>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<1,MONOMIAL>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<1,MONOMIAL>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<1,SCALAR>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<1,SCALAR>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
template void FE<1,BERNSTEIN>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<1,BERNSTEIN>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<1,SZABAB>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<1,SZABAB>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
#endif
template void FE<1,XYZ>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<1,XYZ>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);

template void FE<2,LAGRANGE>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,LAGRANGE>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<2,LAGRANGE>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,HIERARCHIC>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,HIERARCHIC>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<2,HIERARCHIC>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,L2_HIERARCHIC>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,L2_HIERARCHIC>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<2,L2_HIERARCHIC>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,CLOUGH>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,CLOUGH>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<2,CLOUGH>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,HERMITE>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,HERMITE>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<2,HERMITE>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,MONOMIAL>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,MONOMIAL>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<2,MONOMIAL>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,SCALAR>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,SCALAR>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<2,SCALAR>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
template void FE<2,BERNSTEIN>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,BERNSTEIN>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<2,BERNSTEIN>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,SZABAB>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,SZABAB>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<2,SZABAB>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
#endif
template void FE<2,XYZ>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<2,XYZ>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<2,XYZ>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);

// Intel 9.1 complained it needed this in devel mode.
template void FE<2,XYZ>::init_face_shape_functions(const std::vector<Point>&, const Elem*);

template void FE<3,LAGRANGE>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,LAGRANGE>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<3,LAGRANGE>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,HIERARCHIC>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,HIERARCHIC>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<3,HIERARCHIC>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,L2_HIERARCHIC>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,L2_HIERARCHIC>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<3,L2_HIERARCHIC>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,CLOUGH>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,CLOUGH>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<3,CLOUGH>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,HERMITE>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,HERMITE>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<3,HERMITE>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,MONOMIAL>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,MONOMIAL>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<3,MONOMIAL>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,SCALAR>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,SCALAR>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<3,SCALAR>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
template void FE<3,BERNSTEIN>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,BERNSTEIN>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<3,BERNSTEIN>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,SZABAB>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,SZABAB>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<3,SZABAB>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
#endif
template void FE<3,XYZ>::reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);
template void FE<3,XYZ>::side_map(Elem const*, Elem const*, const unsigned int, const std::vector<Point>&, std::vector<Point>&);
template void FE<3,XYZ>::edge_reinit(Elem const*, unsigned int, Real, const std::vector<Point>* const, const std::vector<Real>* const);

// Intel 9.1 complained it needed this in devel mode.
template void FE<3,XYZ>::init_face_shape_functions(const std::vector<Point>&, const Elem*);

} // namespace libMesh
