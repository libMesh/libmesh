// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for std::sqrt


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/tensor_value.h"  // May be necessary if destructors
// get instantiated here

namespace libMesh
{

//-------------------------------------------------------
// Full specializations for useless methods in 0D, 1D
#define REINIT_ERROR(_dim, _type, _func)                        \
  template <>                                                   \
  void FE<_dim,_type>::_func(const Elem *,                      \
                             const unsigned int,                \
                             const Real,                        \
                             const std::vector<Point> * const,  \
                             const std::vector<Real> * const)

#define SIDEMAP_ERROR(_dim, _type, _func)                       \
  template <>                                                   \
  void FE<_dim,_type>::_func(const Elem *,                      \
                             const Elem *,                      \
                             const unsigned int,                \
                             const std::vector<Point> &,        \
                             std::vector<Point> &)

#define FACE_EDGE_SHAPE_ERROR(_dim, _func)                              \
  template <>                                                           \
  void FEMap::_func<_dim>(const std::vector<Point> &,                   \
                          const Elem *)                                 \
  {                                                                     \
    libmesh_error_msg("ERROR: This method makes no sense for low-D elements!"); \
  }

// 0D error instantiations
#define ERRORS_IN_0D(_type) \
REINIT_ERROR(0, _type, reinit)      { libmesh_error_msg("ERROR: Cannot reinit 0D " #_type " elements!"); } \
REINIT_ERROR(0, _type, edge_reinit) { libmesh_error_msg("ERROR: Cannot edge_reinit 0D " #_type " elements!"); } \
SIDEMAP_ERROR(0, _type, side_map)   { libmesh_error_msg("ERROR: Cannot side_map 0D " #_type " elements!"); }

ERRORS_IN_0D(CLOUGH)
ERRORS_IN_0D(HERMITE)
ERRORS_IN_0D(HIERARCHIC)
ERRORS_IN_0D(L2_HIERARCHIC)
ERRORS_IN_0D(LAGRANGE)
ERRORS_IN_0D(L2_LAGRANGE)
ERRORS_IN_0D(LAGRANGE_VEC)
ERRORS_IN_0D(MONOMIAL)
ERRORS_IN_0D(MONOMIAL_VEC)
ERRORS_IN_0D(NEDELEC_ONE)
ERRORS_IN_0D(SCALAR)
ERRORS_IN_0D(XYZ)

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
ERRORS_IN_0D(BERNSTEIN)
ERRORS_IN_0D(SZABAB)
ERRORS_IN_0D(RATIONAL_BERNSTEIN)
#endif

// 1D error instantiations
REINIT_ERROR(1, CLOUGH, edge_reinit)        { libmesh_error_msg("ERROR: Cannot edge_reinit 1D CLOUGH elements!"); }
REINIT_ERROR(1, HERMITE, edge_reinit)       { libmesh_error_msg("ERROR: Cannot edge_reinit 1D HERMITE elements!"); }
REINIT_ERROR(1, HIERARCHIC, edge_reinit)    { libmesh_error_msg("ERROR: Cannot edge_reinit 1D HIERARCHIC elements!"); }
REINIT_ERROR(1, L2_HIERARCHIC, edge_reinit) { libmesh_error_msg("ERROR: Cannot edge_reinit 1D L2_HIERARCHIC elements!"); }
REINIT_ERROR(1, LAGRANGE, edge_reinit)      { libmesh_error_msg("ERROR: Cannot edge_reinit 1D LAGRANGE elements!"); }
REINIT_ERROR(1, LAGRANGE_VEC, edge_reinit)  { libmesh_error_msg("ERROR: Cannot edge_reinit 1D LAGRANGE_VEC elements!"); }
REINIT_ERROR(1, L2_LAGRANGE, edge_reinit)   { libmesh_error_msg("ERROR: Cannot edge_reinit 1D L2_LAGRANGE elements!"); }
REINIT_ERROR(1, XYZ, edge_reinit)           { libmesh_error_msg("ERROR: Cannot edge_reinit 1D XYZ elements!"); }
REINIT_ERROR(1, MONOMIAL, edge_reinit)      { libmesh_error_msg("ERROR: Cannot edge_reinit 1D MONOMIAL elements!"); }
REINIT_ERROR(1, MONOMIAL_VEC, edge_reinit)      { libmesh_error_msg("ERROR: Cannot edge_reinit 1D MONOMIAL_VEC elements!"); }
REINIT_ERROR(1, SCALAR, edge_reinit)        { libmesh_error_msg("ERROR: Cannot edge_reinit 1D SCALAR elements!"); }
REINIT_ERROR(1, NEDELEC_ONE, reinit)        { libmesh_error_msg("ERROR: Cannot reinit 1D NEDELEC_ONE elements!"); }
REINIT_ERROR(1, NEDELEC_ONE, edge_reinit)   { libmesh_error_msg("ERROR: Cannot edge_reinit 1D NEDELEC_ONE elements!"); }
SIDEMAP_ERROR(1, NEDELEC_ONE, side_map)     { libmesh_error_msg("ERROR: Cannot side_map 1D NEDELEC_ONE elements!"); }
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
REINIT_ERROR(1, BERNSTEIN, edge_reinit)     { libmesh_error_msg("ERROR: Cannot edge_reinit 1D BERNSTEIN elements!"); }
REINIT_ERROR(1, SZABAB, edge_reinit)        { libmesh_error_msg("ERROR: Cannot edge_reinit 1D SZABAB elements!"); }
REINIT_ERROR(1, RATIONAL_BERNSTEIN, edge_reinit) { libmesh_error_msg("ERROR: Cannot edge_reinit 1D RATIONAL_BERNSTEIN elements!"); }
#endif


//-------------------------------------------------------
// Methods for 2D, 3D
template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::reinit(const Elem * elem,
                       const unsigned int s,
                       const Real /* tolerance */,
                       const std::vector<Point> * const pts,
                       const std::vector<Real> * const weights)
{
  libmesh_assert(elem);
  libmesh_assert (this->qrule != nullptr || pts != nullptr);
  // We now do this for 1D elements!
  // libmesh_assert_not_equal_to (Dim, 1);

  // If we called this function redundantly (e.g. in an FEMContext
  // that is asked not to do any side calculations) then let's skip the
  // whole inverse_map process that calculates side points
  if (this->calculating_nothing())
    {
      this->calculations_started = true; // Ironic
      return;
    }

  // We're (possibly re-) calculating now!  FIXME - we currently
  // expect to be able to use side_map and JxW later, but we could
  // optimize further here.
  this->_fe_map->add_calculations();
  this->_fe_map->get_JxW();
  this->_fe_map->get_xyz();
  this->determine_calculations();

  // Build the side of interest
  const std::unique_ptr<const Elem> side(elem->build_side_ptr(s));

  // Find the max p_level to select
  // the right quadrature rule for side integration
  unsigned int side_p_level = elem->p_level();
  if (elem->neighbor_ptr(s) != nullptr)
    side_p_level = std::max(side_p_level, elem->neighbor_ptr(s)->p_level());

  // Initialize the shape functions at the user-specified
  // points
  if (pts != nullptr)
    {
      // The shape functions do not correspond to the qrule
      this->shapes_on_quadrature = false;

      // Initialize the face shape functions
      this->_fe_map->template init_face_shape_functions<Dim>(*pts, side.get());

      // Compute the Jacobian*Weight on the face for integration
      if (weights != nullptr)
        {
          this->_fe_map->compute_face_map (Dim, *weights, side.get());
        }
      else
        {
          std::vector<Real> dummy_weights (pts->size(), 1.);
          this->_fe_map->compute_face_map (Dim, dummy_weights, side.get());
        }
    }
  // If there are no user specified points, we use the
  // quadrature rule
  else
    {
      // initialize quadrature rule
      this->qrule->init(side->type(), side_p_level);

      if (this->qrule->shapes_need_reinit())
        this->shapes_on_quadrature = false;

      // FIXME - could this break if the same FE object was used
      // for both volume and face integrals? - RHS
      // We might not need to reinitialize the shape functions
      if ((this->get_type() != elem->type())    ||
          (side->type() != last_side)           ||
          (this->get_p_level() != side_p_level) ||
          this->shapes_need_reinit()            ||
          !this->shapes_on_quadrature)
        {
          // Set the element type and p_level
          this->elem_type = elem->type();

          // Set the last_side
          last_side = side->type();

          // Set the last p level
          this->_p_level = side_p_level;

          // Initialize the face shape functions
          this->_fe_map->template init_face_shape_functions<Dim>(this->qrule->get_points(),  side.get());
        }

      // Compute the Jacobian*Weight on the face for integration
      this->_fe_map->compute_face_map (Dim, this->qrule->get_weights(), side.get());

      // The shape functions correspond to the qrule
      this->shapes_on_quadrature = true;
    }

  // make a copy of the Jacobian for integration
  const std::vector<Real> JxW_int(this->_fe_map->get_JxW());

  // make a copy of shape on quadrature info
  bool shapes_on_quadrature_side = this->shapes_on_quadrature;

  // Find where the integration points are located on the
  // full element.
  const std::vector<Point> * ref_qp;
  if (pts != nullptr)
    ref_qp = pts;
  else
    ref_qp = &this->qrule->get_points();

  std::vector<Point> qp;
  this->side_map(elem, side.get(), s, *ref_qp, qp);

  // compute the shape function and derivative values
  // at the points qp
  this->reinit  (elem, &qp);

  this->shapes_on_quadrature = shapes_on_quadrature_side;

  // copy back old data
  this->_fe_map->get_JxW() = JxW_int;
}



template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::edge_reinit(const Elem * elem,
                            const unsigned int e,
                            const Real tolerance,
                            const std::vector<Point> * const pts,
                            const std::vector<Real> * const weights)
{
  libmesh_assert(elem);
  libmesh_assert (this->qrule != nullptr || pts != nullptr);
  // We don't do this for 1D elements!
  libmesh_assert_not_equal_to (Dim, 1);

  // We're (possibly re-) calculating now!  Time to determine what.
  // FIXME - we currently just assume that we're using JxW and calling
  // edge_map later.
  this->_fe_map->add_calculations();
  this->_fe_map->get_JxW();
  this->_fe_map->get_xyz();
  this->determine_calculations();

  // Build the side of interest
  const std::unique_ptr<const Elem> edge(elem->build_edge_ptr(e));

  // Initialize the shape functions at the user-specified
  // points
  if (pts != nullptr)
    {
      // The shape functions do not correspond to the qrule
      this->shapes_on_quadrature = false;

      // Initialize the edge shape functions
      this->_fe_map->template init_edge_shape_functions<Dim> (*pts, edge.get());

      // Compute the Jacobian*Weight on the face for integration
      if (weights != nullptr)
        {
          this->_fe_map->compute_edge_map (Dim, *weights, edge.get());
        }
      else
        {
          std::vector<Real> dummy_weights (pts->size(), 1.);
          this->_fe_map->compute_edge_map (Dim, dummy_weights, edge.get());
        }
    }
  // If there are no user specified points, we use the
  // quadrature rule
  else
    {
      // initialize quadrature rule
      this->qrule->init(edge->type(), elem->p_level());

      if (this->qrule->shapes_need_reinit())
        this->shapes_on_quadrature = false;

      // We might not need to reinitialize the shape functions
      if ((this->get_type() != elem->type())                   ||
          (edge->type() != static_cast<int>(last_edge))        || // Comparison between enum and unsigned, cast the unsigned to int
          this->shapes_need_reinit()                           ||
          !this->shapes_on_quadrature)
        {
          // Set the element type
          this->elem_type = elem->type();

          // Set the last_edge
          last_edge = edge->type();

          // Initialize the edge shape functions
          this->_fe_map->template init_edge_shape_functions<Dim> (this->qrule->get_points(), edge.get());
        }

      // Compute the Jacobian*Weight on the face for integration
      this->_fe_map->compute_edge_map (Dim, this->qrule->get_weights(), edge.get());

      // The shape functions correspond to the qrule
      this->shapes_on_quadrature = true;
    }

  // make a copy of the Jacobian for integration
  const std::vector<Real> JxW_int(this->_fe_map->get_JxW());

  // Find where the integration points are located on the
  // full element.
  std::vector<Point> qp;
  FEMap::inverse_map (Dim, elem, this->_fe_map->get_xyz(), qp, tolerance);

  // compute the shape function and derivative values
  // at the points qp
  this->reinit  (elem, &qp);

  // copy back old data
  this->_fe_map->get_JxW() = JxW_int;
}

template <unsigned int Dim, FEFamily T>
void FE<Dim,T>::side_map (const Elem * elem,
                          const Elem * side,
                          const unsigned int s,
                          const std::vector<Point> & reference_side_points,
                          std::vector<Point> &       reference_points)
{
  // We're calculating mappings - we need at least first order info
  this->calculate_phi = true;
  this->determine_calculations();

  unsigned int side_p_level = elem->p_level();
  if (elem->neighbor_ptr(s) != nullptr)
    side_p_level = std::max(side_p_level, elem->neighbor_ptr(s)->p_level());

  if (side->type() != last_side ||
      side_p_level != this->_p_level ||
      !this->shapes_on_quadrature)
    {
      // Set the element type
      this->elem_type = elem->type();
      this->_p_level = side_p_level;

      // Set the last_side
      last_side = side->type();

      // Initialize the face shape functions
      this->_fe_map->template init_face_shape_functions<Dim>(reference_side_points, side);
    }

  const unsigned int n_points =
    cast_int<unsigned int>(reference_side_points.size());
  reference_points.resize(n_points);
  for (unsigned int i = 0; i < n_points; i++)
    reference_points[i].zero();

  std::vector<unsigned int> elem_nodes_map;
  elem_nodes_map.resize(side->n_nodes());
  for (auto j : side->node_index_range())
    for (auto i : elem->node_index_range())
      if (side->node_id(j) == elem->node_id(i))
        elem_nodes_map[j] = i;
  std::vector<Point> refspace_nodes;
  this->get_refspace_nodes(elem->type(), refspace_nodes);

  const std::vector<std::vector<Real>> & psi_map = this->_fe_map->get_psi();

  // sum over the nodes
  for (auto i : index_range(psi_map))
    {
      const Point & side_node = refspace_nodes[elem_nodes_map[i]];
      for (unsigned int p=0; p<n_points; p++)
        reference_points[p].add_scaled (side_node, psi_map[i][p]);
    }
}

template<unsigned int Dim>
void FEMap::init_face_shape_functions(const std::vector<Point> & qp,
                                      const Elem * side)
{
  // Start logging the shape function initialization
  LOG_SCOPE("init_face_shape_functions()", "FEMap");

  libmesh_assert(side);

  // We're calculating now!
  this->determine_calculations();

  // The element type and order to use in
  // the map
  const FEFamily mapping_family = FEMap::map_fe_type(*side);
  const ElemType mapping_elem_type (side->type());
  const FEType map_fe_type(side->default_order(), mapping_family);

  // The number of quadrature points.
  const unsigned int n_qp = cast_int<unsigned int>(qp.size());

  // Do not use the p_level(), if any, that is inherited by the side.
  const unsigned int n_mapping_shape_functions =
    FEInterface::n_shape_functions(map_fe_type, /*extra_order=*/0, side);

  // resize the vectors to hold current data
  // Psi are the shape functions used for the FE mapping
  if (calculate_xyz)
    this->psi_map.resize        (n_mapping_shape_functions);

  if (Dim > 1)
    {
      if (calculate_dxyz)
        this->dpsidxi_map.resize    (n_mapping_shape_functions);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      if (calculate_d2xyz)
        this->d2psidxi2_map.resize  (n_mapping_shape_functions);
#endif
    }

  if (Dim == 3)
    {
      if (calculate_dxyz)
        this->dpsideta_map.resize     (n_mapping_shape_functions);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      if (calculate_d2xyz)
        {
          this->d2psidxideta_map.resize (n_mapping_shape_functions);
          this->d2psideta2_map.resize   (n_mapping_shape_functions);
        }
#endif
    }

  FEInterface::shape_ptr shape_ptr =
    FEInterface::shape_function(map_fe_type, side);

  FEInterface::shape_deriv_ptr shape_deriv_ptr =
    FEInterface::shape_deriv_function(map_fe_type, side);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  FEInterface::shape_second_deriv_ptr shape_second_deriv_ptr =
    FEInterface::shape_second_deriv_function(Dim-1, map_fe_type, mapping_elem_type);
#endif

  for (unsigned int i=0; i<n_mapping_shape_functions; i++)
    {
      // Allocate space to store the values of the shape functions
      // and their first and second derivatives at the quadrature points.
      if (calculate_xyz)
        this->psi_map[i].resize        (n_qp);
      if (Dim > 1)
        {
          if (calculate_dxyz)
            this->dpsidxi_map[i].resize    (n_qp);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          if (calculate_d2xyz)
            this->d2psidxi2_map[i].resize  (n_qp);
#endif
        }
      if (Dim == 3)
        {
          if (calculate_dxyz)
            this->dpsideta_map[i].resize     (n_qp);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          if (calculate_d2xyz)
            {
              this->d2psidxideta_map[i].resize (n_qp);
              this->d2psideta2_map[i].resize   (n_qp);
            }
#endif
        }


      // Compute the value of mapping shape function i, and its first
      // and second derivatives at quadrature point p
      for (unsigned int p=0; p<n_qp; p++)
        {
          if (calculate_xyz)
            this->psi_map[i][p]            = shape_ptr             (map_fe_type, side, i,    qp[p], false);
          if (Dim > 1)
            {
              if (calculate_dxyz)
                this->dpsidxi_map[i][p]    = shape_deriv_ptr       (map_fe_type, side, i, 0, qp[p], false);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
              if (calculate_d2xyz)
                this->d2psidxi2_map[i][p]  = shape_second_deriv_ptr(map_fe_type, side, i, 0, qp[p], false);
#endif
            }
          // libMesh::out << "this->d2psidxi2_map["<<i<<"][p]=" << d2psidxi2_map[i][p] << std::endl;

          // If we are in 3D, then our sides are 2D faces.
          // For the second derivatives, we must also compute the cross
          // derivative d^2() / dxi deta
          if (Dim == 3)
            {
              if (calculate_dxyz)
                this->dpsideta_map[i][p]       = shape_deriv_ptr       (map_fe_type, side, i, 1, qp[p], false);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
              if (calculate_d2xyz)
                {
                  this->d2psidxideta_map[i][p] = shape_second_deriv_ptr(map_fe_type, side, i, 1, qp[p], false);
                  this->d2psideta2_map[i][p]   = shape_second_deriv_ptr(map_fe_type, side, i, 2, qp[p], false);
                }
#endif
            }
        }
    }
}

template<unsigned int Dim>
void FEMap::init_edge_shape_functions(const std::vector<Point> & qp,
                                      const Elem * edge)
{
  // Start logging the shape function initialization
  LOG_SCOPE("init_edge_shape_functions()", "FEMap");

  libmesh_assert(edge);

  // We're calculating now!
  this->determine_calculations();

  // The element type and order to use in
  // the map
  const FEFamily mapping_family = FEMap::map_fe_type(*edge);
  const ElemType mapping_elem_type (edge->type());
  const FEType map_fe_type(edge->default_order(), mapping_family);

  // The number of quadrature points.
  const unsigned int n_qp = cast_int<unsigned int>(qp.size());

  // Do not use the p_level(), if any, that is inherited by the side.
  const unsigned int n_mapping_shape_functions =
    FEInterface::n_shape_functions(map_fe_type, /*extra_order=*/0, edge);

  // resize the vectors to hold current data
  // Psi are the shape functions used for the FE mapping
  if (calculate_xyz)
    this->psi_map.resize        (n_mapping_shape_functions);
  if (calculate_dxyz)
    this->dpsidxi_map.resize    (n_mapping_shape_functions);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  if (calculate_d2xyz)
    this->d2psidxi2_map.resize  (n_mapping_shape_functions);
#endif

  FEInterface::shape_ptr shape_ptr =
    FEInterface::shape_function(map_fe_type, edge);

  FEInterface::shape_deriv_ptr shape_deriv_ptr =
    FEInterface::shape_deriv_function(map_fe_type, edge);

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  FEInterface::shape_second_deriv_ptr shape_second_deriv_ptr =
    FEInterface::shape_second_deriv_function(1, map_fe_type, mapping_elem_type);
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

  for (unsigned int i=0; i<n_mapping_shape_functions; i++)
    {
      // Allocate space to store the values of the shape functions
      // and their first and second derivatives at the quadrature points.
      if (calculate_xyz)
        this->psi_map[i].resize        (n_qp);
      if (calculate_dxyz)
        this->dpsidxi_map[i].resize    (n_qp);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      if (calculate_d2xyz)
        this->d2psidxi2_map[i].resize  (n_qp);
#endif

      // Compute the value of mapping shape function i, and its first
      // and second derivatives at quadrature point p
      for (unsigned int p=0; p<n_qp; p++)
        {
          if (calculate_xyz)
            this->psi_map[i][p]        = shape_ptr             (map_fe_type, edge, i,    qp[p], false);
          if (calculate_dxyz)
            this->dpsidxi_map[i][p]    = shape_deriv_ptr       (map_fe_type, edge, i, 0, qp[p], false);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          if (calculate_d2xyz)
            this->d2psidxi2_map[i][p]  = shape_second_deriv_ptr(map_fe_type, edge, i, 0, qp[p], false);
#endif
        }
    }
}



void FEMap::compute_face_map(int dim, const std::vector<Real> & qw,
                             const Elem * side)
{
  libmesh_assert(side);

  // We're calculating now!
  this->determine_calculations();

  LOG_SCOPE("compute_face_map()", "FEMap");

  // The number of quadrature points.
  const unsigned int n_qp = cast_int<unsigned int>(qw.size());

  const FEFamily mapping_family = FEMap::map_fe_type(*side);
  const Order    mapping_order     (side->default_order());
  const FEType map_fe_type(mapping_order, mapping_family);

  // Do not use the p_level(), if any, that is inherited by the side.
  const unsigned int n_mapping_shape_functions =
    FEInterface::n_shape_functions(map_fe_type, /*extra_order=*/0, side);

  switch (dim)
    {
    case 1:
      {
        // A 1D finite element, currently assumed to be in 1D space
        // This means the boundary is a "0D finite element", a
        // NODEELEM.

        // Resize the vectors to hold data at the quadrature points
        {
          if (calculate_xyz)
            this->xyz.resize(n_qp);
          if (calculate_dxyz)
            normals.resize(n_qp);

          if (calculate_dxyz)
            this->JxW.resize(n_qp);
        }

        // If we have no quadrature points, there's nothing else to do
        if (!n_qp)
          break;

        // We need to look back at the full edge to figure out the normal
        // vector
        const Elem * elem = side->parent();
        libmesh_assert (elem);
        if (calculate_dxyz)
          {
            if (side->node_id(0) == elem->node_id(0))
              normals[0] = Point(-1.);
            else
              {
                libmesh_assert_equal_to (side->node_id(0),
                                         elem->node_id(1));
                normals[0] = Point(1.);
              }
          }

        // Calculate x at the point
        if (calculate_xyz)
          libmesh_assert_equal_to (this->psi_map.size(), 1);
        // In the unlikely event we have multiple quadrature
        // points, they'll be in the same place
        for (unsigned int p=0; p<n_qp; p++)
          {
            if (calculate_xyz)
              {
                this->xyz[p].zero();
                this->xyz[p].add_scaled          (side->point(0), this->psi_map[0][p]);
              }
            if (calculate_dxyz)
              {
                normals[p] = normals[0];
                this->JxW[p] = 1.0*qw[p];
              }
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
          if (calculate_xyz)
            this->xyz.resize(n_qp);
          if (calculate_dxyz)
            {
              this->dxyzdxi_map.resize(n_qp);
              this->tangents.resize(n_qp);
              this->normals.resize(n_qp);

              this->JxW.resize(n_qp);
            }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          if (calculate_d2xyz)
            {
              this->d2xyzdxi2_map.resize(n_qp);
              this->curvatures.resize(n_qp);
            }
#endif
        }

        // Clear the entities that will be summed
        // Compute the tangent & normal at the quadrature point
        for (unsigned int p=0; p<n_qp; p++)
          {
            if (calculate_xyz)
              this->xyz[p].zero();
            if (calculate_dxyz)
              {
                this->tangents[p].resize(LIBMESH_DIM-1); // 1 Tangent in 2D, 2 in 3D
                this->dxyzdxi_map[p].zero();
              }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
            if (calculate_d2xyz)
              this->d2xyzdxi2_map[p].zero();
#endif
          }

        // compute x, dxdxi at the quadrature points
        for (unsigned int i=0; i<n_mapping_shape_functions; i++) // sum over the nodes
          {
            const Point & side_point = side->point(i);

            for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
              {
                if (calculate_xyz)
                  this->xyz[p].add_scaled          (side_point, this->psi_map[i][p]);
                if (calculate_dxyz)
                  this->dxyzdxi_map[p].add_scaled  (side_point, this->dpsidxi_map[i][p]);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                if (calculate_d2xyz)
                  this->d2xyzdxi2_map[p].add_scaled(side_point, this->d2psidxi2_map[i][p]);
#endif
              }
          }

        // Compute the tangent & normal at the quadrature point
        if (calculate_dxyz)
          {
            for (unsigned int p=0; p<n_qp; p++)
              {
                // The first tangent comes from just the edge's Jacobian
                this->tangents[p][0] = this->dxyzdxi_map[p].unit();

#if LIBMESH_DIM == 2
                // For a 2D element living in 2D, the normal is given directly
                // from the entries in the edge Jacobian.
                this->normals[p] = (Point(this->dxyzdxi_map[p](1), -this->dxyzdxi_map[p](0), 0.)).unit();

#elif LIBMESH_DIM == 3
                // For a 2D element living in 3D, there is a second tangent.
                // For the second tangent, we need to refer to the full
                // element's (not just the edge's) Jacobian.
                const Elem * elem = side->parent();
                libmesh_assert(elem);

                // Inverse map xyz[p] to a reference point on the parent...
                Point reference_point = FEMap::inverse_map(2, elem, this->xyz[p]);

                // Get dxyz/dxi and dxyz/deta from the parent map.
                Point dx_dxi  = FEMap::map_deriv (2, elem, 0, reference_point);
                Point dx_deta = FEMap::map_deriv (2, elem, 1, reference_point);

                // The second tangent vector is formed by crossing these vectors.
                tangents[p][1] = dx_dxi.cross(dx_deta).unit();

                // Finally, the normal in this case is given by crossing these
                // two tangents.
                normals[p] = tangents[p][0].cross(tangents[p][1]).unit();
#endif


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                // The curvature is computed via the familiar Frenet formula:
                // curvature = [d^2(x) / d (xi)^2] dot [normal]
                // For a reference, see:
                // F.S. Merritt, Mathematics Manual, 1962, McGraw-Hill, p. 310
                //
                // Note: The sign convention here is different from the
                // 3D case.  Concave-upward curves (smiles) have a positive
                // curvature.  Concave-downward curves (frowns) have a
                // negative curvature.  Be sure to take that into account!
                if (calculate_d2xyz)
                  {
                    const Real numerator   = this->d2xyzdxi2_map[p] * this->normals[p];
                    const Real denominator = this->dxyzdxi_map[p].norm_sq();
                    libmesh_assert_not_equal_to (denominator, 0);
                    curvatures[p] = numerator / denominator;
                  }
#endif
              }

            // compute the jacobian at the quadrature points
            for (unsigned int p=0; p<n_qp; p++)
              {
                const Real the_jac = this->dxyzdxi_map[p].norm();

                libmesh_assert_greater (the_jac, 0.);

                this->JxW[p] = the_jac*qw[p];
              }
          }

        // done computing the map
        break;
      }



    case 3:
      {
        // A 3D finite element living in 3D space.
        // Resize the vectors to hold data at the quadrature points
        {
          if (calculate_xyz)
            this->xyz.resize(n_qp);
          if (calculate_dxyz)
            {
              this->dxyzdxi_map.resize(n_qp);
              this->dxyzdeta_map.resize(n_qp);
              this->tangents.resize(n_qp);
              this->normals.resize(n_qp);
              this->JxW.resize(n_qp);
            }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          if (calculate_d2xyz)
            {
              this->d2xyzdxi2_map.resize(n_qp);
              this->d2xyzdxideta_map.resize(n_qp);
              this->d2xyzdeta2_map.resize(n_qp);
              this->curvatures.resize(n_qp);
            }
#endif
        }

        // Clear the entities that will be summed
        for (unsigned int p=0; p<n_qp; p++)
          {
            if (calculate_xyz)
              this->xyz[p].zero();
            if (calculate_dxyz)
              {
                this->tangents[p].resize(LIBMESH_DIM-1); // 1 Tangent in 2D, 2 in 3D
                this->dxyzdxi_map[p].zero();
                this->dxyzdeta_map[p].zero();
              }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
            if (calculate_d2xyz)
              {
                this->d2xyzdxi2_map[p].zero();
                this->d2xyzdxideta_map[p].zero();
                this->d2xyzdeta2_map[p].zero();
              }
#endif
          }

        // compute x, dxdxi at the quadrature points
        for (unsigned int i=0; i<n_mapping_shape_functions; i++) // sum over the nodes
          {
            const Point & side_point = side->point(i);

            for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
              {
                if (calculate_xyz)
                  this->xyz[p].add_scaled         (side_point, this->psi_map[i][p]);
                if (calculate_dxyz)
                  {
                    this->dxyzdxi_map[p].add_scaled (side_point, this->dpsidxi_map[i][p]);
                    this->dxyzdeta_map[p].add_scaled(side_point, this->dpsideta_map[i][p]);
                  }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                if (calculate_d2xyz)
                  {
                    this->d2xyzdxi2_map[p].add_scaled   (side_point, this->d2psidxi2_map[i][p]);
                    this->d2xyzdxideta_map[p].add_scaled(side_point, this->d2psidxideta_map[i][p]);
                    this->d2xyzdeta2_map[p].add_scaled  (side_point, this->d2psideta2_map[i][p]);
                  }
#endif
              }
          }

        // Compute the tangents, normal, and curvature at the quadrature point
        if (calculate_dxyz)
          {
            for (unsigned int p=0; p<n_qp; p++)
              {
                const Point n  = this->dxyzdxi_map[p].cross(this->dxyzdeta_map[p]);
                this->normals[p]     = n.unit();
                this->tangents[p][0] = this->dxyzdxi_map[p].unit();
                this->tangents[p][1] = n.cross(this->dxyzdxi_map[p]).unit();

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
                if (calculate_d2xyz)
                  {
                    // Compute curvature using the typical nomenclature
                    // of the first and second fundamental forms.
                    // For reference, see:
                    // 1) http://mathworld.wolfram.com/MeanCurvature.html
                    //    (note -- they are using inward normal)
                    // 2) F.S. Merritt, Mathematics Manual, 1962, McGraw-Hill
                    const Real L  = -this->d2xyzdxi2_map[p]    * this->normals[p];
                    const Real M  = -this->d2xyzdxideta_map[p] * this->normals[p];
                    const Real N  = -this->d2xyzdeta2_map[p]   * this->normals[p];
                    const Real E  =  this->dxyzdxi_map[p].norm_sq();
                    const Real F  =  this->dxyzdxi_map[p]      * this->dxyzdeta_map[p];
                    const Real G  =  this->dxyzdeta_map[p].norm_sq();

                    const Real numerator   = E*N -2.*F*M + G*L;
                    const Real denominator = E*G - F*F;
                    libmesh_assert_not_equal_to (denominator, 0.);
                    curvatures[p] = 0.5*numerator/denominator;
                  }
#endif
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


                const Real the_jac = std::sqrt(g11*g22 - g12*g21);

                libmesh_assert_greater (the_jac, 0.);

                this->JxW[p] = the_jac*qw[p];
              }
          }

        // done computing the map
        break;
      }


    default:
      libmesh_error_msg("Invalid dimension dim = " << dim);
    }
}




void FEMap::compute_edge_map(int dim,
                             const std::vector<Real> & qw,
                             const Elem * edge)
{
  libmesh_assert(edge);

  if (dim == 2)
    {
      // A 2D finite element living in either 2D or 3D space.
      // The edges here are the sides of the element, so the
      // (misnamed) compute_face_map function does what we want
      this->compute_face_map(dim, qw, edge);
      return;
    }

  libmesh_assert_equal_to (dim, 3);  // 1D is unnecessary and currently unsupported

  LOG_SCOPE("compute_edge_map()", "FEMap");

  // We're calculating now!
  this->determine_calculations();

  // The number of quadrature points.
  const unsigned int n_qp = cast_int<unsigned int>(qw.size());

  // Resize the vectors to hold data at the quadrature points
  if (calculate_xyz)
    this->xyz.resize(n_qp);
  if (calculate_dxyz)
    {
      this->dxyzdxi_map.resize(n_qp);
      this->dxyzdeta_map.resize(n_qp);
      this->tangents.resize(n_qp);
      this->normals.resize(n_qp);
      this->JxW.resize(n_qp);
    }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  if (calculate_d2xyz)
    {
      this->d2xyzdxi2_map.resize(n_qp);
      this->d2xyzdxideta_map.resize(n_qp);
      this->d2xyzdeta2_map.resize(n_qp);
      this->curvatures.resize(n_qp);
    }
#endif

  // Clear the entities that will be summed
  for (unsigned int p=0; p<n_qp; p++)
    {
      if (calculate_xyz)
        this->xyz[p].zero();
      if (calculate_dxyz)
        {
          this->tangents[p].resize(1);
          this->dxyzdxi_map[p].zero();
          this->dxyzdeta_map[p].zero();
        }
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
      if (calculate_d2xyz)
        {
          this->d2xyzdxi2_map[p].zero();
          this->d2xyzdxideta_map[p].zero();
          this->d2xyzdeta2_map[p].zero();
        }
#endif
    }

  // compute x, dxdxi at the quadrature points
  for (unsigned int i=0,
       psi_map_size=cast_int<unsigned int>(psi_map.size());
       i != psi_map_size; i++) // sum over the nodes
    {
      const Point & edge_point = edge->point(i);

      for (unsigned int p=0; p<n_qp; p++) // for each quadrature point...
        {
          if (calculate_xyz)
            this->xyz[p].add_scaled             (edge_point, this->psi_map[i][p]);
          if (calculate_dxyz)
            this->dxyzdxi_map[p].add_scaled     (edge_point, this->dpsidxi_map[i][p]);
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          if (calculate_d2xyz)
            this->d2xyzdxi2_map[p].add_scaled   (edge_point, this->d2psidxi2_map[i][p]);
#endif
        }
    }

  // Compute the tangents at the quadrature point
  // FIXME: normals (plural!) and curvatures are uncalculated
  if (calculate_dxyz)
    for (unsigned int p=0; p<n_qp; p++)
      {
        const Point n  = this->dxyzdxi_map[p].cross(this->dxyzdeta_map[p]);
        this->tangents[p][0] = this->dxyzdxi_map[p].unit();

        // compute the jacobian at the quadrature points
        const Real the_jac = std::sqrt(this->dxdxi_map(p)*this->dxdxi_map(p) +
                                       this->dydxi_map(p)*this->dydxi_map(p) +
                                       this->dzdxi_map(p)*this->dzdxi_map(p));

        libmesh_assert_greater (the_jac, 0.);

        this->JxW[p] = the_jac*qw[p];
      }
}


// Explicit FEMap Instantiations
FACE_EDGE_SHAPE_ERROR(0,init_face_shape_functions)
template void FEMap::init_face_shape_functions<1>(const std::vector<Point> &, const Elem *);
template void FEMap::init_face_shape_functions<2>(const std::vector<Point> &, const Elem *);
template void FEMap::init_face_shape_functions<3>(const std::vector<Point> &, const Elem *);

FACE_EDGE_SHAPE_ERROR(0,init_edge_shape_functions)
template void FEMap::init_edge_shape_functions<1>(const std::vector<Point> &, const Elem *);
template void FEMap::init_edge_shape_functions<2>(const std::vector<Point> &, const Elem *);
template void FEMap::init_edge_shape_functions<3>(const std::vector<Point> &, const Elem *);

//--------------------------------------------------------------
// Explicit FE instantiations
#define REINIT_AND_SIDE_MAPS(_type) \
template void FE<1,_type>::reinit(Elem const *, unsigned int, Real, const std::vector<Point> * const, const std::vector<Real> * const); \
template void FE<1,_type>::side_map(Elem const *, Elem const *, const unsigned int, const std::vector<Point> &, std::vector<Point> &); \
template void FE<2,_type>::reinit(Elem const *, unsigned int, Real, const std::vector<Point> * const, const std::vector<Real> * const); \
template void FE<2,_type>::side_map(Elem const *, Elem const *, const unsigned int, const std::vector<Point> &, std::vector<Point> &); \
template void FE<2,_type>::edge_reinit(Elem const *, unsigned int, Real, const std::vector<Point> * const, const std::vector<Real> * const); \
template void FE<3,_type>::reinit(Elem const *, unsigned int, Real, const std::vector<Point> * const, const std::vector<Real> * const); \
template void FE<3,_type>::side_map(Elem const *, Elem const *, const unsigned int, const std::vector<Point> &, std::vector<Point> &); \
template void FE<3,_type>::edge_reinit(Elem const *, unsigned int, Real, const std::vector<Point> * const, const std::vector<Real> * const)

REINIT_AND_SIDE_MAPS(LAGRANGE);
REINIT_AND_SIDE_MAPS(LAGRANGE_VEC);
REINIT_AND_SIDE_MAPS(L2_LAGRANGE);
REINIT_AND_SIDE_MAPS(HIERARCHIC);
REINIT_AND_SIDE_MAPS(L2_HIERARCHIC);
REINIT_AND_SIDE_MAPS(CLOUGH);
REINIT_AND_SIDE_MAPS(HERMITE);
REINIT_AND_SIDE_MAPS(MONOMIAL);
REINIT_AND_SIDE_MAPS(MONOMIAL_VEC);
REINIT_AND_SIDE_MAPS(SCALAR);
REINIT_AND_SIDE_MAPS(XYZ);

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
REINIT_AND_SIDE_MAPS(BERNSTEIN);
REINIT_AND_SIDE_MAPS(SZABAB);
REINIT_AND_SIDE_MAPS(RATIONAL_BERNSTEIN);
#endif

template void FE<2,SUBDIVISION>::reinit(Elem const *, unsigned int, Real, const std::vector<Point> * const, const std::vector<Real> * const);
template void FE<2,NEDELEC_ONE>::reinit(Elem const *, unsigned int, Real, const std::vector<Point> * const, const std::vector<Real> * const);
template void FE<2,NEDELEC_ONE>::side_map(Elem const *, Elem const *, const unsigned int, const std::vector<Point> &, std::vector<Point> &);
template void FE<2,NEDELEC_ONE>::edge_reinit(Elem const *, unsigned int, Real, const std::vector<Point> * const, const std::vector<Real> * const);

template void FE<3,NEDELEC_ONE>::reinit(Elem const *, unsigned int, Real, const std::vector<Point> * const, const std::vector<Real> * const);
template void FE<3,NEDELEC_ONE>::side_map(Elem const *, Elem const *, const unsigned int, const std::vector<Point> &, std::vector<Point> &);
template void FE<3,NEDELEC_ONE>::edge_reinit(Elem const *, unsigned int, Real, const std::vector<Point> * const, const std::vector<Real> * const);

// Intel 9.1 complained it needed this in devel mode.
//template void FE<2,XYZ>::init_face_shape_functions(const std::vector<Point> &, const Elem *);
//template void FE<3,XYZ>::init_face_shape_functions(const std::vector<Point> &, const Elem *);

} // namespace libMesh
