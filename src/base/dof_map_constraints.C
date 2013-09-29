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

// C++ Includes -------------------------------------
#include <set>
#include <algorithm> // for std::count, std::fill
#include <sstream>
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>

// Local Includes -----------------------------------
#include "libmesh/boundary_info.h" // needed for dirichlet constraints
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/elem_range.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_type.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_inserter_iterator.h"
#include "libmesh/numeric_vector.h" // for enforce_constraints_exactly()
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel_elem.h"
#include "libmesh/parallel_node.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/periodic_boundary_base.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/quadrature.h" // for dirichlet constraints
#include "libmesh/raw_accessor.h"
#include "libmesh/system.h" // needed by enforce_constraints_exactly()
#include "libmesh/threads.h"
#include "libmesh/tensor_tools.h"


// Anonymous namespace to hold helper classes
namespace {

using namespace libMesh;

  class ComputeConstraints
  {
  public:
    ComputeConstraints (DofConstraints &constraints,
			DofMap &dof_map,
#ifdef LIBMESH_ENABLE_PERIODIC
			PeriodicBoundaries &periodic_boundaries,
#endif
			const MeshBase &mesh,
			const unsigned int variable_number) :
      _constraints(constraints),
      _dof_map(dof_map),
#ifdef LIBMESH_ENABLE_PERIODIC
      _periodic_boundaries(periodic_boundaries),
#endif
      _mesh(mesh),
      _variable_number(variable_number)
    {}

    void operator()(const ConstElemRange &range) const
    {
      const Variable &var_description = _dof_map.variable(_variable_number);

#ifdef LIBMESH_ENABLE_PERIODIC
      AutoPtr<PointLocatorBase> point_locator;
      bool have_periodic_boundaries = !_periodic_boundaries.empty();
      if (have_periodic_boundaries)
        point_locator = _mesh.sub_point_locator();
#endif

      for (ConstElemRange::const_iterator it = range.begin(); it!=range.end(); ++it)
	if (var_description.active_on_subdomain((*it)->subdomain_id()))
	  {
#ifdef LIBMESH_ENABLE_AMR
	    FEInterface::compute_constraints (_constraints,
					      _dof_map,
					      _variable_number,
					      *it);
#endif
#ifdef LIBMESH_ENABLE_PERIODIC
	    // FIXME: periodic constraints won't work on a non-serial
	    // mesh unless it's kept ghost elements from opposing
	    // boundaries!
            if (have_periodic_boundaries)
	      FEInterface::compute_periodic_constraints (_constraints,
						         _dof_map,
						         _periodic_boundaries,
						         _mesh,
						         point_locator.get(),
						         _variable_number,
						         *it);
#endif
	  }
    }

  private:
    DofConstraints &_constraints;
    DofMap &_dof_map;
#ifdef LIBMESH_ENABLE_PERIODIC
    PeriodicBoundaries &_periodic_boundaries;
#endif
    const MeshBase &_mesh;
    const unsigned int _variable_number;
  };

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  class ComputeNodeConstraints
  {
  public:
    ComputeNodeConstraints (NodeConstraints &node_constraints,
			    DofMap &dof_map,
#ifdef LIBMESH_ENABLE_PERIODIC
			    PeriodicBoundaries &periodic_boundaries,
#endif
			    const MeshBase &mesh) :
      _node_constraints(node_constraints),
      _dof_map(dof_map),
#ifdef LIBMESH_ENABLE_PERIODIC
      _periodic_boundaries(periodic_boundaries),
#endif
      _mesh(mesh)
    {}

    void operator()(const ConstElemRange &range) const
    {
#ifdef LIBMESH_ENABLE_PERIODIC
      AutoPtr<PointLocatorBase> point_locator;
      bool have_periodic_boundaries = !_periodic_boundaries.empty();
      if (have_periodic_boundaries)
        point_locator = _mesh.sub_point_locator();
#endif

      for (ConstElemRange::const_iterator it = range.begin(); it!=range.end(); ++it)
	{
#ifdef LIBMESH_ENABLE_AMR
	  FEBase::compute_node_constraints (_node_constraints, *it);
#endif
#ifdef LIBMESH_ENABLE_PERIODIC
	  // FIXME: periodic constraints won't work on a non-serial
	  // mesh unless it's kept ghost elements from opposing
	  // boundaries!
          if (have_periodic_boundaries)
	    FEBase::compute_periodic_node_constraints (_node_constraints,
						       _periodic_boundaries,
						       _mesh,
						       point_locator.get(),
						       *it);
#endif
	  }
    }

  private:
    NodeConstraints &_node_constraints;
    DofMap &_dof_map;
#ifdef LIBMESH_ENABLE_PERIODIC
    PeriodicBoundaries &_periodic_boundaries;
#endif
    const MeshBase &_mesh;
  };
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS


#ifdef LIBMESH_ENABLE_DIRICHLET

  /**
   * This class implements turning an arbitrary
   * boundary function into Dirichlet constraints.  It
   * may be executed in parallel on multiple threads.
   */
  class ConstrainDirichlet
  {
  private:
    DofMap                  &dof_map;
    const MeshBase          &mesh;
    const Real               time;
    const DirichletBoundary  dirichlet;

    template<typename OutputType>
    void apply_dirichlet_impl( const ConstElemRange &range, FunctionBase<Number> *f, FunctionBase<Gradient> *g,
			       const std::set<boundary_id_type> &b, DenseMatrix<Real>& Ke,
			       DenseVector<Number>& Fe, DenseVector<Number>& Ue,
			       const unsigned int var, const Variable&variable,
			       const FEType& fe_type ) const
    {
      typedef OutputType                                                      OutputShape;
      typedef typename TensorTools::IncrementRank<OutputShape>::type          OutputGradient;
      //typedef typename TensorTools::IncrementRank<OutputGradient>::type       OutputTensor;
      typedef typename TensorTools::MakeNumber<OutputShape>::type             OutputNumber;
      typedef typename TensorTools::IncrementRank<OutputNumber>::type         OutputNumberGradient;
      //typedef typename TensorTools::IncrementRank<OutputNumberGradient>::type OutputNumberTensor;

      // The dimensionality of the current mesh
      const unsigned int dim = mesh.mesh_dimension();

      // Boundary info for the current mesh
      const BoundaryInfo& boundary_info = *mesh.boundary_info;

      unsigned int n_vec_dim = FEInterface::n_vec_dim(mesh, fe_type);

      const unsigned int var_component =
	variable.first_scalar_number();

      // Get FE objects of the appropriate type
      AutoPtr<FEGenericBase<OutputType> > fe = FEGenericBase<OutputType>::build(dim, fe_type);

      // Prepare variables for projection
      AutoPtr<QBase> qedgerule (fe_type.default_quadrature_rule(1));
      AutoPtr<QBase> qsiderule (fe_type.default_quadrature_rule(dim-1));

      // The values of the shape functions at the quadrature
      // points
      const std::vector<std::vector<OutputShape> >& phi = fe->get_phi();

      // The gradients of the shape functions at the quadrature
      // points on the child element.
      const std::vector<std::vector<OutputGradient> > *dphi = NULL;

      const FEContinuity cont = fe->get_continuity();

      if (cont == C_ONE)
	{
	  // We'll need gradient data for a C1 projection
	  libmesh_assert(g);

	  const std::vector<std::vector<OutputGradient> >&
	    ref_dphi = fe->get_dphi();
	  dphi = &ref_dphi;
	}

      // The Jacobian * quadrature weight at the quadrature points
      const std::vector<Real>& JxW =
	fe->get_JxW();

      // The XYZ locations of the quadrature points
      const std::vector<Point>& xyz_values =
	fe->get_xyz();

      // The global DOF indices
      std::vector<dof_id_type> dof_indices;
      // Side/edge local DOF indices
      std::vector<unsigned int> side_dofs;

      // Iterate over all the elements in the range
      for (ConstElemRange::const_iterator elem_it=range.begin(); elem_it != range.end(); ++elem_it)
	{
	  const Elem* elem = *elem_it;

          // We only calculate Dirichlet constraints on active
          // elements
          if (!elem->active())
            continue;

	  // Per-subdomain variables don't need to be projected on
	  // elements where they're not active
	  if (!variable.active_on_subdomain(elem->subdomain_id()))
	    continue;

	  // Find out which nodes, edges and sides are on a requested
	  // boundary:
	  std::vector<bool> is_boundary_node(elem->n_nodes(), false),
	    is_boundary_edge(elem->n_edges(), false),
	    is_boundary_side(elem->n_sides(), false);
	  for (unsigned char s=0; s != elem->n_sides(); ++s)
	    {
	      // First see if this side has been requested
	      const std::vector<boundary_id_type>& bc_ids =
		boundary_info.boundary_ids (elem, s);
	      bool do_this_side = false;
	      for (std::size_t i=0; i != bc_ids.size(); ++i)
		if (b.count(bc_ids[i]))
		  {
		    do_this_side = true;
		    break;
		  }
	      if (!do_this_side)
		continue;

	      is_boundary_side[s] = true;

	      // Then see what nodes and what edges are on it
	      for (unsigned int n=0; n != elem->n_nodes(); ++n)
		if (elem->is_node_on_side(n,s))
		  is_boundary_node[n] = true;
	      for (unsigned int e=0; e != elem->n_edges(); ++e)
		if (elem->is_edge_on_side(e,s))
		  is_boundary_edge[e] = true;
	    }

	  // Update the DOF indices for this element based on
	  // the current mesh
	  dof_map.dof_indices (elem, dof_indices, var);

	  // The number of DOFs on the element
	  const unsigned int n_dofs =
	    libmesh_cast_int<unsigned int>(dof_indices.size());

	  // Fixed vs. free DoFs on edge/face projections
	  std::vector<char> dof_is_fixed(n_dofs, false); // bools
	  std::vector<int> free_dof(n_dofs, 0);

	  // The element type
	  const ElemType elem_type = elem->type();

	  // The number of nodes on the new element
	  const unsigned int n_nodes = elem->n_nodes();

	  // Zero the interpolated values
	  Ue.resize (n_dofs); Ue.zero();

	  // In general, we need a series of
	  // projections to ensure a unique and continuous
	  // solution.  We start by interpolating boundary nodes, then
	  // hold those fixed and project boundary edges, then hold
	  // those fixed and project boundary faces,

	  // Interpolate node values first
	  unsigned int current_dof = 0;
	  for (unsigned int n=0; n!= n_nodes; ++n)
	    {
	      // FIXME: this should go through the DofMap,
	      // not duplicate dof_indices code badly!
	      const unsigned int nc =
		FEInterface::n_dofs_at_node (dim, fe_type, elem_type,
					     n);
	      if (!elem->is_vertex(n) || !is_boundary_node[n])
		{
		  current_dof += nc;
		  continue;
		}
	      if (cont == DISCONTINUOUS)
		{
		  libmesh_assert_equal_to (nc, 0);
		}
	      // Assume that C_ZERO elements have a single nodal
	      // value shape function
	      else if (cont == C_ZERO)
		{
		  libmesh_assert_equal_to (nc, n_vec_dim);
		  for( unsigned int c = 0; c < n_vec_dim; c++ )
		    {
		      Ue(current_dof+c) =
			f->component(var_component+c,
				     elem->point(n), time);
		      dof_is_fixed[current_dof+c] = true;
		    }
		  current_dof += n_vec_dim;
		}
	      // The hermite element vertex shape functions are weird
	      else if (fe_type.family == HERMITE)
		{
		  Ue(current_dof) =
		    f->component(var_component,
				 elem->point(n), time);
		  dof_is_fixed[current_dof] = true;
		  current_dof++;
		  Gradient grad =
		    g->component(var_component,
				 elem->point(n), time);
		  // x derivative
		  Ue(current_dof) = grad(0);
		  dof_is_fixed[current_dof] = true;
		  current_dof++;
		  if (dim > 1)
		    {
		      // We'll finite difference mixed derivatives
		      Point nxminus = elem->point(n),
			nxplus = elem->point(n);
		      nxminus(0) -= TOLERANCE;
		      nxplus(0) += TOLERANCE;
		      Gradient gxminus =
			g->component(var_component,
				     nxminus, time);
		      Gradient gxplus =
			g->component(var_component,
				     nxplus, time);
		      // y derivative
		      Ue(current_dof) = grad(1);
		      dof_is_fixed[current_dof] = true;
		      current_dof++;
		      // xy derivative
		      Ue(current_dof) = (gxplus(1) - gxminus(1))
			/ 2. / TOLERANCE;
		      dof_is_fixed[current_dof] = true;
		      current_dof++;

		      if (dim > 2)
			{
			  // z derivative
			  Ue(current_dof) = grad(2);
			  dof_is_fixed[current_dof] = true;
			  current_dof++;
			  // xz derivative
			  Ue(current_dof) = (gxplus(2) - gxminus(2))
			    / 2. / TOLERANCE;
			  dof_is_fixed[current_dof] = true;
			  current_dof++;
			  // We need new points for yz
			  Point nyminus = elem->point(n),
			    nyplus = elem->point(n);
			  nyminus(1) -= TOLERANCE;
			  nyplus(1) += TOLERANCE;
			  Gradient gyminus =
			    g->component(var_component,
					 nyminus, time);
			  Gradient gyplus =
			    g->component(var_component,
					 nyplus, time);
			  // xz derivative
			  Ue(current_dof) = (gyplus(2) - gyminus(2))
			    / 2. / TOLERANCE;
			  dof_is_fixed[current_dof] = true;
			  current_dof++;
			  // Getting a 2nd order xyz is more tedious
			  Point nxmym = elem->point(n),
			    nxmyp = elem->point(n),
			    nxpym = elem->point(n),
			    nxpyp = elem->point(n);
			  nxmym(0) -= TOLERANCE;
			  nxmym(1) -= TOLERANCE;
			  nxmyp(0) -= TOLERANCE;
			  nxmyp(1) += TOLERANCE;
			  nxpym(0) += TOLERANCE;
			  nxpym(1) -= TOLERANCE;
			  nxpyp(0) += TOLERANCE;
			  nxpyp(1) += TOLERANCE;
			  Gradient gxmym =
			    g->component(var_component,
					 nxmym, time);
			  Gradient gxmyp =
			    g->component(var_component,
					 nxmyp, time);
			  Gradient gxpym =
			    g->component(var_component,
					 nxpym, time);
			  Gradient gxpyp =
			    g->component(var_component,
					 nxpyp, time);
			  Number gxzplus = (gxpyp(2) - gxmyp(2))
			    / 2. / TOLERANCE;
			  Number gxzminus = (gxpym(2) - gxmym(2))
			    / 2. / TOLERANCE;
			  // xyz derivative
			  Ue(current_dof) = (gxzplus - gxzminus)
			    / 2. / TOLERANCE;
			  dof_is_fixed[current_dof] = true;
			  current_dof++;
			}
		    }
		}
	      // Assume that other C_ONE elements have a single nodal
	      // value shape function and nodal gradient component
	      // shape functions
	      else if (cont == C_ONE)
		{
		  libmesh_assert_equal_to (nc, 1 + dim);
		  Ue(current_dof) =
		    f->component(var_component,
				 elem->point(n), time);
		  dof_is_fixed[current_dof] = true;
		  current_dof++;
		  Gradient grad =
		    g->component(var_component,
				 elem->point(n), time);
		  for (unsigned int i=0; i!= dim; ++i)
		    {
		      Ue(current_dof) = grad(i);
		      dof_is_fixed[current_dof] = true;
		      current_dof++;
		    }
		}
	      else
		libmesh_error();
	    }

	  // In 3D, project any edge values next
	  if (dim > 2 && cont != DISCONTINUOUS)
	    for (unsigned int e=0; e != elem->n_edges(); ++e)
	      {
		if (!is_boundary_edge[e])
		  continue;

		FEInterface::dofs_on_edge(elem, dim, fe_type, e,
					  side_dofs);

		// Some edge dofs are on nodes and already
		// fixed, others are free to calculate
		unsigned int free_dofs = 0;
		for (unsigned int i=0; i != side_dofs.size(); ++i)
		  if (!dof_is_fixed[side_dofs[i]])
		    free_dof[free_dofs++] = i;

		// There may be nothing to project
		if (!free_dofs)
		  continue;

		Ke.resize (free_dofs, free_dofs); Ke.zero();
		Fe.resize (free_dofs); Fe.zero();
		// The new edge coefficients
		DenseVector<Number> Uedge(free_dofs);

		// Initialize FE data on the edge
		fe->attach_quadrature_rule (qedgerule.get());
		fe->edge_reinit (elem, e);
		const unsigned int n_qp = qedgerule->n_points();

		// Loop over the quadrature points
		for (unsigned int qp=0; qp<n_qp; qp++)
		  {
		    // solution at the quadrature point
		    OutputNumber fineval = 0;
		    libMesh::RawAccessor<OutputNumber> f_accessor( fineval, dim );

		    for( unsigned int c = 0; c < n_vec_dim; c++)
		      f_accessor(c) = f->component(var_component+c, xyz_values[qp], time);

		    // solution grad at the quadrature point
		    OutputNumberGradient finegrad;
		    libMesh::RawAccessor<OutputNumberGradient> g_accessor( finegrad, dim );

		    unsigned int g_rank;
		    switch( FEInterface::field_type( fe_type ) )
		      {
		      case TYPE_SCALAR:
			{
			  g_rank = 1;
			  break;
			}
		      case TYPE_VECTOR:
			{
			  g_rank = 2;
			  break;
			}
		      default:
			libmesh_error();
		      }

		    if (cont == C_ONE)
		      for( unsigned int c = 0; c < n_vec_dim; c++)
			for( unsigned int d = 0; d < g_rank; d++ )
			  g_accessor(c + d*dim ) =
			  g->component(var_component,
				       xyz_values[qp], time)(c);

		    // Form edge projection matrix
		    for (unsigned int sidei=0, freei=0;
			 sidei != side_dofs.size(); ++sidei)
		      {
			unsigned int i = side_dofs[sidei];
			// fixed DoFs aren't test functions
			if (dof_is_fixed[i])
			  continue;
			for (unsigned int sidej=0, freej=0;
			     sidej != side_dofs.size(); ++sidej)
			  {
			    unsigned int j = side_dofs[sidej];
			    if (dof_is_fixed[j])
			      Fe(freei) -= phi[i][qp] * phi[j][qp] *
				JxW[qp] * Ue(j);
			    else
			      Ke(freei,freej) += phi[i][qp] *
				phi[j][qp] * JxW[qp];
			    if (cont == C_ONE)
			      {
				if (dof_is_fixed[j])
				  Fe(freei) -= ((*dphi)[i][qp].contract((*dphi)[j][qp]) ) *
				    JxW[qp] * Ue(j);
				else
				  Ke(freei,freej) += ((*dphi)[i][qp].contract((*dphi)[j][qp]))
				    * JxW[qp];
			      }
			    if (!dof_is_fixed[j])
			      freej++;
			  }
			Fe(freei) += phi[i][qp] * fineval * JxW[qp];
			if (cont == C_ONE)
			  Fe(freei) += (finegrad.contract( (*dphi)[i][qp]) ) *
			    JxW[qp];
			freei++;
		      }
		  }

		Ke.cholesky_solve(Fe, Uedge);

		// Transfer new edge solutions to element
		for (unsigned int i=0; i != free_dofs; ++i)
		  {
		    Number &ui = Ue(side_dofs[free_dof[i]]);
		    libmesh_assert(std::abs(ui) < TOLERANCE ||
				   std::abs(ui - Uedge(i)) < TOLERANCE);
		    ui = Uedge(i);
		    dof_is_fixed[side_dofs[free_dof[i]]] = true;
		  }
	      }

	  // Project any side values (edges in 2D, faces in 3D)
	  if (dim > 1 && cont != DISCONTINUOUS)
	    for (unsigned int s=0; s != elem->n_sides(); ++s)
	      {
		if (!is_boundary_side[s])
		  continue;

		FEInterface::dofs_on_side(elem, dim, fe_type, s,
					  side_dofs);

		// Some side dofs are on nodes/edges and already
		// fixed, others are free to calculate
		unsigned int free_dofs = 0;
		for (unsigned int i=0; i != side_dofs.size(); ++i)
		  if (!dof_is_fixed[side_dofs[i]])
		    free_dof[free_dofs++] = i;

		// There may be nothing to project
		if (!free_dofs)
		  continue;

		Ke.resize (free_dofs, free_dofs); Ke.zero();
		Fe.resize (free_dofs); Fe.zero();
		// The new side coefficients
		DenseVector<Number> Uside(free_dofs);

		// Initialize FE data on the side
		fe->attach_quadrature_rule (qsiderule.get());
		fe->reinit (elem, s);
		const unsigned int n_qp = qsiderule->n_points();

		// Loop over the quadrature points
		for (unsigned int qp=0; qp<n_qp; qp++)
		  {
		    // solution at the quadrature point
		    OutputNumber fineval = 0;
		    libMesh::RawAccessor<OutputNumber> f_accessor( fineval, dim );

		    for( unsigned int c = 0; c < n_vec_dim; c++)
		      f_accessor(c) = f->component(var_component+c, xyz_values[qp], time);

		    // solution grad at the quadrature point
		    OutputNumberGradient finegrad;
		    libMesh::RawAccessor<OutputNumberGradient> g_accessor( finegrad, dim );

		    unsigned int g_rank;
		    switch( FEInterface::field_type( fe_type ) )
		      {
		      case TYPE_SCALAR:
			{
			  g_rank = 1;
			  break;
			}
		      case TYPE_VECTOR:
			{
			  g_rank = 2;
			  break;
			}
		      default:
			libmesh_error();
		      }

		    if (cont == C_ONE)
		      for( unsigned int c = 0; c < n_vec_dim; c++)
			for( unsigned int d = 0; d < g_rank; d++ )
			  g_accessor(c + d*dim ) =
			  g->component(var_component,
				       xyz_values[qp], time)(c);

		    // Form side projection matrix
		    for (unsigned int sidei=0, freei=0;
			 sidei != side_dofs.size(); ++sidei)
		      {
			unsigned int i = side_dofs[sidei];
			// fixed DoFs aren't test functions
			if (dof_is_fixed[i])
			  continue;
			for (unsigned int sidej=0, freej=0;
			     sidej != side_dofs.size(); ++sidej)
			  {
			    unsigned int j = side_dofs[sidej];
			    if (dof_is_fixed[j])
			      Fe(freei) -= phi[i][qp] * phi[j][qp] *
				JxW[qp] * Ue(j);
			    else
			      Ke(freei,freej) += phi[i][qp] *
				phi[j][qp] * JxW[qp];
			    if (cont == C_ONE)
			      {
				if (dof_is_fixed[j])
				  Fe(freei) -= ((*dphi)[i][qp].contract((*dphi)[j][qp])) *
				    JxW[qp] * Ue(j);
				else
				  Ke(freei,freej) += ((*dphi)[i][qp].contract((*dphi)[j][qp]))
				    * JxW[qp];
			      }
			    if (!dof_is_fixed[j])
			      freej++;
			  }
			Fe(freei) += (fineval * phi[i][qp]) * JxW[qp];
			if (cont == C_ONE)
			  Fe(freei) += (finegrad.contract((*dphi)[i][qp])) *
			    JxW[qp];
			freei++;
		      }
		  }

		Ke.cholesky_solve(Fe, Uside);

		// Transfer new side solutions to element
		for (unsigned int i=0; i != free_dofs; ++i)
		  {
		    Number &ui = Ue(side_dofs[free_dof[i]]);
		    libmesh_assert(std::abs(ui) < TOLERANCE ||
				   std::abs(ui - Uside(i)) < TOLERANCE);
		    ui = Uside(i);
		    dof_is_fixed[side_dofs[free_dof[i]]] = true;
		  }
	      }

	  // Lock the DofConstraints since it is shared among threads.
	  {
	    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

	    for (unsigned int i = 0; i < n_dofs; i++)
	      {
		DofConstraintRow empty_row;
		if (dof_is_fixed[i] &&
		    !dof_map.is_constrained_dof(dof_indices[i]))
		  dof_map.add_constraint_row
		    (dof_indices[i], empty_row,
		     Ue(i), /* forbid_constraint_overwrite = */ true);
	      }
	  }
	}

    } // apply_dirichlet_impl

  public:
    ConstrainDirichlet (DofMap &dof_map_in,
		        const MeshBase &mesh_in,
		        const Real time_in,
		        const DirichletBoundary &dirichlet_in) :
    dof_map(dof_map_in),
    mesh(mesh_in),
    time(time_in),
    dirichlet(dirichlet_in) { }

    ConstrainDirichlet (const ConstrainDirichlet &in) :
    dof_map(in.dof_map),
    mesh(in.mesh),
    time(in.time),
    dirichlet(in.dirichlet) { }

    void operator()(const ConstElemRange &range) const
    {
      FunctionBase<Number> *f = dirichlet.f.get();
      FunctionBase<Gradient> *g = dirichlet.g.get();
      const std::set<boundary_id_type> &b = dirichlet.b;

      // We need data to project
      libmesh_assert(f);

      /**
       * This method examines an arbitrary boundary solution to calculate
       * corresponding Dirichlet constraints on the current mesh.  The
       * input function \p f gives the arbitrary solution.
       */

      // The element matrix and RHS for projections.
      // Note that Ke is always real-valued, whereas
      // Fe may be complex valued if complex number
      // support is enabled
      DenseMatrix<Real> Ke;
      DenseVector<Number> Fe;
      // The new element coefficients
      DenseVector<Number> Ue;

      // Loop over all the variables we've been requested to project
      for (unsigned int v=0; v!=dirichlet.variables.size(); v++)
	{
	  const unsigned int var = dirichlet.variables[v];

	  const Variable& variable = dof_map.variable(var);

	  const FEType& fe_type = variable.type();

	  if (fe_type.family == SCALAR)
	    continue;

	  switch( FEInterface::field_type( fe_type ) )
	    {
	    case TYPE_SCALAR:
	      {
		this->apply_dirichlet_impl<Real>( range, f, g, b, Ke, Fe, Ue, var, variable, fe_type );
		break;
	      }
	    case TYPE_VECTOR:
	      {
		this->apply_dirichlet_impl<RealGradient>( range, f, g, b, Ke, Fe, Ue, var, variable, fe_type );
		break;
	      }
	    default:
	      libmesh_error();
	    }

	}
    }

  }; // class ConstrainDirichlet


#endif // LIBMESH_ENABLE_DIRICHLET


} // anonymous namespace



namespace libMesh
{

// ------------------------------------------------------------
// DofMap member functions

#ifdef LIBMESH_ENABLE_CONSTRAINTS


dof_id_type DofMap::n_constrained_dofs() const
{
  parallel_object_only();

  dof_id_type nc_dofs = this->n_local_constrained_dofs();
  this->comm().sum(nc_dofs);
  return nc_dofs;
}


dof_id_type DofMap::n_local_constrained_dofs() const
{
  const DofConstraints::const_iterator lower =
    _dof_constraints.lower_bound(this->first_dof()),
                                       upper =
    _dof_constraints.upper_bound(this->end_dof()-1);

  return std::distance(lower, upper);
}


void DofMap::create_dof_constraints(const MeshBase& mesh, Real time)
{
  parallel_object_only();

  START_LOG("create_dof_constraints()", "DofMap");

  libmesh_assert (mesh.is_prepared());

  const unsigned int dim = mesh.mesh_dimension();

  // We might get constraint equations from AMR hanging nodes in 2D/3D
  // or from boundary conditions in any dimension
  const bool possible_local_constraints = false
#ifdef LIBMESH_ENABLE_AMR
    || dim > 1
#endif
#ifdef LIBMESH_ENABLE_PERIODIC
    || !_periodic_boundaries->empty()
#endif
#ifdef LIBMESH_ENABLE_DIRICHLET
    || !_dirichlet_boundaries->empty()
#endif
  ;

  // Even if we don't have constraints, another processor might.
  bool possible_global_constraints = possible_local_constraints;
#if defined(LIBMESH_ENABLE_PERIODIC) || defined(LIBMESH_ENABLE_DIRICHLET) || defined(LIBMESH_ENABLE_AMR)
  libmesh_assert(this->comm().verify(mesh.is_serial()));

  if (!mesh.is_serial())
    this->comm().max(possible_global_constraints);
#endif

  if (!possible_global_constraints)
    {
      // make sure we stop logging though
      STOP_LOG("create_dof_constraints()", "DofMap");
      return;
    }

  // Here we build the hanging node constraints.  This is done
  // by enforcing the condition u_a = u_b along hanging sides.
  // u_a = u_b is collocated at the nodes of side a, which gives
  // one row of the constraint matrix.

  // define the range of elements of interest
  ConstElemRange range;
  if (possible_local_constraints)
    {
      // With SerialMesh or a serial ParallelMesh, every processor
      // computes every constraint
      MeshBase::const_element_iterator
        elem_begin = mesh.elements_begin(),
        elem_end   = mesh.elements_end();

      // With a parallel ParallelMesh, processors compute only
      // their local constraints
      if (!mesh.is_serial())
        {
	  elem_begin = mesh.local_elements_begin();
	  elem_end   = mesh.local_elements_end();
        }

      // set the range to contain the specified elements
      range.reset (elem_begin, elem_end);
    }
  else
    range.reset (mesh.elements_end(), mesh.elements_end());

  // compute_periodic_constraints requires a point_locator() from our
  // Mesh, that point_locator() construction is threaded.  Rather than
  // nest threads within threads we'll make sure it's preconstructed.
#ifdef LIBMESH_ENABLE_PERIODIC
  if (!_periodic_boundaries->empty() && !range.empty())
    mesh.sub_point_locator();
#endif

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  // recalculate node constraints from scratch
  _node_constraints.clear();

  Threads::parallel_for (range,
			 ComputeNodeConstraints (_node_constraints,
					         *this,
#ifdef LIBMESH_ENABLE_PERIODIC
					         *_periodic_boundaries,
#endif
					         mesh));
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS


  // recalculate dof constraints from scratch
  _dof_constraints.clear();

  // Look at all the variables in the system.  Reset the element
  // range at each iteration -- there is no need to reconstruct it.
  for (unsigned int variable_number=0; variable_number<this->n_variables();
       ++variable_number, range.reset())
    Threads::parallel_for (range,
			   ComputeConstraints (_dof_constraints,
					       *this,
#ifdef LIBMESH_ENABLE_PERIODIC
					       *_periodic_boundaries,
#endif
					       mesh,
					       variable_number));

#ifdef LIBMESH_ENABLE_DIRICHLET
  for (DirichletBoundaries::iterator
         i = _dirichlet_boundaries->begin();
       i != _dirichlet_boundaries->end(); ++i, range.reset())
    {
      Threads::parallel_for
	(range,
	 ConstrainDirichlet(*this, mesh, time, **i)
	 );
    }
#endif // LIBMESH_ENABLE_DIRICHLET

  STOP_LOG("create_dof_constraints()", "DofMap");
}



void DofMap::add_constraint_row (const dof_id_type dof_number,
				 const DofConstraintRow& constraint_row,
                                 const Number constraint_rhs,
				 const bool forbid_constraint_overwrite)
{
  // Optionally allow the user to overwrite constraints.  Defaults to false.
  if (forbid_constraint_overwrite)
    if (this->is_constrained_dof(dof_number))
      {
	libMesh::err << "ERROR: DOF " << dof_number << " was already constrained!"
		      << std::endl;
	libmesh_error();
      }

  std::pair<dof_id_type, std::pair<DofConstraintRow,Number> > kv(dof_number, std::make_pair(constraint_row, constraint_rhs));

  _dof_constraints.insert(kv);
}



void DofMap::print_dof_constraints(std::ostream& os,
                                   bool print_nonlocal) const
{
  parallel_object_only();

  std::string local_constraints =
    this->get_local_constraints(print_nonlocal);

  if (this->processor_id())
    {
      this->comm().send(0, local_constraints);
    }
  else
    {
      os << "Processor 0:\n";
      os << local_constraints;

      for (processor_id_type i=1; i<this->n_processors(); ++i)
        {
          this->comm().receive(i, local_constraints);
          os << "Processor " << i << ":\n";
          os << local_constraints;
        }
    }
}

std::string DofMap::get_local_constraints(bool print_nonlocal) const
{
  std::ostringstream os;
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  if (print_nonlocal)
    os << "All ";
  else
    os << "Local ";

  os << "Node Constraints:"
     << std::endl;

  for (NodeConstraints::const_iterator it=_node_constraints.begin();
       it != _node_constraints.end(); ++it)
    {
      const Node *node = it->first;

      // Skip non-local nodes if requested
      if (!print_nonlocal &&
          node->processor_id() != this->processor_id())
        continue;

      const NodeConstraintRow& row = it->second.first;
      const Point& offset = it->second.second;

      os << "Constraints for Node id " << node->id()
	 << ": \t";

      for (NodeConstraintRow::const_iterator pos=row.begin();
	   pos != row.end(); ++pos)
	os << " (" << pos->first->id() << ","
	   << pos->second << ")\t";

      os << "rhs: " << offset;

      os << std::endl;
    }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

  if (print_nonlocal)
    os << "All ";
  else
    os << "Local ";

  os << "DoF Constraints:"
     << std::endl;

  for (DofConstraints::const_iterator it=_dof_constraints.begin();
       it != _dof_constraints.end(); ++it)
    {
      const dof_id_type i = it->first;

      // Skip non-local dofs if requested
      if (!print_nonlocal &&
          ((i < this->first_dof()) ||
           (i >= this->end_dof())))
        continue;

      const DofConstraintRow& row = it->second.first;
      const Number rhs = it->second.second;

      os << "Constraints for DoF " << i
	 << ": \t";

      for (DofConstraintRow::const_iterator pos=row.begin();
	   pos != row.end(); ++pos)
	os << " (" << pos->first << ","
	   << pos->second << ")\t";

      os << "rhs: " << rhs;

      os << std::endl;
    }

  return os.str();
}



void DofMap::constrain_element_matrix (DenseMatrix<Number>& matrix,
				       std::vector<dof_id_type>& elem_dofs,
				       bool asymmetric_constraint_rows) const
{
  libmesh_assert_equal_to (elem_dofs.size(), matrix.m());
  libmesh_assert_equal_to (elem_dofs.size(), matrix.n());

  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // The constrained matrix is built up as C^T K C.
  DenseMatrix<Number> C;


  this->build_constraint_matrix (C, elem_dofs);

  START_LOG("constrain_elem_matrix()", "DofMap");

  // It is possible that the matrix is not constrained at all.
  if ((C.m() == matrix.m()) &&
      (C.n() == elem_dofs.size())) // It the matrix is constrained
    {
      // Compute the matrix-matrix-matrix product C^T K C
      matrix.left_multiply_transpose  (C);
      matrix.right_multiply (C);


      libmesh_assert_equal_to (matrix.m(), matrix.n());
      libmesh_assert_equal_to (matrix.m(), elem_dofs.size());
      libmesh_assert_equal_to (matrix.n(), elem_dofs.size());


      for (unsigned int i=0; i<elem_dofs.size(); i++)
	// If the DOF is constrained
	if (this->is_constrained_dof(elem_dofs[i]))
	  {
	    for (unsigned int j=0; j<matrix.n(); j++)
	      matrix(i,j) = 0.;

	    matrix(i,i) = 1.;

            if (asymmetric_constraint_rows)
              {
	        DofConstraints::const_iterator
	          pos = _dof_constraints.find(elem_dofs[i]);

	        libmesh_assert (pos != _dof_constraints.end());

	        const DofConstraintRow& constraint_row = pos->second.first;

		// This is an overzealous assertion in the presence of
		// heterogenous constraints: we now can constrain "u_i = c"
		// with no other u_j terms involved.
		//
	        // libmesh_assert (!constraint_row.empty());

	        for (DofConstraintRow::const_iterator
		       it=constraint_row.begin(); it != constraint_row.end();
		     ++it)
	          for (unsigned int j=0; j<elem_dofs.size(); j++)
		    if (elem_dofs[j] == it->first)
		      matrix(i,j) = -it->second;
	      }
	  }
    } // end if is constrained...

  STOP_LOG("constrain_elem_matrix()", "DofMap");
}



void DofMap::constrain_element_matrix_and_vector (DenseMatrix<Number>& matrix,
						  DenseVector<Number>& rhs,
						  std::vector<dof_id_type>& elem_dofs,
						  bool asymmetric_constraint_rows) const
{
  libmesh_assert_equal_to (elem_dofs.size(), matrix.m());
  libmesh_assert_equal_to (elem_dofs.size(), matrix.n());
  libmesh_assert_equal_to (elem_dofs.size(), rhs.size());

  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // The constrained matrix is built up as C^T K C.
  // The constrained RHS is built up as C^T F
  DenseMatrix<Number> C;

  this->build_constraint_matrix (C, elem_dofs);

  START_LOG("cnstrn_elem_mat_vec()", "DofMap");

  // It is possible that the matrix is not constrained at all.
  if ((C.m() == matrix.m()) &&
      (C.n() == elem_dofs.size())) // It the matrix is constrained
    {
      // Compute the matrix-matrix-matrix product C^T K C
      matrix.left_multiply_transpose  (C);
      matrix.right_multiply (C);


      libmesh_assert_equal_to (matrix.m(), matrix.n());
      libmesh_assert_equal_to (matrix.m(), elem_dofs.size());
      libmesh_assert_equal_to (matrix.n(), elem_dofs.size());


      for (unsigned int i=0; i<elem_dofs.size(); i++)
	if (this->is_constrained_dof(elem_dofs[i]))
	  {
	    for (unsigned int j=0; j<matrix.n(); j++)
	      matrix(i,j) = 0.;

	    // If the DOF is constrained
	    matrix(i,i) = 1.;

            // This will put a nonsymmetric entry in the constraint
            // row to ensure that the linear system produces the
            // correct value for the constrained DOF.
            if (asymmetric_constraint_rows)
              {
	        DofConstraints::const_iterator
	          pos = _dof_constraints.find(elem_dofs[i]);

	        libmesh_assert (pos != _dof_constraints.end());

	        const DofConstraintRow& constraint_row = pos->second.first;

// p refinement creates empty constraint rows
//	    libmesh_assert (!constraint_row.empty());

	        for (DofConstraintRow::const_iterator
		       it=constraint_row.begin(); it != constraint_row.end();
		     ++it)
	          for (unsigned int j=0; j<elem_dofs.size(); j++)
		    if (elem_dofs[j] == it->first)
		      matrix(i,j) = -it->second;
              }
	  }


      // Compute the matrix-vector product C^T F
      DenseVector<Number> old_rhs(rhs);

      // compute matrix/vector product
      C.vector_mult_transpose(rhs, old_rhs);
    } // end if is constrained...

  STOP_LOG("cnstrn_elem_mat_vec()", "DofMap");
}



void DofMap::heterogenously_constrain_element_matrix_and_vector
  (DenseMatrix<Number>& matrix,
   DenseVector<Number>& rhs,
   std::vector<dof_id_type>& elem_dofs,
   bool asymmetric_constraint_rows) const
{
  libmesh_assert_equal_to (elem_dofs.size(), matrix.m());
  libmesh_assert_equal_to (elem_dofs.size(), matrix.n());
  libmesh_assert_equal_to (elem_dofs.size(), rhs.size());

  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // The constrained matrix is built up as C^T K C.
  // The constrained RHS is built up as C^T (F - K H)
  DenseMatrix<Number> C;
  DenseVector<Number> H;

  this->build_constraint_matrix_and_vector (C, H, elem_dofs);

  START_LOG("hetero_cnstrn_elem_mat_vec()", "DofMap");

  // It is possible that the matrix is not constrained at all.
  if ((C.m() == matrix.m()) &&
      (C.n() == elem_dofs.size())) // It the matrix is constrained
    {
      // Compute matrix/vector product K H
      DenseVector<Number> KH;
      matrix.vector_mult(KH, H);

      // Compute the matrix-vector product C^T (F - KH)
      DenseVector<Number> F_minus_KH(rhs);
      F_minus_KH -= KH;
      C.vector_mult_transpose(rhs, F_minus_KH);

      // Compute the matrix-matrix-matrix product C^T K C
      matrix.left_multiply_transpose  (C);
      matrix.right_multiply (C);

      libmesh_assert_equal_to (matrix.m(), matrix.n());
      libmesh_assert_equal_to (matrix.m(), elem_dofs.size());
      libmesh_assert_equal_to (matrix.n(), elem_dofs.size());

      for (unsigned int i=0; i<elem_dofs.size(); i++)
	if (this->is_constrained_dof(elem_dofs[i]))
	  {
	    for (unsigned int j=0; j<matrix.n(); j++)
	      matrix(i,j) = 0.;

	    // If the DOF is constrained
	    matrix(i,i) = 1.;

            // This will put a nonsymmetric entry in the constraint
            // row to ensure that the linear system produces the
            // correct value for the constrained DOF.
            if (asymmetric_constraint_rows)
              {
	        DofConstraints::const_iterator
	          pos = _dof_constraints.find(elem_dofs[i]);

	        libmesh_assert (pos != _dof_constraints.end());

	        const DofConstraintRow& constraint_row = pos->second.first;

// p refinement creates empty constraint rows
//	    libmesh_assert (!constraint_row.empty());

	        for (DofConstraintRow::const_iterator
		       it=constraint_row.begin(); it != constraint_row.end();
		     ++it)
	          for (unsigned int j=0; j<elem_dofs.size(); j++)
		    if (elem_dofs[j] == it->first)
		      matrix(i,j) = -it->second;

		rhs(i) = pos->second.second;
              }
	    else
              rhs(i) = 0.;
	  }

    } // end if is constrained...

  STOP_LOG("hetero_cnstrn_elem_mat_vec()", "DofMap");
}



void DofMap::constrain_element_matrix (DenseMatrix<Number>& matrix,
				       std::vector<dof_id_type>& row_dofs,
				       std::vector<dof_id_type>& col_dofs,
				       bool asymmetric_constraint_rows) const
{
  libmesh_assert_equal_to (row_dofs.size(), matrix.m());
  libmesh_assert_equal_to (col_dofs.size(), matrix.n());

  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // The constrained matrix is built up as R^T K C.
  DenseMatrix<Number> R;
  DenseMatrix<Number> C;

  // Safeguard against the user passing us the same
  // object for row_dofs and col_dofs.  If that is done
  // the calls to build_matrix would fail
  std::vector<dof_id_type> orig_row_dofs(row_dofs);
  std::vector<dof_id_type> orig_col_dofs(col_dofs);

  this->build_constraint_matrix (R, orig_row_dofs);
  this->build_constraint_matrix (C, orig_col_dofs);

  START_LOG("constrain_elem_matrix()", "DofMap");

  row_dofs = orig_row_dofs;
  col_dofs = orig_col_dofs;


  // It is possible that the matrix is not constrained at all.
  if ((R.m() == matrix.m()) &&
      (R.n() == row_dofs.size()) &&
      (C.m() == matrix.n()) &&
      (C.n() == col_dofs.size())) // If the matrix is constrained
    {
      // K_constrained = R^T K C
      matrix.left_multiply_transpose  (R);
      matrix.right_multiply (C);


      libmesh_assert_equal_to (matrix.m(), row_dofs.size());
      libmesh_assert_equal_to (matrix.n(), col_dofs.size());


      for (unsigned int i=0; i<row_dofs.size(); i++)
	if (this->is_constrained_dof(row_dofs[i]))
	  {
	    for (unsigned int j=0; j<matrix.n(); j++)
            {
              if(row_dofs[i] != col_dofs[j])
                matrix(i,j) = 0.;
              else // If the DOF is constrained
                matrix(i,j) = 1.;
            }

            if (asymmetric_constraint_rows)
              {
	        DofConstraints::const_iterator
	          pos = _dof_constraints.find(row_dofs[i]);

	        libmesh_assert (pos != _dof_constraints.end());

	        const DofConstraintRow& constraint_row = pos->second.first;

	        libmesh_assert (!constraint_row.empty());

	        for (DofConstraintRow::const_iterator
		       it=constraint_row.begin(); it != constraint_row.end();
		     ++it)
	          for (unsigned int j=0; j<col_dofs.size(); j++)
		    if (col_dofs[j] == it->first)
		      matrix(i,j) = -it->second;
              }
	  }
    } // end if is constrained...

  STOP_LOG("constrain_elem_matrix()", "DofMap");
}



void DofMap::constrain_element_vector (DenseVector<Number>&       rhs,
				       std::vector<dof_id_type>& row_dofs,
				       bool) const
{
  libmesh_assert_equal_to (rhs.size(), row_dofs.size());

  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // The constrained RHS is built up as R^T F.
  DenseMatrix<Number> R;

  this->build_constraint_matrix (R, row_dofs);

  START_LOG("constrain_elem_vector()", "DofMap");

  // It is possible that the vector is not constrained at all.
  if ((R.m() == rhs.size()) &&
      (R.n() == row_dofs.size())) // if the RHS is constrained
    {
      // Compute the matrix-vector product
      DenseVector<Number> old_rhs(rhs);
      R.vector_mult_transpose(rhs, old_rhs);

      libmesh_assert_equal_to (row_dofs.size(), rhs.size());

      for (unsigned int i=0; i<row_dofs.size(); i++)
	if (this->is_constrained_dof(row_dofs[i]))
	  {
	    // If the DOF is constrained
            DofConstraints::const_iterator
                  pos = _dof_constraints.find(row_dofs[i]);

            libmesh_assert (pos != _dof_constraints.end());

	    rhs(i) = 0;
	  }
    } // end if the RHS is constrained.

  STOP_LOG("constrain_elem_vector()", "DofMap");
}



void DofMap::constrain_element_dyad_matrix (DenseVector<Number>& v,
					    DenseVector<Number>& w,
					    std::vector<dof_id_type>& row_dofs,
					    bool) const
{
  libmesh_assert_equal_to (v.size(), row_dofs.size());
  libmesh_assert_equal_to (w.size(), row_dofs.size());

  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // The constrained RHS is built up as R^T F.
  DenseMatrix<Number> R;

  this->build_constraint_matrix (R, row_dofs);

  START_LOG("cnstrn_elem_dyad_mat()", "DofMap");

  // It is possible that the vector is not constrained at all.
  if ((R.m() == v.size()) &&
      (R.n() == row_dofs.size())) // if the RHS is constrained
    {
      // Compute the matrix-vector products
      DenseVector<Number> old_v(v);
      DenseVector<Number> old_w(w);

      // compute matrix/vector product
      R.vector_mult_transpose(v, old_v);
      R.vector_mult_transpose(w, old_w);

      libmesh_assert_equal_to (row_dofs.size(), v.size());
      libmesh_assert_equal_to (row_dofs.size(), w.size());

      /* Constrain only v, not w.  */

      for (unsigned int i=0; i<row_dofs.size(); i++)
	if (this->is_constrained_dof(row_dofs[i]))
	  {
	    // If the DOF is constrained
            DofConstraints::const_iterator
                  pos = _dof_constraints.find(row_dofs[i]);

            libmesh_assert (pos != _dof_constraints.end());

	    v(i) = 0;
	  }
    } // end if the RHS is constrained.

  STOP_LOG("cnstrn_elem_dyad_mat()", "DofMap");
}



void DofMap::constrain_nothing (std::vector<dof_id_type>& dofs) const
{
  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // All the work is done by \p build_constraint_matrix.  We just need
  // a dummy matrix.
  DenseMatrix<Number> R;
  this->build_constraint_matrix (R, dofs);
}



void DofMap::enforce_constraints_exactly (const System &system,
                                          NumericVector<Number> *v,
                                          bool homogeneous) const
{
  parallel_object_only();

  if (!this->n_constrained_dofs())
    return;

  START_LOG("enforce_constraints_exactly()","DofMap");

  if (!v)
    v = system.solution.get();

  NumericVector<Number> *v_local  = NULL; // will be initialized below
  NumericVector<Number> *v_global = NULL; // will be initialized below
  AutoPtr<NumericVector<Number> > v_built;
  if (v->type() == SERIAL)
    {
      v_built = NumericVector<Number>::build(this->comm());
      v_built->init(this->n_dofs(), this->n_local_dofs(), true, PARALLEL);
      v_built->close();

      for (dof_id_type i=v_built->first_local_index();
           i<v_built->last_local_index(); i++)
        v_built->set(i, (*v)(i));
      v_built->close();
      v_global = v_built.get();

      v_local = v;
      libmesh_assert (v_local->closed());
    }
  else if (v->type() == PARALLEL)
    {
      v_built = NumericVector<Number>::build(this->comm());
      v_built->init (v->size(), v->size(), true, SERIAL);
      v->localize(*v_built);
      v_built->close();
      v_local = v_built.get();

      v_global = v;
    }
  else if (v->type() == GHOSTED)
    {
      v_local = v;
      v_global = v;
    }
  else // unknown v->type()
    {
      libMesh::err << "ERROR: Unknown v->type() == " << v->type()
		    << std::endl;
      libmesh_error();
    }

  // We should never hit these asserts because we should error-out in
  // else clause above.  Just to be sure we don't try to use v_local
  // and v_global uninitialized...
  libmesh_assert(v_local);
  libmesh_assert(v_global);
  libmesh_assert_equal_to (this, &(system.get_dof_map()));

  DofConstraints::const_iterator c_it = _dof_constraints.begin();
  const DofConstraints::const_iterator c_end = _dof_constraints.end();

  for ( ; c_it != c_end; ++c_it)
    {
      dof_id_type constrained_dof = c_it->first;
      if (constrained_dof < this->first_dof() ||
          constrained_dof >= this->end_dof())
        continue;

      const DofConstraintRow constraint_row = c_it->second.first;

      Number exact_value = homogeneous ?
        0 : c_it->second.second;
      for (DofConstraintRow::const_iterator
	   j=constraint_row.begin(); j != constraint_row.end();
	   ++j)
        exact_value += j->second * (*v_local)(j->first);

      v_global->set(constrained_dof, exact_value);
    }

  // If the old vector was serial, we probably need to send our values
  // to other processors
  if (v->type() == SERIAL)
    {
#ifndef NDEBUG
      v_global->close();
#endif
      v_global->localize (*v);
    }
  v->close();

  STOP_LOG("enforce_constraints_exactly()","DofMap");
}



std::pair<Real, Real>
DofMap::max_constraint_error (const System &system,
                              NumericVector<Number> *v) const
{
  if (!v)
    v = system.solution.get();
  NumericVector<Number> &vec = *v;

  // We'll assume the vector is closed
  libmesh_assert (vec.closed());

  Real max_absolute_error = 0., max_relative_error = 0.;

  const MeshBase &mesh = system.get_mesh();

  libmesh_assert_equal_to (this, &(system.get_dof_map()));

  // indices on each element
  std::vector<dof_id_type> local_dof_indices;

  MeshBase::const_element_iterator       elem_it  =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator elem_end =
    mesh.active_local_elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    {
      const Elem* elem = *elem_it;

      this->dof_indices(elem, local_dof_indices);
      std::vector<dof_id_type> raw_dof_indices = local_dof_indices;

      // Constraint matrix for each element
      DenseMatrix<Number> C;

      this->build_constraint_matrix (C, local_dof_indices);

      // Continue if the element is unconstrained
      if (!C.m())
        continue;

      libmesh_assert_equal_to (C.m(), raw_dof_indices.size());
      libmesh_assert_equal_to (C.n(), local_dof_indices.size());

      for (unsigned int i=0; i!=C.m(); ++i)
        {
          // Recalculate any constrained dof owned by this processor
          dof_id_type global_dof = raw_dof_indices[i];
          if (this->is_constrained_dof(global_dof) &&
              global_dof >= vec.first_local_index() &&
              global_dof < vec.last_local_index())
          {
            DofConstraints::const_iterator
                  pos = _dof_constraints.find(global_dof);

            libmesh_assert (pos != _dof_constraints.end());

            Number exact_value = pos->second.second;
            for (unsigned int j=0; j!=C.n(); ++j)
              {
                if (local_dof_indices[j] != global_dof)
                  exact_value += C(i,j) *
                    vec(local_dof_indices[j]);
              }

            max_absolute_error = std::max(max_absolute_error,
              std::abs(vec(global_dof) - exact_value));
            max_relative_error = std::max(max_relative_error,
              std::abs(vec(global_dof) - exact_value)
              / std::abs(exact_value));
          }
        }
    }

  return std::pair<Real, Real>(max_absolute_error, max_relative_error);
}



void DofMap::build_constraint_matrix (DenseMatrix<Number>& C,
				      std::vector<dof_id_type>& elem_dofs,
				      const bool called_recursively) const
{
  if (!called_recursively) START_LOG("build_constraint_matrix()", "DofMap");

  // Create a set containing the DOFs we already depend on
  typedef std::set<dof_id_type> RCSet;
  RCSet dof_set;

  bool we_have_constraints = false;

  // Next insert any other dofs the current dofs might be constrained
  // in terms of.  Note that in this case we may not be done: Those
  // may in turn depend on others.  So, we need to repeat this process
  // in that case until the system depends only on unconstrained
  // degrees of freedom.
  for (unsigned int i=0; i<elem_dofs.size(); i++)
    if (this->is_constrained_dof(elem_dofs[i]))
      {
        we_have_constraints = true;

	// If the DOF is constrained
	DofConstraints::const_iterator
	  pos = _dof_constraints.find(elem_dofs[i]);

	libmesh_assert (pos != _dof_constraints.end());

	const DofConstraintRow& constraint_row = pos->second.first;

// Constraint rows in p refinement may be empty
//	libmesh_assert (!constraint_row.empty());

	for (DofConstraintRow::const_iterator
	       it=constraint_row.begin(); it != constraint_row.end();
	     ++it)
	  dof_set.insert (it->first);
      }

  // May be safe to return at this point
  // (but remember to stop the perflog)
  if (!we_have_constraints)
    {
      STOP_LOG("build_constraint_matrix()", "DofMap");
      return;
    }

  for (unsigned int i=0; i != elem_dofs.size(); ++i)
    dof_set.erase (elem_dofs[i]);

  // If we added any DOFS then we need to do this recursively.
  // It is possible that we just added a DOF that is also
  // constrained!
  //
  // Also, we need to handle the special case of an element having DOFs
  // constrained in terms of other, local DOFs
  if (!dof_set.empty() ||  // case 1: constrained in terms of other DOFs
      !called_recursively) // case 2: constrained in terms of our own DOFs
    {
      const unsigned int old_size =
        libmesh_cast_int<unsigned int>(elem_dofs.size());

      // Add new dependency dofs to the end of the current dof set
      elem_dofs.insert(elem_dofs.end(),
		       dof_set.begin(), dof_set.end());

      // Now we can build the constraint matrix.
      // Note that resize also zeros for a DenseMatrix<Number>.
      C.resize (old_size,
		libmesh_cast_int<unsigned int>(elem_dofs.size()));

      // Create the C constraint matrix.
      for (unsigned int i=0; i != old_size; i++)
	if (this->is_constrained_dof(elem_dofs[i]))
	  {
	    // If the DOF is constrained
	    DofConstraints::const_iterator
	      pos = _dof_constraints.find(elem_dofs[i]);

	    libmesh_assert (pos != _dof_constraints.end());

	    const DofConstraintRow& constraint_row = pos->second.first;

// p refinement creates empty constraint rows
//	    libmesh_assert (!constraint_row.empty());

	    for (DofConstraintRow::const_iterator
		   it=constraint_row.begin(); it != constraint_row.end();
		 ++it)
	      for (unsigned int j=0; j != elem_dofs.size(); j++)
		if (elem_dofs[j] == it->first)
		  C(i,j) = it->second;
	  }
	else
	  {
            C(i,i) = 1.;
	  }

      // May need to do this recursively.  It is possible
      // that we just replaced a constrained DOF with another
      // constrained DOF.
      DenseMatrix<Number> Cnew;

      this->build_constraint_matrix (Cnew, elem_dofs, true);

      if ((C.n() == Cnew.m()) &&
	  (Cnew.n() == elem_dofs.size())) // If the constraint matrix
	C.right_multiply(Cnew);           // is constrained...

      libmesh_assert_equal_to (C.n(), elem_dofs.size());
    }

  if (!called_recursively) STOP_LOG("build_constraint_matrix()", "DofMap");
}



void DofMap::build_constraint_matrix_and_vector
  (DenseMatrix<Number>& C,
   DenseVector<Number>& H,
   std::vector<dof_id_type>& elem_dofs,
   const bool called_recursively) const
{
  if (!called_recursively)
    START_LOG("build_constraint_matrix_and_vector()", "DofMap");

  // Create a set containing the DOFs we already depend on
  typedef std::set<dof_id_type> RCSet;
  RCSet dof_set;

  bool we_have_constraints = false;

  // Next insert any other dofs the current dofs might be constrained
  // in terms of.  Note that in this case we may not be done: Those
  // may in turn depend on others.  So, we need to repeat this process
  // in that case until the system depends only on unconstrained
  // degrees of freedom.
  for (unsigned int i=0; i<elem_dofs.size(); i++)
    if (this->is_constrained_dof(elem_dofs[i]))
      {
        we_have_constraints = true;

	// If the DOF is constrained
	DofConstraints::const_iterator
	  pos = _dof_constraints.find(elem_dofs[i]);

	libmesh_assert (pos != _dof_constraints.end());

	const DofConstraintRow& constraint_row = pos->second.first;

// Constraint rows in p refinement may be empty
//	libmesh_assert (!constraint_row.empty());

	for (DofConstraintRow::const_iterator
	       it=constraint_row.begin(); it != constraint_row.end();
	     ++it)
	  dof_set.insert (it->first);
      }

  // May be safe to return at this point
  // (but remember to stop the perflog)
  if (!we_have_constraints)
    {
      STOP_LOG("build_constraint_matrix_and_vector()", "DofMap");
      return;
    }

  for (unsigned int i=0; i != elem_dofs.size(); ++i)
    dof_set.erase (elem_dofs[i]);

  // If we added any DOFS then we need to do this recursively.
  // It is possible that we just added a DOF that is also
  // constrained!
  //
  // Also, we need to handle the special case of an element having DOFs
  // constrained in terms of other, local DOFs
  if (!dof_set.empty() ||  // case 1: constrained in terms of other DOFs
      !called_recursively) // case 2: constrained in terms of our own DOFs
    {
      const unsigned int old_size =
        libmesh_cast_int<unsigned int>(elem_dofs.size());

      // Add new dependency dofs to the end of the current dof set
      elem_dofs.insert(elem_dofs.end(),
		       dof_set.begin(), dof_set.end());

      // Now we can build the constraint matrix and vector.
      // Note that resize also zeros for a DenseMatrix and DenseVector
      C.resize (old_size,
		libmesh_cast_int<unsigned int>(elem_dofs.size()));
      H.resize (old_size);

      // Create the C constraint matrix.
      for (unsigned int i=0; i != old_size; i++)
	if (this->is_constrained_dof(elem_dofs[i]))
	  {
	    // If the DOF is constrained
	    DofConstraints::const_iterator
	      pos = _dof_constraints.find(elem_dofs[i]);

	    libmesh_assert (pos != _dof_constraints.end());

	    const DofConstraintRow& constraint_row = pos->second.first;

// p refinement creates empty constraint rows
//	    libmesh_assert (!constraint_row.empty());

	    for (DofConstraintRow::const_iterator
		   it=constraint_row.begin(); it != constraint_row.end();
		 ++it)
	      for (unsigned int j=0; j != elem_dofs.size(); j++)
		if (elem_dofs[j] == it->first)
		  C(i,j) = it->second;

            H(i) = pos->second.second;
	  }
	else
	  {
            C(i,i) = 1.;
	  }

      // May need to do this recursively.  It is possible
      // that we just replaced a constrained DOF with another
      // constrained DOF.
      DenseMatrix<Number> Cnew;
      DenseVector<Number> Hnew;

      this->build_constraint_matrix_and_vector (Cnew, Hnew, elem_dofs, true);

      if ((C.n() == Cnew.m()) &&          // If the constraint matrix
	  (Cnew.n() == elem_dofs.size())) // is constrained...
        {
	  C.right_multiply(Cnew);

          // If x = Cy + h and y = Dz + g
          // Then x = (CD)z + (Cg + h)
          C.vector_mult_add(H, 1, Hnew);
        }

      libmesh_assert_equal_to (C.n(), elem_dofs.size());
    }

  if (!called_recursively)
    STOP_LOG("build_constraint_matrix_and_vector()", "DofMap");
}


void DofMap::allgather_recursive_constraints(MeshBase& mesh)
{
  // This function must be run on all processors at once
  parallel_object_only();

  // Return immediately if there's nothing to gather
  if (this->n_processors() == 1)
    return;

  // We might get to return immediately if none of the processors
  // found any constraints
  unsigned int has_constraints = !_dof_constraints.empty()
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
                                 || !_node_constraints.empty()
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS
                                 ;
  this->comm().max(has_constraints);
  if (!has_constraints)
    return;

  // We might have calculated constraints for constrained dofs
  // which have support on other processors.
  // Push these out first.
  {
  std::vector<std::set<dof_id_type> > pushed_ids(this->n_processors());

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  std::vector<std::set<dof_id_type> > pushed_node_ids(this->n_processors());
#endif

  MeshBase::element_iterator
    foreign_elem_it  = mesh.active_not_local_elements_begin(),
    foreign_elem_end = mesh.active_not_local_elements_end();

  // Collect the constraints to push to each processor
  for (; foreign_elem_it != foreign_elem_end; ++foreign_elem_it)
    {
      Elem *elem = *foreign_elem_it;

      std::vector<dof_id_type> my_dof_indices;
      this->dof_indices (elem, my_dof_indices);

      for (unsigned int i=0; i != my_dof_indices.size(); ++i)
        if (this->is_constrained_dof(my_dof_indices[i]))
          pushed_ids[elem->processor_id()].insert(my_dof_indices[i]);

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
      for (unsigned int n=0; n != elem->n_nodes(); ++n)
        if (this->is_constrained_node(elem->get_node(n)))
          pushed_node_ids[elem->processor_id()].insert(elem->node(n));
#endif
    }

  // Now trade constraint rows
  for (processor_id_type p = 0; p != this->n_processors(); ++p)
    {
      // Push to processor procup while receiving from procdown
      processor_id_type procup = (this->processor_id() + p) %
                                  this->n_processors();
      processor_id_type procdown = (this->n_processors() +
                                    this->processor_id() - p) %
                                    this->n_processors();

      // Pack the dof constraint rows and rhs's to push to procup
      const std::size_t pushed_ids_size = pushed_ids[procup].size();
      std::vector<std::vector<dof_id_type> > pushed_keys(pushed_ids_size);
      std::vector<std::vector<Real> > pushed_vals(pushed_ids_size);
      std::vector<Number> pushed_rhss(pushed_ids_size);
      std::set<dof_id_type>::const_iterator it = pushed_ids[procup].begin();
      for (std::size_t i = 0; it != pushed_ids[procup].end();
           ++i, ++it)
        {
          const dof_id_type pushed_id = *it;
          DofConstraintRow &row = _dof_constraints[pushed_id].first;
	  std::size_t row_size = row.size();
          pushed_keys[i].reserve(row_size);
          pushed_vals[i].reserve(row_size);
          for (DofConstraintRow::iterator j = row.begin();
               j != row.end(); ++j)
            {
              pushed_keys[i].push_back(j->first);
              pushed_vals[i].push_back(j->second);
            }
          pushed_rhss[i] = _dof_constraints[pushed_id].second;
        }

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
      // Pack the node constraint rows to push to procup
      const std::size_t pushed_nodes_size = pushed_node_ids[procup].size();
      std::vector<std::vector<dof_id_type> > pushed_node_keys(pushed_nodes_size);
      std::vector<std::vector<Real> > pushed_node_vals(pushed_nodes_size);
      std::vector<Point> pushed_node_offsets(pushed_nodes_size);
      std::set<dof_id_type>::const_iterator node_it = pushed_node_ids[procup].begin();
      for (std::size_t i = 0; node_it != pushed_node_ids[procup].end();
           ++i, ++node_it)
        {
          const Node *node = mesh.node_ptr(*node_it);
          NodeConstraintRow &row = _node_constraints[node].first;
	  std::size_t row_size = row.size();
          pushed_node_keys[i].reserve(row_size);
          pushed_node_vals[i].reserve(row_size);
          for (NodeConstraintRow::iterator j = row.begin();
               j != row.end(); ++j)
            {
              pushed_node_keys[i].push_back(j->first->id());
              pushed_node_vals[i].push_back(j->second);
            }
          pushed_node_offsets[i] = _node_constraints[node].second;
        }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

      // Trade pushed dof constraint rows
      std::vector<dof_id_type> pushed_ids_from_me
        (pushed_ids[procup].begin(), pushed_ids[procup].end());
      std::vector<dof_id_type> pushed_ids_to_me;
      std::vector<std::vector<dof_id_type> > pushed_keys_to_me;
      std::vector<std::vector<Real> > pushed_vals_to_me;
      std::vector<Number> pushed_rhss_to_me;
      this->comm().send_receive(procup, pushed_ids_from_me,
					procdown, pushed_ids_to_me);
      this->comm().send_receive(procup, pushed_keys,
					procdown, pushed_keys_to_me);
      this->comm().send_receive(procup, pushed_vals,
					procdown, pushed_vals_to_me);
      this->comm().send_receive(procup, pushed_rhss,
					procdown, pushed_rhss_to_me);
      libmesh_assert_equal_to (pushed_ids_to_me.size(), pushed_keys_to_me.size());
      libmesh_assert_equal_to (pushed_ids_to_me.size(), pushed_vals_to_me.size());
      libmesh_assert_equal_to (pushed_ids_to_me.size(), pushed_rhss_to_me.size());

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
      // Trade pushed node constraint rows
      std::vector<dof_id_type> pushed_node_ids_from_me
        (pushed_node_ids[procup].begin(), pushed_node_ids[procup].end());
      std::vector<dof_id_type> pushed_node_ids_to_me;
      std::vector<std::vector<dof_id_type> > pushed_node_keys_to_me;
      std::vector<std::vector<Real> > pushed_node_vals_to_me;
      std::vector<Point> pushed_node_offsets_to_me;
      this->comm().send_receive(procup, pushed_node_ids_from_me,
					procdown, pushed_node_ids_to_me);
      this->comm().send_receive(procup, pushed_node_keys,
					procdown, pushed_node_keys_to_me);
      this->comm().send_receive(procup, pushed_node_vals,
					procdown, pushed_node_vals_to_me);
      this->comm().send_receive(procup, pushed_node_offsets,
					procdown, pushed_node_offsets_to_me);

      // Note that we aren't doing a send_receive on the Nodes
      // themselves.  At this point we should only be pushing out
      // "raw" constraints, and there should be no
      // constrained-by-constrained-by-etc. situations that could
      // involve non-semilocal nodes.

      libmesh_assert_equal_to (pushed_node_ids_to_me.size(), pushed_node_keys_to_me.size());
      libmesh_assert_equal_to (pushed_node_ids_to_me.size(), pushed_node_vals_to_me.size());
      libmesh_assert_equal_to (pushed_node_ids_to_me.size(), pushed_node_offsets_to_me.size());
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

      // Add the dof constraints that I've been sent
      for (std::size_t i = 0; i != pushed_ids_to_me.size(); ++i)
        {
          libmesh_assert_equal_to (pushed_keys_to_me[i].size(), pushed_vals_to_me[i].size());

          dof_id_type constrained = pushed_ids_to_me[i];

          // If we don't already have a constraint for this dof,
          // add the one we were sent
          if (!this->is_constrained_dof(constrained))
            {
              DofConstraintRow &row = _dof_constraints[constrained].first;
              for (std::size_t j = 0; j != pushed_keys_to_me[i].size(); ++j)
                {
                  row[pushed_keys_to_me[i][j]] = pushed_vals_to_me[i][j];
                }
              _dof_constraints[constrained].second = pushed_rhss_to_me[i];
            }
        }

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
      // Add the node constraints that I've been sent
      for (std::size_t i = 0; i != pushed_node_ids_to_me.size(); ++i)
        {
          libmesh_assert_equal_to (pushed_node_keys_to_me[i].size(), pushed_node_vals_to_me[i].size());

          dof_id_type constrained_id = pushed_node_ids_to_me[i];

          // If we don't already have a constraint for this node,
          // add the one we were sent
          const Node *constrained = mesh.node_ptr(constrained_id);
          if (!this->is_constrained_node(constrained))
            {
              NodeConstraintRow &row = _node_constraints[constrained].first;
              for (std::size_t j = 0; j != pushed_node_keys_to_me[i].size(); ++j)
                {
                  const Node *key_node = mesh.node_ptr(pushed_node_keys_to_me[i][j]);
                  libmesh_assert(key_node);
                  row[key_node] = pushed_node_vals_to_me[i][j];
                }
              _node_constraints[constrained].second = pushed_node_offsets_to_me[i];
            }
        }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS
    }
  }

  // Now start checking for any other constraints we need
  // to know about, requesting them recursively.

  // Create sets containing the DOFs and nodes we already depend on
  typedef std::set<dof_id_type> DoF_RCSet;
  DoF_RCSet unexpanded_dofs;

  for (DofConstraints::iterator i = _dof_constraints.begin();
	 i != _dof_constraints.end(); ++i)
    {
      unexpanded_dofs.insert(i->first);
    }

  typedef std::set<const Node *> Node_RCSet;
  Node_RCSet unexpanded_nodes;

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  for (NodeConstraints::iterator i = _node_constraints.begin();
	 i != _node_constraints.end(); ++i)
    {
      unexpanded_nodes.insert(i->first);
    }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

  // We have to keep recursing while the unexpanded set is
  // nonempty on *any* processor
  bool unexpanded_set_nonempty = !unexpanded_dofs.empty() ||
                                 !unexpanded_nodes.empty();
  this->comm().max(unexpanded_set_nonempty);

  while (unexpanded_set_nonempty)
    {
      // Let's make sure we don't lose sync in this loop.
      parallel_object_only();

      // Request sets
      DoF_RCSet   dof_request_set;
      Node_RCSet node_request_set;

      // Request sets to send to each processor
      std::vector<std::vector<dof_id_type> >
        requested_dof_ids(this->n_processors()),
        requested_node_ids(this->n_processors());

      // And the sizes of each
      std::vector<dof_id_type>
        dof_ids_on_proc(this->n_processors(), 0),
        node_ids_on_proc(this->n_processors(), 0);

      // Fill (and thereby sort and uniq!) the main request sets
      for (DoF_RCSet::iterator i = unexpanded_dofs.begin();
           i != unexpanded_dofs.end(); ++i)
        {
          DofConstraintRow &row = _dof_constraints[*i].first;
          for (DofConstraintRow::iterator j = row.begin();
               j != row.end(); ++j)
            {
              const dof_id_type constraining_dof = j->first;

              // If it's non-local and we haven't already got a
              // constraint for it, we might need to ask for one
              if (((constraining_dof < this->first_dof()) ||
                   (constraining_dof >= this->end_dof())) &&
                  !_dof_constraints.count(constraining_dof))
                dof_request_set.insert(constraining_dof);
            }
        }

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
      for (Node_RCSet::iterator i = unexpanded_nodes.begin();
           i != unexpanded_nodes.end(); ++i)
        {
          NodeConstraintRow &row = _node_constraints[*i].first;
          for (NodeConstraintRow::iterator j = row.begin();
               j != row.end(); ++j)
            {
              const Node* const node = j->first;
              libmesh_assert(node);

              // If it's non-local and we haven't already got a
              // constraint for it, we might need to ask for one
              if ((node->processor_id() != this->processor_id()) &&
                  !_node_constraints.count(node))
                node_request_set.insert(node);
            }
        }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

      // Clear the unexpanded constraint sets; we're about to expand
      // them
      unexpanded_dofs.clear();
      unexpanded_nodes.clear();

      // Count requests by processor
      processor_id_type proc_id = 0;
      for (DoF_RCSet::iterator i = dof_request_set.begin();
           i != dof_request_set.end(); ++i)
        {
          while (*i >= _end_df[proc_id])
            proc_id++;
          dof_ids_on_proc[proc_id]++;
        }

      for (Node_RCSet::iterator i = node_request_set.begin();
           i != node_request_set.end(); ++i)
        {
          libmesh_assert(*i);
          libmesh_assert_less ((*i)->processor_id(), this->n_processors());
          node_ids_on_proc[(*i)->processor_id()]++;
        }

      for (processor_id_type p = 0; p != this->n_processors(); ++p)
        {
          requested_dof_ids[p].reserve(dof_ids_on_proc[p]);
          requested_node_ids[p].reserve(node_ids_on_proc[p]);
        }

      // Prepare each processor's request set
      proc_id = 0;
      for (DoF_RCSet::iterator i = dof_request_set.begin();
           i != dof_request_set.end(); ++i)
        {
          while (*i >= _end_df[proc_id])
            proc_id++;
          requested_dof_ids[proc_id].push_back(*i);
        }

      for (Node_RCSet::iterator i = node_request_set.begin();
           i != node_request_set.end(); ++i)
        {
          requested_node_ids[(*i)->processor_id()].push_back((*i)->id());
        }

      // Now request constraint rows from other processors
      for (processor_id_type p=1; p != this->n_processors(); ++p)
        {
          // Trade my requests with processor procup and procdown
          processor_id_type procup = (this->processor_id() + p) %
                                      this->n_processors();
          processor_id_type procdown = (this->n_processors() +
                                        this->processor_id() - p) %
                                        this->n_processors();
          std::vector<dof_id_type> dof_request_to_fill,
                                   node_request_to_fill;
          this->comm().send_receive(procup, requested_dof_ids[procup],
					    procdown, dof_request_to_fill);
          this->comm().send_receive(procup, requested_node_ids[procup],
					    procdown, node_request_to_fill);

          // Fill those requests
          std::vector<std::vector<dof_id_type> > dof_row_keys(dof_request_to_fill.size()),
                                                 node_row_keys(node_request_to_fill.size());

          // FIXME - this could be an unordered set, given a
          // hash<pointers> specialization
          std::set<const Node*> nodes_requested;

          std::vector<std::vector<Real> > dof_row_vals(dof_request_to_fill.size()),
                                          node_row_vals(node_request_to_fill.size());
          std::vector<Number> dof_row_rhss(dof_request_to_fill.size());
          std::vector<Point>  node_row_rhss(node_request_to_fill.size());
          for (std::size_t i=0; i != dof_request_to_fill.size(); ++i)
            {
              dof_id_type constrained = dof_request_to_fill[i];
              if (_dof_constraints.count(constrained))
                {
                  DofConstraintRow &row = _dof_constraints[constrained].first;
		  std::size_t row_size = row.size();
                  dof_row_keys[i].reserve(row_size);
                  dof_row_vals[i].reserve(row_size);
                  for (DofConstraintRow::iterator j = row.begin();
                       j != row.end(); ++j)
                    {
                      dof_row_keys[i].push_back(j->first);
                      dof_row_vals[i].push_back(j->second);

                      // We should never have a 0 constraint
                      // coefficient; that's implicit via sparse
                      // constraint storage
                      libmesh_assert(j->second);
                    }
                  dof_row_rhss[i] = _dof_constraints[constrained].second;
                }
            }

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
          for (std::size_t i=0; i != node_request_to_fill.size(); ++i)
            {
              dof_id_type constrained_id = node_request_to_fill[i];
              const Node *constrained_node = mesh.node_ptr(constrained_id);
              if (_node_constraints.count(constrained_node))
                {
                  const NodeConstraintRow &row = _node_constraints[constrained_node].first;
		  std::size_t row_size = row.size();
                  node_row_keys[i].reserve(row_size);
                  node_row_vals[i].reserve(row_size);
                  for (NodeConstraintRow::const_iterator j = row.begin();
                       j != row.end(); ++j)
                    {
                      const Node* node = j->first;
                      node_row_keys[i].push_back(node->id());
                      node_row_vals[i].push_back(j->second);

                      // If we're not sure whether our send
		      // destination already has this node, let's give
		      // it a copy.
		      if (node->processor_id() != procdown)
                        nodes_requested.insert(node);

                      // We can have 0 nodal constraint
                      // coefficients, where no Lagrange constrant
                      // exists but non-Lagrange basis constraints
                      // might.
                      // libmesh_assert(j->second);
                    }
                  node_row_rhss[i] = _node_constraints[constrained_node].second;
                }
            }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

          // Trade back the results
          std::vector<std::vector<dof_id_type> > dof_filled_keys,
                                                 node_filled_keys;
          std::vector<std::vector<Real> > dof_filled_vals,
                                          node_filled_vals;
          std::vector<Number> dof_filled_rhss;
          std::vector<Point> node_filled_rhss;
          this->comm().send_receive(procdown, dof_row_keys,
					    procup, dof_filled_keys);
          this->comm().send_receive(procdown, dof_row_vals,
					    procup, dof_filled_vals);
          this->comm().send_receive(procdown, dof_row_rhss,
					    procup, dof_filled_rhss);
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
          this->comm().send_receive(procdown, node_row_keys,
					    procup, node_filled_keys);
          this->comm().send_receive(procdown, node_row_vals,
					    procup, node_filled_vals);
          this->comm().send_receive(procdown, node_row_rhss,
					    procup, node_filled_rhss);

          // Constraining nodes might not even exist on our subset of
          // a distributed mesh, so let's make them exist.
          this->comm().send_receive_packed_range
            (procdown, &mesh, nodes_requested.begin(), nodes_requested.end(),
             procup,   &mesh, mesh_inserter_iterator<Node>(mesh));

#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS
          libmesh_assert_equal_to (dof_filled_keys.size(), requested_dof_ids[procup].size());
          libmesh_assert_equal_to (dof_filled_vals.size(), requested_dof_ids[procup].size());
          libmesh_assert_equal_to (dof_filled_rhss.size(), requested_dof_ids[procup].size());
          libmesh_assert_equal_to (node_filled_keys.size(), requested_node_ids[procup].size());
          libmesh_assert_equal_to (node_filled_vals.size(), requested_node_ids[procup].size());
          libmesh_assert_equal_to (node_filled_rhss.size(), requested_node_ids[procup].size());

          // Add any new constraint rows we've found
          for (std::size_t i=0; i != requested_dof_ids[procup].size(); ++i)
            {
              libmesh_assert_equal_to (dof_filled_keys[i].size(), dof_filled_vals[i].size());
              // FIXME - what about empty p constraints!?
              if (!dof_filled_keys[i].empty())
                {
                  dof_id_type constrained = requested_dof_ids[procup][i];
                  DofConstraintRow &row = _dof_constraints[constrained].first;
                  for (std::size_t j = 0; j != dof_filled_keys[i].size(); ++j)
                    row[dof_filled_keys[i][j]] = dof_filled_vals[i][j];
                  _dof_constraints[constrained].second = dof_filled_rhss[i];

                  // And prepare to check for more recursive constraints
                  unexpanded_dofs.insert(constrained);
                }
            }

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
          for (std::size_t i=0; i != requested_node_ids[procup].size(); ++i)
            {
              libmesh_assert_equal_to (node_filled_keys[i].size(), node_filled_vals[i].size());
              if (!node_filled_keys[i].empty())
                {
                  dof_id_type constrained_id = requested_node_ids[procup][i];
                  const Node* constrained_node = mesh.node_ptr(constrained_id);
                  NodeConstraintRow &row = _node_constraints[constrained_node].first;
                  for (std::size_t j = 0; j != node_filled_keys[i].size(); ++j)
                    {
                      const Node* key_node =
                        mesh.node_ptr(node_filled_keys[i][j]);
                      libmesh_assert(key_node);
                      row[key_node] = node_filled_vals[i][j];
                    }
                  _node_constraints[constrained_node].second = node_filled_rhss[i];

                  // And prepare to check for more recursive constraints
                  unexpanded_nodes.insert(constrained_node);
                }
            }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS
        }

      // We have to keep recursing while the unexpanded set is
      // nonempty on *any* processor
      unexpanded_set_nonempty = !unexpanded_dofs.empty() ||
                                !unexpanded_nodes.empty();
      this->comm().max(unexpanded_set_nonempty);
    }
}



void DofMap::process_constraints (MeshBase& mesh)
{
  // With a parallelized Mesh, we've computed our local constraints,
  // but they may depend on non-local constraints that we'll need to
  // take into account.
  if (!mesh.is_serial())
    this->allgather_recursive_constraints(mesh);

  // Create a set containing the DOFs we already depend on
  typedef std::set<dof_id_type> RCSet;
  RCSet unexpanded_set;

  for (DofConstraints::iterator i = _dof_constraints.begin();
	 i != _dof_constraints.end(); ++i)
    unexpanded_set.insert(i->first);

  while (!unexpanded_set.empty())
    for (RCSet::iterator i = unexpanded_set.begin();
	 i != unexpanded_set.end(); /* nothing */)
      {
	// If the DOF is constrained
	DofConstraints::iterator
	  pos = _dof_constraints.find(*i);

	libmesh_assert (pos != _dof_constraints.end());

	DofConstraintRow& constraint_row = pos->second.first;

	Number& constraint_rhs = pos->second.second;

	std::vector<dof_id_type> constraints_to_expand;

	for (DofConstraintRow::const_iterator
	       it=constraint_row.begin(); it != constraint_row.end();
	     ++it)
	  if (it->first != *i && this->is_constrained_dof(it->first))
            {
              unexpanded_set.insert(it->first);
	      constraints_to_expand.push_back(it->first);
	    }

	for (std::size_t j = 0; j != constraints_to_expand.size();
	     ++j)
	  {
            dof_id_type expandable = constraints_to_expand[j];

            const Real this_coef = constraint_row[expandable];

	    DofConstraints::const_iterator
	      subpos = _dof_constraints.find(expandable);

	    libmesh_assert (subpos != _dof_constraints.end());

	    const DofConstraintRow& subconstraint_row = subpos->second.first;

	    for (DofConstraintRow::const_iterator
	           it=subconstraint_row.begin();
		   it != subconstraint_row.end(); ++it)
              {
		constraint_row[it->first] += it->second * this_coef;
	      }
            constraint_rhs += subpos->second.second * this_coef;

	    constraint_row.erase(expandable);
          }

	if (constraints_to_expand.empty())
	  unexpanded_set.erase(i++);
	else
	  i++;
      }

  // In parallel we can't guarantee that nodes/dofs which constrain
  // others are on processors which are aware of that constraint, yet
  // we need such awareness for sparsity pattern generation.  So send
  // other processors any constraints they might need to know about.
  if (!mesh.is_serial())
    this->scatter_constraints(mesh);

  // Now that we have our root constraint dependencies sorted out, add
  // them to the send_list
  this->add_constraints_to_send_list();
}


void DofMap::scatter_constraints(MeshBase& mesh)
{
  // At this point each processor with a constrained node knows
  // the corresponding constraint row, but we also need each processor
  // with a constrainer node to know the corresponding row(s).

  // This function must be run on all processors at once
  parallel_object_only();

  // Return immediately if there's nothing to gather
  if (this->n_processors() == 1)
    return;

  // We might get to return immediately if none of the processors
  // found any constraints
  unsigned int has_constraints = !_dof_constraints.empty()
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
                                 || !_node_constraints.empty()
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS
                                 ;
  this->comm().max(has_constraints);
  if (!has_constraints)
    return;

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  std::vector<std::set<dof_id_type> > pushed_node_ids(this->n_processors());
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

  std::vector<std::set<dof_id_type> > pushed_ids(this->n_processors());

  // Collect the dof constraints I need to push to each processor
  dof_id_type constrained_proc_id = 0;
  for (DofConstraints::iterator i = _dof_constraints.begin();
	 i != _dof_constraints.end(); ++i)
    {
      const dof_id_type constrained = i->first;
      while (constrained >= _end_df[constrained_proc_id])
        constrained_proc_id++;

      if (constrained_proc_id != this->processor_id())
        continue;

      DofConstraintRow &row = i->second.first;
      for (DofConstraintRow::iterator j = row.begin();
           j != row.end(); ++j)
        {
          const dof_id_type constraining = j->first;

          processor_id_type constraining_proc_id = 0;
          while (constraining >= _end_df[constraining_proc_id])
            constraining_proc_id++;

          if (constraining_proc_id != this->processor_id() &&
              constraining_proc_id != constrained_proc_id)
            pushed_ids[constraining_proc_id].insert(constrained);
        }
    }

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  // Collect the node constraints to push to each processor
  for (NodeConstraints::iterator i = _node_constraints.begin();
	 i != _node_constraints.end(); ++i)
    {
      const Node *constrained = i->first;

      if (constrained->processor_id() != this->processor_id())
        continue;

      NodeConstraintRow &row = i->second.first;
      for (NodeConstraintRow::iterator j = row.begin();
           j != row.end(); ++j)
        {
          const Node *constraining = j->first;

          if (constraining->processor_id() != this->processor_id() &&
              constraining->processor_id() != constrained->processor_id())
            pushed_node_ids[constraining->processor_id()].insert(constrained->id());
        }
    }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

  // Now trade constraint rows
  for (processor_id_type p = 0; p != this->n_processors(); ++p)
    {
      // Push to processor procup while receiving from procdown
      processor_id_type procup = (this->processor_id() + p) %
                                  this->n_processors();
      processor_id_type procdown = (this->n_processors() +
                                    this->processor_id() - p) %
                                    this->n_processors();

      // Pack the dof constraint rows and rhs's to push to procup
      const std::size_t pushed_ids_size = pushed_ids[procup].size();
      std::vector<std::vector<dof_id_type> > pushed_keys(pushed_ids_size);
      std::vector<std::vector<Real> > pushed_vals(pushed_ids_size);
      std::vector<Number> pushed_rhss(pushed_ids_size);

      std::set<dof_id_type>::const_iterator it;
      std::size_t push_i;
      for (push_i = 0, it = pushed_ids[procup].begin();
           it != pushed_ids[procup].end(); ++push_i, ++it)
        {
          const dof_id_type constrained = *it;
          DofConstraintRow &row = _dof_constraints[constrained].first;
	  std::size_t row_size = row.size();
          pushed_keys[push_i].reserve(row_size);
          pushed_vals[push_i].reserve(row_size);
          for (DofConstraintRow::iterator j = row.begin();
               j != row.end(); ++j)
            {
              pushed_keys[push_i].push_back(j->first);
              pushed_vals[push_i].push_back(j->second);
            }
          pushed_rhss[push_i] = _dof_constraints[constrained].second;
        }

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
      // Pack the node constraint rows to push to procup
      const std::size_t pushed_node_ids_size = pushed_node_ids[procup].size();
      std::vector<std::vector<dof_id_type> > pushed_node_keys(pushed_node_ids_size);
      std::vector<std::vector<Real> > pushed_node_vals(pushed_node_ids_size);
      std::vector<Point> pushed_node_offsets(pushed_node_ids_size);
      std::set<const Node*> pushed_nodes;

      for (push_i = 0, it = pushed_node_ids[procup].begin();
           it != pushed_node_ids[procup].end(); ++push_i, ++it)
        {
          const Node *constrained = mesh.node_ptr(*it);

          if (constrained->processor_id() != procdown)
            pushed_nodes.insert(constrained);

          NodeConstraintRow &row = _node_constraints[constrained].first;
	  std::size_t row_size = row.size();
          pushed_node_keys[push_i].reserve(row_size);
          pushed_node_vals[push_i].reserve(row_size);
          for (NodeConstraintRow::iterator j = row.begin();
               j != row.end(); ++j)
            {
              const Node* constraining = j->first;

              pushed_node_keys[push_i].push_back(constraining->id());
              pushed_node_vals[push_i].push_back(j->second);

              if (constraining->processor_id() != procup)
                pushed_nodes.insert(constraining);
            }
          pushed_node_offsets[push_i] = _node_constraints[constrained].second;
        }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

      // Trade pushed dof constraint rows
      std::vector<dof_id_type> pushed_ids_from_me
        (pushed_ids[procup].begin(), pushed_ids[procup].end());
      std::vector<dof_id_type> pushed_ids_to_me;
      std::vector<std::vector<dof_id_type> > pushed_keys_to_me;
      std::vector<std::vector<Real> > pushed_vals_to_me;
      std::vector<Number> pushed_rhss_to_me;
      this->comm().send_receive(procup, pushed_ids_from_me,
					procdown, pushed_ids_to_me);
      this->comm().send_receive(procup, pushed_keys,
					procdown, pushed_keys_to_me);
      this->comm().send_receive(procup, pushed_vals,
					procdown, pushed_vals_to_me);
      this->comm().send_receive(procup, pushed_rhss,
					procdown, pushed_rhss_to_me);
      libmesh_assert_equal_to (pushed_ids_to_me.size(), pushed_keys_to_me.size());
      libmesh_assert_equal_to (pushed_ids_to_me.size(), pushed_vals_to_me.size());
      libmesh_assert_equal_to (pushed_ids_to_me.size(), pushed_rhss_to_me.size());

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
      // Trade pushed node constraint rows
      std::vector<dof_id_type> pushed_node_ids_from_me
        (pushed_node_ids[procup].begin(), pushed_node_ids[procup].end());
      std::vector<dof_id_type> pushed_node_ids_to_me;
      std::vector<std::vector<dof_id_type> > pushed_node_keys_to_me;
      std::vector<std::vector<Real> > pushed_node_vals_to_me;
      std::vector<Point> pushed_node_offsets_to_me;
      this->comm().send_receive(procup, pushed_node_ids_from_me,
					procdown, pushed_node_ids_to_me);
      this->comm().send_receive(procup, pushed_node_keys,
					procdown, pushed_node_keys_to_me);
      this->comm().send_receive(procup, pushed_node_vals,
					procdown, pushed_node_vals_to_me);
      this->comm().send_receive(procup, pushed_node_offsets,
					procdown, pushed_node_offsets_to_me);

      // Constraining nodes might not even exist on our subset of
      // a distributed mesh, so let's make them exist.
      this->comm().send_receive_packed_range
        (procup, &mesh, pushed_nodes.begin(), pushed_nodes.end(),
         procdown, &mesh, mesh_inserter_iterator<Node>(mesh));

      libmesh_assert_equal_to (pushed_node_ids_to_me.size(), pushed_node_keys_to_me.size());
      libmesh_assert_equal_to (pushed_node_ids_to_me.size(), pushed_node_vals_to_me.size());
      libmesh_assert_equal_to (pushed_node_ids_to_me.size(), pushed_node_offsets_to_me.size());
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

      // Add the dof constraints that I've been sent
      for (std::size_t i = 0; i != pushed_ids_to_me.size(); ++i)
        {
          libmesh_assert_equal_to (pushed_keys_to_me[i].size(), pushed_vals_to_me[i].size());

          dof_id_type constrained = pushed_ids_to_me[i];

          // If we don't already have a constraint for this dof,
          // add the one we were sent
          if (!this->is_constrained_dof(constrained))
            {
              DofConstraintRow &row = _dof_constraints[constrained].first;
              for (std::size_t j = 0; j != pushed_keys_to_me[i].size(); ++j)
                {
                  row[pushed_keys_to_me[i][j]] = pushed_vals_to_me[i][j];
                }
              _dof_constraints[constrained].second = pushed_rhss_to_me[i];
            }
        }

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
      // Add the node constraints that I've been sent
      for (std::size_t i = 0; i != pushed_node_ids_to_me.size(); ++i)
        {
          libmesh_assert_equal_to (pushed_node_keys_to_me[i].size(), pushed_node_vals_to_me[i].size());

          dof_id_type constrained_id = pushed_node_ids_to_me[i];

          // If we don't already have a constraint for this node,
          // add the one we were sent
          const Node *constrained = mesh.node_ptr(constrained_id);
          if (!this->is_constrained_node(constrained))
            {
              NodeConstraintRow &row = _node_constraints[constrained].first;
              for (std::size_t j = 0; j != pushed_node_keys_to_me[i].size(); ++j)
                {
                  const Node *key_node = mesh.node_ptr(pushed_node_keys_to_me[i][j]);
                  libmesh_assert(key_node);
                  row[key_node] = pushed_node_vals_to_me[i][j];
                }
              _node_constraints[constrained].second = pushed_node_offsets_to_me[i];
            }
        }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS
    }

  // Next we need to push constraints to processors which don't own
  // the constrained dof, don't own the constraining dof, but own an
  // element supporting the constraining dof.
  //
  // We need to be able to quickly look up constrained dof ids by what
  // constrains them, so that we can handle the case where we see a
  // foreign element containing one of our constraining DoF ids and we
  // need to push that constraint.
  //
  // Getting distributed adaptive sparsity patterns right is hard.

  typedef std::map<dof_id_type, std::set<dof_id_type> > DofConstrainsMap;
  DofConstrainsMap dof_id_constrains;

  for (DofConstraints::iterator i = _dof_constraints.begin();
	 i != _dof_constraints.end(); ++i)
    {
      const dof_id_type constrained = i->first;
      DofConstraintRow &row = i->second.first;
      for (DofConstraintRow::iterator j = row.begin();
           j != row.end(); ++j)
        {
          const dof_id_type constraining = j->first;

          dof_id_type constraining_proc_id = 0;
          while (constraining >= _end_df[constraining_proc_id])
            constraining_proc_id++;

          if (constraining_proc_id == this->processor_id())
            dof_id_constrains[constraining].insert(constrained);
        }
    }

  // Loop over all foreign elements, find any supporting our
  // constrained dof indices.
  pushed_ids.clear();
  pushed_ids.resize(this->n_processors());

  MeshBase::const_element_iterator it = mesh.active_not_local_elements_begin(),
                                  end = mesh.active_not_local_elements_end();
  for (; it != end; ++it)
    {
      const Elem *elem = *it;

      std::vector<dof_id_type> my_dof_indices;
      this->dof_indices (elem, my_dof_indices);

      for (std::size_t i=0; i != my_dof_indices.size(); ++i)
        {
	  DofConstrainsMap::const_iterator dcmi =
            dof_id_constrains.find(my_dof_indices[i]);
          if (dcmi != dof_id_constrains.end())
            {
              for (DofConstrainsMap::mapped_type::const_iterator mti =
		     dcmi->second.begin();
                   mti != dcmi->second.end(); ++mti)
                {
                  const dof_id_type constrained = *mti;

                  dof_id_type the_constrained_proc_id = 0;
                  while (constrained >= _end_df[the_constrained_proc_id])
                    the_constrained_proc_id++;

                  const dof_id_type elemproc = elem->processor_id();
                  if (elemproc != the_constrained_proc_id)
                    pushed_ids[elemproc].insert(constrained);
                }
            }
        }
    }

  // One last trade of constraint rows
  for (processor_id_type p = 0; p != this->n_processors(); ++p)
    {
      // Push to processor procup while receiving from procdown
      processor_id_type procup = (this->processor_id() + p) %
                                  this->n_processors();
      processor_id_type procdown = (this->n_processors() +
                                    this->processor_id() - p) %
                                    this->n_processors();

      // Pack the dof constraint rows and rhs's to push to procup
      const std::size_t pushed_ids_size = pushed_ids[procup].size();
      std::vector<std::vector<dof_id_type> > pushed_keys(pushed_ids_size);
      std::vector<std::vector<Real> > pushed_vals(pushed_ids_size);
      std::vector<Number> pushed_rhss(pushed_ids_size);

      // As long as we're declaring them outside the loop, let's initialize them too!
      std::set<dof_id_type>::const_iterator pushed_ids_iter = pushed_ids[procup].begin();
      std::size_t push_i = 0;
      for ( ; pushed_ids_iter != pushed_ids[procup].end(); ++push_i, ++pushed_ids_iter)
        {
          const dof_id_type constrained = *pushed_ids_iter;
          DofConstraintRow &row = _dof_constraints[constrained].first;
	  std::size_t row_size = row.size();
          pushed_keys[push_i].reserve(row_size);
          pushed_vals[push_i].reserve(row_size);
          for (DofConstraintRow::iterator j = row.begin();
               j != row.end(); ++j)
            {
              pushed_keys[push_i].push_back(j->first);
              pushed_vals[push_i].push_back(j->second);
            }
          pushed_rhss[push_i] = _dof_constraints[constrained].second;
        }

      // Trade pushed dof constraint rows
      std::vector<dof_id_type> pushed_ids_from_me
        (pushed_ids[procup].begin(), pushed_ids[procup].end());
      std::vector<dof_id_type> pushed_ids_to_me;
      std::vector<std::vector<dof_id_type> > pushed_keys_to_me;
      std::vector<std::vector<Real> > pushed_vals_to_me;
      std::vector<Number> pushed_rhss_to_me;
      this->comm().send_receive(procup, pushed_ids_from_me,
					procdown, pushed_ids_to_me);
      this->comm().send_receive(procup, pushed_keys,
					procdown, pushed_keys_to_me);
      this->comm().send_receive(procup, pushed_vals,
					procdown, pushed_vals_to_me);
      this->comm().send_receive(procup, pushed_rhss,
					procdown, pushed_rhss_to_me);
      libmesh_assert_equal_to (pushed_ids_to_me.size(), pushed_keys_to_me.size());
      libmesh_assert_equal_to (pushed_ids_to_me.size(), pushed_vals_to_me.size());
      libmesh_assert_equal_to (pushed_ids_to_me.size(), pushed_rhss_to_me.size());

      // Add the dof constraints that I've been sent
      for (std::size_t i = 0; i != pushed_ids_to_me.size(); ++i)
        {
          libmesh_assert_equal_to (pushed_keys_to_me[i].size(), pushed_vals_to_me[i].size());

          dof_id_type constrained = pushed_ids_to_me[i];

          // If we don't already have a constraint for this dof,
          // add the one we were sent
          if (!this->is_constrained_dof(constrained))
            {
              DofConstraintRow &row = _dof_constraints[constrained].first;
              for (std::size_t j = 0; j != pushed_keys_to_me[i].size(); ++j)
                {
                  row[pushed_keys_to_me[i][j]] = pushed_vals_to_me[i][j];
                }
              _dof_constraints[constrained].second = pushed_rhss_to_me[i];
            }
        }
    }
}


void DofMap::add_constraints_to_send_list()
{
  // This function must be run on all processors at once
  parallel_object_only();

  // Return immediately if there's nothing to gather
  if (this->n_processors() == 1)
    return;

  // We might get to return immediately if none of the processors
  // found any constraints
  unsigned int has_constraints = !_dof_constraints.empty();
  this->comm().max(has_constraints);
  if (!has_constraints)
    return;

  for (DofConstraints::iterator i = _dof_constraints.begin();
	 i != _dof_constraints.end(); ++i)
    {
      dof_id_type constrained_dof = i->first;

      // We only need the dependencies of our own constrained dofs
      if (constrained_dof < this->first_dof() ||
          constrained_dof >= this->end_dof())
        continue;

      DofConstraintRow& constraint_row = i->second.first;
      for (DofConstraintRow::const_iterator
	   j=constraint_row.begin(); j != constraint_row.end();
	   ++j)
        {
          dof_id_type constraint_dependency = j->first;

          // No point in adding one of our own dofs to the send_list
          if (constraint_dependency >= this->first_dof() &&
              constraint_dependency < this->end_dof())
            continue;

          _send_list.push_back(constraint_dependency);
        }
    }
}



#endif // LIBMESH_ENABLE_CONSTRAINTS


#ifdef LIBMESH_ENABLE_AMR

void DofMap::constrain_p_dofs (unsigned int var,
                               const Elem *elem,
                               unsigned int s,
                               unsigned int p)
{
  // We're constraining dofs on elem which correspond to p refinement
  // levels above p - this only makes sense if elem's p refinement
  // level is above p.
  libmesh_assert_greater (elem->p_level(), p);
  libmesh_assert_less (s, elem->n_sides());

  const unsigned int sys_num = this->sys_number();
  const unsigned int dim = elem->dim();
  ElemType type = elem->type();
  FEType low_p_fe_type = this->variable_type(var);
  FEType high_p_fe_type = this->variable_type(var);
  low_p_fe_type.order = static_cast<Order>(low_p_fe_type.order + p);
  high_p_fe_type.order = static_cast<Order>(high_p_fe_type.order +
                                            elem->p_level());

  const unsigned int n_nodes = elem->n_nodes();
  for (unsigned int n = 0; n != n_nodes; ++n)
    if (elem->is_node_on_side(n, s))
      {
        const Node * const node = elem->get_node(n);
        const unsigned int low_nc =
	  FEInterface::n_dofs_at_node (dim, low_p_fe_type, type, n);
        const unsigned int high_nc =
	  FEInterface::n_dofs_at_node (dim, high_p_fe_type, type, n);

	// since we may be running this method concurretly
	// on multiple threads we need to acquire a lock
	// before modifying the _dof_constraints object.
	Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

        if (elem->is_vertex(n))
          {
	    // Add "this is zero" constraint rows for high p vertex
            // dofs
            for (unsigned int i = low_nc; i != high_nc; ++i)
              {
                _dof_constraints[node->dof_number(sys_num,var,i)].first.clear();
                _dof_constraints[node->dof_number(sys_num,var,i)].second = 0.;
              }
          }
        else
          {
            const unsigned int total_dofs = node->n_comp(sys_num, var);
            libmesh_assert_greater_equal (total_dofs, high_nc);
	    // Add "this is zero" constraint rows for high p
            // non-vertex dofs, which are numbered in reverse
            for (unsigned int j = low_nc; j != high_nc; ++j)
              {
                const unsigned int i = total_dofs - j - 1;
                _dof_constraints[node->dof_number(sys_num,var,i)].first.clear();
                _dof_constraints[node->dof_number(sys_num,var,i)].second = 0.;
              }
          }
      }
}

#endif // LIBMESH_ENABLE_AMR


#ifdef LIBMESH_ENABLE_DIRICHLET
void DofMap::add_dirichlet_boundary (const DirichletBoundary& dirichlet_boundary)
{
  _dirichlet_boundaries->push_back(new DirichletBoundary(dirichlet_boundary));
}


void DofMap::remove_dirichlet_boundary (const DirichletBoundary& boundary_to_remove)
{
  // Find a boundary condition matching the one to be removed
  std::vector<DirichletBoundary *>::iterator it = _dirichlet_boundaries->begin();
  std::vector<DirichletBoundary *>::iterator end = _dirichlet_boundaries->end();
  for (; it != end; ++it)
    {
      DirichletBoundary *bdy = *it;

      if ((bdy->b == boundary_to_remove.b) &&
	   bdy->variables == boundary_to_remove.variables)
	break;
    }

  // Delete it and remove it
  libmesh_assert (it != end);
  delete *it;
  _dirichlet_boundaries->erase(it);
}


DirichletBoundaries::~DirichletBoundaries()
{
  for (std::vector<DirichletBoundary *>::iterator it = begin(); it != end(); ++it)
    delete *it;
}

#endif // LIBMESH_ENABLE_DIRICHLET


#ifdef LIBMESH_ENABLE_PERIODIC

void DofMap::add_periodic_boundary (const PeriodicBoundaryBase& periodic_boundary)
{
  // See if we already have a periodic boundary associated myboundary...
  PeriodicBoundaryBase* existing_boundary = _periodic_boundaries->boundary(periodic_boundary.myboundary);

  if ( existing_boundary == NULL )
  {
    // ...if not, clone the input (and its inverse) and add them to the PeriodicBoundaries object
    PeriodicBoundaryBase* boundary = periodic_boundary.clone().release();
    PeriodicBoundaryBase* inverse_boundary = periodic_boundary.clone(PeriodicBoundaryBase::INVERSE).release();

    // _periodic_boundaries takes ownership of the pointers
    _periodic_boundaries->insert(std::make_pair(boundary->myboundary, boundary));
    _periodic_boundaries->insert(std::make_pair(inverse_boundary->myboundary, inverse_boundary));
  }
  else
  {
    // ...otherwise, merge this object's variable IDs with the existing boundary object's.
    existing_boundary->merge(periodic_boundary);

    // Do the same merging process for the inverse boundary.  Note: the inverse better already exist!
    PeriodicBoundaryBase* inverse_boundary = _periodic_boundaries->boundary(periodic_boundary.pairedboundary);
    libmesh_assert(inverse_boundary);
    inverse_boundary->merge(periodic_boundary);
  }
}




void DofMap::add_periodic_boundary (const PeriodicBoundaryBase& boundary, const PeriodicBoundaryBase& inverse_boundary)
{
  libmesh_assert_equal_to (boundary.myboundary, inverse_boundary.pairedboundary);
  libmesh_assert_equal_to (boundary.pairedboundary, inverse_boundary.myboundary);

  // Allocate copies on the heap.  The _periodic_boundaries object will manage this memory.
  // Note: this also means that the copy constructor for the PeriodicBoundary (or user class
  // derived therefrom) must be implemented!
  // PeriodicBoundary* p_boundary = new PeriodicBoundary(boundary);
  // PeriodicBoundary* p_inverse_boundary = new PeriodicBoundary(inverse_boundary);

  // We can't use normal copy construction since this leads to slicing with derived classes.
  // Use clone() "virtual constructor" instead.  But, this *requires* user to override the clone()
  // method.  Note also that clone() allocates memory.  In this case, the _periodic_boundaries object
  // takes responsibility for cleanup.
  PeriodicBoundaryBase* p_boundary = boundary.clone().release();
  PeriodicBoundaryBase* p_inverse_boundary = inverse_boundary.clone().release();

  // Add the periodic boundary and its inverse to the PeriodicBoundaries data structure.  The
  // PeriodicBoundaries data structure takes ownership of the pointers.
  _periodic_boundaries->insert(std::make_pair(p_boundary->myboundary, p_boundary));
  _periodic_boundaries->insert(std::make_pair(p_inverse_boundary->myboundary, p_inverse_boundary));
}


#endif


} // namespace libMesh
