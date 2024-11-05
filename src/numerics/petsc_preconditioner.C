// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC

// Local Includes
#include "libmesh/petsc_preconditioner.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/enum_preconditioner_type.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"

namespace libMesh
{

template <typename T>
PetscPreconditioner<T>::PetscPreconditioner (const libMesh::Parallel::Communicator & comm_in) :
  Preconditioner<T>(comm_in)
{}



template <typename T>
void PetscPreconditioner<T>::apply(const NumericVector<T> & x, NumericVector<T> & y)
{
  PetscVector<T> & x_pvec = cast_ref<PetscVector<T> &>(const_cast<NumericVector<T> &>(x));
  PetscVector<T> & y_pvec = cast_ref<PetscVector<T> &>(const_cast<NumericVector<T> &>(y));

  Vec x_vec = x_pvec.vec();
  Vec y_vec = y_pvec.vec();

  LibmeshPetscCall(PCApply(_pc, x_vec, y_vec));
}




template <typename T>
void PetscPreconditioner<T>::init ()
{
  libmesh_error_msg_if(!this->_matrix, "ERROR: No matrix set for PetscPreconditioner, but init() called");

  // Clear the preconditioner in case it has been created in the past
  if (!this->_is_initialized)
    {
      // Should probably use PCReset(), but it's not working at the moment so we'll destroy instead
      if (_pc)
        _pc.destroy();

      LibmeshPetscCall(PCCreate(this->comm().get(), _pc.get()));

      auto pmatrix = cast_ptr<PetscMatrixBase<T> *>(this->_matrix);
      _mat = pmatrix->mat();
    }

  LibmeshPetscCall(PCSetOperators(_pc, _mat, _mat));

  // Set the PCType.  Note: this used to be done *before* the call to
  // PCSetOperators(), and only when !_is_initialized, but
  // 1.) Some preconditioners (those employing sub-preconditioners,
  // for example) have to call PCSetUp(), and can only do this after
  // the operators have been set.
  // 2.) It should be safe to call set_petsc_preconditioner_type()
  // multiple times.
  set_petsc_preconditioner_type(this->_preconditioner_type, *_pc);

  this->_is_initialized = true;
}



template <typename T>
void PetscPreconditioner<T>::clear()
{
  // Calls custom deleter
  _pc.destroy();
}



template <typename T>
PC PetscPreconditioner<T>::pc()
{
  return _pc;
}



template <typename T>
void PetscPreconditioner<T>::set_petsc_preconditioner_type (const PreconditionerType & preconditioner_type, PC & pc)
{
  // get the communicator from the PETSc object
  Parallel::communicator comm;
  PetscErrorCode ierr = PetscObjectGetComm((PetscObject)pc, & comm);
  if (ierr != LIBMESH_PETSC_SUCCESS)
    libmesh_error_msg("Error retrieving communicator");

  #define CasePCSetType(PreconditionerType, PCType)                       \
  case PreconditionerType:                                                \
    LibmeshPetscCallA(comm, PCSetType (pc, const_cast<KSPType>(PCType))); \
    break;

  switch (preconditioner_type)
    {
    CasePCSetType(IDENTITY_PRECOND,     PCNONE)
    CasePCSetType(CHOLESKY_PRECOND,     PCCHOLESKY)
    CasePCSetType(ICC_PRECOND,          PCICC)
    CasePCSetType(ILU_PRECOND,          PCILU)
    CasePCSetType(LU_PRECOND,           PCLU)
    CasePCSetType(ASM_PRECOND,          PCASM)
    CasePCSetType(JACOBI_PRECOND,       PCJACOBI)
    CasePCSetType(BLOCK_JACOBI_PRECOND, PCBJACOBI)
    CasePCSetType(SOR_PRECOND,          PCSOR)
    CasePCSetType(EISENSTAT_PRECOND,    PCEISENSTAT)
    CasePCSetType(AMG_PRECOND,          PCHYPRE)
    CasePCSetType(SVD_PRECOND,          PCSVD)
    CasePCSetType(USER_PRECOND,         PCMAT)
    CasePCSetType(SHELL_PRECOND,        PCSHELL)

    default:
      libMesh::err << "ERROR:  Unsupported PETSC Preconditioner: "
                   << Utility::enum_to_string(preconditioner_type) << std::endl
                   << "Continuing with PETSC defaults" << std::endl;
    }

  // Set additional options if we are doing AMG and
  // HYPRE is available
#ifdef LIBMESH_HAVE_PETSC_HYPRE
  if (preconditioner_type == AMG_PRECOND)
    LibmeshPetscCallA(comm, PCHYPRESetType(pc, "boomeramg"));
#endif

  // Let the commandline override stuff
  LibmeshPetscCallA(comm, PCSetFromOptions(pc));
}



#ifdef LIBMESH_HAVE_PETSC_HYPRE
template <typename T>
void PetscPreconditioner<T>::set_petsc_aux_data(PC & pc, System & sys, const unsigned v)
{
  // Get the communicator from the PETSc object
  Parallel::communicator comm;
  PetscErrorCode ierr = PetscObjectGetComm((PetscObject)pc, &comm);
  libmesh_error_msg_if(ierr != LIBMESH_PETSC_SUCCESS,
                       "Error retrieving communicator");

  // Make sure the preconditioner options are set
  LibmeshPetscCallA(comm, PCSetFromOptions(pc));

  // Get the type of preconditioner we are using
  PCType pc_type = nullptr;
  LibmeshPetscCallA(comm, PCGetType(pc, &pc_type));

  // Check if hypre ams/ads, otherwise we quit with nothing to do
  if (pc_type && std::string(pc_type) == PCHYPRE)
  {
    // Get the hypre preconditioner we are using
    PCType hypre_type = nullptr;
    LibmeshPetscCallA(comm, PCHYPREGetType(pc, &hypre_type));

    // If not ams/ads, we quit with nothing to do
    if (std::string(hypre_type) == "ams")
    {
      // If multiple variables, we error out as senseless
      libmesh_error_msg_if(sys.n_vars() > 1,
                           "Error applying hypre AMS to a system with multiple "
                           "variables");
      // If not a 1st order Nédélec or a 2d 1st order Raviart-Thomas system, we
      // error out as we do not support anything else at the moment
      libmesh_error_msg_if(sys.variable(v).type() != FEType(1, NEDELEC_ONE) &&
                          (sys.variable(v).type() != FEType(1, RAVIART_THOMAS) ||
                           sys.get_mesh().mesh_dimension() != 2),
                           "Error applying hypre AMS to a system "
                           "whose variable is not 1st order Nedelec or 1st "
                           "order Raviart-Thomas on a 2d mesh");
      set_hypre_ams_data(pc, sys, v);
    }
    else if (std::string(hypre_type) == "ads")
    {
      // If multiple variables, we error out as senseless
      libmesh_error_msg_if(sys.n_vars() > 1,
                           "Error applying hypre ADS to a system with multiple "
                           "variables");
      // If not a 3d 1st order Raviart-Thomas system, we error out as we do not
      // support anything else at the moment
      libmesh_error_msg_if(sys.variable(v).type() != FEType(1, RAVIART_THOMAS) ||
                           sys.get_mesh().mesh_dimension() != 3,
                           "Error applying hypre ADS to a system "
                           "whose variable is not 1st "
                           "order Raviart-Thomas on a 3d mesh");
      set_hypre_ads_data(pc, sys, v);
    }
  }
}



template <typename T>
void PetscPreconditioner<T>::set_hypre_ams_data(PC & pc, System & sys, const unsigned v)
{
  // Get the communicator from the PETSc object
  Parallel::communicator comm;
  PetscErrorCode ierr = PetscObjectGetComm((PetscObject)pc, &comm);
  libmesh_error_msg_if(ierr != LIBMESH_PETSC_SUCCESS,
                       "Error retrieving communicator");
  Parallel::Communicator Comm(comm);

  // Dummy Lagrange system defined over the same mesh so we can enumerate the vertices
  System & lagrange_sys = sys.get_equation_systems().add_system<System>("__hypre_ams_vertices");
  lagrange_sys.hide_output() = true;
  lagrange_sys.add_variable("__lagrange");
  lagrange_sys.reinit_mesh();

  // Global (i.e. total) and local (i.e. to this processor) number of edges and vertices
  const std::vector<dof_id_type> n_all_edges = sys.get_dof_map().n_dofs_per_processor(v);
  const dof_id_type n_glb_edges = std::reduce(n_all_edges.begin(), n_all_edges.end());
  const dof_id_type n_loc_edges = n_all_edges[global_processor_id()];
  const dof_id_type n_glb_verts = lagrange_sys.n_dofs();
  const dof_id_type n_loc_verts = lagrange_sys.n_local_dofs();

  // We require contiguous indexing through the edges: if the system has multiple
  // variables/splits, we therefore effectively need to find local dof numbers

  // The number of edges in all preceding processors
  dof_id_signed_type edge_offset = std::reduce(n_all_edges.begin(),
                                               n_all_edges.begin() + global_processor_id());

  // Whether dofs are in variable-major or node-major order
  bool var_major = !libMesh::on_command_line("--node-major-dofs");

  // If variable-major order, we need only subtract all the preceding dofs; but
  // if node-major order, we need to enumerate the dofs on this processor
  std::vector<dof_id_type> idx_edges;
  if (var_major)
  {
    edge_offset -= sys.get_dof_map().first_dof();
    for (auto i : make_range(v))
      edge_offset -= sys.get_dof_map().n_local_dofs(i);
  }
  else
    sys.get_dof_map().local_variable_indices(idx_edges, sys.get_mesh(), v);

  // Create the discrete grandient matrix, representing the edges in terms of its vertices
  // Preallocate 2 diagonal + 2 off-diagonal nonzeros as the vertices could fall on either
  PetscMatrix<Real> G(Comm, n_glb_edges, n_glb_verts, n_loc_edges, n_loc_verts, 2, 2);

  // Create vectors for the coordinates of the vertices
  PetscVector<Real> x(Comm, n_glb_verts, n_loc_verts);
  PetscVector<Real> y(Comm, n_glb_verts, n_loc_verts);
  PetscVector<Real> z(Comm, n_glb_verts, n_loc_verts);

  // Create vectors for the mat-vec products, representing const vector fields in the Nédélec basis
  PetscVector<Real> Gx(Comm, n_glb_edges, n_loc_edges);
  PetscVector<Real> Gy(Comm, n_glb_edges, n_loc_edges);
  PetscVector<Real> Gz(Comm, n_glb_edges, n_loc_edges);

  // Populate the discrete gradient matrix and the coordinate vectors
  for (const auto & elem : sys.get_mesh().active_local_element_ptr_range())
    for (auto edge : make_range(elem->n_edges()))
    {
      // The edge's first vertex: if owned, populate coordinate vectors
      const Node & vert_node = elem->node_ref(elem->local_edge_node(edge, 0));
      const dof_id_type vert_dof = vert_node.dof_number(lagrange_sys.number(), 0, 0);

      if (vert_node.processor_id() == global_processor_id())
      {
        x.set(vert_dof, vert_node(0));
        y.set(vert_dof, vert_node(1));
        z.set(vert_dof, vert_node(2));
      }

      // The edge's second vertex: if owned, populate coordinate vectors
      const Node & wert_node = elem->node_ref(elem->local_edge_node(edge, 1));
      const dof_id_type wert_dof = wert_node.dof_number(lagrange_sys.number(), 0, 0);

      if (wert_node.processor_id() == global_processor_id())
      {
        x.set(wert_dof, wert_node(0));
        y.set(wert_dof, wert_node(1));
        z.set(wert_dof, wert_node(2));
      }

      // The edge's (middle) node: if owned, populate discrete gradient matrix
      const Node & edge_node = elem->node_ref(elem->local_edge_node(edge, 2));
      const dof_id_type edge_dof = edge_node.dof_number(sys.number(), v, 0);

      if (edge_node.processor_id() == global_processor_id())
      {
        const dof_id_type seq_edge_dof = edge_offset +
          (var_major ? edge_dof
                     : std::distance(idx_edges.begin(),
                       std::find(idx_edges.begin(), idx_edges.end(), edge_dof)));

        const Real sign = elem->positive_edge_orientation(edge) ? 1 : -1;

        G.set(seq_edge_dof, vert_dof,  sign);
        G.set(seq_edge_dof, wert_dof, -sign);
      }
    }

  // Assemble the discrete gradient matrix and the coordinate vectors
  G.close();
  x.close();
  y.close();
  z.close();

  // Compute the matrix-vector products
  G.vector_mult(Gx, x);
  G.vector_mult(Gy, y);
  G.vector_mult(Gz, z);

  // Hand over the matrix and vectors
  LibmeshPetscCallA(comm, PCHYPRESetDiscreteGradient(pc, G.mat()));
  LibmeshPetscCallA(comm, PCHYPRESetEdgeConstantVectors(pc, Gx.vec(), Gy.vec(), Gz.vec()));
}



template <typename T>
void PetscPreconditioner<T>::set_hypre_ads_data(PC & pc, System & sys, const unsigned v)
{
}
#endif


//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT PetscPreconditioner<Number>;

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_PETSC
