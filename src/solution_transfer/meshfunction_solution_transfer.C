#include "libmesh/meshfunction_solution_transfer.h"

#include "libmesh/system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/mesh_function.h"

namespace libMesh {

MeshFunctionSolutionTransfer::MeshFunctionSolutionTransfer()
{}

MeshFunctionSolutionTransfer::~MeshFunctionSolutionTransfer()
{}

void
MeshFunctionSolutionTransfer::transfer(const Variable & from_var, const Variable & to_var)
{
  // This only works when transferring to a Lagrange variable
  libmesh_assert(to_var.type().family == LAGRANGE);

  unsigned int to_var_num = to_var.number();

  System * from_sys = from_var.system();
  System * to_sys = to_var.system();

  // Only works with a serialized mesh to transfer from!
  libmesh_assert(from_sys->get_mesh().is_serial());

  unsigned int to_sys_num = to_sys->number();

  EquationSystems & from_es = from_sys->get_equation_systems();

  //Create a serialized version of the solution vector
  NumericVector<Number> * serialized_solution = NumericVector<Number>::build().release();
  serialized_solution->init(from_sys->n_dofs(), false, SERIAL);

  // Need to pull down a full copy of this vector on every processor so we can get values in parallel
  from_sys->solution->localize(*serialized_solution);

  MeshFunction from_func(from_es, *serialized_solution, from_sys->get_dof_map(), to_var_num);
  from_func.init();
  
  MeshBase::const_node_iterator nd     = to_sys->get_mesh().local_nodes_begin();
  MeshBase::const_node_iterator nd_end = to_sys->get_mesh().local_nodes_end();
  
  // Now loop over the nodes of the 'To' mesh setting values for each variable.
  for(;nd != nd_end; ++nd)
    // 0 is for the value component
    to_sys->solution->set((*nd)->dof_number(to_sys_num, to_var_num, 0), from_func(**nd));

  to_sys->solution->close();
  to_sys->update();

  delete serialized_solution;
}

} // namespace libMesh
