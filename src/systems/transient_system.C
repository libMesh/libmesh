// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/rb_construction.h"
#include "libmesh/eigen_system.h"

namespace libMesh
{


// ------------------------------------------------------------
// TransientSystem implementation
template <class Base>
TransientSystem<Base>::TransientSystem (EquationSystems & es,
                                        const std::string & name_in,
                                        const unsigned int number_in) :

  Base                 (es, name_in, number_in)
{
  this->add_old_vectors();
}



template <class Base>
TransientSystem<Base>::~TransientSystem ()
{
  this->clear();

  // We still have std::unique_ptrs for API compatibility, but
  // now that we're System::add_vector()ing these, we can trust
  // the base class to handle memory management
  old_local_solution.release();
  older_local_solution.release();
}



template <class Base>
void TransientSystem<Base>::clear ()
{
  // clear the parent data
  Base::clear();

  // the old & older local solutions
  // are now deleted by System!
  // old_local_solution->clear();
  // older_local_solution->clear();

  // FIXME: This preserves maximum backwards compatibility,
  // but is probably grossly unnecessary:
  old_local_solution.release();
  older_local_solution.release();

  // Restore us to a "basic" state
  this->add_old_vectors();
}



template <class Base>
void TransientSystem<Base>::re_update ()
{
  // re_update the parent system
  Base::re_update ();

  const std::vector<dof_id_type> & send_list = this->get_dof_map().get_send_list ();

  const dof_id_type first_local_dof = Base::get_dof_map().first_dof();
  const dof_id_type end_local_dof  = Base::get_dof_map().end_dof();

  // Check sizes
  libmesh_assert_greater_equal (end_local_dof, first_local_dof);
  libmesh_assert_greater_equal (older_local_solution->size(), send_list.size());
  libmesh_assert_greater_equal (old_local_solution->size(), send_list.size());

  // Even if we don't have to do anything ourselves, localize() may
  // use parallel_only tools
  // if (first_local_dof == end_local_dof)
  //   return;

  // Update the old & older solutions with the send_list,
  // which may have changed since their last update.
  older_local_solution->localize (first_local_dof,
                                  end_local_dof-1,
                                  send_list);

  old_local_solution->localize (first_local_dof,
                                end_local_dof-1,
                                send_list);
}




template <class Base>
Number TransientSystem<Base>::old_solution (const dof_id_type global_dof_number) const
{
  // Check the sizes
  libmesh_assert_less (global_dof_number, this->get_dof_map().n_dofs());
  libmesh_assert_less (global_dof_number, old_local_solution->size());

  return (*old_local_solution)(global_dof_number);
}



template <class Base>
Number TransientSystem<Base>::older_solution (const dof_id_type global_dof_number) const
{
  // Check the sizes
  libmesh_assert_less (global_dof_number, this->get_dof_map().n_dofs());
  libmesh_assert_less (global_dof_number, older_local_solution->size());

  return (*older_local_solution)(global_dof_number);
}



template <class Base>
void TransientSystem<Base>::add_old_vectors()
{
#ifdef LIBMESH_ENABLE_GHOSTED
  old_local_solution =
    std::unique_ptr<NumericVector<Number>>
    (&(this->add_vector("_transient_old_local_solution", true, GHOSTED)));
  older_local_solution =
    std::unique_ptr<NumericVector<Number>>
    (&(this->add_vector("_transient_older_local_solution", true, GHOSTED)));
#else
  old_local_solution =
    std::unique_ptr<NumericVector<Number>>
    (&(this->add_vector("_transient_old_local_solution", true, SERIAL)));
  older_local_solution =
    std::unique_ptr<NumericVector<Number>>
    (&(this->add_vector("_transient_older_local_solution", true, SERIAL)));
#endif
}

// ------------------------------------------------------------
// TransientSystem instantiations
template class TransientSystem<LinearImplicitSystem>;
template class TransientSystem<NonlinearImplicitSystem>;
template class TransientSystem<ExplicitSystem>;
template class TransientSystem<System>;
template class TransientSystem<RBConstruction>;
#ifdef LIBMESH_HAVE_SLEPC
template class TransientSystem<EigenSystem>;
#endif

} // namespace libMesh
