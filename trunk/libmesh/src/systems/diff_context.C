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


#include "libmesh/diff_context.h"
#include "libmesh/diff_system.h"

namespace libMesh
{



DiffContext::DiffContext (const System& sys) :
  time(sys.time),
  system_time(sys.time),
  elem_solution_derivative(1.),
  fixed_solution_derivative(0.),
  dof_indices_var(sys.n_vars()),
  _deltat(NULL),
  _system(sys),
  _is_adjoint(false)
{
  // Finally initialize solution/residual/jacobian data structures
  unsigned int n_vars = sys.n_vars();

  elem_subsolutions.reserve(n_vars);
  elem_subresiduals.reserve(n_vars);
  elem_subjacobians.resize(n_vars);
  if (sys.use_fixed_solution)
    elem_fixed_subsolutions.reserve(n_vars);

  // If the user resizes sys.qoi, it will invalidate us
  unsigned int n_qoi = sys.qoi.size();
  elem_qoi.resize(n_qoi);
  elem_qoi_derivative.resize(n_qoi);
  elem_qoi_subderivatives.resize(n_qoi);
  for (unsigned int q=0; q != n_qoi; ++q)
    elem_qoi_subderivatives[q].reserve(n_vars);

  for (unsigned int i=0; i != n_vars; ++i)
    {
      elem_subsolutions.push_back(new DenseSubVector<Number>(elem_solution));
      elem_subresiduals.push_back(new DenseSubVector<Number>(elem_residual));
      for (unsigned int q=0; q != n_qoi; ++q)
        elem_qoi_subderivatives[q].push_back(new DenseSubVector<Number>(elem_qoi_derivative[q]));
      elem_subjacobians[i].reserve(n_vars);

      if (sys.use_fixed_solution)
        elem_fixed_subsolutions.push_back
	  (new DenseSubVector<Number>(elem_fixed_solution));

      for (unsigned int j=0; j != n_vars; ++j)
        {
          elem_subjacobians[i].push_back
            (new DenseSubMatrix<Number>(elem_jacobian));
        }
    }
}



DiffContext::~DiffContext ()
{
  for (unsigned int i=0; i != elem_subsolutions.size(); ++i)
    {
      delete elem_subsolutions[i];
      delete elem_subresiduals[i];
      for (unsigned int q=0; q != elem_qoi_subderivatives.size(); ++q)
        delete elem_qoi_subderivatives[q][i];
      if (!elem_fixed_subsolutions.empty())
        delete elem_fixed_subsolutions[i];

      for (unsigned int j=0; j != elem_subjacobians[i].size(); ++j)
        delete elem_subjacobians[i][j];
    }

  // We also need to delete all the DenseSubVectors from the localized_vectors map
  // localized_vectors iterators
  std::map<const NumericVector<Number>*, std::pair<DenseVector<Number>, std::vector<DenseSubVector<Number>*> > >::iterator localized_vectors_it = localized_vectors.begin();
  std::map<const NumericVector<Number>*, std::pair<DenseVector<Number>, std::vector<DenseSubVector<Number>*> > >::iterator localized_vectors_end = localized_vectors.end();

  // Loop over every localized_vector
  for(; localized_vectors_it != localized_vectors_end; ++localized_vectors_it)
    {
      // Grab the DenseSubVector to be deleted
      std::vector<DenseSubVector<Number>* >&  localized_vector_dsv = localized_vectors_it->second.second;

      // Loop over that vector and delete each entry
      for(unsigned int i=0; i != localized_vector_dsv.size(); ++i)
	delete localized_vector_dsv[i];
	      
    }
      
}


  void DiffContext::set_deltat_pointer(Real* dt)
  {
    // We may actually want to be able to set this pointer to NULL, so
    // don't report an error for that.
    _deltat = dt;
  }

  Real DiffContext::get_deltat_value()
  {
    libmesh_assert(_deltat);

    return *_deltat;
  }

  void DiffContext::add_localized_vector (NumericVector<Number> & _localized_vector, const System & _sys)
  {
    // Make an empty pair keyed with a reference to this _localized_vector
    localized_vectors[&_localized_vector] = std::make_pair(DenseVector<Number>(), std::vector<DenseSubVector<Number>*>());

    unsigned int n_vars = _sys.n_vars();

    localized_vectors[&_localized_vector].second.reserve(n_vars);
    
    // Fill the DenseSubVector with n_vars copies of DenseVector
    for(unsigned int i=0; i != n_vars; ++i)
      localized_vectors[&_localized_vector].second.push_back(new DenseSubVector<Number>(localized_vectors[&_localized_vector].first));

  }

DenseVector<Number>& DiffContext::get_localized_vector (const NumericVector<Number> & _localized_vector)
{
  return localized_vectors[&_localized_vector].first;
}

const DenseVector<Number>& DiffContext::get_localized_vector (const NumericVector<Number> & _localized_vector) const
{
  std::map<const NumericVector<Number>*, std::pair<DenseVector<Number>, std::vector<DenseSubVector<Number>*> > >::const_iterator
    localized_vectors_it = localized_vectors.find(&_localized_vector);
  libmesh_assert(localized_vectors_it != localized_vectors.end());
  return localized_vectors_it->second.first;
}

DenseSubVector<Number>& DiffContext::get_localized_subvector (const NumericVector<Number> & _localized_vector, unsigned int _var)
{
  return *localized_vectors[&_localized_vector].second[_var];
}

const DenseSubVector<Number>& DiffContext::get_localized_subvector (const NumericVector<Number> & _localized_vector, unsigned int _var) const
{
  std::map<const NumericVector<Number>*, std::pair<DenseVector<Number>, std::vector<DenseSubVector<Number>*> > >::const_iterator
    localized_vectors_it = localized_vectors.find(&_localized_vector);
  libmesh_assert(localized_vectors_it != localized_vectors.end());
  return *localized_vectors_it->second.second[_var];
}

} // namespace libMesh



