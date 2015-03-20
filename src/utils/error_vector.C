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
#include <limits>

// Local includes
#include "libmesh/elem.h"
#include "libmesh/error_vector.h"
#include "libmesh/libmesh_logging.h"

#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/gmv_io.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/exodusII_io.h"

namespace libMesh
{



// ------------------------------------------------------------
// ErrorVector class member functions
ErrorVectorReal ErrorVector::minimum() const
{
  START_LOG ("minimum()", "ErrorVector");

  const dof_id_type n = cast_int<dof_id_type>(this->size());
  ErrorVectorReal min = std::numeric_limits<ErrorVectorReal>::max();

  for (dof_id_type i=0; i<n; i++)
    {
      // Only positive (or zero) values in the error vector
      libmesh_assert_greater_equal ((*this)[i], 0.);
      if (this->is_active_elem(i))
        min = std::min (min, (*this)[i]);
    }
  STOP_LOG ("minimum()", "ErrorVector");

  // ErrorVectors are for positive values
  libmesh_assert_greater_equal (min, 0.);

  return min;
}



Real ErrorVector::mean() const
{
  START_LOG ("mean()", "ErrorVector");

  const dof_id_type n = cast_int<dof_id_type>(this->size());

  Real the_mean  = 0;
  dof_id_type nnz = 0;

  for (dof_id_type i=0; i<n; i++)
    if (this->is_active_elem(i))
      {
        the_mean += ( static_cast<Real>((*this)[i]) - the_mean ) / (nnz + 1);

        nnz++;
      }

  STOP_LOG ("mean()", "ErrorVector");

  return the_mean;
}




Real ErrorVector::median()
{
  const dof_id_type n = cast_int<dof_id_type>(this->size());

  if (n == 0)
    return 0.;


  // Build a StatisticsVector<ErrorVectorReal> containing
  // only our active entries and take its mean
  StatisticsVector<ErrorVectorReal> sv;

  sv.reserve (n);

  for (dof_id_type i=0; i<n; i++)
    if(this->is_active_elem(i))
      sv.push_back((*this)[i]);

  return sv.median();
}




Real ErrorVector::median() const
{
  ErrorVector ev = (*this);

  return ev.median();
}




Real ErrorVector::variance(const Real mean_in) const
{
  const dof_id_type n = cast_int<dof_id_type>(this->size());

  START_LOG ("variance()", "ErrorVector");

  Real the_variance = 0;
  dof_id_type nnz = 0;

  for (dof_id_type i=0; i<n; i++)
    if (this->is_active_elem(i))
      {
        const Real delta = ( static_cast<Real>((*this)[i]) - mean_in );
        the_variance += (delta * delta - the_variance) / (nnz + 1);

        nnz++;
      }

  STOP_LOG ("variance()", "ErrorVector");

  return the_variance;
}




std::vector<dof_id_type> ErrorVector::cut_below(Real cut) const
{
  START_LOG ("cut_below()", "ErrorVector");

  const dof_id_type n = cast_int<dof_id_type>(this->size());

  std::vector<dof_id_type> cut_indices;
  cut_indices.reserve(n/2);  // Arbitrary

  for (dof_id_type i=0; i<n; i++)
    if (this->is_active_elem(i))
      {
        if ((*this)[i] < cut)
          {
            cut_indices.push_back(i);
          }
      }

  STOP_LOG ("cut_below()", "ErrorVector");

  return cut_indices;
}




std::vector<dof_id_type> ErrorVector::cut_above(Real cut) const
{
  START_LOG ("cut_above()", "ErrorVector");

  const dof_id_type n = cast_int<dof_id_type>(this->size());

  std::vector<dof_id_type> cut_indices;
  cut_indices.reserve(n/2);  // Arbitrary

  for (dof_id_type i=0; i<n; i++)
    if (this->is_active_elem(i))
      {
        if ((*this)[i] > cut)
          {
            cut_indices.push_back(i);
          }
      }

  STOP_LOG ("cut_above()", "ErrorVector");

  return cut_indices;
}



bool ErrorVector::is_active_elem (dof_id_type i) const
{
  libmesh_assert_less (i, this->size());

  if (_mesh)
    {
      libmesh_assert(_mesh->elem(i));
      return _mesh->elem(i)->active();
    }
  else
    return ((*this)[i] != 0.);
}


void ErrorVector::plot_error(const std::string& filename,
                             const MeshBase& oldmesh) const
{
  UniquePtr<MeshBase> meshptr = oldmesh.clone();
  MeshBase &mesh = *meshptr;

  // The all_first_order routine requires that renumbering be allowed
  mesh.allow_renumbering(true);

  mesh.all_first_order();
  EquationSystems temp_es (mesh);
  ExplicitSystem& error_system
    = temp_es.add_system<ExplicitSystem> ("Error");
  error_system.add_variable("error", CONSTANT, MONOMIAL);
  temp_es.init();

  const DofMap& error_dof_map = error_system.get_dof_map();

  MeshBase::const_element_iterator       el     =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();
  std::vector<dof_id_type> dof_indices;

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      error_dof_map.dof_indices(elem, dof_indices);

      const dof_id_type elem_id = elem->id();

      //0 for the monomial basis
      const dof_id_type solution_index = dof_indices[0];

      // libMesh::out << "elem_number=" << elem_number << std::endl;
      libmesh_assert_less (elem_id, (*this).size());

      // We may have zero error values in special circumstances
      // libmesh_assert_greater ((*this)[elem_id], 0.);
      error_system.solution->set(solution_index, (*this)[elem_id]);
    }

  // We may have to renumber if the original numbering was not
  // contiguous.  Since this is just a temporary mesh, that's probably
  // fine.
  if (mesh.max_elem_id() != mesh.n_elem() ||
      mesh.max_node_id() != mesh.n_nodes())
    {
      mesh.allow_renumbering(true);
      mesh.renumber_nodes_and_elements();
    }

  if (filename.rfind(".gmv") < filename.size())
    {
      GMVIO(mesh).write_discontinuous_gmv(filename,
                                          temp_es, false);
    }
  else if (filename.rfind(".plt") < filename.size())
    {
      TecplotIO (mesh).write_equation_systems
        (filename, temp_es);
    }
#ifdef LIBMESH_HAVE_EXODUS_API
  else if( (filename.rfind(".exo") < filename.size()) ||
           (filename.rfind(".e") < filename.size()) )
    {
      ExodusII_IO io(mesh);
      io.write(filename);
      io.write_element_data(temp_es);
    }
#endif
  else
    {
      libmesh_here();
      libMesh::err << "Warning: ErrorVector::plot_error currently only"
                   << " supports .gmv and .plt and .exo/.e (if enabled) output;" << std::endl;
      libMesh::err << "Could not recognize filename: " << filename
                   << std::endl;
    }
}

} // namespace libMesh
