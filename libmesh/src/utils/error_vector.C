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


// C++ includes
#include <limits>

// Local includes
#include "elem.h"
#include "error_vector.h"
#include "libmesh_logging.h"

#include "dof_map.h"
#include "equation_systems.h"
#include "explicit_system.h"
#include "mesh_base.h"
#include "numeric_vector.h"
#include "gmv_io.h"
#include "tecplot_io.h"

namespace libMesh
{



// ------------------------------------------------------------
// ErrorVector class member functions
ErrorVectorReal ErrorVector::minimum() const
{
  START_LOG ("minimum()", "ErrorVector");

  const unsigned int n = this->size();
  ErrorVectorReal min = std::numeric_limits<ErrorVectorReal>::max();

  for (unsigned int i=0; i<n; i++)
    {
      // Only positive (or zero) values in the error vector
      libmesh_assert((*this)[i] >= 0.);
      if (this->is_active_elem(i))
        min = std::min (min, (*this)[i]);
    }
  STOP_LOG ("minimum()", "ErrorVector");

  // ErrorVectors are for positive values
  libmesh_assert (min >= 0.);

  return min;
}



Real ErrorVector::mean() const
{
  START_LOG ("mean()", "ErrorVector");

  const unsigned int n = this->size();

  Real mean  = 0;
  unsigned int nnz = 0;

  for (unsigned int i=0; i<n; i++)
    if (this->is_active_elem(i))
      {
	mean += ( static_cast<Real>((*this)[i]) - mean ) / (nnz + 1);

	nnz++;
      }

  STOP_LOG ("mean()", "ErrorVector");

  return mean;
}




Real ErrorVector::median()
{
  const unsigned int n   = this->size();

  if (n == 0)
    return 0.;


  // Build a StatisticsVector<ErrorVectorReal> containing
  // only our active entries and take its mean
  StatisticsVector<ErrorVectorReal> sv;

  sv.reserve (n);

  for (unsigned int i=0; i<n; i++)
    if(this->is_active_elem(i))
      sv.push_back((*this)[i]);

  return sv.median();
}




Real ErrorVector::median() const
{
  ErrorVector ev = (*this);

  return ev.median();
}




Real ErrorVector::variance(const Real mean) const
{
  const unsigned int n   = this->size();

  START_LOG ("variance()", "ErrorVector");

  Real variance = 0;
  unsigned int nnz = 0;

  for (unsigned int i=0; i<n; i++)
    if (this->is_active_elem(i))
      {
	const Real delta = ( static_cast<Real>((*this)[i]) - mean );
	variance += (delta * delta - variance) / (nnz + 1);

	nnz++;
      }

  STOP_LOG ("variance()", "ErrorVector");

  return variance;
}




std::vector<unsigned int> ErrorVector::cut_below(Real cut) const
{
  START_LOG ("cut_below()", "ErrorVector");

  const unsigned int n = this->size();

  std::vector<unsigned int> cut_indices;
  cut_indices.reserve(n/2);  // Arbitrary

  for (unsigned int i=0; i<n; i++)
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




std::vector<unsigned int> ErrorVector::cut_above(Real cut) const
{
  START_LOG ("cut_above()", "ErrorVector");

  const unsigned int n   = this->size();

  std::vector<unsigned int> cut_indices;
  cut_indices.reserve(n/2);  // Arbitrary

  for (unsigned int i=0; i<n; i++)
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



bool ErrorVector::is_active_elem (unsigned int i) const
{
  libmesh_assert (i < this->size());

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
  AutoPtr<MeshBase> meshptr = oldmesh.clone();
  MeshBase &mesh = *meshptr;
  mesh.all_first_order();
  EquationSystems temp_es (mesh);
  ExplicitSystem& error_system
    = temp_es.add_system<ExplicitSystem> ("Error");
  error_system.add_variable("e", CONSTANT, MONOMIAL);
  temp_es.init();

  const DofMap& error_dof_map = error_system.get_dof_map();

  MeshBase::const_element_iterator       el     =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();
  std::vector<unsigned int> dof_indices;

  for ( ; el != end_el; ++el)
  {
    Elem* elem = *el;

    error_dof_map.dof_indices(elem, dof_indices);

    const unsigned int elem_id = elem->id();

    //0 for the monomial basis
    const unsigned int solution_index = dof_indices[0];
    // libMesh::out << "elem_number=" << elem_number << std::endl;
    libmesh_assert (elem_id < (*this).size());
//  We may have zero error values in special circumstances
//    libmesh_assert ((*this)[elem_id] > 0.);
    error_system.solution->set(solution_index, (*this)[elem_id]);
  }

  if (filename.rfind(".gmv") < filename.size())
    {
//      GMVIO(mesh).write_equation_systems(filename,
//                                         temp_es);
// write_discontinuous_gmv doesn't work on second order meshes??
// [RHS]
      GMVIO(mesh).write_discontinuous_gmv(filename,
                                          temp_es, false);
    }
  else if (filename.rfind(".plt") < filename.size())
    {
      TecplotIO (mesh).write_equation_systems
        (filename, temp_es);
    }
  else
    {
      libmesh_here();
      libMesh::err << "Warning: ErrorVector::plot_error currently only"
                    << " supports .gmv and .plt output;" << std::endl;
      libMesh::err << "Could not recognize filename: " << filename
                    << std::endl;
    }
}

} // namespace libMesh

