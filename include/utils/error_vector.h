// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_ERROR_VECTOR_H
#define LIBMESH_ERROR_VECTOR_H

// Local Includes
#include "libmesh/statistics.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// Forward Declarations
template <typename> class MeshBaseTempl;
typedef MeshBaseTempl<Real> MeshBase;
class Mesh;

/**
 * The \p ErrorVector is a specialization of the
 * \p StatisticsVector for error data computed on a finite element
 * mesh.  In general, when computing the error on a mesh only the
 * active elements are considered, but the \p ErrorVector is sized
 * according to the total number of elements in the mesh.  The
 * \p ErrorVector is thus padded with zeros for all the inactive
 * elements, and this must be taken into account when calculating
 * the statistics.  Since the error is a positive quantity this class
 * assumes it contains positive data (i.e. min_val >= 0.).
 *
 * \author Benjamin S. Kirk
 * \date 2003
 */
class ErrorVector : public StatisticsVector<ErrorVectorReal>
{
public:

  /**
   * ErrorVector constructor; sets initial length to \p i.
   *
   * If mesh is not null, MeshBase::elem() and Elem::is_active() will
   * be used to distinguish active and inactive elements.  If mesh is null,
   * ErrorVector will assume that all 0.0 error values correspond to inactive
   * elements and all non-zero error values correspond to active elements.
   */
  ErrorVector(dof_id_type i=0, MeshBase * mesh = nullptr) :
    StatisticsVector<ErrorVectorReal> (i),
    _mesh(mesh)
  {}

  /**
   * ErrorVector constructor; sets initial length to \p i and initial values to \p val.
   *
   * If mesh is not null, MeshBase::elem() and Elem::is_active() will
   * be used to distinguish active and inactive elements.  If mesh is null,
   * ErrorVector will assume that all 0.0 error values correspond to inactive
   * elements and all non-zero error values correspond to active elements.
   */
  ErrorVector(dof_id_type i, ErrorVectorReal val) :
    StatisticsVector<ErrorVectorReal> (i,val) {}

  /**
   * \returns The minimum nonzero value in the data set.
   */
  virtual ErrorVectorReal minimum() const override;

  /**
   * \returns The mean value of the data set. Ignores
   * zero values.
   */
  virtual Real mean() const override;

  /**
   * \returns The median (e.g. the middle) value of the data set,
   * ignoring inactive elements.
   *
   * This function modifies the original data by sorting, so it can't
   * be called on const objects.  Source: GNU Scientific Library
   */
  virtual Real median() override;

  /**
   * A const version of the median function.
   * Requires twice the memory of original
   * data set but does not change the original.
   */
  virtual Real median() const override;

  /**
   * \returns The variance of the data set ignoring inactive elements.
   *
   * Uses a recurrence relation to prevent data overflow for large
   * sums.
   *
   * \note The variance is equal to the standard deviation squared.
   * The variance is normalized by N in this case.  Source: GNU
   * Scientific Library.
   */
  virtual Real variance() const override
  { return this->variance(this->mean()); }

  /**
   * \returns The variance of the data set ignoring inactive elements
   * and given the \p mean.
   *
   * This is useful for efficiency when you have already calculated
   * the mean. Uses a recurrence relation to prevent data overflow for
   * large sums.
   *
   * \note The variance is equal to the standard deviation squared.
   * Source: GNU Scientific Library.
   */
  virtual Real variance(const Real mean) const override;

  /**
   * \returns A vector of dof_id_types which correspond
   * to the indices of every member of the data set
   * below the cutoff value cut ignoring inactive elements.
   */
  virtual std::vector<dof_id_type> cut_below(Real cut) const override;

  /**
   * \returns A vector of dof_id_types which correspond
   * to the indices of every member of the data set
   * above the cutoff value cut ignoring inactive elements.
   */
  virtual std::vector<dof_id_type> cut_above(Real cut) const override;

  /**
   * Plots a data file, of a type determined by looking at
   * the file extension in \p filename, of the error values on
   * the active elements of \p mesh.
   */
  void plot_error(const std::string & filename,
                  const MeshBase & mesh) const;

protected:
  /**
   * Utility function to decide whether element i is active
   */
  bool is_active_elem (dof_id_type i) const;

  /**
   * Pointer to the mesh, which may be used to decide which
   * elements are active
   */
  MeshBase * _mesh;
};

} // namespace libMesh

#endif // LIBMESH_ERROR_VECTOR_H
