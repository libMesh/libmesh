// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_SENSITIVITY_DATA_H
#define LIBMESH_SENSITIVITY_DATA_H


// Local Includes -----------------------------------
#include "libmesh/libmesh_common.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/system.h"

// C++ Includes   -----------------------------------
#include <vector>

namespace libMesh
{

// Forward declaractions
class QoISet;

/**
 * Data structure for holding completed parameter sensitivity
 * calculations.
 */
class SensitivityData
{
public:
  class Row
  {
  public:
    Row(SensitivityData & sd, unsigned int qoi) : _sd(sd), _qoi(qoi) {}

    Number & operator[] (unsigned int parameter) { return _sd.derivative(_qoi, parameter); }
  private:
    SensitivityData & _sd;
    unsigned int _qoi;
  };

  class ConstRow
  {
  public:
    ConstRow(const SensitivityData & sd, unsigned int qoi) : _sd(sd), _qoi(qoi) {}

    const Number & operator[] (unsigned int parameter) { return _sd.derivative(_qoi, parameter); }
  private:
    const SensitivityData & _sd;
    unsigned int _qoi;
  };

  /**
   * Default constructor: empty data set
   */
  SensitivityData() {}

  /**
   * Constructor from QoISet and ParameterVector: allocates space
   * for all required sensitivities
   */
  SensitivityData(const QoISet & qoi_indices,
                  const System & sys,
                  const ParameterVector & parameter_vector);

  /**
   * Clears and deallocates all data
   */
  void clear() { _grad_data.clear(); }

  /**
   * Given QoISet and ParameterVector, allocates space
   * for all required first derivative data
   */
  void allocate_data(const QoISet & qoi_indices,
                     const System & sys,
                     const ParameterVector & parameter_vector);

  /**
   * Given QoISet and ParameterVector, allocates space
   * for all required second derivative data
   */
  void allocate_hessian_data(const QoISet & qoi_indices,
                             const System & sys,
                             const ParameterVector & parameter_vector);

  /**
   * Returns the parameter sensitivity derivative for the specified
   * quantity of interest for the specified parameter
   */
  const Number & derivative (unsigned int qoi_index,
                             unsigned int parameter_index) const;

  /**
   * Returns the parameter sensitivity derivative for the specified
   * quantity of interest for the specified pair of parameters
   */
  const Number & second_derivative (unsigned int qoi_index,
                                    unsigned int parameter_index1,
                                    unsigned int parameter_index2) const;

  /**
   * Gets/sets the parameter sensitivity derivative for the specified
   * quantity of interest for the specified parameter
   */
  Number & derivative (unsigned int qoi_index,
                       unsigned int parameter_index);

  /**
   * Gets/sets the parameter sensitivity second derivative for the
   * specified quantity of interest for the specified pair of
   * parameters
   */
  Number & second_derivative (unsigned int qoi_index,
                              unsigned int parameter_index1,
                              unsigned int parameter_index2);

  /**
   * Vector address type operator: sd[q][p] is an alias for
   * sd.derivative(q,p)
   */
  ConstRow operator[] (unsigned int qoi) const { return ConstRow(*this, qoi); }

  Row operator[] (unsigned int qoi) { return Row(*this, qoi); }

private:
  /**
   * Data storage; currently pretty trivial
   */
  std::vector<std::vector<Number> > _grad_data;
  std::vector<std::vector<std::vector<Number> > > _hess_data;
};



// ------------------------------------------------------------
// SensitivityData inline methods



inline
SensitivityData::SensitivityData(const QoISet & qoi_indices,
                                 const System & sys,
                                 const ParameterVector & parameter_vector)
{
  this->allocate_data(qoi_indices, sys, parameter_vector);
}



inline
void SensitivityData::allocate_data(const QoISet & qoi_indices,
                                    const System & sys,
                                    const ParameterVector & parameter_vector)
{
  const std::size_t Np = parameter_vector.size();
  const unsigned int Nq =
    cast_int<unsigned int>(sys.qoi.size());

  if (_grad_data.size() < Nq)
    _grad_data.resize(Nq);

  for (unsigned int i=0; i != Nq; ++i)
    if (qoi_indices.has_index(i))
      {
        _grad_data[i].clear();
        _grad_data[i].resize(Np);
      }
}



inline
void SensitivityData::allocate_hessian_data(const QoISet & qoi_indices,
                                            const System & sys,
                                            const ParameterVector & parameter_vector)
{
  const std::size_t Np = parameter_vector.size();
  const unsigned int Nq =
    cast_int<unsigned int>(sys.qoi.size());

  if (_hess_data.size() < Nq)
    _hess_data.resize(Nq);

  for (unsigned int i=0; i != Nq; ++i)
    if (qoi_indices.has_index(i))
      {
        _hess_data[i].clear();
        _hess_data[i].resize(Np);
        for (std::size_t j=0; j != Np; ++j)
          _hess_data[i][j].resize(Np);
      }
}



inline
const Number & SensitivityData::derivative(unsigned int qoi_index,
                                           unsigned int parameter_index) const
{
  libmesh_assert_less (qoi_index, _grad_data.size());
  libmesh_assert_less (parameter_index, _grad_data[qoi_index].size());

  return _grad_data[qoi_index][parameter_index];
}



inline
Number & SensitivityData::derivative(unsigned int qoi_index,
                                     unsigned int parameter_index)
{
  libmesh_assert_less (qoi_index, _grad_data.size());
  libmesh_assert_less (parameter_index, _grad_data[qoi_index].size());

  return _grad_data[qoi_index][parameter_index];
}



inline
const Number & SensitivityData::second_derivative(unsigned int qoi_index,
                                                  unsigned int parameter_index1,
                                                  unsigned int parameter_index2) const
{
  libmesh_assert_less (qoi_index, _hess_data.size());
  libmesh_assert_less (parameter_index1, _hess_data[qoi_index].size());
  libmesh_assert_less (parameter_index2, _hess_data[qoi_index][parameter_index1].size());

  return _hess_data[qoi_index][parameter_index1][parameter_index2];
}



inline
Number & SensitivityData::second_derivative(unsigned int qoi_index,
                                            unsigned int parameter_index1,
                                            unsigned int parameter_index2)
{
  libmesh_assert_less (qoi_index, _hess_data.size());
  libmesh_assert_less (parameter_index1, _hess_data[qoi_index].size());
  libmesh_assert_less (parameter_index2, _hess_data[qoi_index][parameter_index1].size());

  return _hess_data[qoi_index][parameter_index1][parameter_index2];
}

} // namespace libMesh

#endif // LIBMESH_SENSITIVITY_DATA_H
