

// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __sensitivity_data_h__
#define __sensitivity_data_h__


// C++ Includes   -----------------------------------
#include <vector>

// Local Includes -----------------------------------
#include "libmesh_common.h"

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
    Row(SensitivityData &sd, unsigned int qoi) : _sd(sd), _qoi(qoi) {}

    Number& operator[] (unsigned int parameter) { return _sd.value(_qoi, parameter); }
  private:
    SensitivityData &_sd;
    unsigned int _qoi;
  };

  class ConstRow
  {
  public:
    ConstRow(const SensitivityData &sd, unsigned int qoi) : _sd(sd), _qoi(qoi) {}

    const Number& operator[] (unsigned int parameter) { return _sd.value(_qoi, parameter); }
  private:
    const SensitivityData &_sd;
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
  SensitivityData(const QoISet& qoi_indices,
                  const System& sys,
                  const ParameterVector& parameter_vector);

  /**
   * Clears and deallocates all data
   */
  void clear() { _data.clear(); }

  /**
   * Given QoISet and ParameterVector, allocates space
   * for all required sensitivities
   */
  void allocate_data(const QoISet &qoi_indices,
                     const System& sys,
                     const ParameterVector& parameter_vector);

  /**
   * Returns the parameter sensitivity derivative for the specified
   * quantity of interest for the specified parameter
   */
  const Number& value(unsigned int qoi_index,
                      unsigned int parameter_index) const;

  /**
   * Gets/sets the parameter sensitivity derivative for the specified
   * quantity of interest for the specified parameter
   */
  Number& value(unsigned int qoi_index,
                unsigned int parameter_index);

  ConstRow operator[] (unsigned int qoi) const { return ConstRow(*this, qoi); }

  Row operator[] (unsigned int qoi) { return Row(*this, qoi); }

private: 
  /**
   * Data storage; currently pretty trivial
   */
  std::vector<std::vector<Number> > _data;
};



// ------------------------------------------------------------
// SensitivityData inline methods



inline
SensitivityData::SensitivityData(const QoISet &qoi_indices,
                                 const System& sys,
                                 const ParameterVector& parameter_vector)
{
  this->allocate_data(qoi_indices, sys, parameter_vector);
}



inline
void SensitivityData::allocate_data(const QoISet &qoi_indices,
                                    const System& sys,
                                    const ParameterVector& parameter_vector)
{
  const unsigned int Np = parameter_vector.size();
  const unsigned int Nq = sys.qoi.size();

  if (_data.size() < Nq)
    _data.resize(Nq);

  for (unsigned int i=0; i != Nq; ++i)
    if (qoi_indices.has_index(i))
      {
        _data[i].clear();
        _data[i].resize(Np);
      }
}



inline
const Number& SensitivityData::value(unsigned int qoi_index,
                                     unsigned int parameter_index) const
{
  libmesh_assert(qoi_index < _data.size());
  libmesh_assert(parameter_index < _data[qoi_index].size());

  return _data[qoi_index][parameter_index];
}



inline
Number& SensitivityData::value(unsigned int qoi_index,
                               unsigned int parameter_index)
{
  libmesh_assert(qoi_index < _data.size());
  libmesh_assert(parameter_index < _data[qoi_index].size());

  return _data[qoi_index][parameter_index];
}

#endif // #define __sensitivity_data_h__
