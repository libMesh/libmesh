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



#ifndef LIBMESH_SYSTEM_NORM_H
#define LIBMESH_SYSTEM_NORM_H

// Local includes
#include "libmesh/libmesh_common.h" // for Real
#include "libmesh/enum_norm_type.h"

// C++ includes
#include <vector>

namespace libMesh
{

/**
 * This class defines a norm/seminorm to be applied to a NumericVector which
 * contains coefficients in a finite element space.
 *
 * Discrete vector norms and weighted l2 combinations of Sobolev norms and
 * seminorms are representable.
 *
 * \author Roy H. Stogner
 * \date 2008
 */
class SystemNorm
{
public:

  /**
   * Constructor, defaults to DISCRETE_L2
   */
  SystemNorm();

  /**
   * Constructor, for discrete vector norms, systems with one variable,
   * and systems for which the same norm type should be used with a
   * weight of one on each variable.
   *
   * This is deliberately an implicit constructor; we want user code
   * to be able to include lines like "error_norm = L2"
   */
  SystemNorm(const FEMNormType & t);

  /**
   * Constructor, for unweighted sobolev norms on systems with multiple
   * variables.
   *
   * For a system with n variables, the final norm will be the l2 norm of the
   * n-vector of the norms in each variable.
   */
  explicit
  SystemNorm(const std::vector<FEMNormType> & norms);

  /**
   * Constructor, for weighted sobolev norms on systems with multiple
   * variables.
   *
   * For a system with n variables, the final norm will be the l2 norm of the
   * n-vector of the norms in each variable, each multiplied by weight.
   */
  SystemNorm(const std::vector<FEMNormType> & norms,
             std::vector<Real> & weights);

  /**
   * Constructor, for weighted sobolev norms on systems with multiple
   * variables and their adjoints
   *
   * For a system with n variables, the final norm computed will be of the form
   * norm_u^T*R*norm_z where R is a scaling matrix
   */
  SystemNorm(const std::vector<FEMNormType> & norms,
             std::vector<std::vector<Real> > & weights);

  /**
   * Copy Constructor
   */
  SystemNorm(const SystemNorm & s);

  /**
   * Returns true if this is purely a discrete norm
   */
  bool is_discrete() const;

  /**
   * Returns the weighted norm v^T*W*v where W represents our
   * weights matrix or weights vector times identity matrix.
   */
  Real calculate_norm(const std::vector<Real> & v);

  /**
   * Returns the weighted inner product v1^T*W*v2 where R is our weights
   */
  Real calculate_norm(const std::vector<Real> & v1,
                      const std::vector<Real> & v2);

  /**
   * Returns true if no weight matrix W is specified or an identiy matrix is specified, otherwise returns false
   */
  bool is_identity();

  /**
   * Returns the type of the norm in variable \p var
   */
  FEMNormType type(unsigned int var) const;

  /**
   * Sets the type of the norm in variable \p var
   */
  void set_type(unsigned int var, const FEMNormType & t);

  /**
   * Returns the weight corresponding to the norm in variable \p var
   */
  Real weight(unsigned int var) const;

  /**
   * Sets the weight corresponding to the norm in variable \p var
   */
  void set_weight(unsigned int var, Real w);

  /**
   * Sets the weight corresponding to the norm from the variable pair v1(var1) coming from v2(var2). See calculate_norm
   */
  void set_off_diagonal_weight(unsigned int i, unsigned int j, Real w);

  /**
   * Returns the squared weight corresponding to the norm in variable
   * \p var.  We cache that at construction time to save a few flops.
   */
  Real weight_sq(unsigned int var) const;



private:
  std::vector<FEMNormType> _norms;

  std::vector<Real> _weights;
  std::vector<Real> _weights_sq;

  /**
   * One more data structure needed to store the off diagonal
   * components for the generalize SystemNorm case
   */
  std::vector<std::vector<Real> > _off_diagonal_weights;
};



// ------------------------------------------------------------
// SystemNorm inline methods

inline
SystemNorm::SystemNorm() :
  _norms(1, DISCRETE_L2), _weights(1, 1.0), _weights_sq(1, 1.0)
{
}


inline
SystemNorm::SystemNorm(const FEMNormType & t) :
  _norms(1, t), _weights(1, 1.0), _weights_sq(1, 1.0)
{
}


inline
SystemNorm::SystemNorm(const std::vector<FEMNormType> & norms) :
  _norms(norms), _weights(1, 1.0), _weights_sq(1, 1.0)
{
  if (_norms.empty())
    _norms.push_back(DISCRETE_L2);
}


inline
SystemNorm::SystemNorm(const std::vector<FEMNormType> & norms,
                       std::vector<Real> & weights) :
  _norms(norms), _weights(weights), _weights_sq(_weights.size(), 0.0)
{
  if (_norms.empty())
    _norms.push_back(DISCRETE_L2);

  if (_weights.empty())
    {
      _weights.push_back(1.0);
      _weights_sq.push_back(1.0);
    }
  else
    for (std::size_t i=0; i != _weights.size(); ++i)
      _weights_sq[i] = _weights[i] * _weights[i];
}

inline
SystemNorm::SystemNorm(const std::vector<FEMNormType> & norms,
                       std::vector<std::vector<Real> > & weights):
  _norms(norms),
  _weights(weights.size()),
  _weights_sq(weights.size()),
  _off_diagonal_weights(weights)
{
  if(_norms.empty())
    _norms.push_back(DISCRETE_L2);

  if (_weights.empty())
    {
      _weights.push_back(1.0);
      _weights_sq.push_back(1.0);
    }
  else
    {
      // Loop over the entries of the user provided matrix and store its entries in
      // the _off_diagonal_weights or _diagonal_weights
      for(std::size_t i=0; i!=_off_diagonal_weights.size(); ++i)
        {
          if(_off_diagonal_weights[i].size() > i)
            {
              _weights[i] = _off_diagonal_weights[i][i];
              _off_diagonal_weights[i][i] = 0;
            }
          else
            _weights[i] = 1.0;
        }
      for (std::size_t i=0; i != _weights.size(); ++i)
        _weights_sq[i] = _weights[i] * _weights[i];
    }
}

inline
SystemNorm::SystemNorm(const SystemNorm & s) :
  _norms(s._norms), _weights(s._weights), _weights_sq(s._weights_sq)
{
}


inline
bool SystemNorm::is_discrete() const
{
  libmesh_assert (!_norms.empty());

  if (_norms[0] == DISCRETE_L1 ||
      _norms[0] == DISCRETE_L2 ||
      _norms[0] == DISCRETE_L_INF)
    return true;

  return false;
}


inline
FEMNormType SystemNorm::type(unsigned int var) const
{
  libmesh_assert (!_norms.empty());

  std::size_t i = (var < _norms.size()) ? var : _norms.size() - 1;

  return _norms[i];
}



inline
void SystemNorm::set_type(unsigned int var, const FEMNormType & t)
{
  libmesh_assert (!_norms.empty());

  if (var >= _norms.size())
    _norms.resize(var+1, t);

  _norms[var] = t;
}


inline
Real SystemNorm::weight(unsigned int var) const
{
  libmesh_assert (!_weights.empty());

  return (var < _weights.size()) ? _weights[var] : 1.0;
}


inline
void SystemNorm::set_weight(unsigned int var, Real w)
{
  libmesh_assert (!_weights.empty());

  if (var >= _weights.size())
    {
      _weights.resize(var+1, 1.0);
      _weights_sq.resize(var+1, 1.0);
    }

  _weights[var] = w;
  _weights_sq[var] = w*w;
}

inline
void SystemNorm::set_off_diagonal_weight(unsigned int i,
                                         unsigned int j,
                                         Real w)
{
  libmesh_assert (!_weights.empty());

  if (i >= _off_diagonal_weights.size())
    {
      _off_diagonal_weights.resize(i+1);
    }

  if (j >= _off_diagonal_weights[i].size())
    {
      _off_diagonal_weights[i].resize(j+1, 0.);
    }

  _off_diagonal_weights[i][j] = w;

}


inline
Real SystemNorm::weight_sq(unsigned int var) const
{
  libmesh_assert (!_weights_sq.empty());

  return (var < _weights_sq.size()) ? _weights_sq[var] : 1.0;
}


inline
Real SystemNorm::calculate_norm(const std::vector<Real> & v1,
                                const std::vector<Real> & v2)
{
  // The vectors are assumed to both be vectors of the (same number
  // of) components
  std::size_t vsize = v1.size();
  libmesh_assert_equal_to (vsize, v2.size());

  // We'll support implicitly defining weights, but if the user sets
  // more weights than he uses then something's probably wrong
  std::size_t diagsize = this->_weights.size();
  libmesh_assert_greater_equal (vsize, diagsize);

  // Initialize the variable val
  Real val = 0.;

  // Loop over all the components of the system with explicit
  // weights
  for(std::size_t i = 0; i != diagsize; i++)
    {
      val += this->_weights[i] * v1[i] * v2[i];
    }
  // Loop over all the components of the system with implicit
  // weights
  for(std::size_t i = diagsize; i < vsize; i++)
    {
      val += v1[i] * v2[i];
    }

  // Loop over the components of the system
  std::size_t nrows = this->_off_diagonal_weights.size();
  libmesh_assert_less_equal (vsize, nrows);

  for(std::size_t i = 0; i != nrows; i++)
    {
      std::size_t ncols = this->_off_diagonal_weights[i].size();
      for(std::size_t j=0; j != ncols; j++)
        {
          // Note that the diagonal weights here were set to zero
          // in the constructor
          val += this->_off_diagonal_weights[i][j] * v1[i] * v2[j];
        }
    }

  return(val);
}

inline
Real SystemNorm::calculate_norm(const std::vector<Real> & v1)
{
  return this->calculate_norm(v1,v1);
}

inline
bool SystemNorm::is_identity()
{
  std::size_t nrows = this->_off_diagonal_weights.size();

  // If any of the off-diagonal elements is not 0, then we are in the non-identity case
  for(std::size_t i = 0; i != nrows; i++)
    {
      std::size_t ncols = this->_off_diagonal_weights[i].size();
      for(std::size_t j = 0; j != ncols; j++)
        {
          if(_off_diagonal_weights[i][j] != 0)
            {
              return(false);
            }
        }
    }

  // If any of the diagonal elements is not 1, then we are in the non-identity case
  nrows = this->_weights.size();
  for(std::size_t i = 0; i != nrows; i++)
    if(_weights[i] != 1)
      return(false);

  // If all the off-diagonals elements are 0, and diagonal elements 1, then we are in an identity case
  return(true);
}

} // namespace libMesh

#endif // LIBMESH_SYSTEM_NORM_H
