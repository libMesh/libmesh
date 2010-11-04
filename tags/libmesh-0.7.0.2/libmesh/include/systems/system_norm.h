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



#ifndef __system_norm_h__
#define __system_norm_h__

// C++ includes
#include <vector>

// Local includes
#include "libmesh_common.h" // for Real
#include "enum_norm_type.h"

namespace libMesh
{

// Forward Declarations

/**
 * This class defines a norm/seminorm to be applied to a NumericVector which
 * contains coefficients in a finite element space.
 *
 * Discrete vector norms and weighted l2 combinations of Sobolev norms and
 * seminorms are representable.
 *
 * @author Roy H. Stogner 2008
 */

// ------------------------------------------------------------
// SystemNorm class definition
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
   */
  SystemNorm(const FEMNormType &t);
  
  /**
   * Constructor, for unweighted sobolev norms on systems with multiple
   * variables.
   *
   * For a system with n variables, the final norm will be the l2 norm of the
   * n-vector of the norms in each variable.
   */
  SystemNorm(const std::vector<FEMNormType> &norms);

  /**
   * Constructor, for weighted sobolev norms on systems with multiple
   * variables.
   *
   * For a system with n variables, the final norm will be the l2 norm of the
   * n-vector of the norms in each variable, each multiplied by weight.
   */
  SystemNorm(const std::vector<FEMNormType> &norms, std::vector<Real> &weights);

  /**
   * Copy Constructor
   */
  SystemNorm(const SystemNorm &s);
  
  /**
   * Returns true if this is purely a discrete norm
   */
  bool is_discrete() const;

  /**
   * Returns the type of the norm in variable \p var
   */
  FEMNormType type(unsigned int var) const;

  /**
   * Sets the type of the norm in variable \p var
   */
  void set_type(unsigned int var, const FEMNormType& t);

  /**
   * Returns the weight corresponding to the norm in variable \p var
   */
  Real weight(unsigned int var) const;

  /**
   * Sets the weight corresponding to the norm in variable \p var
   */
  void set_weight(unsigned int var, Real w);

  /**
   * Returns the squared weight corresponding to the norm in variable
   * \p var.  We cache that at construction time to save a few flops.
   */
  Real weight_sq(unsigned int var) const;

private:
  std::vector<FEMNormType> _norms;

  std::vector<Real> _weights;
  std::vector<Real> _weights_sq;
};



// ------------------------------------------------------------
// SystemNorm inline methods

inline
SystemNorm::SystemNorm() :
    _norms(1, DISCRETE_L2), _weights(1, 1.0), _weights_sq(1, 1.0)
{ 
}


inline
SystemNorm::SystemNorm(const FEMNormType &t) :
    _norms(1, t), _weights(1, 1.0), _weights_sq(1, 1.0)
{ 
}


inline
SystemNorm::SystemNorm(const std::vector<FEMNormType> &norms) :
    _norms(norms), _weights(1, 1.0), _weights_sq(1, 1.0)
{
  if (_norms.empty())
    _norms.push_back(DISCRETE_L2);
}


inline
SystemNorm::SystemNorm(const std::vector<FEMNormType> &norms,
		       std::vector<Real> &weights) :
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
    for (unsigned int i=0; i != _weights.size(); ++i)
      _weights_sq[i] = _weights[i] * _weights[i];
}


inline
SystemNorm::SystemNorm(const SystemNorm &s) :
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

  unsigned int i = (var < _norms.size()) ? var : _norms.size() - 1;

  return _norms[i];
}



inline
void SystemNorm::set_type(unsigned int var, const FEMNormType &t)
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
Real SystemNorm::weight_sq(unsigned int var) const
{
  libmesh_assert (!_weights_sq.empty());
  
  return (var < _weights_sq.size()) ? _weights_sq[var] : 1.0;
}


} // namespace libMesh

#endif // #define __system_norm_h__
