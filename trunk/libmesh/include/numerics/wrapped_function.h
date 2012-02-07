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



#ifndef __wrapped_function_h__
#define __wrapped_function_h__

// C++ includes



// Local Includes
#include "function_base.h"
#include "libmesh_common.h"

namespace libMesh
{



// Forward Declarations
template <typename T> class DenseVector;
class Point;


/**
 */

// ------------------------------------------------------------
// WrappedFunction class definition
template <typename Output=Number>
class WrappedFunction : public FunctionBase<Output>
{
public:

  /**
   * Constructor to wrap scalar-valued function pointers.
   */
  WrappedFunction (const System &sys,
                   Output fptr(const Point& p,
                               const Parameters& parameters,
                               const std::string& sys_name,
                               const std::string& unknown_name) = NULL,
                   const Parameters* parameters = NULL,
                   unsigned int varnum=0)
    : _sys(sys),
      _parameters(parameters),
      _fptr(fptr),
      _varnum(varnum)
  {
    this->_initialized = true;
    if (!parameters)
      _parameters = &sys.get_equation_systems().parameters;
  }

  virtual void init () {}

  virtual void clear () {}

  /**
   * @returns the scalar value of variable varnum at coordinate \p p
   * and time \p time.
   */
  virtual Output operator() (const Point& p,
			     const Real time = 0.);

  /**
   * Return function for vectors.
   * Returns in \p output the values of all system variables at the
   * coordinate \p p and for time \p time.
   */
  virtual void operator() (const Point& p,
			   const Real time,
			   DenseVector<Output>& output);

  /**
   * @returns the vector component \p i at coordinate
   * \p p and time \p time.
   */
  virtual Output operator() (const unsigned int i,
                             const Point& p,
			     const Real time=0.);

protected:

  const System& _sys;

  const Parameters* _parameters;

  Output (*_fptr)(const Point& p,
                  const Parameters& parameters,
                  const std::string& sys_name,
                  const std::string& unknown_name);

  unsigned int _varnum;
};


// ------------------------------------------------------------
// WrappedFunction inline methods


template <typename Output>
inline
Output WrappedFunction<Output>::operator() (const Point& p,
			                    const Real time)
{
  libmesh_assert(_fptr);
  libmesh_assert(_parameters);
  return _fptr(p,
               *_parameters,
               _sys.name(),
               _sys.variable_name(_varnum));
}

/**
 * Return function for vectors.
 * Returns in \p output the values of all system variables at the
 * coordinate \p p and for time \p time.
 */
template <typename Output>
inline
void WrappedFunction<Output>::operator() (const Point& p,
			                  const Real time,
			                  DenseVector<Output>& output)
{
  libmesh_assert(_fptr);
  libmesh_assert(_parameters);

  // We fill each entry of output with a single scalar component of
  // the data in our System
  const unsigned int size = output.size();
  libmesh_assert(size == _sys.n_components());

  // Loop over variables, then over each component in
  // vector-valued variables, evaluating each.
  const unsigned int n_vars = _sys.n_vars();
  for (unsigned int v = 0; v != n_vars; ++v)
    {
      const unsigned int n_components =
        _sys.variable(v).n_components();
      if (n_components == 1)
        output(_sys.variable_scalar_number(v,0)) =
	   _fptr(p, *_parameters, _sys.name(), _sys.variable_name(v));
      else
        {
          // Right now our only non-scalar variable type is the
          // SCALAR variables.  The irony is priceless.
          libmesh_assert(_sys.variable(v).type().family == SCALAR);

	  // We pass the point (j,0,0) to an old-style fptr function
	  // pointer to distinguish the different scalars within the
	  // SCALAR variable.
          for (unsigned int j=1; j != n_components; ++j)
            output(_sys.variable_scalar_number(v,j)) =
	       _fptr(Point(j,0,0), *_parameters,
                     _sys.name(), _sys.variable_name(v));
        }
    }
}


/**
 * @returns the vector component \p i at coordinate
 * \p p and time \p time.
 */
template <typename Output>
inline
Output WrappedFunction<Output>::operator() (const unsigned int i,
                                            const Point& p,
			                    const Real time)
{
  libmesh_assert(_fptr);
  libmesh_assert(_parameters);

  // Loop over variables, then over each component in
  // vector-valued variables.
  const unsigned int n_vars = _sys.n_vars();
  for (unsigned int v = 0; v != n_vars; ++v)
    {
      const unsigned int n_components =
        _sys.variable(v).n_components();
      if (n_components == 1 &&
          i == _sys.variable_scalar_number(v,0))
        return _fptr(p, *_parameters, _sys.name(), _sys.variable_name(v));
      else if (i >= _sys.variable_scalar_number(v,0) &&
               i <= _sys.variable_scalar_number(v,n_components-1))
        {
          // Right now our only non-scalar variable type is the
          // SCALAR variables.  The irony is priceless.
          libmesh_assert(_sys.variable(i).type().family == SCALAR);

	  // We pass the point (j,0,0) to an old-style fptr function
	  // pointer to distinguish the different scalars within the
	  // SCALAR variable.
          for (unsigned int j=1; j != n_components; ++j)
            if (i == _sys.variable_scalar_number(v,j))
              return _fptr(Point(j,0,0), *_parameters,
			   _sys.name(), _sys.variable_name(j));
        }
    }
  libMesh::err << "Component index " << i <<
		  " not found in system " << _sys.name() << std::endl;
  libmesh_error();
  return Output();
}



} // namespace libMesh

#endif // __wrapped_function__
