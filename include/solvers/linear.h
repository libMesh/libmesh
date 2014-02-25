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



#ifndef LIBMESH_LINEAR_H
#define LIBMESH_LINEAR_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/solver.h"

// C++ includes

namespace libMesh
{


/**
 * This is a generic class that defines a linear to be used in a
 * simulation.  A user can define a linear by deriving from this
 * class and implementing certain functions.
 *
 * @author Benjamin S. Kirk, 2003-2004.
 */

// ------------------------------------------------------------
// Linear class definition

template <class T = Solver>
class Linear : public T
{
public:

  /**
   * Constructor. Requires a reference to a system to be solved.
   */
  explicit
  Linear (EquationSystems& es);

  /**
   * Constructor.  Requires a referece to the \p EquationSystems object.
   */
  Linear (EquationSystems& es,
          const std::string& name,
          const unsigned int number);

  /**
   * Destructor.
   */
  ~Linear ();
};



// ------------------------------------------------------------
// Linear inline members
template <class T>
inline
Linear<T>::Linear(EquationSystems& es) :
  T (es)
{
  libmesh_deprecated();
}



template <class T>
inline
Linear<T>::Linear (EquationSystems& es,
                   const std::string& name,
                   const unsigned int number) :
  T (es, name, number)
{
  libmesh_deprecated();
}



template <class T>
inline
Linear<T>::~Linear ()
{
}


} // namespace libMesh


#endif // LIBMESH_LINEAR_H
