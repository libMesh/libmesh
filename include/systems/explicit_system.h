// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_EXPLICIT_SYSTEM_H
#define LIBMESH_EXPLICIT_SYSTEM_H

// Local Includes
#include "libmesh/system.h"

namespace libMesh
{

/**
 * \brief Manages consistently variables, degrees of freedom, and coefficient
 * vectors for explicit systems.
 *
 * The ExplicitSystem provides only "right hand side" storage, which
 * should be sufficient for solving most types of explicit problems, i.e.,
 * problems that do not need to store a sparse matrix.
 *
 * The ExplicitSystem class is meant for problems where a given right hand
 * side is directly applied to the (differential) operator. In this case the
 * assembly routine can directly compute the product A*x without constructing
 * a sparse matrix first.
 *
 * \note Additional vectors/matrices can be added via parent class
 * interfaces.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 */
class ExplicitSystem : public System
{
public:

  /**
   * Constructor.
   */
  ExplicitSystem (EquationSystems & es,
                  const std::string & name,
                  const unsigned int number);

  /**
   * Special functions.
   * - This class has the same restrictions/defaults as its base class.
   * - The destructor is defaulted out-of-line.
   */
  ExplicitSystem (const ExplicitSystem &) = delete;
  ExplicitSystem & operator= (const ExplicitSystem &) = delete;
  ExplicitSystem (ExplicitSystem &&) = default;
  ExplicitSystem & operator= (ExplicitSystem &&) = delete;
  virtual ~ExplicitSystem ();

  /**
   * The type of system.
   */
  typedef ExplicitSystem sys_type;

  /**
   * The type of the parent
   */
  typedef System Parent;

  /**
   * \returns A reference to *this.
   */
  sys_type & system () { return *this; }

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear () override;

  /**
   * Prepares \p qoi for quantity of interest assembly, then calls
   * user qoi function.
   * Can be overridden in derived classes.
   */
  virtual void assemble_qoi (const QoISet & qoi_indices = QoISet()) override;

  /**
   * Prepares \p adjoint_rhs for quantity of interest derivative assembly,
   * then calls user qoi derivative function.
   * Can be overridden in derived classes.
   */
  virtual void assemble_qoi_derivative (const QoISet & qoi_indices = QoISet(),
                                        bool include_liftfunc = true,
                                        bool apply_constraints = true) override;

  /**
   * For explicit systems, just assemble the system which should directly
   * compute A*x.
   */
  virtual void solve () override;

  /**
   * \returns \p "Explicit".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const override { return "Explicit"; }

  /**
   * The system matrix.  Implicit systems are characterized by
   * the need to solve the linear system Ax=b.  This is the
   * right-hand-side vector b.
   */
  NumericVector<Number> * rhs;

private:

  /**
   * Add the system right-hand-side vector to the \p _vectors data structure.
   * Useful in initialization.
   */
  void add_system_rhs ();
};

} // namespace libMesh

#endif // LIBMESH_EXPLICIT_SYSTEM_H
