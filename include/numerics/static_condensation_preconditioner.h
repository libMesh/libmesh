// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_STATIC_CONDENSATION_PRECONDITIONER_H
#define LIBMESH_STATIC_CONDENSATION_PRECONDITIONER_H

#include "libmesh/preconditioner.h"
#include "libmesh/static_condensation.h"

namespace libMesh
{
class StaticCondensationPreconditioner : public Preconditioner<Number>
{
public:
  StaticCondensationPreconditioner(StaticCondensation & sc);

  virtual bool initialized() const override { return _sc.initialized(); }

  virtual void init() override { _sc.init(); }

  virtual void setup() override { _sc.setup(); }

  virtual void apply(const NumericVector<Number> & full_rhs,
                     NumericVector<Number> & full_sol) override;

  virtual void clear() override { _sc.clear(); }

  virtual void zero() override { _sc.zero(); }

private:
  StaticCondensation & _sc;
};

inline StaticCondensationPreconditioner::StaticCondensationPreconditioner(StaticCondensation & sc)
  : Preconditioner<Number>(sc.comm()), _sc(sc)
{
}

inline void
StaticCondensationPreconditioner::apply(const NumericVector<Number> & full_rhs,
                                        NumericVector<Number> & full_sol)
{
  _sc.apply(full_rhs, full_sol);
}

} // namespace libMesh

#endif // LIBMESH_STATIC_CONDENSATION_PRECONDITIONER_H
