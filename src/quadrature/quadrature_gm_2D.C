// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local includes
#include "libmesh/quadrature_gm.h"

namespace libMesh
{



void QGrundmann_Moller::init_2D(const ElemType type_in,
                                unsigned int p)
{
  // Nearly all GM rules contain negative weights, so if you are not
  // allowing rules with negative weights, we cannot continue!
  if (!allow_rules_with_negative_weights)
    libmesh_error_msg("You requested a Grundmann-Moller rule but\n"        \
                      << "are not allowing rules with negative weights!\n" \
                      << "Either select a different quadrature class or\n" \
                      << "set allow_rules_with_negative_weights==true.");

  switch (type_in)
    {
    case TRI3:
    case TRI6:
      {
        switch(_order + 2*p)
          {

          default:
            {
              // Untested above _order=23 but should work...
              gm_rule((_order + 2*p)/2, /*dim=*/2);
              return;
            }
          } // end switch (order)
      } // end case TRI3, TRI6

    default:
      libmesh_error_msg("ERROR: Unsupported element type: " << type_in);
    } // end switch (type_in)
}

} // namespace libMesh
