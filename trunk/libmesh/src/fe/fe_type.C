
// $Id: fe_type.C,v 1.3 2006-04-05 16:14:28 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include "fe_type.h"
#include "quadrature_clough.h"
#include "quadrature_gauss.h"

// ---------------------------------------
// FEType class members

AutoPtr<QBase>
FEType::default_quadrature_rule (const unsigned int dim,
                                 const int extraorder) const
{
  // Clough elements have at least piecewise cubic functions
  if (family == CLOUGH)
    return AutoPtr<QBase>
      (new QClough(dim,
        static_cast<Order>
          (std::max(static_cast<unsigned int>
            (this->default_quadrature_order()),7 + extraorder))));
  
  return AutoPtr<QBase>
    (new QGauss(dim, static_cast<Order>(this->default_quadrature_order()
                                        + extraorder)));
}
