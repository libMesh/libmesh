// $Id: quadrature_rules.h,v 1.1 2003-02-06 17:58:34 ddreyer Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __quadrature_rules_h__
#define __quadrature_rules_h__

// C++ includes
#include <string>

// Local includes
#include "enum_quadrature_type.h"


/**
 * A namespace for utility functions for
 * quadrature rules.
 */
namespace QuadratureRules
{

  /**
   * The number of quadrature rules that
   * are defined (INVALD_Q_RULE excluded).  
   * You might have to update this if you
   * add a new one!
   */
  const unsigned int num_rules = 5;

  /**
   * Returns a standard string representation
   * for the specific quadrature rule.
   */
  std::string name (const QuadratureType t);
}

#endif




