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
   * The types of quadrature rules that may 
   * be used for numerical integration
   * over geometric entities.
   */
  const QuadratureType valid_elem_rules[] = {QGAUSS,
					     QSIMPSON,
					     QTRAP};

  /**
   * The number of valid quadrature rules for
   * numerical integration over geometric entities.
   */
  const unsigned int num_valid_elem_rules = 3;


  /**
   * Returns a standard string representation
   * for the specific quadrature rule.
   */
  std::string name (const QuadratureType t);
}

#endif




