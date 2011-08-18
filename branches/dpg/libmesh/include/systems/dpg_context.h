
// $Id: fem_context.h 4784 2011-08-03 15:29:54Z trumanellis $

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



#ifndef __dpg_context_h__
#define __dpg_context_h__

// C++ includes
#include <map>

// Local Includes
#include "fem_context.h"
#include "vector_value.h"
#include "fe_type.h"

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
#include "tensor_value.h"
#endif

namespace libMesh
{

  // Forward Declarations
  class Elem;
  class FEBase;
  class QBase;
  class Point;
  template <typename T> class NumericVector;

/**
 * This class provides all data required for a physics package
 * (e.g. an DPGSystem subclass) to perform local element residual
 * and jacobian integrations.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Truman E. Ellis 2011 
 */

// ------------------------------------------------------------
// DPGContext class definition

class DPGContext : public FEMContext
{
public:

  /**
   * Constructor.  Allocates some but fills no data structures.
   */
  DPGContext (const System &sys);

  /**
   * Destructor.
   */
  virtual ~DPGContext ();

  /**
   * Finite element objects for each test functions's interior, sides 
   * and edges.
   */
  std::map<FEType, FEBase *> element_test;
  std::map<FEType, FEBase *> side_test;
  std::map<FEType, FEBase *> edge_test;

  /**
   * Pointers to the same finite element objects, but indexed
   * by variable number
   */
  std::vector<FEBase *> element_test_var;
  std::vector<FEBase *> side_test_var;
  std::vector<FEBase *> edge_test_var;

  /**
   * Define data structure to hold the enriched element matrix to solve
   * for the optimal basis functions from the enriched space.
   */
  DenseMatrix<Number> K_test;

  /**
   * Define data structure to hold the right-hand-side vectors for
   * each field and flux DOF. External vector cooresponds to each 
   * field or flux variable, internal vector handles the respective
   * DOF.
   */
  std::vector< std::vector<DenseVector<Number> > > F_tests;
};



// ------------------------------------------------------------
// DPGContext inline methods



} // namespace libMesh

#endif
