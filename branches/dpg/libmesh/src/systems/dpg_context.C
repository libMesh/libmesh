// $Id: fem_context.C 4784 2011-08-03 15:29:54Z trumanellis $

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



#include "dof_map.h"
#include "elem.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "dpg_context.h"
#include "libmesh_logging.h"
#include "mesh_base.h"
#include "quadrature.h"
#include "system.h"
#include "time_solver.h"
#include "unsteady_solver.h" // For euler_residual

namespace libMesh
{





DPGContext::DPGContext (const System &sys)
  : FEMContext(sys)
{
}


DPGContext::~DPGContext()
{
  // We don't want to store AutoPtrs in STL containers, but we don't
  // want to leak memory either
  for (std::map<FEType, FEBase *>::iterator i = element_test.begin();
       i != element_test.end(); ++i)
    delete i->second;
  element_test.clear();

  for (std::map<FEType, FEBase *>::iterator i = side_test.begin();
       i != side_test.end(); ++i)
    delete i->second;
  side_test.clear();

  for (std::map<FEType, FEBase *>::iterator i = edge_test.begin();
       i != edge_test.end(); ++i)
    delete i->second;
  edge_test.clear();
}


} // namespace libMesh
