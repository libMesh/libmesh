// $Id: hp_singular.C,v 1.1 2006-10-05 20:50:15 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2006  Benjamin S. Kirk, John W. Peterson
  
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


// C++ includes

// Local Includes
#include "elem.h"
#include "hp_singular.h"
#include "libmesh_logging.h"
#include "mesh.h"
#include "system.h"

//-----------------------------------------------------------------
// HPSingularity implementations


void HPSingularity::select_refinement (System &system)
{
  START_LOG("select_refinement()", "HPSingularity");

  // The current mesh
  const MeshBase& mesh = system.get_mesh();

  MeshBase::const_element_iterator       elem_it  =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator elem_end =
    mesh.active_local_elements_end(); 

  for (; elem_it != elem_end; ++elem_it)
    {
      const Elem* elem = *elem_it;

      // We're only checking elements that are already flagged for h
      // refinement
      if (elem->refinement_flag() != Elem::REFINE)
        continue;
    }

  STOP_LOG("select_refinement()", "HPCoarsenTest");
}
