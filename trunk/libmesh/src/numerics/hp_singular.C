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

// C++ includes

// Local Includes
#include "elem.h"
#include "hp_singular.h"
#include "libmesh_logging.h"
#include "mesh.h"
#include "system.h"

#ifdef ENABLE_AMR

//-----------------------------------------------------------------
// HPSingularity implementations


void HPSingularity::select_refinement (System &system)
{
  START_LOG("select_refinement()", "HPSingularity");

  // The current mesh
  MeshBase& mesh = system.get_mesh();

  MeshBase::element_iterator       elem_it  =
    mesh.active_local_elements_begin();
  const MeshBase::element_iterator elem_end =
    mesh.active_local_elements_end(); 

  for (; elem_it != elem_end; ++elem_it)
    {
      Elem* elem = *elem_it;

      // We're only checking elements that are already flagged for h
      // refinement
      if (elem->refinement_flag() != Elem::REFINE)
        continue;

      elem->set_p_refinement_flag(Elem::REFINE);
      elem->set_refinement_flag(Elem::DO_NOTHING);

      for (std::list<Point>::iterator ppoint =
             singular_points.begin();
           ppoint != singular_points.end(); ++ppoint)
        {
          if (elem->contains_point(*ppoint))
            {
              elem->set_p_refinement_flag(Elem::DO_NOTHING);
              elem->set_refinement_flag(Elem::REFINE);
              break;
            }
        }
    }

  STOP_LOG("select_refinement()", "HPSingularity");
}

#endif // #ifdef ENABLE_AMR
