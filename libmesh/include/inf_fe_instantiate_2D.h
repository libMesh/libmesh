// $Id: inf_fe_instantiate_2D.h,v 1.1 2003-01-24 17:24:38 jwpeterson Exp $

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


// Note: _No_ include guards!



/**
 * Sanity check, _without_ prior inclusion of mesh_config.h.
 *
 * This file is no typical header file. It is only to be 
 * included at the _end_ of an implementation file, so that 
 * the proper variations of the InfFE class are instantiated.
 */
#ifdef ENABLE_INFINITE_ELEMENTS



// Local includes
#include "inf_fe_macro.h"


/** 
 * Collect all 2D explicit instantiations for class InfFE
 */
SAY_THAT_INF_FE_IS(template,2,CARTESIAN) 

/* SAY_THAT_INF_FE_IS(template,2,SPHERICAL) */

/* SAY_THAT_INF_FE_IS(template,2,ELLIPSOIDAL) */



#endif
