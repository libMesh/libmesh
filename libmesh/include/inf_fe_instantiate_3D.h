// $Id: inf_fe_instantiate_3D.h,v 1.5 2003-09-25 21:46:55 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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
 * Sanity check, _without_ prior inclusion of libmesh_config.h.
 *
 * This file is no typical header file. It is only to be 
 * included at the _end_ of an implementation file, so that 
 * the proper variations of the InfFE class are instantiated.
 */
#ifdef ENABLE_INFINITE_ELEMENTS



// Local includes
#include "inf_fe_macro.h"


/** 
 * Collect all 3D explicit instantiations for class InfFE
 */
INSTANTIATE_INF_FE(3,CARTESIAN);

/* INSTANTIATE_INF_FE(3,SPHERICAL); */

/* INSTANTIATE_INF_FE(3,ELLIPSOIDAL); */



#endif
