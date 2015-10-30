// $Id: equation_systems_macro.h,v 1.1 2003-04-09 19:26:57 ddreyer Exp $

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


#ifndef __equation_systems_macro_h__
#define __equation_systems_macro_h__


// Local includes
#include "mesh_config.h"


/**
 * This macro helps in instantiating specific versions
 * of the \p EquationSystems<T_sys> class.  Once you add
 * another child of \p SystemBase, simply add it here, too.
 * Then \p EquationSystems<T_sys> will get instantiated
 * properly.
 */
#if defined(USE_COMPLEX_NUMBERS) 

# define INSTANTIATE_EQUATION_SYSTEMS template class EquationSystems<GeneralSystem>;  \
                                      template class EquationSystems<ThinSystem>;     \
                                      template class EquationSystems<NewmarkSystem>;  \
                                      template class EquationSystems<FrequencySystem>

#elif defined(USE_REAL_NUMBERS)

# define INSTANTIATE_EQUATION_SYSTEMS template class EquationSystems<GeneralSystem>;  \
                                      template class EquationSystems<ThinSystem>;     \
                                      template class EquationSystems<NewmarkSystem>

#else

 CHOKE THIS WHILE COMPILING...

#endif


#endif // ifndef __equation_systems_macro_h__
