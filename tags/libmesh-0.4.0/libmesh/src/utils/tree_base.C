// $Id: tree_base.C,v 1.2 2003-05-15 23:34:36 benkirk Exp $

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



// C++ includes

// Local includes
#include "tree_base.h"
#include "tree.h"


// ------------------------------------------------------------
// TreeBase class methods
AutoPtr<TreeBase> TreeBase::build (const unsigned int N,
				   const MeshBase& m,
				   const unsigned int level)
{
  switch (N)
    {
    case 4:
      {
	AutoPtr<TreeBase> ap(new Tree<4>(m, level));
	return ap;
      }

    case 8:
      {
	AutoPtr<TreeBase> ap(new Tree<8>(m, level));
	return ap;
      }

    default:
      {
	std::cerr << "ERROR: Bad N = " << N << std::endl;
	error();
      }
    }

  error();
  AutoPtr<TreeBase> ap(NULL);
  return ap;
}







