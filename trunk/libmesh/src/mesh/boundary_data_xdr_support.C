// $Id: boundary_data_xdr_support.C,v 1.1 2003-05-14 11:54:37 ddreyer Exp $

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
#include <fstream>

// Local includes
#include "boundary_data.h"
#include "mesh_base.h"



//------------------------------------------------------
// BoundaryData functions
void BoundaryData::read_xdr (const std::string& /* name */)
{
  std::cerr << "ERROR: not yet implemented." << std::endl;
  error();
}




void BoundaryData::read_xdr_binary (const std::string& /* name */)
{
  std::cerr << "ERROR: not yet implemented." << std::endl;
  error();
}





void BoundaryData::write_xdr (const std::string& /* name */)
{
  std::cerr << "ERROR: not yet implemented." << std::endl;
  error();
}




void BoundaryData::write_xdr_binary (const std::string& /* name */)
{
  std::cerr << "ERROR: not yet implemented." << std::endl;
  error();
}


