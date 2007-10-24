// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __enum_xdr_mode_h__
#define __enum_xdr_mode_h__

/*
 * The \p libMeshEnums namespace is the namespace all \p enum definitions
 * should be put into.
 */

// ------------------------------------------------------------
// enum XdrMode definition
namespace libMeshEnums {
  
  /**
   * Defines an \p enum for read/write mode in Xdr format.
   * \p READ, \p WRITE perform reading and writing in ASCII format,
   * and \p DECODE, \p ENCODE do the same in binary format.
   */
  enum XdrMODE
    {
      UNKNOWN = -1, ENCODE=0, DECODE, WRITE, READ
    };
}

using namespace libMeshEnums;

#endif // #define __enum_xdr_mode_h__




