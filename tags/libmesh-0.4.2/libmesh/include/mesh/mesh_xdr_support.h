// $Id: mesh_xdr_support.h,v 1.2 2004-01-03 15:37:42 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __mesh_xdr_support_h__
#define __mesh_xdr_support_h__

// Local includes
#include "libmesh_common.h"
#include "xdrIO.h" // for XdrIO::XdrIO_TYPE
#include <string>
#include <vector>

// Forward declarations
class Node;
class Mesh;
class Elem;
class BoundaryInfo;

/**
 * This class allows us to
 * reuse code, while still
 * maintaining a uniform
 * interface to XDR in mesh_base.h
 * There are two major functions
 * here:
 *
 * @begin{itemize}
 * @item mesh_interface
 * @item soln_interface
 * @end{itemize}
 *
 * When someone calls read_xdr(filename)
 * The program calls XdrInterface::mesh_interface(filename, DECODE, dim)
 *
 * @author John W. Peterson, 2002
 */

class XdrInterface
{
 public:
  /**
   * Constructor.
   */
  XdrInterface() {}
  
  /**
   * Destructor.
   */
  ~XdrInterface() {}

  /**
   * Read OR Write Mesh files
   */
  void mesh_interface(const std::string& name,
		      const XdrIO::XdrIO_TYPE access,
		      std::vector<Node*>& nodes,
		      std::vector<Elem*>& elements,
		      BoundaryInfo& boundary_info,
		      Mesh& mesh);
  
  /**
   * Read OR Write Soln files
   */
  void soln_interface(const std::string& name,
		      const XdrIO::XdrIO_TYPE access,
		      std::vector<Number>& soln,
		      std::vector<std::string>& var_names,
		      Mesh& mesh);
  
  /*
  void soln_interface(const std::string& name,
		      const XdrIO::XdrIO_TYPE access);
  */

 private:

  /**
   * Actual implementation of reading OR writing Soln files.
   */
  void soln_interface_impl(const std::string& name,
			   const XdrIO::XdrIO_TYPE access,
			   std::vector<Real>& soln,
			   std::vector<std::string>& var_names,
			   Mesh& mesh);
  
};

#endif
