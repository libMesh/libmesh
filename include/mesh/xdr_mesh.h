// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_XDR_MESH_H
#define LIBMESH_XDR_MESH_H

// Local Includes
#include "libmesh/xdr_mgf.h"

// C++ includes

namespace libMesh
{

// forward declarations
class XdrMHEAD;

/**
 * The \p XdrMESH class.
 * This class is responsible
 * for reading/writing
 * information about the mesh
 * to \p xdr style binary files.
 *
 * \author Bill Barth
 * \author Robert McLay
 * \date 2000
 */
class XdrMESH: public XdrMGF
{
public:

  /**
   * Constructor.  Initializes
   * \p m_dim to -1.
   */
  XdrMESH() : m_dim(-1) {}

  /**
   * Calls the \p init method
   * in the parent class, \p XdrMGF
   * with the appropriate parameters.
   *
   * \param type One of: \p UNKNOWN, \p ENCODE, \p DECODE
   * \param fn const char pointer which points to the filename
   * \param icnt Number to be appended to file e.g. \p name.mesh.0000
   * \param dim Problem dimension (always three in MGF)
   */
  void init(XdrIO_TYPE type, const char * fn, int icnt, int dim=3)
  { XdrMGF::init(type, fn, "mesh", icnt); m_dim = dim;}

  /**
   * Destructor.
   */
  ~XdrMESH() {}

  /**
   * Read/Write the mesh_base.header.
   * Uses \p xdr_int found
   * in \p rpc/rpc.h.
   *
   * \param hd Pointer to an \p xdr mesh_base.header object
   * @return 1 on success
   */
  int header(XdrMHEAD * hd);

  /**
   * Read/Write an integer connectivity array
   *
   * \param array Pointer to an array of \p ints
   * \param numvar Total number of variables to be read/written
   * \param num Basically a dummy parameter
   * @return numvar*num
   */
  int Icon(int * array, int numvar, int num)   { return dataBlk(array, numvar, num);}

  /**
   * Read/Write a coord of appropriate size.
   *
   * \param array Pointer to an array of \p Reals
   * \param dim The size of "array"
   * \param size Size of \p array (number of elements)
   * @return dim*size
   */
  int coord(Real * array, int dim, int size)  { return dataBlk(array, dim, size);}

  /**
   * Read/Write a BC of appropriate size
   *
   * \param array Pointer to an array of \p Reals
   * \param size Size of \p array (number of elements)
   * @return 3*size
   */
  int BC(int * array, int size)                { return dataBlk(array, 3, size);}


private:

  /**
   * Dimension of the mesh
   */
  int m_dim;
};

} // namespace libMesh

#endif // LIBMESH_XDR_MESH_H
