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

#ifndef LIBMESH_XDR_SOLN_H
#define LIBMESH_XDR_SOLN_H

// Local includes
#include "libmesh/xdr_mgf.h"

// C++ includes

namespace libMesh
{

// forward declarations
class XdrSHEAD;

/**
 * The \p XdrSOLN class.
 * This class is responsible
 * for reading/writing
 * information about the solution
 * to \p xdr style binary files.
 *
 * \author Bill Barth
 * \author Robert McLay
 * \date 2000
 */
class XdrSOLN: public XdrMGF
{
public:
  /**
   * Constructor.
   * Initializes \p m_wrtVar to -1.
   */
  XdrSOLN() : m_wrtVar(-1) {}

  /**
   * Calls the \p init method
   * in the parent class, \p XdrMGF
   * with the appropriate parameters.
   *
   * \param type One of: \p UNKNOWN, \p ENCODE, \p DECODE
   * \param fn const char pointer to a file name
   * \param icnt Number to be appended to file e.g. \p name.soln.0000
   */
  void init(XdrIO_TYPE type, const char * fn, int icnt)
  {XdrMGF::init (type, fn, "soln",icnt);}

  /**
   * Destructor.
   */
  ~XdrSOLN() {}

  /**
   * Read/Write the solution header.
   * Uses \p xdr_int found
   * in \p rpc/rpc.h.
   *
   * \param hd Pointer to an \p xdr solution header object
   * @return 1 on success
   */
  int header(XdrSHEAD * hd);

  /**
   * Read/Write solution values.
   *
   * \param array Pointer to array of \p Reals to be read/written
   * \param size Size of individual variables to be written
   * @return m_wrtVar*size
   */
  int values(Real * array, int size) { return dataBlk(array, m_wrtVar, size);}

private:
  int m_wrtVar;
};


} // namespace libMesh

#endif // LIBMESH_XDR_SOLN_H
