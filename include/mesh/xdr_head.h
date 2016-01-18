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

#ifndef LIBMESH_XDR_HEAD_H
#define LIBMESH_XDR_HEAD_H

// Local includes
#include "libmesh/xdr_mgf.h"

namespace libMesh
{

/**
 * The \p XdrHEAD class.  This is a base class for deriving either
 * solution (\p XdrSHEAD) or mesh (\p XdrMHEAD) header interface
 * classes.
 *
 * \author Bill Barth
 * \author Robert McLay
 * \date 2000
 */
class XdrHEAD
{
public:
  /**
   * Constructor.
   */
  XdrHEAD();

  /**
   * Destructor.
   */
  virtual ~XdrHEAD();

  /**
   * Set the mesh/solution file id.
   */
  void setId(const char * id)           { delete [] mp_id; mp_id = cpyString(id); }

  /**
   * Get the mesh/solution file id.
   */
  const char * getId() const            { return mp_id; }

  /**
   * Set the mesh/solution file title.
   */
  void setTitle(const char * title)     { delete [] mp_title; mp_title = cpyString(title); }

  /**
   * Get the mesh/solution file title.
   */
  const char * getTitle() const         { return mp_title; }

  /**
   * Set the total number of
   * nodes in the mesh/solution file.
   */
  void setNumNodes(int numNodes)       { m_numNodes = numNodes; }

  /**
   * Get the total number of
   * nodes in the mesh/solution file.
   */
  int  getNumNodes() const             { return m_numNodes; }

  /**
   * Set the number of
   * boundary conditions in the
   * mesh/solution file.
   */
  void setNumBCs(int numBCs)           { m_numBCs = numBCs; }

  /**
   * Get the number of
   * boundary conditions in
   * them mesh/solution file.
   */
  int  getNumBCs() const               { return m_numBCs; }

  /**
   * Set the string size of the
   * mesh/solution file. (?)
   */
  void setStrSize(int strSize)         { m_strSize = strSize; }

  //     /**
  //      * Set the string size of the
  //      * mesh /solutionfile. (?)
  //      */
  //     int  getStrSize() const              { return m_strSize; }

protected:

  /**
   * Number of variables written
   * to output, e.g. u,v,w,p,T = 5
   */
  int m_wrtVar;

  /**
   * Total number of variables,
   * may differ from the total
   * number of variables actually
   * written.
   */
  int m_numvar;

  /**
   * The mesh file number
   * which corresponds to a given
   * solution file.
   */
  int m_meshCnt;

  /**
   * The internal solution number.
   */
  int m_kstep;

  /**
   * Number of elemetns in the
   * solution/mesh.
   */
  int m_numel;

  /**
   * Number of nodes in the
   * solution/mesh.
   */
  int m_numNodes;

  /**
   * Total mesh weighting i.e.
   * How many nodes are there
   * and where are they?
   */
  int m_sumWghts;

  /**
   * Number of boundary
   * conditions in the solution/mesh.
   */
  int m_numBCs;

  /**
   * String size (Not sure of what?)
   */
  int m_strSize;

  /**
   * An ID string for the file.
   */
  char * mp_id;

  /**
   * A title string for the file.
   */
  char * mp_title;

  /**
   * User's simulation title
   */
  char * mp_userTitle;

  /**
   * List of null-separated variable names.
   */
  char * mp_varTitle;

  /**
   * Current solution time.
   */
  xdr_Real m_time;

  /**
   * Uses std::memcpy to create an exact
   * copy of \p src, then returns
   * that copy.  Note: I don't know
   * where the memory allocated
   * for this copy gets deleted!
   *
   * @return Copy of \p src
   */
  char * cpyString(const char * src, int len = -1);

private:
  XdrHEAD(const XdrHEAD &);
  const XdrHEAD & operator=(const XdrHEAD &);
};


} // namespace libMesh


#endif // LIBMESH_XDR_HEAD_H
