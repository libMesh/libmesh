// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_MATLAB_IO_H
#define LIBMESH_MATLAB_IO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_input.h"

// C++ includes

namespace libMesh
{

// Forward declarations
template <typename> class MeshBaseTempl;
typedef MeshBaseTempl<Real> MeshBase;

/**
 * This class implements reading meshes in the Matlab PDE toolkit
 * in a proprietary  format.
 *
 * A VALID INPUT FILE for this type of mesh should be
 * generated in Matlab with the following steps:
 * 1.) Draw the domain and triangulate it in the GUI
 * 2.) Export the mesh to matlab using Mesh->Export Mesh
 * 3.) Create a file with this script:
 * \code
 * fid = fopen('filename', 'w');
 * fprintf(fid, '%d %d \n', length(p), length(t));
 * fprintf(fid, '%f %f \n', p);
 * fprintf(fid, '%d %d %d %d \n', t);
 * fclose(fid);
 * \endcode
 *
 * What's going on here?
 * There is no standard for exporting PDE toolkit meshes
 * to files in Matlab.  When you choose "export mesh" in the GUI,
 * it returns three matrices that it likes to call
 * p, e, and t.  All meshes (as far as I can tell) that
 * come from the PDE toolkit are 2D triangle meshes.
 *
 * p is the point matrix...
 * Row 1: x coordinate
 * Row 2: y coordinate
 *
 * e is the edge matrix ...
 * Row 1: starting point number          (dummy)
 * Row 2: ending point number            (dummy)
 * Row 3: starting parameter value (?)   (dummy)
 * Row 4: ending parameter value (?)     (dummy)
 * Row 5: boundary segment number (?)    (dummy)
 * Row 6: left-hand subdomain number     (dummy)
 * Row 7: right-hand subdomain number    (dummy)
 *
 * t is the triangle matrix ...
 * Row 1: Node number 1
 * Row 2: Node number 2
 * Row 3: Node number 3
 * Row 4: subdomain number               (dummy)
 *
 * There are some important things to notice here:
 * o The implied ordering of the p matrix is 1..N
 * o The e matrix is entirely irrelevant in this code
 * o All of the matrices are row based
 *
 * \author John W. Peterson
 * \date 2004
 */
class MatlabIO : public MeshInput<MeshBase>
{
public:
  /**
   *  Constructor.  Takes a non-const Mesh reference which it
   * will fill up with elements.
   */
  explicit
  MatlabIO (MeshBase &);

  /**
   * Reads in a matlab data file based on the string
   * you pass it.
   */
  virtual void read (const std::string & name) override;

private:
  /**
   * Implementation of the read() function.  This function
   * is called by the public interface function and implements
   * reading the file.
   */
  void read_stream (std::istream & in);
};


// ------------------------------------------------------------
// MatlabIO inline members
inline
MatlabIO::MatlabIO (MeshBase & mesh_in) :
  MeshInput<MeshBase>  (mesh_in)
{}

} // namespace libMesh


#endif // LIBMESH_MATLAB_IO_H
