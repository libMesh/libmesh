// $Id: libmesh_base.h,v 1.2 2003-09-02 18:02:38 benkirk Exp $

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



#ifndef __libmesh_base_h__
#define __libmesh_base_h__


// C++ includes

// Local includes





/**
 * The \p libMeshBase is the base class for \p libMesh, however it is
 * separate for a reason.  \p libMeshBase provides the local processor id
 * and the number of processors in the simulation, which is in turn used
 * in the \p here() macro (and many other places).  However, we need to
 * split this off from \p libMesh, otherwise circular dependencies would arise
 * because \p libMesh depends on \p here(), get it?
 */
class libMeshBase
{

protected:

  /**
   * This class only contains static members and should never be
   * instantiated, so it has a protected default constructor.
   */
  libMeshBase() {}
  
public:
  
  /**
   * @returns the number of processors used in the current simulation.
   */
  static unsigned int n_processors();

  /**
   * @returns the index of the local processor.
   */
  static unsigned int processor_id();

  
protected:

  /**
   * Total number of processors used.
   */
  static int _n_processors;

  /**
   * The local processor id.
   */
  static int _processor_id;
};



// ------------------------------------------------------------
// libMeshBase inline member functions
inline
unsigned int libMeshBase::n_processors()
{
  return static_cast<unsigned int>(_n_processors);
}



inline
unsigned int libMeshBase::processor_id()
{
  return static_cast<unsigned int>(_processor_id);
}



#endif // #define __libmesh_base_h__
