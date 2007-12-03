// $Id: legacy_xdr_io.C 2560 2007-12-03 17:52:20Z benkirk $

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



// C++ includes
#include <iostream>
#include <iomanip>

#include <vector>
#include <string>

// Local includes
#include "xdr_io.h"
#include "legacy_xdr_io.h"
#include "xdr_cxx.h"
#include "enum_xdr_mode.h"
#include "mesh_base.h"







// ------------------------------------------------------------
// XdrIO members
XdrIO::XdrIO (MeshBase& mesh, const bool binary) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase>(mesh),
  _binary (binary)
{
}




XdrIO::XdrIO (const MeshBase& mesh, const bool binary) :
  MeshOutput<MeshBase>(mesh),
  _binary (binary)
{
}




XdrIO::~XdrIO ()
{
}




bool & XdrIO::binary ()
{
  return _binary;
}




bool XdrIO::binary () const
{
  return _binary;
}


void XdrIO::read (const std::string& name)
{
  Xdr io (name, this->binary() ? DECODE : READ);

  // get the version string
  std::string version;
  io.data(version);
  
  std::cout << "version=" << version << std::endl;

  // Check for a legacy version format.
  if (!(version.find("libMesh") < version.size()))
    {
      io.close();
      LegacyXdrIO(MeshInput<MeshBase>::mesh(), this->binary()).read(name);
      return;
    }
  
  
}



void XdrIO::write (const std::string& name)
{

  if (true)
    {
      LegacyXdrIO(MeshOutput<MeshBase>::mesh(), this->binary()).write(name);
      return;
    }

  
  Xdr io (name, this->binary() ? ENCODE : WRITE);

  // set the version string
  std::string version("libMesh-0.7.0+");
  io.data(version, "# File Format Identifier");
  
  
}
