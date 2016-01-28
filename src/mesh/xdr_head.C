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

// Local includes
#include "libmesh/xdr_head.h"

namespace libMesh
{

// ------------------------------------------------------------
// XdrHEAD members
XdrHEAD::XdrHEAD()
{
  m_wrtVar = 0;
  m_numvar = 0;

  m_meshCnt = 0;
  m_kstep = 0;

  m_numel = 0;
  m_numNodes = 0;
  m_sumWghts = 0;
  m_numBCs = 0;
  m_strSize = 0;
  mp_id = 0;
  mp_title = 0;
  mp_userTitle = 0;
  mp_varTitle = 0;

  m_time = 0;
}



XdrHEAD::~XdrHEAD()
{
  delete [] mp_id;
  delete [] mp_title;
  delete [] mp_userTitle;
  delete [] mp_varTitle;
}



char * XdrHEAD::cpyString(const char * src, int len)
{
  char * temp = libmesh_nullptr;
  int myLen = len;
  if(src)
    {
      if (myLen == -1)
        myLen = cast_int<int>(std::strlen(src))+1;
      temp = new char[myLen];
      temp = (char *) std::memcpy(temp, src, (myLen)*sizeof(char));
    }
  return temp;
}

} // namespace libMesh
