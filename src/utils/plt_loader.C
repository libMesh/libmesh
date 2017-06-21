// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA



#include "libmesh/plt_loader.h"

namespace libMesh
{



//---------------------------------------------------------
// PltLoader static data
const unsigned int PltLoader::NNodes[4] = {3, 4, 4, 8};



//-----------------------------------------------------------------------------
// PltLoader members
void PltLoader::clear ()
{
  // clear vectors & strings.  Using .erase() for strings instead of .clear()
  // since GCC 2.95.3 does not support .clear().
  _version.erase();
  _title.erase();

  _var_names.clear();
  _var_types.clear();
  _zone_types.clear();
  _zone_names.clear();
  _zone_pack.clear();
  _imax.clear();
  _jmax.clear();
  _kmax.clear();
  _data.clear();
  _conn.clear();

  // reinitialize
  _is_foreign = false;
  _n_vars     = 0;
  _n_zones    = 0;
}



void PltLoader::set_n_vars (const unsigned int nv)
{
  _n_vars = nv;

  _var_types.resize (this->n_vars());
  _var_names.resize (this->n_vars());

  // Default to float data
  std::fill (_var_types.begin(), _var_types.end(), 1);

  // If the number of zones is set, resize the data.
  if (this->n_zones())
    {
      _data.resize  (this->n_zones());

      for (unsigned int z=0; z<this->n_zones(); z++)
        _data[z].resize  (this->n_vars());
    }
}



void PltLoader::set_n_zones (const unsigned int nz)
{
  _n_zones = nz;

  _zone_types.resize (this->n_zones());
  _zone_names.resize (this->n_zones());
  _zone_pack.resize  (this->n_zones());

  _imax.resize (this->n_zones());
  _jmax.resize (this->n_zones());
  _kmax.resize (this->n_zones());

  _data.resize (this->n_zones());
  _conn.resize        (this->n_zones());

  // If the number of variables are set, resize the data.
  if (this->n_vars())
    for (unsigned int z=0; z<this->n_zones(); z++)
      _data[z].resize (this->n_vars());
}

} // namespace libMesh
