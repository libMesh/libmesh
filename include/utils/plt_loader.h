// Copyright (C) 2002-2007  Benjamin S. Kirk

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



#ifndef LIBMESH_PLT_LOADER_H
#define LIBMESH_PLT_LOADER_H

// Local includes
#include "libmesh/libmesh_common.h"

// C++ includes
#include <string>
#include <vector>

namespace libMesh
{



/**
 * This class will read a binary \p .plt file.  These types of files
 * are for use with Amtec's <a href="http://www.tecplot.com">Tecplot</a>
 * visualization package.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 */
class PltLoader
{
public:

  /**
   * Constructor.  Initializes data.
   */
  PltLoader  (const bool v=false);

  /**
   * Constructor.  Reads the file specified by \p name.
   */
  PltLoader  (const std::string & name, const bool v=false);

  /**
   * Destructor.
   */
  ~PltLoader ();

  /**
   * Clear all data and return to a pristine state.
   */
  void clear ();

  /**
   * Returns the verbosity.
   */
  bool verbose () const { return _verbose; }

  /**
   * Reads the .plt file specified by \p name.
   */
  void read (const std::string & name);

  /**
   * Writes an ASCII Tecplot file.  The optional parameter \p version
   * specifies the version format to write.
   */
  void write_dat (const std::string & name,
                  const unsigned int version=10) const;

  // BSK - this functionality requires FORTRAN subrouitine calls,
  //       and there is no need to "dirty up" \p libMesh with FORTRAN
  //       just to enable these methods.
  //   /**
  //    * Writes a plot3d files.  The grid will be in basename.g and
  //    * the solution will be in basename.q.  It is assumed that the
  //    * first three variables from the .plt file are the (x,y,z)
  //    * locations of the grid points.  The optional parameter \p reverse
  //    * specifies if the output file will have reversed byte ordering.
  //    */
  //   void write_plot3d (const std::string & basename,
  //      const bool reverse=false,
  //      const bool gridonly=false) const;

  //   /**
  //    * Writes a Cart3D .tri component file.  The number of components
  //    * will be the number of zones in the .plt file.
  //    */
  //   void write_tri (const std::string & name,
  //   const bool reverse=false,
  //   const bool gridonly=false) const;



  //--------------------------------------------------------------
  // Data access



  /**
   * Enum defining the zone type in the Tecplot binary file,
   * for use with the old .plt format.
   */
  enum OldZoneType { BLOCK=0,
                     POINT,
                     FEBLOCK,
                     FEPOINT };

  /**
   * Enum defining the zone type in the Tecplot binary file,
   * for use with the new .plt format.
   */
  enum NewZoneType { ORDERED=0,
                     FELINESEG,
                     FETRIANGLE,
                     FEQUADRILATERAL,
                     FETETRAHEDRON,
                     FEBRICK };

  /**
   * Enum defining the data type of each variable.
   */
  enum DataType { FLOAT=1,
                  DOUBLE,
                  LONGINT,
                  SHORTINT,
                  BYTE,
                  BIT};

  /**
   * Enum defining the finite element types
   */
  enum FEType { TRI=0,
                QUAD,
                TET,
                HEX };



  //--------------------------------------------------------------
  // Public Data Access



  /**
   * @returns the Tecplot version number string.  This identifies the
   * version of Tecplot (or preplot) that wrote the binary file.  Currently,
   * PltLoader understands versions "#!TDV7X " and "#!TDV1XX"
   */
  const std::string & version () const { return _version; }

  /**
   * @returns \p true if the binary type of the file is different than the
   * machine that is reading it.  If this is the case we must perform an
   * endian-swap on all input data.
   */
  bool is_foreign () const { return _is_foreign; }

  /**
   * @returns the data set title
   */
  const std::string & title () const { return _title; }

  /**
   * @returns the number of variables in the data set.
   */
  unsigned int n_vars () const { return _n_vars; }

  /**
   * @returns the name of variable \p v.
   */
  const std::string & var_name (const unsigned int v) const;

  /**
   * @returns the type of variable \p v
   */
  unsigned int var_type (const unsigned int v) const;

  /**
   * @returns the number of zones.
   */
  unsigned int n_zones () const { return _n_zones; }

  /**
   * @returns the type of zone \p z
   */
  unsigned int zone_type (const unsigned int z) const;

  /**
   * @returns the name of zone \p z.
   */
  const std::string & zone_name (const unsigned int z) const;

  /**
   * @returns the data packing flag for zone \p z.
   */
  unsigned int zone_pack (const unsigned int z) const;

  /**
   * @returns \p imax for zone \p z.
   */
  unsigned int imax (const unsigned int z) const;

  /**
   * @returns \p jmax for zone \p z.
   */
  unsigned int jmax (const unsigned int z) const;

  /**
   * @returns \p kmax for zone \p z.
   */
  unsigned int kmax (const unsigned int z) const;

  /**
   * @returns the number of nodes in the mesh (for unstructured meshes).
   */
  unsigned int n_nodes (const unsigned int z) const;

  /**
   * @returns the number of elements in the mesh (for unstructured meshes).
   */
  unsigned int n_elem (const unsigned int z) const;

  /**
   * @returns the element type for the \p zth zone (for unstructured meshes).
   */
  FEType elem_type (const unsigned int z) const;

  /**
   * @returns a reference to the data read from the file
   */
  const std::vector<std::vector<std::vector<float> > > & get_data () const;

  /**
   * Enum defining the number of nodes for each element type.
   */
  static const unsigned int NNodes[4];


private:


  /**
   * Read the header of the binary file.
   */
  void read_header (std::istream & in);

  /**
   * Read data from the binary file.
   */
  void read_data (std::istream & in);

  /**
   * Read data for the zth zone in BLOCK structured format.
   */
  void read_block_data (std::istream & in, const unsigned int zn);

  /**
   * Read data for the zth zone in POINT structured format.
   */
  void read_point_data (std::istream & in, const unsigned int zn);

  /**
   * Read data for the zth zone in FEBLOCK unstructured format.
   */
  void read_feblock_data (std::istream & in, const unsigned int zn);

  /**
   * Read data for the zth zone in FEPOINT unstructured format.
   */
  void read_fepoint_data (std::istream & in, const unsigned int zn);


  //--------------------------------------------------------------
  // Private Data Access



  /**
   * @returns the Tecplot version number string.
   */
  std::string & version () { return _version; }

  /**
   * @returns \p true if the binary type of the file is different than the
   * machine that is reading it.  If this is the case we must perform an
   * endian-swap on all input data.
   */
  bool & is_foreign () { return _is_foreign; }

  /**
   * @returns the data set title
   */
  std::string & title () { return _title; }

  /**
   * @returns the number of variables in the data set.
   */
  void set_n_vars (const unsigned int nv);

  /**
   * @returns the name of variable \p v.
   */
  std::string & var_name (const unsigned int v);

  /**
   * @returns the type of variable \p v
   */
  unsigned int & var_type (const unsigned int v);

  /**
   * @returns the number of zones.
   */
  void set_n_zones (const unsigned int nz);

  /**
   * @returns the type of zone \p z
   */
  unsigned int & zone_type (const unsigned int z);

  /**
   * @returns the name of zone \p z.
   */
  std::string & zone_name (const unsigned int z);

  /**
   * @returns the data pack flag for zone \p z.
   */
  unsigned int & zone_pack (const unsigned int z);

  /**
   * @returns \p imax for zone \p z.
   */
  unsigned int & imax (const unsigned int z);

  /**
   * @returns \p jmax for zone \p z.
   */
  unsigned int & jmax (const unsigned int z);

  /**
   * @returns \p kmax for zone \p z.
   */
  unsigned int & kmax (const unsigned int z);


  //--------------------------------------------------------------
  // Private Data


  /**
   * Verbosity
   */
  const bool _verbose;

  /**
   * The Tecplot Version number string.
   */
  std::string _version;

  /**
   * Is the data foreign?
   */
  bool _is_foreign;

  /**
   * The Tecplot data set title.
   */
  std::string _title;

  /**
   * The number of variables in the data set.
   */
  unsigned int _n_vars;

  /**
   * The name for each variable.
   */
  std::vector<std::string> _var_names;

  /**
   * The type of each variable.  Must be one of the
   * enumerated \p DataType types.
   */
  std::vector<unsigned int> _var_types;

  /**
   * The number of zones.
   */
  unsigned int _n_zones;

  /**
   * The type of each zone.
   */
  std::vector<unsigned int> _zone_types;

  /**
   * The name of each zone.
   */
  std::vector<std::string> _zone_names;

  /**
   * The data packing for each zone (new version only)
   */
  std::vector<unsigned int> _zone_pack;

  /**
   * The (imax,jmax,kmax) value for each zone.
   */
  std::vector<unsigned int> _imax;
  std::vector<unsigned int> _jmax;
  std::vector<unsigned int> _kmax;

  /**
   * Vector to hold the data.
   */
  std::vector<std::vector<std::vector<float> > >  _data;

  /**
   * Vectors to hold the connectivity for each zone
   * (only for unstructured files).
   */
  std::vector<std::vector<int> > _conn;

  /**
   * Scratch data & relevant sizes.
   */
  mutable char buf[512];
};



//---------------------------------------------------------
// PltLoader inline members
inline
PltLoader::PltLoader (const bool v) :
  _verbose      (v),
  _is_foreign   (false),
  _n_vars       (0),
  _n_zones      (0)
{
}



inline
PltLoader::PltLoader (const std::string & name, const bool v) :
  _verbose      (v),
  _is_foreign   (false),
  _n_vars       (0),
  _n_zones      (0)
{
  this->read (name);
}



inline
PltLoader::~PltLoader()
{
}



inline
const std::string & PltLoader::var_name (const unsigned int v) const
{
  libmesh_assert_less (v, this->n_vars());
  libmesh_assert_less (v, _var_names.size());
  libmesh_assert_equal_to (this->n_vars(), _var_names.size());

  return _var_names[v];
}



inline
std::string & PltLoader::var_name (const unsigned int v)
{
  libmesh_assert_less (v, this->n_vars());
  libmesh_assert_less (v, _var_names.size());
  libmesh_assert_equal_to (this->n_vars(), _var_names.size());

  return _var_names[v];
}



inline
unsigned int PltLoader::var_type (const unsigned int v) const
{
  libmesh_assert_less (v, this->n_vars());
  libmesh_assert_less (v, _var_types.size());
  libmesh_assert_equal_to (this->n_vars(), _var_types.size());

  return _var_types[v];
}



inline
unsigned int & PltLoader::var_type (const unsigned int v)
{
  libmesh_assert_less (v, this->n_vars());
  libmesh_assert_less (v, _var_types.size());
  libmesh_assert_equal_to (this->n_vars(), _var_types.size());

  return _var_types[v];
}



inline
unsigned int PltLoader::zone_type (const unsigned int z) const
{
  libmesh_assert_less (z, this->n_zones());
  libmesh_assert_less (z, _zone_types.size());
  libmesh_assert_equal_to (this->n_zones(), _zone_types.size());

  return _zone_types[z];
}



inline
unsigned int & PltLoader::zone_type (const unsigned int z)
{
  libmesh_assert_less (z, this->n_zones());
  libmesh_assert_less (z, _zone_types.size());
  libmesh_assert_equal_to (this->n_zones(), _zone_types.size());

  return _zone_types[z];
}



inline
const std::string & PltLoader::zone_name (const unsigned int z) const
{
  libmesh_assert_less (z, this->n_zones());
  libmesh_assert_less (z, _zone_names.size());
  libmesh_assert_equal_to (this->n_zones(), _zone_names.size());

  return _zone_names[z];
}



inline
std::string & PltLoader::zone_name (const unsigned int z)
{
  libmesh_assert_less (z, this->n_zones());
  libmesh_assert_less (z, _zone_names.size());
  libmesh_assert_equal_to (this->n_zones(), _zone_names.size());

  return _zone_names[z];
}



inline
unsigned int PltLoader::zone_pack (const unsigned int z) const
{
  libmesh_assert_less (z, this->n_zones());
  libmesh_assert_less (z, _zone_pack.size());
  libmesh_assert_equal_to (this->n_zones(), _zone_pack.size());

  return _zone_pack[z];
}



inline
unsigned int & PltLoader::zone_pack (const unsigned int z)
{
  libmesh_assert_less (z, this->n_zones());
  libmesh_assert_less (z, _zone_pack.size());
  libmesh_assert_equal_to (this->n_zones(), _zone_pack.size());

  return _zone_pack[z];
}



inline
unsigned int PltLoader::imax (const unsigned int z) const
{
  libmesh_assert_less (z, this->n_zones());
  libmesh_assert_equal_to (_imax.size(), this->n_zones());

  return _imax[z];
}



inline
unsigned int & PltLoader::imax (const unsigned int z)
{
  libmesh_assert_less (z, this->n_zones());
  libmesh_assert_equal_to (_imax.size(), this->n_zones());

  return _imax[z];
}



inline
unsigned int PltLoader::jmax (const unsigned int z) const
{
  libmesh_assert_less (z, this->n_zones());
  libmesh_assert_equal_to (_jmax.size(), this->n_zones());

  return _jmax[z];
}



inline
unsigned int & PltLoader::jmax (const unsigned int z)
{
  libmesh_assert_less (z, this->n_zones());
  libmesh_assert_equal_to (_jmax.size(), this->n_zones());

  return _jmax[z];
}



inline
unsigned int PltLoader::kmax (const unsigned int z) const
{
  libmesh_assert_less (z, this->n_zones());
  libmesh_assert_equal_to (_kmax.size(), this->n_zones());

  return _kmax[z];
}



inline
unsigned int & PltLoader::kmax (const unsigned int z)
{
  libmesh_assert_less (z, this->n_zones());
  libmesh_assert_equal_to (_kmax.size(), this->n_zones());

  return _kmax[z];
}



inline
unsigned int PltLoader::n_nodes (const unsigned int z) const
{
  libmesh_assert_less (z, this->n_zones());

  // Only for unstructured zones!
  libmesh_assert_greater (this->zone_type(z), 1);

  return this->imax(z);
}



inline
unsigned int PltLoader::n_elem (const unsigned int z) const
{
  libmesh_assert_less (z, this->n_zones());

  // Only for unstructured zones!
  libmesh_assert_greater (this->zone_type(z), 1);

  return this->jmax(z);
}



inline
PltLoader::FEType PltLoader::elem_type (const unsigned int z) const
{
  libmesh_assert_less (z, this->n_zones());

  // Only for unstructured zones!
  libmesh_assert_greater (this->zone_type(z), 1);

  return static_cast<FEType>(this->kmax(z));
}


inline
const std::vector<std::vector<std::vector<float> > > &
PltLoader::get_data () const
{
  return _data;
}




} // namespace libMesh


#endif // LIBMESH_PLT_LOADER_H
