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



// C++ includes
#include <fstream>

// Local includes
#include "libmesh/plt_loader.h"

namespace libMesh
{



// extern "C" {

//   void open_   (const int *, const char *, int *);
//   void idata_  (const int *, int *, int *);
//   void fdata_  (const int *, float*, int *);
//   void ddata_  (const int *, double *, int *);
//   void close_  (const int *);
// }



//-----------------------------------------------------------------------------
// PltLoader write members
// void PltLoader::write_plot3d (const std::string & basename,
//       const bool reverse,
//       const bool gridonly) const
// {
//   const std::string gname = basename + ".g";
//   const std::string qname = basename + ".q";

//   // FORTAN file unit numbers
//   const int gunit = 25;
//   const int qunit = 26;

//   // Tell the user what files we are creating
//   if (this->verbose())
//     {
//       libMesh::out << "Plot3D output will be written to " << gname;

//       if (!gridonly)
// libMesh::out << " and " << qname;

//       libMesh::out << std::endl;
//     }


//   // Open the FORTRAN unformatted file
//   {
//     libmesh_assert_equal_to (gname.size(), qname.size());
//     int len = gname.size();

//     open_ (&gunit, gname.c_str(), &len);

//     if (!gridonly)
//       open_ (&qunit, qname.c_str(), &len);
//   }

//   // Write the headers
//   {
//     std::vector<int> ints;
//     ints.reserve (3*this->n_zones());

//     for (unsigned int zn=0; zn<this->n_zones(); zn++)
//       {
// ints.push_back(this->imax(zn));
// ints.push_back(this->jmax(zn));
// ints.push_back(this->kmax(zn));
//       }

//     int nb  = this->n_zones();
//     int one = 1;
//     int len = ints.size();

//     libmesh_assert_equal_to (static_cast<unsigned int>(len), 3*this->n_zones());

//     idata_ (&gunit, &nb, &one);
//     idata_ (&gunit, &ints[0], &len);

//     if (!gridonly)
//       {
// idata_ (&qunit, &nb, &one);
// idata_ (&qunit, &ints[0], &len);
//       }
//   }


//   // Variables to write to plot3D file
//   std::vector<unsigned int> write_vars;
//   write_vars.reserve (5);

//   std::fill (write_vars.begin(), write_vars.end(), 0);


//   //------------------------------------------------------------------------
//   // Ask the user which variables to write
//   if (!gridonly)
//     {
//       libMesh::out << "Variables:" << std::endl;

//       for (unsigned int v=0; v<this->n_vars(); v++)
// libMesh::out << " " << v << ") \"" << this->var_name(v) << "\""
//   << std::endl;
//       libMesh::out << std::endl;

//       int n_write_vars = 0;

//       while (true)
// {
//   libMesh::out << "How many variables to write to the Plot3D file? 1<=n<=" << this->n_vars()
//     << " "
//     << std::endl
//     << "(-1 writes them all): ";

//   std::cin >> n_write_vars;

//   if (n_write_vars == -1)
//     break;

//   if ((n_write_vars >= 1) &&
//       (n_write_vars <= static_cast<int>(this->n_vars())))
//     break;
// };


//       // The user wants all the variables
//       if ((n_write_vars == -1) ||
//   (n_write_vars == static_cast<int>(this->n_vars())))
// {
//   for (unsigned int wv=0; wv<this->n_vars(); wv++)
//     write_vars.push_back (wv);
// }

//       // The user wants a subset of the variables
//       else
// {
//   libmesh_assert_greater_equal (n_write_vars, 1);
//   libmesh_assert_less (n_write_vars, static_cast<int>(this->n_vars()));

//   libMesh::out << "Select the " << n_write_vars << " variables to write to the Plot3D file: "
//     << std::endl;

//   for (int wv=0; wv<n_write_vars; wv++)
//     {
//       int num=0;

//       std::cin >> num;

//       libmesh_assert_less (num, static_cast<int>(this->n_vars()));

//       write_vars.push_back (num);
//     }
// }

//       libMesh::out << std::endl;
//     } // if (!gridonly)



//   //------------------------------------------------------------------------
//   // Write the coordinates & data for each block
//   for (unsigned int zn=0; zn<this->n_zones(); zn++)
//     {
//       // Write the coordinates
//       {
// std::vector<float> coords;   // the nodal coordinates
// coords.reserve (3*this->imax(zn)*this->jmax(zn)*this->kmax(zn));

// for (unsigned int v=0; v<3; v++)
//   {
//     unsigned int l=0;

//     for (unsigned int k=0; k<this->kmax(zn); k++)
//       for (unsigned int j=0; j<this->jmax(zn); j++)
// for (unsigned int i=0; i<this->imax(zn); i++)
//   {
//     libmesh_assert_less (l, _data[zn][v].size());
//     coords.push_back (_data[zn][v][l++]);
//   }
//   }

// // Write to the grid file
// {
//   int len = coords.size();
//   libmesh_assert_equal_to (static_cast<unsigned int>(len),
//   3*this->imax(zn)*this->jmax(zn)*this->kmax(zn));

//   fdata_ (&gunit, &coords[0], &len);
// }
//       }


//       //------------------------------------------------------------------------
//       // Write the data
//       if (!gridonly)
// {
//   std::vector<float> data;     // arbitrary data
//   std::vector<float> conds(4); // plot3D conditions [FSMACH, ALPHA, RE, TIME]
//   data.reserve (write_vars.size()*this->imax(zn)*this->jmax(zn)*this->kmax(zn));
//   std::fill    (conds.begin(), conds.end(), 0.);

//   if (zn == 0)
//     libMesh::out << "  Writing ";

//   for (std::size_t i=0; i<write_vars.size(); i++)
//     {
//       // Number of the variable to write
//       const unsigned int v = write_vars[i];

//       libmesh_assert_less (v, this->n_vars());

//       // Tell the user what variable we are writing, but only
//       // once per file.
//       if (zn == 0)
// libMesh::out << "\"" << this->var_name(v) << "\" ";

//       unsigned int l=0;

//       for (unsigned int k=0; k<this->kmax(zn); k++)
// for (unsigned int j=0; j<this->jmax(zn); j++)
//   for (unsigned int i=0; i<this->imax(zn); i++)
//     {
//       libmesh_assert_less (l, _data[zn][v].size());
//       data.push_back ((v < this->n_vars()) ?
//       _data[zn][v][l++] : 0.);
//     }
//     }

//   if (zn == 0)
//     libMesh::out << "to " << qname << std::endl;

//   // Write to the solution file
//   {
//     int len = conds.size();

//     fdata_ (&qunit, &conds[0], &len);
//   }

//   // Write to the solution file
//   {
//     int len = data.size();
//     libmesh_assert_equal_to (static_cast<unsigned int>(len),
//     write_vars.size()*this->imax(zn)*this->jmax(zn)*this->kmax(zn));

//     fdata_ (&qunit, &data[0], &len);
//   }
// }
//     }

//   // Close the FORTAN files
//   close_ (&gunit);

//   if (!gridonly)
//     close_ (&qunit);

//   // Possibly reverse the orders
//   if (reverse)
//     {
//       if (this->verbose())
// libMesh::out << "Reversing byte-ordering for output files."
//   << std::endl;

//       Utility::reverse_endian (gname);

//       if (!gridonly)
// Utility::reverse_endian (qname);
//     }
// }



// void PltLoader::write_tri (const std::string & name,
//    const bool reverse,
//    const bool gridonly) const
//   {
//   // Check out
//   // http://people.nas.nasa.gov/~aftosmis/cart3d/cart3dTriangulations.html
//   // for the .tri, .triq format

//   // FORTRAN file unit numbers
//   const int gunit = 25;

//   if (this->verbose())
//     {
//       libMesh::out << "Writing unformatted .tri file " << name
// << std::endl;

//       if (gridonly)
// libMesh::out << "Only writing the grid to " << name
//   << std::endl;
//     }


//   // Open the FORTRAN unformatted file
//   {
//     int len = name.size();

//     open_ (&gunit, name.c_str(), &len);
//   }


//   // Write the header
//   unsigned int n_nodes =0;
//   unsigned int n_tri   =0;
//   unsigned int n_scalar=this->n_vars()-3;

//   {
//     std::vector<int> ints;

//     for (unsigned int zone=0; zone<this->n_zones(); zone++)
//       {
// libmesh_assert_equal_to (this->elem_type(zone), TRI);
// n_nodes += this->n_nodes(zone);
// n_tri   += this->n_elem(zone);
//       }

//     ints.push_back (n_nodes);
//     ints.push_back (n_tri);

//     if (!gridonly)
//       if (this->n_vars() > 3)
// ints.push_back(n_scalar);

//     int len = ints.size();
//     idata_ (&gunit, &ints[0], &len);
//   }

//   // Write the nodal values.
//   {
//     std::vector<float> coords;
//     coords.reserve (3*n_nodes);

//     for (unsigned int zone=0; zone<this->n_zones(); zone++)
//       for (unsigned int n=0; n<this->n_nodes(zone); n++)
// {
//   coords.push_back (_data[zone][0][n]);
//   coords.push_back (_data[zone][1][n]);
//   coords.push_back (_data[zone][2][n]);
// }
//     // Got all the nodes for all the zones


//     int len = coords.size();
//     fdata_ (&gunit, &coords[0], &len);
//   }

//   // Write the connectivity
//   {
//     std::vector<int> conn;
//     conn.reserve (3*n_tri);

//     for (unsigned int zone=0; zone<this->n_zones(); zone++)
//       {
// // The connectivity for this zone
// const std::vector<int> & zconn = _conn[zone];

// libmesh_assert (!zconn.empty());

// // Append the connectivity for this zone to the connectivity
// // array
// conn.insert (conn.end(), zconn.begin(), zconn.end());
//       }

//     int len = conn.size();
//     libmesh_assert_equal_to (static_cast<unsigned int>(len), 3*n_tri);
//     idata_ (&gunit, &conn[0], &len);
//   }


//   // Write the component index for each triangle
//   {
//     std::vector<int> comp;
//     comp.reserve (n_tri);

//     for (unsigned int zone=0; zone<this->n_zones(); zone++)
//       comp.insert (comp.end(), this->n_elem(zone), zone+1);

//     int len = comp.size();
//     libmesh_assert_equal_to (static_cast<unsigned int>(len), n_tri);
//     idata_ (&gunit, &comp[0], &len);
//   }


//   // Possibly write additional values for each node
//   if (!gridonly)
//     if (this->n_vars() > 3)
//       {
// if (this->verbose())
//   {
//     libMesh::out << "Writing variables ";

//     for (unsigned int v=3; v<this->n_vars(); v++)
//       libMesh::out << "\"" << this->var_name(v) << "\" ";

//     libMesh::out << "to the output file " << name
//       << std::endl;
//   }

// std::vector<float> data;

// data.reserve (n_nodes*(this->n_vars()-3));

// for (unsigned int zone=0; zone<this->n_zones(); zone++)
//   for (unsigned int n=0; n<this->n_nodes(zone); n++)
//     for (unsigned int v=3; v<this->n_vars(); v++)
//       data.push_back (_data[zone][v][n]);

// int len = data.size();
// libmesh_assert_equal_to (static_cast<unsigned int>(len),
// n_nodes*(this->n_vars()-3));
// fdata_ (&gunit, &data[0], &len);
//       }


//   // Close the FORTRAN file
//   close_ (&gunit);


//   // Possibly reverse the orders
//   if (reverse)
//     {
//       if (this->verbose())
// libMesh::out << "Reversing byte-ordering for output files."
//   << std::endl;

//       Utility::reverse_endian (name);
//     }
// }



void PltLoader::write_dat (const std::string & name,
                           const unsigned int version_in) const
{
  std::ofstream out_stream (name.c_str());

  out_stream << "TITLE=\""
             << this->title()
             << "\""
             << '\n';

  out_stream << "VARIABLES = ";

  for (unsigned int v=0; v<this->n_vars(); v++)
    out_stream << "\"" << this->var_name(v) << "\"\n";

  for (unsigned int z=0; z<this->n_zones(); z++)
    {
      out_stream << "ZONE T=\"" << this->zone_name(z) << "\"\n";
      out_stream << " I="  << this->imax(z)
                 << ", J=" << this->jmax(z)
                 << ", K=" << this->kmax(z);

      // Write BLOCK data for this zone
      if (this->zone_type(z) == BLOCK)
        {
          if (version_in < 10)
            {
              out_stream << ", F=BLOCK\n";
            }
          else
            {
              out_stream << ", ZONETYPE=Ordered\n"
                         << "DATAPACKING=BLOCK\n";
            }

          out_stream << "DT=(";
          for (unsigned int v=0; v<this->n_vars(); v++)
            out_stream << "SINGLE ";
          out_stream << ")\n";

          out_stream.precision(9);

          for (unsigned int v=0; v<this->n_vars(); v++)
            {
              unsigned int l=0;

              for (unsigned int k=0; k<this->kmax(z); k++)
                for (unsigned int j=0; j<this->jmax(z); j++)
                  for (unsigned int i=0; i<this->imax(z); i++)
                    {
                      // GCC 2.95.3 has scientific in the ios class instead
                      // of in namespace std::
#ifndef LIBMESH_BROKEN_IOSTREAM
                      out_stream << std::scientific
                                 << _data[z][v][l++] << " ";
#else
                      out_stream << std::ios::scientific
                                 << _data[z][v][l++] << " ";
#endif
                      // Throw in a newline every 5 entries to
                      // avoid really long lines.
                      if (l%5 == 0)
                        out_stream << '\n';
                    }

              if (l%5 != 0)
                out_stream << '\n';
            }
        } // end if (this->zone_type(z) == BLOCK)

      // Write POINT data for this zone
      else if (this->zone_type(z) == POINT)
        {
          if (version_in < 10)
            {
              out_stream << ", F=POINT\n";
            }
          else
            {
              out_stream << ", ZONETYPE=Ordered\n"
                         << "DATAPACKING=POINT\n";
            }

          out_stream << "DT=(";
          for (unsigned int v=0; v<this->n_vars(); v++)
            out_stream << "SINGLE ";
          out_stream << ")\n";

          out_stream.precision(9);

          {
            unsigned int l=0;

            for (unsigned int k=0; k<this->kmax(z); k++)
              for (unsigned int j=0; j<this->jmax(z); j++)
                for (unsigned int i=0; i<this->imax(z); i++)
                  {
                    for (unsigned int v=0; v<this->n_vars(); v++)

                      // GCC 2.95.3 has scientific in the ios class instead
                      // of in namespace std::
#ifndef LIBMESH_BROKEN_IOSTREAM
                      out_stream << std::scientific
                                 << _data[z][v][l] << " ";
#else
                    out_stream << std::ios::scientific
                               << _data[z][v][l] << " ";
#endif
                    out_stream << '\n';

                    l++;
                  }
          }
        } // end else if (this->zone_type(z) == POINT)

      // Otherwise, unrecognized zone type
      else
        libmesh_error_msg("Unrecognized zone type: this->zone_type(z)==" << this->zone_type(z));
    }
}

} // namespace libMesh
