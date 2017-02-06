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
#include <iostream>
#include <fstream>
#include <cstring>

// Local includes
#include "libmesh/utility.h"
#include "libmesh/plt_loader.h"

namespace libMesh
{



//-----------------------------------------------------------------------------
// PltLoader reading members
void PltLoader::read (const std::string & name)
{
  std::ifstream in (name.c_str(), std::ios::in|std::ios::binary);

  if (!in.good())
    libmesh_error_msg("Error reading input file " << name);


  if (this->verbose())
    libMesh::out << std::endl
                 << "Reading input file " << name
                 << std::endl
                 << "-------------------------------------------------------------------------"
                 << std::endl;

  this->read_header (in);
  this->read_data   (in);

  if (this->verbose())
    libMesh::out << std::endl
                 << "-------------------------------------------------------------------------"
                 << std::endl;

}



void PltLoader::read_header (std::istream & in)
{
  libmesh_assert (in.good());

  //----------------------------------------------------
  // Read the TECPLOT header

  // Read the version number
  {
    in.read (buf, 8);

    // Using erase for GCC 2.95.3
    this->version().erase();

    for (unsigned int i=0; i<8; i++)
      this->version() += buf[i];

    if (this->verbose())
      libMesh::out << "Tecplot Version: "
                   << this->version()
                   << std::endl;
  }


  //----------------------------------------------------
  // Read plt files written by older versions of Tecplot
  if (this->version().rfind("V7") < this->version().size())
    {
      if (this->verbose())
        libMesh::out << "Reading legacy .plt format (<= v9) ..."
                     << std::endl;

      // Read the value of 1 to determine byte ordering
      {
        int one = 0;
        in.read (buf, LIBMESH_SIZEOF_INT);
        std::memcpy  (&one, buf, LIBMESH_SIZEOF_INT);

        if (one != 1)
          {
            if (this->verbose())
              libMesh::out << "Tecplot data is Foreign!"
                           << std::endl;
            this->is_foreign() = true;

            // Make sure one reversed is one
            Utility::ReverseBytes rb(this->is_foreign());
            libmesh_assert_equal_to (rb(one), 1);
          }
      }

      // A byte-reverser in case the data is foreign
      Utility::ReverseBytes rb(this->is_foreign());

      // Read the title
      {
        int i=0;

        // Using erase for GCC 2.95.3
        this->title().erase();

        do
          {
            in.read (buf, LIBMESH_SIZEOF_INT);
            std::memcpy  (&i, buf, LIBMESH_SIZEOF_INT);
            rb(i);

            // Don't add trailing \0
            if (i)
              this->title() += static_cast<char>(i);
          }
        while (i);
      }

      // Read the number of variables in the data set
      {
        int nv;
        in.read (buf, LIBMESH_SIZEOF_INT);
        std::memcpy  (&nv, buf, LIBMESH_SIZEOF_INT);
        rb(nv);

        this->set_n_vars (nv);
      }

      // Read the variable names
      for (unsigned int v=0; v<this->n_vars(); v++)
        {
          int i=0;

          // Using erase for GCC 2.95.3
          this->var_name(v).erase();

          do
            {
              in.read (buf, LIBMESH_SIZEOF_INT);
              std::memcpy  (&i, buf, LIBMESH_SIZEOF_INT);
              rb(i);

              // Don't add trailing \0
              if (i)
                this->var_name(v) += static_cast<char>(i);
            }
          while (i);
        }



      // Read zones from the header.
      // Continue reading until the end-of-header
      // marker (357.) is found.
      int nz=0;
      std::vector<std::string> zname;
      std::vector<int>         ztype, zimax, zjmax, zkmax;

      {
        float f=0.;

        do
          {
            // find the next Zone marker
            do
              {
                f = 0.;
                in.read (buf, LIBMESH_SIZEOF_FLOAT);
                std::memcpy  (&f, buf, LIBMESH_SIZEOF_FLOAT);
                rb(f);
              }
            while ((f != 299.) &&
                   (f != 357.) &&
                   in.good());


            // Did we overrun the file?
            if (!in.good())
              libmesh_error_msg("ERROR: Unexpected end-of-file!");

            // Found a Zone marker
            else if (f == 299.)
              {
                // Incriment the Zone counter
                nz++;

                // Read the zone name
                {
                  int i=0;
                  std::string name;

                  do
                    {
                      in.read (buf, LIBMESH_SIZEOF_INT);
                      std::memcpy  (&i, buf, LIBMESH_SIZEOF_INT);
                      rb(i);

                      // Don't add trailing \0
                      if (i)
                        name += static_cast<char>(i);
                    }
                  while (i);

                  zname.push_back(name);
                }

                // Read the zone format
                {
                  int zt;
                  in.read (buf, LIBMESH_SIZEOF_INT);
                  std::memcpy  (&zt, buf, LIBMESH_SIZEOF_INT);
                  rb(zt);

                  ztype.push_back(zt);
                  //libMesh::out << "zone type=" << ztype.back() << std::endl;
                }

                // Read the zone color
                {
                  int zc=0;

                  in.read (buf, LIBMESH_SIZEOF_INT);
                  std::memcpy  (&zc, buf, LIBMESH_SIZEOF_INT);
                  rb(zc);

                  //libMesh::out << "zone color=" << zc << std::endl;
                }

                // Read in the block dimensions
                {
                  int
                    i_max=0,
                    j_max=0,
                    k_max=0;

                  in.read (buf, LIBMESH_SIZEOF_INT);
                  std::memcpy  (&i_max, buf, LIBMESH_SIZEOF_INT);
                  rb(i_max);

                  in.read (buf, LIBMESH_SIZEOF_INT);
                  std::memcpy  (&j_max, buf, LIBMESH_SIZEOF_INT);
                  rb(j_max);

                  in.read (buf, LIBMESH_SIZEOF_INT);
                  std::memcpy  (&k_max, buf, LIBMESH_SIZEOF_INT);
                  rb(k_max);

                  zimax.push_back (i_max);
                  zjmax.push_back (j_max);
                  zkmax.push_back (k_max);
                }
              } // else if (f == 299.)
          }
        while ((f != 357.) && in.good());
      }

      // Set the header data
      this->set_n_zones (nz);

      for (unsigned int z=0; z<this->n_zones(); z++)
        {
          this->zone_type(z) = ztype[z];
          this->zone_name(z) = zname[z];
          this->imax(z)      = zimax[z];
          this->jmax(z)      = zjmax[z];
          this->kmax(z)      = zkmax[z];
        }
    }


  //----------------------------------------------------
  // Read plt files written by newer versions of Tecplot
  else if (this->version().rfind("V1") < this->version().size())
    {
      if (this->verbose())
        libMesh::out << "Reading new .plt format (>= v10)..."
                     << std::endl;

      // Read the value of 1 to determine byte ordering
      {
        int one = 0;

        in.read (buf, LIBMESH_SIZEOF_INT);
        std::memcpy  (&one, buf, LIBMESH_SIZEOF_INT);

        if (one != 1)
          {
            if (this->verbose())
              libMesh::err << "Tecplot data is Foreign!"
                           << std::endl;
            this->is_foreign() = true;

            // Make sure one reversed is one
            Utility::ReverseBytes rb(this->is_foreign());
            libmesh_assert_equal_to (rb(one), 1);
          }
      }

      // A byte-reverser in case the data is foreign
      Utility::ReverseBytes rb(this->is_foreign());

      // Read the title
      {
        int i=0;

        // Using erase() for GCC 2.95.3
        this->title().erase();
        do
          {
            in.read (buf, LIBMESH_SIZEOF_INT);
            std::memcpy  (&i, buf, LIBMESH_SIZEOF_INT);
            rb(i);

            // Don't add trailing \0
            if (i)
              this->title() += static_cast<char>(i);
          }
        while (i);
      }

      // Read the number of variables in the data set
      {
        int nv;
        in.read (buf, LIBMESH_SIZEOF_INT);
        std::memcpy  (&nv, buf, LIBMESH_SIZEOF_INT);
        rb(nv);

        this->set_n_vars (nv);
      }

      // Read the variable names
      for (unsigned int v=0; v<this->n_vars(); v++)
        {
          int i=0;

          // Using erase() for GCC 2.95.3
          this->var_name(v).erase();

          do
            {
              in.read (buf, LIBMESH_SIZEOF_INT);
              std::memcpy  (&i, buf, LIBMESH_SIZEOF_INT);
              rb(i);

              // Don't add trailing \0
              if (i)
                this->var_name(v) += static_cast<char>(i);
            }
          while (i);
        }



      // Read zones from the header.
      // Continue reading until the end-of-header
      // marker (357.) is found.
      int nz=0;
      std::vector<std::string> zname;
      std::vector<int>         zpack, ztype, zimax, zjmax, zkmax;

      {
        float f=0.;

        do
          {
            // find the next Zone marker
            do
              {
                f = 0.;
                in.read (buf, LIBMESH_SIZEOF_FLOAT);
                std::memcpy  (&f, buf, LIBMESH_SIZEOF_FLOAT);
                rb(f);
              }
            while ((f != 299.) &&
                   (f != 357.) &&
                   in.good());


            // Did we overrun the file?
            if (!in.good())
              libmesh_error_msg("ERROR: Unexpected end-of-file!");

            // Found a Zone marker
            else if (f == 299.)
              {
                // Incriment the Zone counter
                nz++;

                // Read the zone name
                {
                  int i=0;
                  std::string name;

                  do
                    {
                      in.read (buf, LIBMESH_SIZEOF_INT);
                      std::memcpy  (&i, buf, LIBMESH_SIZEOF_INT);
                      rb(i);

                      // Don't add trailing \0
                      if (i)
                        name += static_cast<char>(i);
                    }
                  while (i);

                  zname.push_back(name);
                }

                // Read the zone color
                {
                  int zc=0;
                  in.read (buf, LIBMESH_SIZEOF_INT);
                  std::memcpy  (&zc, buf, LIBMESH_SIZEOF_INT);
                  rb(zc);
                }

                // Read the zone format
                {
                  int zt;
                  in.read (buf, LIBMESH_SIZEOF_INT);
                  std::memcpy  (&zt, buf, LIBMESH_SIZEOF_INT);
                  rb(zt);

                  ztype.push_back(zt);
                }

                // Read the data packing flag
                {
                  int dp=0;
                  in.read (buf, LIBMESH_SIZEOF_INT);
                  std::memcpy (&dp, buf, LIBMESH_SIZEOF_INT);
                  rb(dp);

                  zpack.push_back (dp);
                }

                // Will we specify the variable location?
                {
                  int svl=0;
                  int  vl=0;
                  in.read (buf, LIBMESH_SIZEOF_INT);
                  std::memcpy  (&svl, buf, LIBMESH_SIZEOF_INT);
                  rb(svl);

                  if (svl)
                    for (unsigned int v=0; v<this->n_vars(); v++)
                      {
                        in.read (buf, LIBMESH_SIZEOF_INT);
                        std::memcpy  (&vl, buf, LIBMESH_SIZEOF_INT);
                        rb(vl);
                        libmesh_assert_equal_to (vl, 0); // Only know about node-based data
                        // right now
                      }

                }

                // Get the number of user-defined face-neighbors
                {
                  int fn=0;
                  in.read (buf, LIBMESH_SIZEOF_INT);
                  std::memcpy  (&fn, buf, LIBMESH_SIZEOF_INT);
                  rb(fn);
                }

                // Read in the block dimensions
                {
                  if (ztype.back() != ORDERED)
                    {
                      int np=0, ne=0;

                      in.read (buf, LIBMESH_SIZEOF_INT);
                      std::memcpy  (&np, buf, LIBMESH_SIZEOF_INT);
                      rb(np);

                      in.read (buf, LIBMESH_SIZEOF_INT);
                      std::memcpy  (&ne, buf, LIBMESH_SIZEOF_INT);
                      rb(ne);

                      zimax.push_back (np);
                      zjmax.push_back (ne);
                      zjmax.push_back (0);
                    }

                  int
                    i_max=0,
                    j_max=0,
                    k_max=0;

                  in.read (buf, LIBMESH_SIZEOF_INT);
                  std::memcpy  (&i_max, buf, LIBMESH_SIZEOF_INT);
                  rb(i_max);

                  in.read (buf, LIBMESH_SIZEOF_INT);
                  std::memcpy  (&j_max, buf, LIBMESH_SIZEOF_INT);
                  rb(j_max);

                  in.read (buf, LIBMESH_SIZEOF_INT);
                  std::memcpy  (&k_max, buf, LIBMESH_SIZEOF_INT);
                  rb(k_max);

                  // These are only useful for orderd data.  Otherwise
                  // we grabbed the relevant values above.
                  if (ztype.back() != ORDERED)
                    {
                      zimax.push_back (i_max);
                      zjmax.push_back (j_max);
                      zkmax.push_back (k_max);
                    }
                }
              } // else if (f == 299.)
          }
        while ((f != 357.) && in.good());
      }

      // Set the header data
      this->set_n_zones (nz);

      for (unsigned int z=0; z<this->n_zones(); z++)
        {
          this->zone_type(z) = ztype[z];
          this->zone_name(z) = zname[z];
          this->zone_pack(z) = zpack[z];
          this->imax(z)      = zimax[z];
          this->jmax(z)      = zjmax[z];
          this->kmax(z)      = zkmax[z];
        }
    }



  //----------------------------------------------------
  // Unrecognized Tecplot Version!
  else
    libmesh_error_msg("ERROR:  This plot file was written by an unrecognized version of Tecplot!:\n" << this->version());








  // Print the data to the screen.
  if (this->verbose())
    {
      libMesh::out << "Tecplot Header: "
                   << this->title() << std::endl;

      libMesh::out << "Variables: ";
      for (unsigned int v=0; v<this->n_vars(); v++)
        libMesh::out << "\"" << this->var_name (v) << "\"" << " ";
      libMesh::out << std::endl;

      libMesh::out << "Variable Types: ";
      for (unsigned int v=0; v<this->n_vars(); v++)
        libMesh::out << this->var_type (v) << " ";
      libMesh::out << std::endl;

      libMesh::out << "Zone Names: ";
      for (unsigned int z=0; z<this->n_zones(); z++)
        libMesh::out << "\"" << this->zone_name (z) << "\"" << " ";
      libMesh::out << std::endl;

      libMesh::out << "Zone Types: ";
      for (unsigned int z=0; z<this->n_zones(); z++)
        {
          libMesh::out << this->zone_type (z) << " ";

          if (this->zone_type (z) != ORDERED)
            libMesh::out << "(" << this->n_nodes(z) << "," << this->n_elem(z) << ") ";
        }
      libMesh::out << std::endl;

      libMesh::out << "Zone Dimensions: " << std::endl;
      for (unsigned int z=0; z<this->n_zones(); z++)
        libMesh::out << "  ("
                     << this->imax(z) << ","
                     << this->jmax(z) << ","
                     << this->kmax(z) << ")"
                     << std::endl;
    }
}



void PltLoader::read_data (std::istream & in)
{
  libmesh_assert (in.good());

  // A byte-reverser in case the data is foreign
  Utility::ReverseBytes rb(this->is_foreign());

  //----------------------------------------------------
  // Read the TECPLOT data for each zone
  if (this->verbose())
    {
      libMesh::out << "Reading Zones";
      libMesh::out.flush();
    }


  for (unsigned int zone=0; zone<this->n_zones(); zone++)
    {
      if (this->verbose())
        {
          libMesh::out << ".";
          libMesh::out.flush();
        }


      //----------------------------------------------------
      // Read plt files written by older versions of Tecplot
      if (this->version().rfind("V7") < this->version().size())
        {
          float f = 0.;

          // Find the next Zone marker.
          do
            {
              f = 0.;
              in.read (buf, LIBMESH_SIZEOF_FLOAT);
              std::memcpy  (&f, buf, LIBMESH_SIZEOF_FLOAT);
              rb(f);
            }
          while ((f != 299.) && in.good());

          // Did we overrun the file?
          if (!in.good())
            libmesh_error_msg("ERROR: Unexpected end-of-file!");

          // Get the number of repeated vars.
          unsigned int n_rep_vars=0;
          std::vector<int> rep_vars;

          {
            in.read (buf, LIBMESH_SIZEOF_INT);
            std::memcpy  (&n_rep_vars, buf, LIBMESH_SIZEOF_INT);
            rb(n_rep_vars);

            rep_vars.resize (n_rep_vars);

            // Get the repeated variables number.
            for (unsigned int v=0; v<n_rep_vars; v++)
              {
                libmesh_error_msg("ERROR:  I don't understand repeated variables yet!");

                in.read (buf, LIBMESH_SIZEOF_INT);
                std::memcpy  (&rep_vars[v], buf, LIBMESH_SIZEOF_INT);
                rb(rep_vars[v]);
              }
          }

          // Get the variable data type
          //libMesh::out << "var_types=";
          for (unsigned int v=0; v<this->n_vars(); v++)
            {
              in.read (buf, LIBMESH_SIZEOF_INT);
              std::memcpy  (&this->var_type(v), buf, LIBMESH_SIZEOF_INT);
              rb(this->var_type(v));

              //libMesh::out << this->var_type(v) << " ";
            }
          //libMesh::out << std::endl;



          // Read the data.
          switch (this->zone_type(zone) )
            {
              // Block-based data.  Structured meshes.
            case BLOCK:
              {
                this->read_block_data (in, zone);
                break;
              }

              // Point-based data.  Structured meshes.
            case POINT:
              {
                this->read_point_data (in, zone);
                break;
              }

              // FE block data.  Unstructured meshes.
            case FEBLOCK:
              {
                this->read_feblock_data (in, zone);

                if (this->verbose())

                  libMesh::out << "Zone " << zone << ":" << std::endl
                               << "  nnodes   =" << this->imax(zone) << std::endl
                               << "  nelem    =" << this->jmax(zone) << std::endl
                               << "  elem_type=" << this->kmax(zone) << std::endl
                               << std::endl;
                break;
              }

              // FE point data.  Unstructured meshes.
            case FEPOINT:
              {
                this->read_fepoint_data (in, zone);
                break;
              }

            default:
              libmesh_error_msg("ERROR: Unsupported Zone type: " << this->zone_type(zone));
            } // end switch on zone type
        }


      //----------------------------------------------------
      // Read plt files written by newer versions of Tecplot
      else if (this->version().rfind("V1") < this->version().size())
        {
          float f = 0.;

          // Find the next Zone marker.
          do
            {
              f = 0.;
              in.read (buf, LIBMESH_SIZEOF_FLOAT);
              std::memcpy  (&f, buf, LIBMESH_SIZEOF_FLOAT);
              rb(f);
            }
          while ((f != 299.) && in.good());

          // Did we overrun the file?
          if (!in.good())
            libmesh_error_msg("ERROR: Unexpected end-of-file!");

          // Get the variable data type
          for (unsigned int v=0; v<this->n_vars(); v++)
            {
              in.read (buf, LIBMESH_SIZEOF_INT);
              std::memcpy  (&this->var_type(v), buf, LIBMESH_SIZEOF_INT);
              rb(this->var_type(v));

              //libMesh::out << this->var_type(v) << " ";
            }

          // Get the variable sharing flag
          {
            int vs=0;
            int sv=0;

            in.read (buf, LIBMESH_SIZEOF_INT);
            std::memcpy  (&vs, buf, LIBMESH_SIZEOF_INT);
            rb(vs);

            if (vs)
              {
                for (unsigned int v=0; v<this->n_vars(); v++)
                  {
                    in.read (buf, LIBMESH_SIZEOF_INT);
                    std::memcpy  (&sv, buf, LIBMESH_SIZEOF_INT);
                    rb(sv);

                    if (sv != -1)
                      libmesh_error_msg("ERROR:  I don't understand variable sharing!");
                  }
              }
          }

          // Get zone to share connectivity with
          {
            int sc=0;
            in.read (buf, LIBMESH_SIZEOF_INT);
            std::memcpy  (&sc, buf, LIBMESH_SIZEOF_INT);
            rb(sc);

            libmesh_assert_equal_to (sc, -1);
          }


          // Read the data.
          if (this->zone_type(zone) == ORDERED)
            {
              // Block-based data.  Structured meshes.
              if (this->zone_pack(zone) == 0)
                this->read_block_data (in, zone);

              // Point-based data.  Structured meshes.
              else if (this->zone_pack(zone) == 1)
                this->read_point_data (in, zone);

              else
                libmesh_error_msg("Unrecognized zone_pack(zone) = " << this->zone_pack(zone));
            }
          else
            {
              // Block-based data.  Unstructured meshes.
              if (this->zone_pack(zone) == 0)
                this->read_feblock_data (in, zone);

              // Point-based data.  Unstructured meshes.
              else if (this->zone_pack(zone) == 1)
                this->read_fepoint_data (in, zone);

              else
                libmesh_error_msg("Unrecognized zone_pack(zone) = " << this->zone_pack(zone));
            }
        }



      //----------------------------------------------------
      // Unrecognized Tecplot Version!
      else
        libmesh_error_msg("ERROR:  This plot file was written by an unrecognized version of Tecplot!:\n" << this->version());

    } // end loop on zones
}



void PltLoader::read_block_data (std::istream & in, const unsigned int zone)
{
  libmesh_assert (in.good());


  // A byte-reverser in case the data is foreign
  Utility::ReverseBytes rb(this->is_foreign());


  for (unsigned int var=0; var<this->n_vars(); var++)
    {

      switch (this->var_type(var))
        {

          // Read a single-precision variable
        case FLOAT:
          {
            std::vector<float> & data = _data[zone][var];

            data.clear();
            data.resize (this->imax(zone)*
                         this->jmax(zone)*
                         this->kmax(zone));

            in.read ((char *) &data[0], LIBMESH_SIZEOF_FLOAT*data.size());

            for (std::size_t i=0; i<data.size(); i++)
              rb(data[i]);

            break;
          }

          // Read a double-precision variable
        case DOUBLE:
          {
            std::vector<double> ddata;
            std::vector<float> & data = _data[zone][var];

            data.clear();
            data.resize (this->imax(zone)*
                         this->jmax(zone)*
                         this->kmax(zone));

            ddata.resize (this->imax(zone)*
                          this->jmax(zone)*
                          this->kmax(zone));

            in.read ((char *) &ddata[0], LIBMESH_SIZEOF_DOUBLE*ddata.size());

            for (std::size_t i=0; i<data.size(); i++)
              data[i] = rb(ddata[i]);

            break;
          }

        default:
          libmesh_error_msg("ERROR: Unsupported data type: " << this->var_type(var));
        }
    }
}



void PltLoader::read_point_data (std::istream & in, const unsigned int zone)
{
  libmesh_assert (in.good());

  // A byte-reverser in case the data is foreign
  Utility::ReverseBytes rb(this->is_foreign());

  // First allocate space
  for (unsigned int var=0; var<this->n_vars(); var++)
    {
      std::vector<float> & data = _data[zone][var];

      data.clear();
      data.reserve (this->imax(zone)*
                    this->jmax(zone)*
                    this->kmax(zone));
    }


  for (unsigned int k=0; k<this->kmax(zone); k++)
    for (unsigned int j=0; j<this->jmax(zone); j++)
      for (unsigned int i=0; i<this->imax(zone); i++)
        for (unsigned int var=0; var<this->n_vars(); var++)
          if (this->var_type(var) == FLOAT)
            {
              float f = 0.;

              libmesh_assert (in.good());

              in.read (buf, LIBMESH_SIZEOF_FLOAT);
              std::memcpy  (&f, buf, LIBMESH_SIZEOF_FLOAT);
              rb(f);

              _data[zone][var].push_back(f);
            }
          else if (this->var_type(var) == DOUBLE)
            {
              double d = 0.;

              libmesh_assert (in.good());

              in.read (buf, LIBMESH_SIZEOF_DOUBLE);
              std::memcpy  (&d, buf, LIBMESH_SIZEOF_DOUBLE);
              rb(d);

              _data[zone][var].push_back(d);
            }
          else
            libmesh_error_msg("ERROR: unsupported data type: " << this->var_type(var));
}



void PltLoader::read_feblock_data (std::istream & in, const unsigned int zone)
{
  libmesh_assert (in.good());

  // A byte-reverser in case the data is foreign
  Utility::ReverseBytes rb(this->is_foreign());

  // Read the variable values at each node.
  for (unsigned int var=0; var<this->n_vars(); var++)
    {
      switch (this->var_type(var))
        {

          // Read a single-precision variable
        case FLOAT:
          {
            std::vector<float> & data = _data[zone][var];

            data.clear();
            data.resize (this->imax(zone));

            in.read ((char *) &data[0], LIBMESH_SIZEOF_FLOAT*data.size());

            for (std::size_t i=0; i<data.size(); i++)
              rb(data[i]);

            break;
          }

          // Read a double-precision variable
        case DOUBLE:
          {
            std::vector<double> ddata;
            std::vector<float> & data = _data[zone][var];

            data.clear();
            data.resize (this->imax(zone));
            ddata.resize (this->imax(zone));

            in.read ((char *) &ddata[0], LIBMESH_SIZEOF_DOUBLE*ddata.size());

            for (std::size_t i=0; i<data.size(); i++)
              data[i] = rb(ddata[i]);

            break;
          }

        default:
          libmesh_error_msg("ERROR: Unsupported data type: " << this->var_type(var));
        }
    }

  // Read the connectivity
  {
    // Get the connectivity repetition flag
    int rep=0;
    in.read ((char *) &rep, LIBMESH_SIZEOF_INT);
    rb(rep);

    if (rep == 1 && this->n_zones() > 1)
      libmesh_error_msg("ERROR:  Repeated connectivity not supported!");

    // Read the connectivity
    else
      {
        libmesh_assert_less (zone, _conn.size());
        libmesh_assert_less (this->kmax(zone), 4);

        _conn[zone].resize (this->jmax(zone)*NNodes[this->kmax(zone)]);

        in.read ((char *) &_conn[zone][0], LIBMESH_SIZEOF_INT*_conn[zone].size());

        for (std::size_t i=0; i<_conn[zone].size(); i++)
          rb(_conn[zone][i]);
      }
  }
}



void PltLoader::read_fepoint_data (std::istream & in, const unsigned int zone)
{
  libmesh_assert (in.good());

  // A byte-reverser in case the data is foreign
  Utility::ReverseBytes rb(this->is_foreign());

  // First allocate space
  for (unsigned int var=0; var<this->n_vars(); var++)
    {
      std::vector<float> & data = _data[zone][var];

      data.clear();
      data.reserve (this->imax(zone));
    }


  for (unsigned int i=0; i<this->imax(zone); i++)
    for (unsigned int var=0; var<this->n_vars(); var++)
      if (this->var_type(var) == FLOAT)
        {
          float f = 0.;

          libmesh_assert (in.good());

          in.read (buf, LIBMESH_SIZEOF_FLOAT);
          std::memcpy  (&f, buf, LIBMESH_SIZEOF_FLOAT);
          rb(f);

          _data[zone][var].push_back(f);
        }
      else if (this->var_type(var) == DOUBLE)
        {
          double d = 0.;

          libmesh_assert (in.good());

          in.read (buf, LIBMESH_SIZEOF_DOUBLE);
          std::memcpy  (&d, buf, LIBMESH_SIZEOF_DOUBLE);
          rb(d);

          _data[zone][var].push_back(d);
        }
      else
        libmesh_error_msg("ERROR: unsupported data type: " << this->var_type(var));

  // Read the connectivity
  {
    // Get the connectivity repetition flag
    int rep=0;

    in.read ((char *) &rep, LIBMESH_SIZEOF_INT);
    rb(rep);

    if (rep == 1)
      libmesh_error_msg("ERROR:  Repeated connectivity not supported!");

    // Read the connectivity
    else
      {
        libmesh_assert_less (zone, _conn.size());
        libmesh_assert_less (this->kmax(zone), 4);

        _conn[zone].resize (this->jmax(zone)*NNodes[this->kmax(zone)]);

        in.read ((char *) &_conn[zone][0], LIBMESH_SIZEOF_INT*_conn[zone].size());

        for (std::size_t i=0; i<_conn[zone].size(); i++)
          rb(_conn[zone][i]);
      }
  }
}

} // namespace libMesh
